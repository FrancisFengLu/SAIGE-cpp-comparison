// loco_engine.cpp
// ---------------------------------------------------------------------------
// Leave-one-chromosome-out (LOCO) for the step-1 null model.
//
// LOCO does NOT re-fit the variance component. tau/theta is estimated once on
// the full-genome GRM by AI-REML; then, holding tau FIXED, only the fixed
// effects (and hence eta / mu / res / V and everything derived from V) are
// re-solved per chromosome.
//
// R reference:
//   R/SAIGE_fitGLMM_fast.R:185-231   Get_Coef_LOCO   (PQL Newton loop)
//   R/SAIGE_fitGLMM_fast.R:714-783   binary/survival driver loop
//   R/SAIGE_fitGLMM_fast.R:996-1033  quantitative driver loop
//   R/Util.R:29-64                   updateChrStartEndIndexVec
//
// The chromosome is excluded SUBTRACTIVELY inside the kernel:
// parallelCrossProd_LOCO computes the full-genome crossproduct, subtracts the
// [startIndex,endIndex] block and renormalizes by
// numberMarker_full - Msub_mafge1perc (SAIGE_step1_fast.cpp).
//
// Output contract: LOCO_FORMAT.md.
// ---------------------------------------------------------------------------

#include "loco_engine.hpp"
#include "null_model_engine.hpp"
#include "saige_ai.hpp"
#include "score.hpp"
#include "SAIGE_step1_fast.hpp"

#include <armadillo>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace saige {

namespace {

inline void loco_log(const std::string& s) {
  std::cout << "[loco] " << s << std::endl;
}

// X is stored row-major in Design (d.X[i*p + j]).
arma::fmat row_major_to_fmat(const std::vector<double>& buf, int n, int p) {
  arma::fmat M(n, std::max(p, 0));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      M(i, j) = static_cast<float>(buf[static_cast<size_t>(i) * p + j]);
  return M;
}

arma::fvec to_fvec(const std::vector<double>& v) {
  arma::fvec o(v.size());
  for (size_t i = 0; i < v.size(); ++i) o(i) = static_cast<float>(v[i]);
  return o;
}

// R: family$linkinv / mu.eta / variance for binomial(logit), matching
// glmm.cpp::irls_binary_build so the LOCO loop stays bit-consistent with the
// full-genome fit's IRLS.
void irls_build_binary(const arma::fvec& eta, const arma::fvec& y,
                       const arma::fvec& offset,
                       arma::fvec& mu, arma::fvec& W, arma::fvec& Y) {
  mu = 1.0f / (1.0f + arma::exp(-eta));
  arma::fvec mu_eta = mu % (1.0f - mu);
  arma::fvec varmu  = mu % (1.0f - mu);
  arma::fvec sqrtW  = mu_eta / arma::sqrt(varmu + 1e-20f);
  W = sqrtW % sqrtW;
  Y = eta - offset + (y - mu) / (mu_eta + 1e-20f);
}

void irls_build_gaussian(const arma::fvec& eta, const arma::fvec& y,
                         const arma::fvec& offset,
                         arma::fvec& mu, arma::fvec& W, arma::fvec& Y) {
  // gaussian(identity): mu = eta, mu.eta = 1, variance = 1
  //   Y = eta - offset + (y - mu)/mu.eta  ==  y - offset
  //   W = (mu.eta/sqrt(var))^2            ==  1
  mu = eta;
  W.set_size(mu.n_elem); W.fill(1.0f);
  Y = eta - offset + (y - mu);
}

// R's Get_Coef convergence test:
//   max(|alpha - alpha0| / (|alpha| + |alpha0| + tol.coef)) < tol.coef
bool coef_converged(const arma::fvec& a, const arma::fvec& a0, float tol_coef) {
  if (a.n_elem == 0) return true;
  float worst = 0.0f;
  for (arma::uword i = 0; i < a.n_elem; ++i) {
    const float d = std::fabs(a(i) - a0(i)) /
                    (std::fabs(a(i)) + std::fabs(a0(i)) + tol_coef);
    worst = std::max(worst, d);
  }
  return worst < tol_coef;
}

void save_pack_files(const std::string& dir,
                     const ScoreNull& sn,
                     const arma::fvec& mu,
                     const std::vector<double>& offset)
{
  fs::create_directories(dir);

  arma::vec mu_d          = arma::conv_to<arma::vec>::from(mu);
  arma::vec res_d         = arma::conv_to<arma::vec>::from(sn.res);
  arma::vec V_d           = arma::conv_to<arma::vec>::from(sn.V);
  arma::vec S_a_d         = arma::conv_to<arma::vec>::from(sn.S_a);
  arma::mat XV_d          = arma::conv_to<arma::mat>::from(sn.XV);
  arma::mat XVX_d         = arma::conv_to<arma::mat>::from(sn.XVX);
  arma::mat XVX_inv_d     = arma::conv_to<arma::mat>::from(sn.XVX_inv);
  arma::mat XXVX_inv_d    = arma::conv_to<arma::mat>::from(sn.XXVX_inv);
  arma::mat XVX_inv_XV_d  = arma::conv_to<arma::mat>::from(sn.XVX_inv_XV);

  mu_d        .save(dir + "/mu.arma",          arma::arma_binary);
  res_d       .save(dir + "/res.arma",         arma::arma_binary);
  V_d         .save(dir + "/V.arma",           arma::arma_binary);
  S_a_d       .save(dir + "/S_a.arma",         arma::arma_binary);
  XV_d        .save(dir + "/XV.arma",          arma::arma_binary);
  XVX_d       .save(dir + "/XVX.arma",         arma::arma_binary);
  XVX_inv_d   .save(dir + "/XVX_inv.arma",     arma::arma_binary);
  XXVX_inv_d  .save(dir + "/XXVX_inv.arma",    arma::arma_binary);
  XVX_inv_XV_d.save(dir + "/XVX_inv_XV.arma",  arma::arma_binary);

  // offset is chromosome-invariant, but LOCO_FORMAT.md lists it in the
  // per-chromosome set so step 2 can point its loader at chr<j>/ wholesale.
  arma::vec off_d = arma::conv_to<arma::vec>::from(offset);
  off_d.save(dir + "/offset.arma", arma::arma_binary);

  // X.arma / y.arma are deliberately NOT written here: they are
  // chromosome-invariant and duplicating X would cost ~570 MB at UKB scale.
}

} // namespace

bool run_loco_batch(const Paths& paths,
                    const FitNullConfig& cfg,
                    const LocoRanges& chr,
                    const Design& design,
                    const std::vector<double>& theta,
                    const std::vector<double>& alpha_fit,
                    const std::vector<double>& offset,
                    const std::vector<double>& eta_full,
                    LocoBatchOut& out)
{
  out.chroms.clear();
  out.obj_noK.clear();

  if (!chr.enabled || chr.start.size() != 22 || chr.end.size() != 22) {
    loco_log("chromosome ranges unavailable -> skipping LOCO");
    return false;
  }
  // R disables LOCO whenever a sparse GRM is used to fit the null model; the
  // LOCO kernels only know how to subtract a marker block out of the dense
  // GRM crossproduct.
  if (get_isUseSparseSigmaforModelFitting()) {
    loco_log("sparse GRM was used to fit the null model -> LOCO not available, skipping");
    return false;
  }

  const int n = design.n;
  if (static_cast<int>(eta_full.size()) != n) {
    loco_log("converged eta unavailable (size " + std::to_string(eta_full.size()) +
             " != n=" + std::to_string(n) + ") -> skipping LOCO");
    return false;
  }
  if (theta.size() < 2) {
    loco_log("theta has <2 components -> skipping LOCO");
    return false;
  }

  const bool is_binary = (cfg.trait == "binary" || cfg.trait == "survival");
  const bool is_quant  = (cfg.trait == "quantitative");
  if (!is_binary && !is_quant) {
    loco_log("unsupported trait '" + cfg.trait + "' -> skipping LOCO");
    return false;
  }
  if (cfg.trait == "survival") {
    // R's survival path threads Lambda0/inC through Get_Coef_LOCO; that branch
    // is not ported. Fail loud rather than emit a wrong chr<j>/ set.
    loco_log("survival LOCO (Lambda0/inC path) is not implemented -> skipping LOCO");
    return false;
  }

  // --- fit-space design (what the GLMM solved on) ---
  const int p_fit = design.p;
  arma::fmat X_fit = (p_fit > 0) ? row_major_to_fmat(design.X, n, p_fit)
                                 : arma::fmat(n, 0);
  // --- export-space design (R's `Xorig` when isCovariateOffset) ---
  const bool have_full = (design.p_full > 0 &&
                          design.X_full.size() ==
                              static_cast<size_t>(n) * static_cast<size_t>(design.p_full));
  const int p_exp = have_full ? design.p_full : p_fit;
  arma::fmat X_exp = have_full ? row_major_to_fmat(design.X_full, n, design.p_full)
                               : X_fit;

  arma::fvec y      = to_fvec(design.y);
  arma::fvec off    = to_fvec(offset.size() == static_cast<size_t>(n)
                                  ? offset : std::vector<double>(n, 0.0));
  arma::fvec tau(2);
  tau(0) = static_cast<float>(theta[0]);
  tau(1) = static_cast<float>(theta[1]);

  const int   maxiter    = std::max(5, cfg.maxiter);
  const int   maxiterPCG = cfg.maxiterPCG > 0 ? cfg.maxiterPCG : 500;
  const float tolPCG     = cfg.tolPCG > 0.0 ? static_cast<float>(cfg.tolPCG) : 1e-5f;
  // R: binary/survival pass tol.coef = tol (the outer AI-REML tol, default
  // 0.02); the quantitative driver omits it, taking Get_Coef_LOCO's default 0.1.
  const float tol_coef   = is_quant ? 0.1f : static_cast<float>(std::max(1e-6, cfg.tol));

  // Carried across chromosomes exactly as R does: `alpha` and `eta` are
  // overwritten each iteration of the j-loop and used to seed the next one.
  arma::fvec alpha = (p_fit > 0 && alpha_fit.size() == static_cast<size_t>(p_fit))
                         ? to_fvec(alpha_fit)
                         : arma::fvec(std::max(p_fit, 0), arma::fill::zeros);
  arma::fvec eta   = to_fvec(eta_full);

  // Populate geno.startIndexVec / endIndexVec and precompute the per-chromosome
  // diagonal of the standardized-genotype crossproduct (R: set_Diagof_StdGeno_LOCO).
  {
    arma::ivec sv(22), ev(22);
    for (int j = 0; j < 22; ++j) { sv(j) = chr.start[j]; ev(j) = chr.end[j]; }
    setStartEndIndexVec(sv, ev);
    set_Diagof_StdGeno_LOCO();
  }

  const std::string model_dir = paths.out_prefix;
  fs::create_directories(model_dir);

  for (int j = 0; j < 22; ++j) {
    const int startIndex = chr.start[j];
    const int endIndex   = chr.end[j];
    if (startIndex < 0 || endIndex < 0) continue;

    loco_log("leave chromosome " + std::to_string(j + 1) + " out  [" +
             std::to_string(startIndex) + "," + std::to_string(endIndex) + "]");

    setStartEndIndex(startIndex, endIndex, j);

    // ---- Get_Coef_LOCO: PQL Newton loop with tau held FIXED ----
    arma::fvec mu, W, Y;
    if (is_binary) irls_build_binary(eta, y, off, mu, W, Y);
    else           irls_build_gaussian(eta, y, off, mu, W, Y);

    arma::fvec alpha0 = alpha;
    int iters = 0;
    for (int it = 0; it < maxiter; ++it) {
      ++iters;
      CoefficientsOut re = getCoefficients_LOCO_cpp(Y, X_fit, W, tau,
                                                    maxiterPCG, tolPCG);
      alpha = re.alpha;
      eta   = re.eta + off;

      if (is_binary) irls_build_binary(eta, y, off, mu, W, Y);
      else           irls_build_gaussian(eta, y, off, mu, W, Y);

      if (coef_converged(alpha, alpha0, tol_coef)) break;
      alpha0 = alpha;
    }
    loco_log("  chr" + std::to_string(j + 1) + ": Get_Coef_LOCO iterations = " +
             std::to_string(iters) + "  mean(mu)=" +
             std::to_string(arma::mean(mu)));

    // ---- obj.noK on the (export-space) X and the new mu ----
    ScoreNull sn;
    if (is_binary) {
      sn = build_score_null_binary(X_exp, y, mu);
    } else {
      sn = build_score_null_quant(X_exp, y, mu, 1.0f / tau(0));
    }

    const std::string dir = model_dir + "/chr" + std::to_string(j + 1);
    save_pack_files(dir, sn, mu, offset);

    out.chroms.push_back(j + 1);
    out.obj_noK.push_back(to_pack(sn, X_exp, y, mu, cfg.trait));
    (void)p_exp;
  }

  if (out.chroms.empty()) {
    loco_log("no autosome produced LOCO output");
    return false;
  }
  loco_log("wrote chr<j>/ for " + std::to_string(out.chroms.size()) + " autosomes");
  return true;
}

void register_default_loco_batch() {
  register_loco_batch(&run_loco_batch);
}

} // namespace saige
