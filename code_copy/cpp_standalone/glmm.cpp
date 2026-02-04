#include "glmm.hpp"
#include "saige_ai.hpp"           // Armadillo/PCG wrappers (pure C++ structs)
#include "score.hpp"   // build_score_null_binary/quant
#include "SAIGE_step1_fast.hpp"
#include <armadillo>
#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <string>

namespace saige {

// Checkpoint output directory
static const std::string CP_DIR = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/checkpoints";

// Helper to save vector to CSV
static void save_vec_csv(const arma::fvec& v, const std::string& path) {
  std::ofstream f(path);
  for (arma::uword i = 0; i < v.n_elem; ++i) f << v[i] << "\n";
}

// Helper to save matrix to CSV
static void save_mat_csv(const arma::fmat& m, const std::string& path) {
  std::ofstream f(path);
  f << std::setprecision(10);
  for (arma::uword i = 0; i < m.n_rows; ++i) {
    for (arma::uword j = 0; j < m.n_cols; ++j) {
      if (j > 0) f << ",";
      f << m(i, j);
    }
    f << "\n";
  }
}

// ---------- small helpers ----------

static inline arma::fvec to_fvec(const std::vector<double>& v) {
  arma::fvec out(v.size());
  for (size_t i = 0; i < v.size(); ++i) out[static_cast<arma::uword>(i)] = static_cast<float>(v[i]);
  return out;
}

static inline arma::fvec map_y(const Design& d) { return to_fvec(d.y); }

static inline arma::fvec map_offset(const Design& d, const std::vector<double>& offset_in) {
  if (!offset_in.empty()) return to_fvec(offset_in);
  arma::fvec out(d.n, arma::fill::zeros);
  return out;
}

static inline arma::fmat map_X_row_major_to_fmat(const Design& d) {
  const int n = d.n, p = d.p;
  arma::fmat X(n, p);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      X(i, j) = static_cast<float>(d.X[static_cast<size_t>(i)*p + j]);
  return X;
}

static inline arma::fvec map_beta_init(const Design& d, const std::vector<double>& init) {
  if (init.empty() || d.p == 0) return arma::fvec(d.p, arma::fill::zeros);
  arma::fvec out(d.p);
  for (int j = 0; j < d.p; ++j) out[j] = static_cast<float>(init[static_cast<size_t>(j)]);
  return out;
}

// relative change helper (original)
template <typename Vec>
static inline double rel_change_inf(const Vec& a, const Vec& b, double eps) {
  double num = 0.0, den = 0.0;
  for (size_t i = 0; i < a.n_elem; ++i) {
    num = std::max(num, std::fabs(static_cast<double>(a[i] - b[i])));
    den = std::max(den, std::fabs(static_cast<double>(a[i])) + std::fabs(static_cast<double>(b[i])) + eps);
  }
  return num / den;
}

// R-style relative change: max(abs(a - b) / (abs(a) + abs(b) + tol))
// This matches R SAIGE line 387: max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol
template <typename Vec>
static inline double rel_change_R_style(const Vec& a, const Vec& b, double tol) {
  double max_rc = 0.0;
  for (size_t i = 0; i < a.n_elem; ++i) {
    double ai = static_cast<double>(a[i]);
    double bi = static_cast<double>(b[i]);
    double rc = std::fabs(ai - bi) / (std::fabs(ai) + std::fabs(bi) + tol);
    max_rc = std::max(max_rc, rc);
  }
  return max_rc;
}

static inline double clamp_nonneg(double x) { return x < 0.0 ? 0.0 : x; }

// ---------- score-null packing into FitNullResult ----------

static inline std::vector<double> flatten_rowmajor(const arma::fmat& M) {
  std::vector<double> v(M.n_rows * M.n_cols);
  for (arma::uword i = 0; i < M.n_rows; ++i)
    for (arma::uword j = 0; j < M.n_cols; ++j)
      v[i * M.n_cols + j] = static_cast<double>(M(i, j));
  return v;
}

static void stash_score_null_into(FitNullResult& out,
                                  const saige::ScoreNull& s,
                                  int n, int p)
{
  out.obj_noK.n   = n;
  out.obj_noK.p   = p;
  out.obj_noK.V   = std::vector<double>(s.V.begin(),   s.V.end());
  out.obj_noK.S_a = std::vector<double>(s.S_a.begin(), s.S_a.end());
  out.obj_noK.XVX     = flatten_rowmajor(s.XVX);
  out.obj_noK.XVX_inv = flatten_rowmajor(s.XVX_inv);
}

// ---------- OPTIONAL: JSON exporter that supports LOCO as well ----------
// Call this from main() AFTER fit_null(...) if you want files on disk.
// It writes baseline obj_noK and, if populated, a "loco" array with per-chr packs.
static void export_score_null_json(const Paths& paths, const FitNullResult& res) {
  const std::string out_path = paths.out_prefix + ".obj_noK.json";
  std::ofstream os(out_path, std::ios::out | std::ios::trunc);
  if (!os) return;

  auto dump_pack = [&](const ScoreNullPack& p, std::ostream& s) {
    s << "{";
    s << "\"n\":" << p.n << ",\"p\":" << p.p << ",";
    s << "\"V\":[";
    for (size_t i=0;i<p.V.size();++i){ if(i) s<<","; s<<p.V[i]; }
    s << "],\"S_a\":[";
    for (size_t i=0;i<p.S_a.size();++i){ if(i) s<<","; s<<p.S_a[i]; }
    s << "],\"XVX\":[";
    for (size_t i=0;i<p.XVX.size();++i){ if(i) s<<","; s<<p.XVX[i]; }
    s << "],\"XVX_inv\":[";
    for (size_t i=0;i<p.XVX_inv.size();++i){ if(i) s<<","; s<<p.XVX_inv[i]; }
    s << "]}";
  };

  os << std::setprecision(10);
  os << "{\n  \"baseline\": ";
  dump_pack(res.obj_noK, os);

  // LOCO section (present if loco_obj_noK has any non-empty entries)
  bool any_loco = false;
  for (const auto& pk : res.loco_obj_noK) if (pk.n > 0) { any_loco = true; break; }
  if (any_loco) {
    os << ",\n  \"loco\": [\n";
    for (size_t c = 0; c < res.loco_obj_noK.size(); ++c) {
      if (res.loco_obj_noK[c].n == 0) continue;
      os << "    {\"chrom\":" << (c+1) << ",\"pack\":";
      dump_pack(res.loco_obj_noK[c], os);
      os << "}";
      // comma handling: emit comma if next non-empty exists
      size_t k = c + 1;
      while (k < res.loco_obj_noK.size() && res.loco_obj_noK[k].n == 0) ++k;
      if (k < res.loco_obj_noK.size()) os << ",";
      os << "\n";
    }
    os << "  ]\n";
  } else {
    os << "\n";
  }
  os << "}\n";
}

// ---------- family pieces (IRLS scaffolding) ----------

static inline void irls_binary_build(const arma::fvec& eta,
                                     const arma::fvec& y,
                                     const arma::fvec& offset,
                                     arma::fvec& mu,
                                     arma::fvec& mu_eta,
                                     arma::fvec& W,
                                     arma::fvec& Y) {
  mu = 1.0f / (1.0f + arma::exp(-eta));
  mu_eta = mu % (1.0f - mu);
  arma::fvec varmu = mu % (1.0f - mu);
  arma::fvec sqrtW = mu_eta / arma::sqrt(varmu + 1e-20f);
  W = sqrtW % sqrtW;                // == varmu (logistic)
  Y = eta - offset + (y - mu) / (mu_eta + 1e-20f);
}

static inline void irls_gaussian_build(const arma::fvec& eta,
                                       const arma::fvec& y,
                                       const arma::fvec& offset,
                                       arma::fvec& mu,
                                       arma::fvec& mu_eta,
                                       arma::fvec& W,
                                       arma::fvec& Y) {
  mu = eta;
  mu_eta.set_size(mu.n_elem); mu_eta.fill(1.0f);
  W.set_size(mu.n_elem);      W.fill(1.0f);
  Y = eta - offset + (y - mu) / mu_eta;   // = y - offset
}

// ---------- Binary solver (AI-REML on tau[1]; tau[0] fixed=1) ----------

inline arma::fvec sigmoid_stable_f(const arma::fvec& eta_f) {
  arma::vec eta = arma::conv_to<arma::vec>::from(eta_f);
  // clamp to avoid exp overflow/underflow in double
  eta = arma::clamp(eta, -40.0, 40.0);
  arma::vec mu = 1.0 / (1.0 + arma::exp(-eta));
  return arma::conv_to<arma::fvec>::from(mu);
}

// robust relative change: max_i |a-b| / (|a|+|b|+tol)
inline double rel_change_tau_Rstyle(const arma::fvec& a,
                                    const arma::fvec& b,
                                    float tol) {
  arma::fvec num = arma::abs(a - b);
  arma::fvec den = arma::abs(a) + arma::abs(b) + tol;
  return static_cast<double>(num.max() / den.max());
}


// FitNullResult binary_glmm_solver(const Paths& paths,
//                                  const FitNullConfig& cfg,
//                                  const Design& d,
//                                  const std::vector<double>& offset_in,
//                                  const std::vector<double>& beta_init_in)
// {
//   const int n = d.n, p = d.p;

//   // Map inputs
//   arma::fmat X = (p > 0) ? map_X_row_major_to_fmat(d) : arma::fmat(n, 0);
//   arma::fvec y = map_y(d);
//   arma::fvec offset = map_offset(d, offset_in);
//   arma::fvec beta_init = map_beta_init(d, beta_init_in);

//   // Basic shape assertions (fail fast if Design/inputs are inconsistent)
//   if ((int)X.n_rows != n) throw std::runtime_error("binary_glmm_solver: X.n_rows != n");
//   if ((int)X.n_cols != p) throw std::runtime_error("binary_glmm_solver: X.n_cols != p");
//   if ((int)y.n_elem != n) throw std::runtime_error("binary_glmm_solver: y.n_elem != n");
//   if ((int)offset.n_elem != n) throw std::runtime_error("binary_glmm_solver: offset.n_elem != n");
//   if (p > 0 && (int)beta_init.n_elem != p)
//     throw std::runtime_error("binary_glmm_solver: beta_init.n_elem != p");

//   // Initial linear predictor
//   arma::fvec eta = (p > 0) ? (X * beta_init + offset) : offset;

//   // Tuning/limits
//   const int   maxiter    = std::max(5, cfg.maxiter);
//   const float tol_coef   = static_cast<float>(std::max(1e-4, cfg.tol));
//   const int   maxiterPCG = cfg.maxiterPCG > 0 ? cfg.maxiterPCG : 500;
//   const float tolPCG     = cfg.tolPCG > 0.0 ? static_cast<float>(cfg.tolPCG) : 1e-5f;
//   const int   nrun       = cfg.nrun > 0 ? cfg.nrun : 30;
//   const float trace_cut  = cfg.traceCVcutoff > 0.0 ? static_cast<float>(cfg.traceCVcutoff) : 0.1f;

//   // Important numerical floors/caps
//   const float W_FLOOR    = 1e-6f;    // floor for IRLS weights mu*(1-mu)
//   const double AI_FLOOR  = 1e-10;    // floor for AI denominator
//   const double TAU1_CAP  = 1e6;      // conservative cap to avoid runaways

//   // Variance components (tau[0] is residual/dispersion; tau[1] random-effect VC)
//   arma::fvec tau(2); tau[0] = 1.0f; tau[1] = 0.5f;

//   // Work vectors
//   arma::fvec mu(n, arma::fill::zeros);
//   arma::fvec W(n,  arma::fill::zeros);
//   arma::fvec Y(n,  arma::fill::zeros);

//   arma::fvec alpha_prev(p, arma::fill::zeros);
//   arma::fvec tau_prev   = tau;

//   for (int it = 0; it < maxiter; ++it) {
//     // -------- (A) Update fixed effects and IRLS quantities for CURRENT tau --------
//     // We call the coefficient/PCG routine first so Sigma_iY/X/cov reflect current tau.
//     // getCoefficients_cpp should solve (X' Σ^-1 X) alpha = X' Σ^-1 Y (or equivalent)
//     // and return alpha, Sigma_iY (n), Sigma_iX (n×p), cov (p×p).
//     auto coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);
//     // ^ Ensure that getCoefficients_cpp fills:
//     //   - coef.alpha (length p, possibly p==0)
//     //   - coef.Sigma_iY (length n)
//     //   - coef.Sigma_iX (n×p)
//     //   - coef.cov     (p×p)

//     if (p > 0 && (int)coef.alpha.n_elem != p)
//       throw std::runtime_error("getCoefficients_cpp: alpha.n_elem != p");
//     if ((int)coef.Sigma_iY.n_elem != n)
//       throw std::runtime_error("getCoefficients_cpp: Sigma_iY.n_elem != n");
//     if ((int)coef.Sigma_iX.n_rows != n || (int)coef.Sigma_iX.n_cols != p)
//       throw std::runtime_error("getCoefficients_cpp: Sigma_iX not n×p");
//     if ((int)coef.cov.n_rows != p || (int)coef.cov.n_cols != p)
//       throw std::runtime_error("getCoefficients_cpp: cov not p×p");

//     arma::fvec alpha = (p > 0) ? coef.alpha : arma::fvec(); // allow p==0
//     eta = (p > 0) ? (X * alpha + offset) : offset;

//     // Stable logistic and floored weights
//     mu = sigmoid_stable_f(eta);
//     W  = mu % (1.0f - mu);
//     W  = arma::clamp(W, W_FLOOR, std::numeric_limits<float>::infinity());

//     // Working response Y = eta - offset + (y - mu) / mu.eta; for logit, mu.eta = mu*(1-mu)
//     // Use the same floor in the denominator to avoid blow-ups.
//     arma::fvec mu_eta = W; // already floored mu*(1-mu)
//     Y = eta - offset + (y - mu) / mu_eta;

//     // -------- (B) Variance-component (AI/PCG) update using FRESH Y/W/etc. --------
//     auto ai = getAIScore_cpp(Y, X, W, tau,
//                              coef.Sigma_iY, coef.Sigma_iX, coef.cov,
//                              nrun, maxiterPCG, tolPCG, trace_cut);

//     // Accumulate in double and apply damping to avoid giant Newton steps
//     const double AI    = std::max(AI_FLOOR, static_cast<double>(ai.AI));
//     const double score = static_cast<double>(ai.YPAPY - ai.Trace);
//     double delta = score / AI;

//     // Simple, effective damping: shrink big steps to (0,1] scale
//     double step  = 1.0 / (1.0 + std::fabs(delta));
//     double tau1_candidate = static_cast<double>(tau[1]) + step * delta;

//     // Enforce nonnegativity + cap
//     double tau1_new = std::min(std::max(0.0, tau1_candidate), TAU1_CAP);

//     tau_prev = tau;
//     tau[1]   = static_cast<float>(tau1_new);

//     // Convergence check (use R-style metric to avoid division by ~0)
//     double rc_tau   = rel_change_tau_Rstyle(tau, tau_prev, tol_coef);
//     double rc_alpha = (p > 0) ? rel_change_tau_Rstyle(alpha, alpha_prev, tol_coef) : 0.0;

//     if (std::max(rc_tau, rc_alpha) < tol_coef) {
//       // finalize and stash score-null
//       arma::fvec mu_final = sigmoid_stable_f(eta); // use stable sigmoid for final mu
//       auto sn = saige::build_score_null_binary(X, y, mu_final);

//       FitNullResult out;
//       out.alpha  = (p > 0)
//                    ? std::vector<double>(alpha.begin(), alpha.end())
//                    : std::vector<double>();
//       out.theta  = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
//       out.offset = offset_in;

//       stash_score_null_into(out, sn, n, p);
//       // export_score_null_json(paths, out);  // optional
//       return out;
//     }

//     alpha_prev = (p > 0) ? alpha : arma::fvec(); // keep previous for next iteration
//   }

//   // -------- fallthrough: not converged within maxiter; finalize safely --------
//   {
//     // Recompute with last eta
//     mu = sigmoid_stable_f(eta);
//     W  = arma::clamp(mu % (1.0f - mu), W_FLOOR, std::numeric_limits<float>::infinity());
//     arma::fvec mu_eta = W;
//     Y = eta - offset + (y - mu) / mu_eta;

//     auto coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);

//     arma::fvec mu_final = sigmoid_stable_f(eta);
//     auto sn = saige::build_score_null_binary(X, y, mu_final);

//     FitNullResult out;
//     out.alpha  = (p > 0)
//                  ? std::vector<double>(coef.alpha.begin(), coef.alpha.end())
//                  : std::vector<double>();
//     out.theta  = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
//     out.offset = offset_in;

//     stash_score_null_into(out, sn, n, p);
//     export_score_null_json(paths, out);
//     return out;
//   }
// }
FitNullResult binary_glmm_solver(const Paths& paths,
                                 const FitNullConfig& cfg,
                                 const Design& d,
                                 const std::vector<double>& offset_in,
                                 const std::vector<double>& beta_init_in)
{
  const int n = d.n, p = d.p;

  arma::fmat X = (p > 0) ? map_X_row_major_to_fmat(d) : arma::fmat(n, 0);
  arma::fvec y = map_y(d);
  arma::fvec offset = map_offset(d, offset_in);
  arma::fvec beta_init = map_beta_init(d, beta_init_in);
  arma::fvec eta = (p > 0) ? (X * beta_init + offset) : offset;

// std::cout << "shape = " << X.n_rows << " x " << X.n_cols << "\n";
// std::cout << "shape = " << y.n_rows << " y " << y.n_cols << "\n";
// std::cout << "shape = " << offset.n_rows << " offset " << offset.n_cols << "\n";
// std::cout << "shape = " << eta.n_rows << " eta " << eta.n_cols << "\n";
// std::cout << "shape = " << beta_init.n_rows << " beta_init " << beta_init.n_cols << std::endl;;

  const int   maxiter    = std::max(5, cfg.maxiter);
  const float tol_coef   = static_cast<float>(std::max(1e-6, cfg.tol));
  const int   maxiterPCG = cfg.maxiterPCG > 0 ? cfg.maxiterPCG : 500;
  const float tolPCG     = cfg.tolPCG > 0.0 ? static_cast<float>(cfg.tolPCG) : 1e-5f;
  const int   nrun       = cfg.nrun > 0 ? cfg.nrun : 30;
  const float trace_cut  = cfg.traceCVcutoff > 0.0 ? static_cast<float>(cfg.traceCVcutoff) : 0.1f;

  // Using tau=0.0 to match R's initial value
  arma::fvec tau(2); tau[0] = 1.0f; tau[1] = 0.0f;

  arma::fvec mu, mu_eta, W, Y;
  arma::fvec alpha_prev(p, arma::fill::zeros);
  arma::fvec tau_prev   = tau;

auto check_dims = [&](const char* where){
  if ((int)X.n_rows != n) throw std::runtime_error(std::string(where)+": X.n_rows!=n");
  if ((int)X.n_cols != p) throw std::runtime_error(std::string(where)+": X.n_cols!=p");
  if ((int)y.n_elem != n) throw std::runtime_error(std::string(where)+": y.n_elem!=n");
  if ((int)offset.n_elem != n) throw std::runtime_error(std::string(where)+": offset.n_elem!=n");
  if (p>0 && (int)beta_init.n_elem != p) throw std::runtime_error(std::string(where)+": beta_init.n_elem!=p");
};
check_dims("inputs");

  // ===== DEBUG: Print initial eta before GLMM loop =====
  std::cout << "\n===== C++ GLMM Solver Initial Values =====" << std::endl;
  std::cout << "beta_init: [";
  for (size_t i = 0; i < beta_init.n_elem; ++i) {
    std::cout << beta_init[i];
    if (i < beta_init.n_elem - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "eta[0:5]: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;
  std::cout << "y[0:5]: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << std::endl;
  std::cout << "offset[0:5]: " << offset[0] << " " << offset[1] << " " << offset[2] << " " << offset[3] << " " << offset[4] << std::endl;
  std::cout << "==========================================\n" << std::endl;

  // ============ CHECKPOINT 2: Initial values before GLMM loop ============
  std::cout << "\n=== C++ CHECKPOINT 2: Initial values ===" << std::endl;
  std::cout << "CP2: tau = " << tau[0] << " " << tau[1] << std::endl;
  std::cout << "CP2: beta_init[0:3] = ";
  for (arma::uword i = 0; i < std::min((arma::uword)3, beta_init.n_elem); ++i) std::cout << beta_init[i] << " ";
  std::cout << std::endl;
  std::cout << "CP2: eta[0:10] = ";
  for (int i = 0; i < std::min(10, n); ++i) std::cout << eta[i] << " ";
  std::cout << std::endl;
  std::cout << "CP2: n = " << n << ", p = " << p << std::endl;
  // Save to files
  {
    std::ofstream f(CP_DIR + "/CPP_CP2_initial.txt");
    f << "tau: " << tau[0] << " " << tau[1] << "\n";
    f << "n: " << n << "\n";
    f << "p: " << p << "\n";
    f << "beta_init: ";
    for (size_t i = 0; i < beta_init.n_elem; ++i) f << beta_init[i] << " ";
    f << "\n";
  }
  save_vec_csv(eta, CP_DIR + "/CPP_CP2_eta.csv");
  save_vec_csv(y, CP_DIR + "/CPP_CP2_y.csv");
  std::cout << "CP2 saved to: " << CP_DIR << "/CPP_CP2_*.csv/txt" << std::endl;
  // ============ END CHECKPOINT 2 ============

  // ============ DEBUG MODE: Stop after first real tau update ============
  const bool STOP_AFTER_FIRST_TAU_UPDATE = false;
  int real_tau_updates = 0;

  // Track alpha across outer iterations (R's alpha0 in Get_Coef)
  // Initialize to beta_init (the GLM alpha) to match R's behavior
  arma::fvec alpha_outer_prev = arma::conv_to<arma::fvec>::from(beta_init);
  CoefficientsOut coef;

  for (int it = 0; it < maxiter; ++it) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== C++ OUTER ITERATION " << it << " START ===" << std::endl;
    std::cout << "tau at start: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "eta[0:5] at start: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    // ======================================================================
    // Inner IRLS loop to converge alpha (matching R's Get_Coef behavior)
    // R compares alpha against alpha0 from BEFORE the inner loop starts
    // ======================================================================
    arma::fvec alpha0_outer = alpha_outer_prev;  // Compare against previous outer iteration's alpha

    for (int inner_it = 0; inner_it < maxiter; ++inner_it) {
      std::cout << "\n--- Inner IRLS iteration " << inner_it << " ---" << std::endl;

      // Build working response Y and weights W from current eta
      irls_binary_build(eta, y, offset, mu, mu_eta, W, Y);

      std::cout << "After irls_binary_build:" << std::endl;
      std::cout << "  mu[0:5]: " << mu[0] << " " << mu[1] << " " << mu[2] << " " << mu[3] << " " << mu[4] << std::endl;
      std::cout << "  W[0:5]: " << W[0] << " " << W[1] << " " << W[2] << " " << W[3] << " " << W[4] << std::endl;
      std::cout << "  Y[0:5]: " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << std::endl;
      std::cout << "  |Y|: " << arma::norm(Y) << ", |W|: " << arma::norm(W) << std::endl;

      // Solve for alpha: (X' Sigma^-1 X)^-1 X' Sigma^-1 Y
      coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);

      std::cout << "After getCoefficients_cpp:" << std::endl;
      std::cout << "  alpha: [" << coef.alpha[0] << ", " << coef.alpha[1] << ", " << coef.alpha[2] << "]" << std::endl;
      std::cout << "  coef.eta[0:5]: " << coef.eta[0] << " " << coef.eta[1] << " " << coef.eta[2] << " " << coef.eta[3] << " " << coef.eta[4] << std::endl;
      std::cout << "  Sigma_iY[0:5]: " << coef.Sigma_iY[0] << " " << coef.Sigma_iY[1] << " " << coef.Sigma_iY[2] << " " << coef.Sigma_iY[3] << " " << coef.Sigma_iY[4] << std::endl;

      // Update eta from the coefficient solve
      eta = coef.eta + offset;
      std::cout << "  eta after update [0:5]: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;

      // Check convergence of alpha against alpha0_outer (R's behavior)
      // R compares new alpha against the alpha from BEFORE the inner loop started
      double rc_alpha_inner = (p > 0) ? rel_change_R_style(coef.alpha, alpha0_outer, tol_coef) : 0.0;
      std::cout << "  rc_alpha (vs alpha0_outer): " << rc_alpha_inner << " (tol=" << tol_coef << ")" << std::endl;

      if (rc_alpha_inner < tol_coef) {
        std::cout << "  [Inner IRLS] CONVERGED after " << (inner_it + 1) << " iterations" << std::endl;
        break;
      }

      // Update alpha0_outer for next inner iteration (R also does alpha0 = alpha)
      alpha0_outer = coef.alpha;
    }

    // Save alpha for next outer iteration's comparison
    alpha_outer_prev = coef.alpha;

    // Final IRLS build with converged eta for this tau iteration
    std::cout << "\n--- Final IRLS build after inner loop convergence ---" << std::endl;
    irls_binary_build(eta, y, offset, mu, mu_eta, W, Y);
    std::cout << "Final W[0:5]: " << W[0] << " " << W[1] << " " << W[2] << " " << W[3] << " " << W[4] << std::endl;
    std::cout << "Final Y[0:5]: " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << std::endl;
    std::cout << "Final |Y|: " << arma::norm(Y) << ", |W|: " << arma::norm(W) << std::endl;

  auto check_coef = [&](const char* where){
  if (p>0 && (int)coef.alpha.n_elem != p) throw std::runtime_error(std::string(where)+": alpha.n_elem!=p");
  if ((int)coef.Sigma_iY.n_elem != n) throw std::runtime_error(std::string(where)+": Sigma_iY.n_elem!=n");
  if ((int)coef.Sigma_iX.n_rows != n) throw std::runtime_error(std::string(where)+": Sigma_iX.n_rows!=n");
  if ((int)coef.Sigma_iX.n_cols != p) throw std::runtime_error(std::string(where)+": Sigma_iX.n_cols!=p");
  if ((int)coef.cov.n_rows != p || (int)coef.cov.n_cols != p) throw std::runtime_error(std::string(where)+": cov not p×p");
};
check_coef("getCoefficients_cpp");

    std::cout << "Finished getCoefficients_cpp"  << std::endl;

    // ============ CHECKPOINT 3: After first getCoefficients (it==0 only) ============
    if (it == 0) {
      std::cout << "\n=== C++ CHECKPOINT 3: After first getCoefficients ===" << std::endl;
      std::cout << "CP3: alpha = ";
      for (size_t i = 0; i < coef.alpha.n_elem; ++i) std::cout << coef.alpha[i] << " ";
      std::cout << std::endl;
      std::cout << "CP3: Sigma_iY[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << coef.Sigma_iY[i] << " ";
      std::cout << std::endl;
      std::cout << "CP3: eta[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << coef.eta[i] << " ";
      std::cout << std::endl;
      std::cout << "CP3: Y[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << Y[i] << " ";
      std::cout << std::endl;
      std::cout << "CP3: W[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << W[i] << " ";
      std::cout << std::endl;
      std::cout << "CP3: cov diag = ";
      for (int i = 0; i < p; ++i) std::cout << coef.cov(i,i) << " ";
      std::cout << std::endl;
      // Save to files
      {
        std::ofstream f(CP_DIR + "/CPP_CP3_coef.txt");
        f << "alpha: ";
        for (size_t i = 0; i < coef.alpha.n_elem; ++i) f << coef.alpha[i] << " ";
        f << "\n";
        f << "cov_diag: ";
        for (int i = 0; i < p; ++i) f << coef.cov(i,i) << " ";
        f << "\n";
      }
      save_vec_csv(coef.Sigma_iY, CP_DIR + "/CPP_CP3_Sigma_iY.csv");
      save_mat_csv(coef.Sigma_iX, CP_DIR + "/CPP_CP3_Sigma_iX.csv");
      save_vec_csv(coef.eta, CP_DIR + "/CPP_CP3_eta.csv");
      save_vec_csv(Y, CP_DIR + "/CPP_CP3_Y.csv");
      save_vec_csv(W, CP_DIR + "/CPP_CP3_W.csv");
      std::cout << "CP3 saved to: " << CP_DIR << "/CPP_CP3_*.csv/txt" << std::endl;
    }
    // ============ END CHECKPOINT 3 ============

    auto ai = getAIScore_cpp(Y, X, W, tau, coef.Sigma_iY, coef.Sigma_iX, coef.cov,
                             nrun, maxiterPCG, tolPCG, trace_cut);

if ((int)Y.n_elem != n) throw std::runtime_error("Y.n_elem!=n before AI");
if ((int)W.n_elem != n) throw std::runtime_error("W.n_elem!=n before AI");

    std::cout << "Finished getAI"  << std::endl;

    // ============ CHECKPOINT 4: After first getAIScore (it==0 only) ============
    if (it == 0) {
      std::cout << "\n=== C++ CHECKPOINT 4: After first getAIScore ===" << std::endl;
      std::cout << "CP4: YPAPY = " << ai.YPAPY << std::endl;
      std::cout << "CP4: Trace = " << ai.Trace << std::endl;
      std::cout << "CP4: AI = " << ai.AI << std::endl;
      std::cout << "CP4: PY[0:10] = ";
      for (int i = 0; i < std::min(10, (int)ai.PY.n_elem); ++i) std::cout << ai.PY[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4: Score (YPAPY-Trace) = " << (ai.YPAPY - ai.Trace) << std::endl;
      // Save to files
      {
        std::ofstream f(CP_DIR + "/CPP_CP4_AIScore.txt");
        f << std::setprecision(10);
        f << "YPAPY: " << ai.YPAPY << "\n";
        f << "Trace: " << ai.Trace << "\n";
        f << "AI: " << ai.AI << "\n";
        f << "Score: " << (ai.YPAPY - ai.Trace) << "\n";
      }
      save_vec_csv(ai.PY, CP_DIR + "/CPP_CP4_PY.csv");
      save_vec_csv(ai.APY, CP_DIR + "/CPP_CP4_APY.csv");
      std::cout << "CP4 saved to: " << CP_DIR << "/CPP_CP4_*.csv/txt" << std::endl;
    }
    // ============ END CHECKPOINT 4 ============

    // ============ CHECKPOINT 4b: Iteration 3 (it==2, R's i==2) intermediate values ============
    // Note: C++ it==2 corresponds to R's i==2 because:
    // - C++ it==0: conservative update (tau stays at 0)
    // - C++ it==1: first standard update (tau becomes 0.281)
    // - C++ it==2: second standard update (R's iteration 2 when tau starts at 0.281)
    if (it == 2) {
      std::cout << "\n=== C++ CHECKPOINT 4b: Iteration 2 (it==1) intermediate values ===" << std::endl;
      std::cout << "CP4b: tau_at_start = [" << tau[0] << ", " << tau[1] << "]" << std::endl;
      std::cout << "CP4b: W[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << W[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: Y[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << Y[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: mu[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << mu[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: eta[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << eta[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: Sigma_iY[0:10] = ";
      for (int i = 0; i < std::min(10, n); ++i) std::cout << coef.Sigma_iY[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: PY[0:10] = ";
      for (int i = 0; i < std::min(10, (int)ai.PY.n_elem); ++i) std::cout << ai.PY[i] << " ";
      std::cout << std::endl;
      std::cout << "CP4b: alpha = [";
      for (size_t i = 0; i < coef.alpha.n_elem; ++i) {
        std::cout << coef.alpha[i];
        if (i < coef.alpha.n_elem - 1) std::cout << ", ";
      }
      std::cout << "]" << std::endl;
      std::cout << "CP4b: YPAPY = " << ai.YPAPY << std::endl;
      std::cout << "CP4b: Trace = " << ai.Trace << std::endl;
      std::cout << "CP4b: AI = " << ai.AI << std::endl;
      std::cout << "CP4b: Score = " << (ai.YPAPY - ai.Trace) << std::endl;

      // Save to files for detailed comparison
      {
        std::ofstream f(CP_DIR + "/CPP_CP4b_iter2.txt");
        f << std::setprecision(12);
        f << "tau: " << tau[0] << " " << tau[1] << "\n";
        f << "YPAPY: " << ai.YPAPY << "\n";
        f << "Trace: " << ai.Trace << "\n";
        f << "AI: " << ai.AI << "\n";
        f << "Score: " << (ai.YPAPY - ai.Trace) << "\n";
        f << "alpha: ";
        for (size_t i = 0; i < coef.alpha.n_elem; ++i) f << coef.alpha[i] << " ";
        f << "\n";
      }
      save_vec_csv(W, CP_DIR + "/CPP_CP4b_W.csv");
      save_vec_csv(Y, CP_DIR + "/CPP_CP4b_Y.csv");
      save_vec_csv(mu, CP_DIR + "/CPP_CP4b_mu.csv");
      save_vec_csv(eta, CP_DIR + "/CPP_CP4b_eta.csv");
      save_vec_csv(coef.Sigma_iY, CP_DIR + "/CPP_CP4b_Sigma_iY.csv");
      save_vec_csv(ai.PY, CP_DIR + "/CPP_CP4b_PY.csv");
      save_vec_csv(ai.APY, CP_DIR + "/CPP_CP4b_APY.csv");
      std::cout << "CP4b saved to: " << CP_DIR << "/CPP_CP4b_*.csv/txt" << std::endl;
    }
    // ============ END CHECKPOINT 4b ============

    double score = static_cast<double>(ai.YPAPY - ai.Trace);
    double AI    = std::max(1e-12, static_cast<double>(ai.AI));

    // R SAIGE uses different update formulas:
    // - First iteration (it==0): conservative update tau_new = tau_old + tau_old^2 * score / n
    // - Subsequent iterations: standard AI-REML tau_new = tau_old + score / AI
    double tau0_val = static_cast<double>(tau[1]);
    double Dtau;
    double tau1_new;

    if (it == 0) {
      // Conservative update for first iteration (R SAIGE line 345)
      Dtau = tau0_val * tau0_val * score / static_cast<double>(n);
      tau1_new = tau0_val + Dtau;
      std::cout << "[CONSERVATIVE UPDATE] tau^2 * score / n = " << tau0_val << "^2 * " << score << " / " << n << " = " << Dtau << std::endl;
    } else {
      // Standard AI-REML update (R SAIGE fitglmmaiRPCG)
      Dtau = score / AI;
      tau1_new = tau0_val + Dtau;
      std::cout << "[STANDARD AI-REML] score / AI = " << score << " / " << AI << " = " << Dtau << std::endl;
    }

    // Step halving when tau goes negative (R SAIGE behavior)
    double step = 1.0;
    while (tau1_new < 0.0 && step > 1e-10) {
        step *= 0.5;
        tau1_new = tau0_val + step * Dtau;
    }
    // Final clamp to ensure non-negative
    tau1_new = std::max(0.0, tau1_new);

    // ===== DETAILED DEBUG OUTPUT FOR ITERATION COMPARISON =====
    std::cout << "\n========== C++ ITERATION " << it << " ==========" << std::endl;
    std::cout << "tau_before_update: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "YPAPY: " << ai.YPAPY << std::endl;
    std::cout << "Trace: " << ai.Trace << std::endl;
    std::cout << "score (YPAPY - Trace): " << score << std::endl;
    std::cout << "AI: " << AI << std::endl;
    std::cout << "Dtau (score/AI): " << Dtau << std::endl;
    std::cout << "step: " << step << std::endl;
    std::cout << "tau1_new: " << tau1_new << std::endl;

    arma::fvec alpha = coef.alpha;
    std::cout << "alpha: [";
    for (size_t i = 0; i < alpha.n_elem; ++i) {
        std::cout << alpha[i];
        if (i < alpha.n_elem - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // NOTE: Do NOT recompute eta here! The correct eta was already set in the inner IRLS loop
    // at line 493 (eta = coef.eta + offset). Line below was a bug that overwrote the correct
    // mixed-model eta with simple X*alpha, causing divergence from R at iteration 3.
    // eta = (p > 0) ? (X * alpha + offset) : offset;  // BUG - REMOVED

    tau_prev = tau;  tau[1] = static_cast<float>(tau1_new);
    // Use R-style convergence: max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol
    double rc_tau   = rel_change_R_style(tau, tau_prev, tol_coef);
    // R only checks tau for convergence in binary case (alpha already converged in Get_Coef)
    // double rc_alpha = (p > 0) ? rel_change_R_style(alpha, alpha_prev, tol_coef) : 0.0;

    std::cout << "tau_after_update: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "rc_tau: " << rc_tau << " (converge if < " << tol_coef << ")" << std::endl;
    std::cout << "==========================================\n" << std::endl;

    // ============ DETAILED OUTPUT FOR EVERY ITERATION ============
    std::cout << "\n===== ITERATION " << it << " DETAILED VALUES =====" << std::endl;
    std::cout << "tau_before: [" << tau_prev[0] << ", " << tau_prev[1] << "]" << std::endl;
    std::cout << "tau_after: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "YPAPY: " << ai.YPAPY << std::endl;
    std::cout << "Trace: " << ai.Trace << std::endl;
    std::cout << "Score: " << score << std::endl;
    std::cout << "AI: " << AI << std::endl;
    std::cout << "Dtau: " << Dtau << std::endl;
    std::cout << "alpha: [" << coef.alpha[0] << ", " << coef.alpha[1] << ", " << coef.alpha[2] << "]" << std::endl;
    std::cout << "First 5 eta: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;
    std::cout << "First 5 W: " << W[0] << " " << W[1] << " " << W[2] << " " << W[3] << " " << W[4] << std::endl;
    std::cout << "First 5 Y: " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << std::endl;
    std::cout << "First 5 Sigma_iY: " << coef.Sigma_iY[0] << " " << coef.Sigma_iY[1] << " " << coef.Sigma_iY[2] << " " << coef.Sigma_iY[3] << " " << coef.Sigma_iY[4] << std::endl;
    std::cout << "First 5 PY: " << ai.PY[0] << " " << ai.PY[1] << " " << ai.PY[2] << " " << ai.PY[3] << " " << ai.PY[4] << std::endl;
    std::cout << "First 5 APY: " << ai.APY[0] << " " << ai.APY[1] << " " << ai.APY[2] << " " << ai.APY[3] << " " << ai.APY[4] << std::endl;
    std::cout << "Norms: |eta|=" << arma::norm(eta) << " |W|=" << arma::norm(W) << " |Y|=" << arma::norm(Y) << std::endl;
    std::cout << "       |Sigma_iY|=" << arma::norm(coef.Sigma_iY) << " |PY|=" << arma::norm(ai.PY) << " |APY|=" << arma::norm(ai.APY) << std::endl;
    std::cout << "============================================\n" << std::endl;

    // Save iteration summary to CSV
    {
      std::ofstream f(CP_DIR + "/CPP_iteration_summary.csv",
                      it == 0 ? std::ios::trunc : std::ios::app);
      if (it == 0) {
        f << "iteration,tau_before,tau_after,Dtau,YPAPY,Trace,Score,AI,alpha0,alpha1,alpha2,eta_norm,W_norm,Y_norm,Sigma_iY_norm,PY_norm,APY_norm\n";
      }
      f << std::setprecision(10);
      f << it << "," << tau_prev[1] << "," << tau[1] << "," << Dtau << ","
        << ai.YPAPY << "," << ai.Trace << "," << score << "," << AI << ","
        << coef.alpha[0] << "," << coef.alpha[1] << "," << coef.alpha[2] << ","
        << arma::norm(eta) << "," << arma::norm(W) << "," << arma::norm(Y) << ","
        << arma::norm(coef.Sigma_iY) << "," << arma::norm(ai.PY) << "," << arma::norm(ai.APY) << "\n";
    }
    // ============ END DETAILED OUTPUT FOR EVERY ITERATION ============

    // ============ EARLY EXIT FOR DEBUGGING: Stop after 4th real tau update ============
    if (STOP_AFTER_FIRST_TAU_UPDATE && std::fabs(Dtau) > 1e-10) {
      real_tau_updates++;
      if (real_tau_updates >= 4) {
        std::cout << "\n" << std::string(70, '*') << std::endl;
        std::cout << "*** STOPPING AFTER 4TH REAL TAU UPDATE (Dtau=" << Dtau << ") ***" << std::endl;
        std::cout << "*** Set STOP_AFTER_FIRST_TAU_UPDATE=false to continue ***" << std::endl;
        std::cout << std::string(70, '*') << std::endl;

        // Save detailed checkpoint for comparison
        std::ofstream f(CP_DIR + "/CPP_FIRST_TAU_UPDATE.txt");
        f << std::setprecision(12);
        f << "iteration: " << it << "\n";
        f << "tau_before: " << tau_prev[0] << " " << tau_prev[1] << "\n";
        f << "tau_after: " << tau[0] << " " << tau[1] << "\n";
        f << "YPAPY: " << ai.YPAPY << "\n";
        f << "Trace: " << ai.Trace << "\n";
        f << "Score: " << score << "\n";
        f << "AI: " << AI << "\n";
        f << "Dtau: " << Dtau << "\n";
        f << "alpha: " << coef.alpha[0] << " " << coef.alpha[1] << " " << coef.alpha[2] << "\n";
        f << "eta_norm: " << arma::norm(eta) << "\n";
        f << "W_norm: " << arma::norm(W) << "\n";
        f << "Y_norm: " << arma::norm(Y) << "\n";
        f << "PY_norm: " << arma::norm(ai.PY) << "\n";
        f << "APY_norm: " << arma::norm(ai.APY) << "\n";
        f.close();

        save_vec_csv(eta, CP_DIR + "/CPP_FIRST_TAU_eta.csv");
        save_vec_csv(W, CP_DIR + "/CPP_FIRST_TAU_W.csv");
        save_vec_csv(Y, CP_DIR + "/CPP_FIRST_TAU_Y.csv");
        save_vec_csv(ai.PY, CP_DIR + "/CPP_FIRST_TAU_PY.csv");
        save_vec_csv(ai.APY, CP_DIR + "/CPP_FIRST_TAU_APY.csv");
        save_vec_csv(coef.Sigma_iY, CP_DIR + "/CPP_FIRST_TAU_Sigma_iY.csv");

        // Return early with current results
        arma::fvec mu_final = 1.0f / (1.0f + arma::exp(-eta));
        auto sn = saige::build_score_null_binary(X, y, mu_final);
        FitNullResult out;
        out.alpha = std::vector<double>(coef.alpha.begin(), coef.alpha.end());
        out.theta = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
        out.offset = offset_in;
        stash_score_null_into(out, sn, n, p);
        return out;
      }
    }
    // ============ END EARLY EXIT ============

    // ============ CHECKPOINT 5: After iteration 1 (it == 0) ============
    if (it == 0) {
      std::cout << "\n=== C++ CHECKPOINT 5: After iteration 1 ===" << std::endl;
      std::cout << "tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
      std::cout << "Score: " << score << std::endl;
      std::cout << "Dtau: " << Dtau << std::endl;

      // Save to files
      std::ofstream f_cp5(CP_DIR + "/CPP_CP5_iteration1.csv");
      f_cp5 << "tau0," << tau[0] << "\n";
      f_cp5 << "tau1," << tau[1] << "\n";
      f_cp5 << "Score," << score << "\n";
      f_cp5 << "Dtau," << Dtau << "\n";
      f_cp5 << "YPAPY," << ai.YPAPY << "\n";
      f_cp5 << "Trace," << ai.Trace << "\n";
      f_cp5 << "AI," << AI << "\n";
      f_cp5.close();
      std::cout << "Saved: " << CP_DIR << "/CPP_CP5_iteration1.csv" << std::endl;
    }

    // R SAIGE: first iteration (conservative update) does NOT check convergence
    // Only check convergence from iteration 1 onwards (after standard AI-REML)
    if (it > 0 && rc_tau < tol_coef) {
      // finalize + stash score-null
      arma::fvec mu_final = 1.0f / (1.0f + arma::exp(-eta));
      auto sn = saige::build_score_null_binary(X, y, mu_final);

      FitNullResult out;
      out.alpha  = std::vector<double>(alpha.begin(), alpha.end());
      out.theta  = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
      out.offset = offset_in;

      stash_score_null_into(out, sn, n, p);
      export_score_null_json(paths, out);
      std::cout << "=== CONVERGED at iteration " << it << " ===" << std::endl;

      // ============ CHECKPOINT 6: Final results (converged) ============
      std::cout << "\n=== C++ CHECKPOINT 6: Final results (CONVERGED) ===" << std::endl;
      std::cout << "Final tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
      std::cout << "Final alpha: [";
      for (size_t i = 0; i < alpha.n_elem; ++i) {
        std::cout << alpha[i];
        if (i < alpha.n_elem - 1) std::cout << ", ";
      }
      std::cout << "]" << std::endl;
      std::cout << "Converged: true" << std::endl;
      std::cout << "Iterations: " << (it + 1) << std::endl;

      // Save to files
      std::ofstream f_cp6(CP_DIR + "/CPP_CP6_final.csv");
      f_cp6 << "tau0," << tau[0] << "\n";
      f_cp6 << "tau1," << tau[1] << "\n";
      f_cp6 << "converged,1\n";
      f_cp6 << "iterations," << (it + 1) << "\n";
      for (size_t i = 0; i < alpha.n_elem; ++i) {
        f_cp6 << "alpha" << i << "," << alpha[i] << "\n";
      }
      f_cp6.close();
      std::cout << "Saved: " << CP_DIR << "/CPP_CP6_final.csv" << std::endl;

      return out;
    }

    alpha_prev = alpha;
  }

  // fallthrough
  irls_binary_build(eta, y, offset, mu, mu_eta, W, Y);
  coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);

  arma::fvec mu_final = 1.0f / (1.0f + arma::exp(-eta));
  auto sn = saige::build_score_null_binary(X, y, mu_final);

  FitNullResult out;
  out.alpha  = std::vector<double>(coef.alpha.begin(), coef.alpha.end());
  out.theta  = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
  out.offset = offset_in;

  stash_score_null_into(out, sn, n, p);
  // export_score_null_json(paths, out);

  // ============ CHECKPOINT 6: Final results (max iterations reached) ============
  std::cout << "\n=== C++ CHECKPOINT 6: Final results (MAX ITERATIONS) ===" << std::endl;
  std::cout << "Final tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
  std::cout << "Final alpha: [";
  for (size_t i = 0; i < coef.alpha.n_elem; ++i) {
    std::cout << coef.alpha[i];
    if (i < coef.alpha.n_elem - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "Converged: false" << std::endl;
  std::cout << "Iterations: " << maxiter << std::endl;

  // Save to files
  std::ofstream f_cp6(CP_DIR + "/CPP_CP6_final.csv");
  f_cp6 << "tau0," << tau[0] << "\n";
  f_cp6 << "tau1," << tau[1] << "\n";
  f_cp6 << "converged,0\n";
  f_cp6 << "iterations," << maxiter << "\n";
  for (size_t i = 0; i < coef.alpha.n_elem; ++i) {
    f_cp6 << "alpha" << i << "," << coef.alpha[i] << "\n";
  }
  f_cp6.close();
  std::cout << "Saved: " << CP_DIR << "/CPP_CP6_final.csv" << std::endl;

  return out;
}

// ---------- Quantitative solver (2×2 AI-REML on [tau0, tau1]) ----------
// Ported from binary_glmm_solver with adaptations:
// - tau starts at [1, 0] (R's default: residual=1, genetic=0)
// - Both tau[0] and tau[1] are free (binary fixes tau[0]=1)
// - Uses irls_gaussian_build (W=1, identity link)
// - 2x2 AI matrix, 2 scores, 2 traces
// - Conservative first iteration, standard AI-REML after
// - Step halving when any tau < 0

FitNullResult quant_glmm_solver(const Paths& paths,
                                const FitNullConfig& cfg,
                                const Design& d,
                                const std::vector<double>& offset_in,
                                const std::vector<double>& beta_init_in)
{
  const int n = d.n, p = d.p;

  arma::fmat X = (p > 0) ? map_X_row_major_to_fmat(d) : arma::fmat(n, 0);
  arma::fvec y = map_y(d);
  arma::fvec offset = map_offset(d, offset_in);
  arma::fvec beta_init = map_beta_init(d, beta_init_in);
  arma::fvec eta = (p > 0) ? (X * beta_init + offset) : offset;

  const int   maxiter    = std::max(5, cfg.maxiter);
  const float tol_coef   = static_cast<float>(std::max(1e-6, cfg.tol));
  const int   maxiterPCG = cfg.maxiterPCG > 0 ? cfg.maxiterPCG : 500;
  const float tolPCG     = cfg.tolPCG > 0.0 ? static_cast<float>(cfg.tolPCG) : 1e-5f;
  const int   nrun       = cfg.nrun > 0 ? cfg.nrun : 30;
  const float trace_cut  = cfg.traceCVcutoff > 0.0 ? static_cast<float>(cfg.traceCVcutoff) : 0.1f;

  // R quantitative initializes tau = [1, 0] (residual=1, genetic=0)
  arma::fvec tau(2); tau[0] = 1.0f; tau[1] = 0.0f;

  arma::fvec mu, mu_eta, W, Y;
  arma::fvec alpha_prev(p, arma::fill::zeros);
  arma::fvec tau_prev = tau;

  // Track alpha across outer iterations (matching R's Get_Coef alpha0)
  arma::fvec alpha_outer_prev = arma::conv_to<arma::fvec>::from(beta_init);
  CoefficientsOut coef;

  // ===== Debug: Initial values =====
  std::cout << "\n===== C++ Quantitative GLMM Solver Initial Values =====" << std::endl;
  std::cout << "tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
  std::cout << "beta_init: [";
  for (size_t i = 0; i < beta_init.n_elem; ++i) {
    std::cout << beta_init[i];
    if (i < beta_init.n_elem - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "eta[0:5]: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;
  std::cout << "y[0:5]: " << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << std::endl;
  std::cout << "n=" << n << ", p=" << p << std::endl;
  std::cout << "==========================================\n" << std::endl;

  // Checkpoint 2
  save_vec_csv(eta, CP_DIR + "/CPP_CP2_eta.csv");
  save_vec_csv(y, CP_DIR + "/CPP_CP2_y.csv");
  {
    std::ofstream f(CP_DIR + "/CPP_CP2_initial.txt");
    f << "tau: " << tau[0] << " " << tau[1] << "\n";
    f << "n: " << n << "\np: " << p << "\n";
    f << "beta_init: ";
    for (size_t i = 0; i < beta_init.n_elem; ++i) f << beta_init[i] << " ";
    f << "\n";
  }

  for (int it = 0; it < maxiter; ++it) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "=== C++ QUANT OUTER ITERATION " << it << " START ===" << std::endl;
    std::cout << "tau at start: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    // ======================================================================
    // Inner IRLS loop (for gaussian W=1, converges in 1 step)
    // ======================================================================
    arma::fvec alpha0_outer = alpha_outer_prev;

    for (int inner_it = 0; inner_it < maxiter; ++inner_it) {
      std::cout << "\n--- Inner IRLS iteration " << inner_it << " ---" << std::endl;

      irls_gaussian_build(eta, y, offset, mu, mu_eta, W, Y);

      std::cout << "After irls_gaussian_build:" << std::endl;
      std::cout << "  W[0:5]: " << W[0] << " " << W[1] << " " << W[2] << " " << W[3] << " " << W[4] << std::endl;
      std::cout << "  Y[0:5]: " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << std::endl;
      std::cout << "  |Y|: " << arma::norm(Y) << ", |W|: " << arma::norm(W) << std::endl;

      coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);

      std::cout << "After getCoefficients_cpp:" << std::endl;
      std::cout << "  alpha: [";
      for (size_t i = 0; i < coef.alpha.n_elem; ++i) {
        std::cout << coef.alpha[i];
        if (i < coef.alpha.n_elem - 1) std::cout << ", ";
      }
      std::cout << "]" << std::endl;
      std::cout << "  coef.eta[0:5]: " << coef.eta[0] << " " << coef.eta[1] << " " << coef.eta[2] << " " << coef.eta[3] << " " << coef.eta[4] << std::endl;

      // Update eta from coefficient solve (NOT X*alpha -- same fix as binary)
      eta = coef.eta + offset;
      std::cout << "  eta after update [0:5]: " << eta[0] << " " << eta[1] << " " << eta[2] << " " << eta[3] << " " << eta[4] << std::endl;

      // Check convergence of alpha
      double rc_alpha_inner = (p > 0) ? rel_change_R_style(coef.alpha, alpha0_outer, tol_coef) : 0.0;
      std::cout << "  rc_alpha (vs alpha0_outer): " << rc_alpha_inner << " (tol=" << tol_coef << ")" << std::endl;

      if (rc_alpha_inner < tol_coef) {
        std::cout << "  [Inner IRLS] CONVERGED after " << (inner_it + 1) << " iterations" << std::endl;
        break;
      }
      alpha0_outer = coef.alpha;
    }

    alpha_outer_prev = coef.alpha;

    // Final IRLS build with converged eta
    std::cout << "\n--- Final IRLS build after inner loop ---" << std::endl;
    irls_gaussian_build(eta, y, offset, mu, mu_eta, W, Y);
    std::cout << "Final W[0:5]: " << W[0] << " " << W[1] << " " << W[2] << " " << W[3] << " " << W[4] << std::endl;
    std::cout << "Final Y[0:5]: " << Y[0] << " " << Y[1] << " " << Y[2] << " " << Y[3] << " " << Y[4] << std::endl;
    std::cout << "Final |Y|: " << arma::norm(Y) << ", |W|: " << arma::norm(W) << std::endl;

    // ===== AI-REML step: 2D score and AI =====
    auto aiq = getAIScore_q_cpp(Y, X, W, tau, coef.Sigma_iY, coef.Sigma_iX,
                                coef.cov, nrun, maxiterPCG, tolPCG, trace_cut);

    // 2D scores
    double score0 = static_cast<double>(aiq.YPA0PY - aiq.Trace[0]);
    double score1 = static_cast<double>(aiq.YPAPY  - aiq.Trace[1]);

    // Debug output
    std::cout << "\n========== C++ QUANT ITERATION " << it << " ==========" << std::endl;
    std::cout << "tau_before: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "YPA0PY: " << aiq.YPA0PY << std::endl;
    std::cout << "YPAPY: " << aiq.YPAPY << std::endl;
    std::cout << "Trace: [" << aiq.Trace[0] << ", " << aiq.Trace[1] << "]" << std::endl;
    std::cout << "Score: [" << score0 << ", " << score1 << "]" << std::endl;
    std::cout << "AI:" << std::endl;
    std::cout << "  [" << aiq.AI(0,0) << ", " << aiq.AI(0,1) << "]" << std::endl;
    std::cout << "  [" << aiq.AI(1,0) << ", " << aiq.AI(1,1) << "]" << std::endl;

    arma::fvec tau_new(2);

    if (it == 0) {
      // Conservative update (R SAIGE lines 891-892):
      // tau[i] = max(0, tau[i] + tau[i]^2 * score[i] / n)
      double tau0_d = static_cast<double>(tau[0]);
      double tau1_d = static_cast<double>(tau[1]);
      double Dtau0 = tau0_d * tau0_d * score0 / static_cast<double>(n);
      double Dtau1 = tau1_d * tau1_d * score1 / static_cast<double>(n);
      tau_new[0] = static_cast<float>(std::max(0.0, tau0_d + Dtau0));
      tau_new[1] = static_cast<float>(std::max(0.0, tau1_d + Dtau1));
      std::cout << "[CONSERVATIVE] Dtau: [" << Dtau0 << ", " << Dtau1 << "]" << std::endl;
    } else {
      // Standard AI-REML: delta = solve(AI_2x2, score_2)
      arma::fvec s(2);
      s[0] = static_cast<float>(score0);
      s[1] = static_cast<float>(score1);

      arma::fmat AI = aiq.AI;
      if (!AI.is_sympd()) { AI = 0.5f * (AI + AI.t()); }
      arma::fvec delta = arma::solve(AI, s, arma::solve_opts::likely_sympd + arma::solve_opts::fast);

      tau_new = tau + delta;

      // Step halving if any component negative (R SAIGE behavior)
      double step = 1.0;
      while ((tau_new[0] < 0.0f || tau_new[1] < 0.0f) && step > 1e-10) {
        step *= 0.5;
        tau_new = tau + static_cast<float>(step) * delta;
      }
      // Final clamp
      tau_new[0] = static_cast<float>(std::max(0.0, static_cast<double>(tau_new[0])));
      tau_new[1] = static_cast<float>(std::max(0.0, static_cast<double>(tau_new[1])));

      std::cout << "[STANDARD AI-REML] delta: [" << delta[0] << ", " << delta[1] << "]" << std::endl;
      std::cout << "[STANDARD AI-REML] step: " << step << std::endl;
    }

    tau_prev = tau;
    tau = tau_new;

    // R-style convergence: max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol
    double rc_tau = rel_change_R_style(tau, tau_prev, tol_coef);

    std::cout << "tau_after: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
    std::cout << "rc_tau: " << rc_tau << " (converge if < " << tol_coef << ")" << std::endl;
    std::cout << "alpha: [";
    for (size_t i = 0; i < coef.alpha.n_elem; ++i) {
      std::cout << coef.alpha[i];
      if (i < coef.alpha.n_elem - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "==========================================\n" << std::endl;

    // Save iteration to CSV
    {
      std::ofstream f(CP_DIR + "/CPP_iteration_summary.csv",
                      it == 0 ? std::ios::trunc : std::ios::app);
      if (it == 0) {
        f << "iteration,tau0_before,tau1_before,tau0_after,tau1_after,score0,score1,YPA0PY,YPAPY,Trace0,Trace1\n";
      }
      f << std::setprecision(10);
      f << it << "," << tau_prev[0] << "," << tau_prev[1] << "," << tau[0] << "," << tau[1]
        << "," << score0 << "," << score1 << "," << aiq.YPA0PY << "," << aiq.YPAPY
        << "," << aiq.Trace[0] << "," << aiq.Trace[1] << "\n";
    }

    // Convergence check (only after iteration 1, matching R)
    if (it > 0 && rc_tau < tol_coef) {
      arma::fvec mu_final = eta;  // identity link
      float tau0_inv = (tau[0] > 0.0f) ? 1.0f / tau[0] : 0.0f;
      auto sn = saige::build_score_null_quant(X, y, mu_final, tau0_inv);

      FitNullResult out;
      out.alpha = std::vector<double>(coef.alpha.begin(), coef.alpha.end());
      out.theta = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
      out.offset = offset_in;

      stash_score_null_into(out, sn, n, p);
      export_score_null_json(paths, out);

      std::cout << "=== CONVERGED at iteration " << it << " ===" << std::endl;
      std::cout << "Final tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
      std::cout << "Iterations: " << (it + 1) << std::endl;

      // Save final checkpoint
      {
        std::ofstream f(CP_DIR + "/CPP_CP6_final.csv");
        f << "tau0," << tau[0] << "\ntau1," << tau[1] << "\n";
        f << "converged,1\niterations," << (it + 1) << "\n";
        for (size_t i = 0; i < coef.alpha.n_elem; ++i)
          f << "alpha" << i << "," << coef.alpha[i] << "\n";
      }

      return out;
    }

    alpha_prev = coef.alpha;
  }

  // Fallthrough: max iterations reached
  irls_gaussian_build(eta, y, offset, mu, mu_eta, W, Y);
  coef = getCoefficients_cpp(Y, X, W, tau, maxiterPCG, tolPCG);

  arma::fvec mu_final = eta;
  float tau0_inv = (tau[0] > 0.0f) ? 1.0f / tau[0] : 0.0f;
  auto sn = saige::build_score_null_quant(X, y, mu_final, tau0_inv);

  FitNullResult out;
  out.alpha = std::vector<double>(coef.alpha.begin(), coef.alpha.end());
  out.theta = {static_cast<double>(tau[0]), static_cast<double>(tau[1])};
  out.offset = offset_in;

  stash_score_null_into(out, sn, n, p);

  std::cout << "\n=== MAX ITERATIONS REACHED ===" << std::endl;
  std::cout << "Final tau: [" << tau[0] << ", " << tau[1] << "]" << std::endl;
  std::cout << "Iterations: " << maxiter << std::endl;

  {
    std::ofstream f(CP_DIR + "/CPP_CP6_final.csv");
    f << "tau0," << tau[0] << "\ntau1," << tau[1] << "\n";
    f << "converged,0\niterations," << maxiter << "\n";
    for (size_t i = 0; i < coef.alpha.n_elem; ++i)
      f << "alpha" << i << "," << coef.alpha[i] << "\n";
  }

  return out;
}

// ---------- registration ----------

void register_default_solvers() {
  register_binary_solver(&binary_glmm_solver);
  register_quant_solver (&quant_glmm_solver);
}

} // namespace saige
