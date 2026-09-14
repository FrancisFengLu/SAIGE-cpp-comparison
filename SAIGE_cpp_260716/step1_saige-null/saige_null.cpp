// saige_null.cpp
// ------------------------------------------------------------------
// High-level orchestration for fitting the null model and (optionally)
// estimating variance ratios. This implements saige::fit_null().
// ------------------------------------------------------------------

#include "saige_null.hpp"
#include "preprocess_engine.hpp"
#include <chrono>
#include "null_model_engine.hpp"
#include "variance_ratio_engine.hpp"
#include "variance_ratio_compute.hpp"
#include "glmm.hpp"
#include "SAIGE_step1_fast.hpp"
#include <filesystem>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace saige {

static inline void ensure_parent_dir_(const std::string& path) {
  if (path.empty()) return;
  fs::path p(path);
  auto dir = p.parent_path();
  if (!dir.empty()) fs::create_directories(dir);
}

static inline Paths sanitize_paths_(Paths paths) {
  // If out_prefix_vr is empty, default to out_prefix
  if (paths.out_prefix_vr.empty()) {
    paths.out_prefix_vr = paths.out_prefix;
  }
  return paths;
}

static inline void register_real_vr() {
  saige::register_vr_runner(&saige::compute_variance_ratio);
}

FitNullResult fit_null(const FitNullConfig& cfg_in,
                       const Paths& paths_in,
                       const Design& design_in)
{
  // --- Sanitize/normalize inputs (paths may be adjusted) ---
  Paths paths = sanitize_paths_(paths_in);
  FitNullConfig cfg = cfg_in;

  // Ensure parent dirs exist for outputs we know about
  ensure_parent_dir_(paths.out_prefix + ".dummy");      // we don't know exact extension here
  ensure_parent_dir_(paths.out_prefix_vr + ".dummy");   // same for VR

  // --- Preprocess: align samples, inv-norm / survival binning, LOCO ranges ---
  PreprocessEngine pre(paths, cfg);
  PreOut prep = pre.run(design_in);

  // --- Fit the null model (and LOCO offsets if hook is registered) ---
  NullModelEngine nme(paths, prep.cfg, prep.chr);

  register_real_vr();
  auto T_glmm = std::chrono::steady_clock::now();
  FitNullResult fit = nme.run(prep.design);
  {
    auto t = std::chrono::steady_clock::now();
    printf("[TIMER-MAIN] %-40s %8.2fs\n", "GLMM null model fitting",
           std::chrono::duration<double>(t - T_glmm).count());
  }
  std::cout << "NullModelEngine Called" << std::endl;

  // --- Variance ratio (optional) ---
  const bool have_plink =
      !paths.bed.empty() && !paths.bim.empty() && !paths.fam.empty();
  const bool do_vr =
      (prep.cfg.num_markers_for_vr > 0) && have_plink;

  if (do_vr) {
    auto T_vr = std::chrono::steady_clock::now();
    VarianceRatioEngine vre(paths, prep.cfg, prep.chr);
    std::cout << "VarianceRatioEngine Called" << std::endl;
    fit = vre.run(fit, prep.design);
    {
      auto t = std::chrono::steady_clock::now();
      printf("[TIMER-MAIN] %-40s %8.2fs\n", "Variance ratio estimation",
             std::chrono::duration<double>(t - T_vr).count());
    }
  }

  return fit;
}

// ---------------------------------------------------------------------------
// Tier-2 lockstep multi-phenotype fit.
//
// Same three phases as fit_null(), but the middle one is shared:
//   prep   : per trait (preprocess + NullModelEngine::prep). No GRM touched.
//   solve  : ONE call to multi_glmm_solver for all P traits.
//   export : per trait (NullModelEngine::export_result + variance ratio).
// ---------------------------------------------------------------------------
std::vector<FitNullResult> fit_null_multi(const FitNullConfig& cfg_in,
                                          const std::vector<Paths>& paths_in,
                                          const std::vector<Design>& designs_in)
{
  const size_t P = designs_in.size();
  if (P == 0) throw std::runtime_error("fit_null_multi: no phenotypes");
  if (paths_in.size() != P)
    throw std::runtime_error("fit_null_multi: paths/designs length mismatch");

  std::vector<Paths>         paths(P);
  std::vector<PreOut>        prep(P);
  std::vector<NullModelEngine> engines;
  std::vector<std::shared_ptr<NullPrep>> nprep(P);
  engines.reserve(P);

  register_real_vr();

  auto T_prep = std::chrono::steady_clock::now();
  for (size_t i = 0; i < P; ++i) {
    paths[i] = sanitize_paths_(paths_in[i]);
    ensure_parent_dir_(paths[i].out_prefix + ".dummy");
    ensure_parent_dir_(paths[i].out_prefix_vr + ".dummy");

    PreprocessEngine pre(paths[i], cfg_in);
    prep[i] = pre.run(designs_in[i]);
    engines.emplace_back(paths[i], prep[i].cfg, prep[i].chr);
  }
  for (size_t i = 0; i < P; ++i)
    nprep[i] = engines[i].prep(prep[i].design);
  {
    auto t = std::chrono::steady_clock::now();
    printf("[TIMER-MAIN] %-40s %8.2fs\n", "lockstep prep (all traits)",
           std::chrono::duration<double>(t - T_prep).count());
  }

  // One driver fits one family, and every trait shares the solver-tuning
  // constants (they come from one YAML). Refuse rather than silently fit some
  // trait with somebody else's tolerance.
  const bool is_binary = nullprep_is_binary(*nprep[0]);
  for (size_t i = 1; i < P; ++i) {
    if (nullprep_is_binary(*nprep[i]) != is_binary)
      throw std::runtime_error("fit_null_multi: cannot lockstep binary and "
                               "quantitative traits in one batch");
    const FitNullConfig& a = prep[0].cfg;
    const FitNullConfig& b = prep[i].cfg;
    if (a.maxiter != b.maxiter || a.tol != b.tol || a.tolPCG != b.tolPCG ||
        a.maxiterPCG != b.maxiterPCG || a.nrun != b.nrun ||
        a.traceCVcutoff != b.traceCVcutoff)
      throw std::runtime_error("fit_null_multi: traits disagree on the AI-REML "
                               "tuning constants");
  }

  std::vector<const Design*>               d_ptr(P);
  std::vector<const std::vector<double>*>  off_ptr(P);
  std::vector<const std::vector<double>*>  bi_ptr(P);
  for (size_t i = 0; i < P; ++i) {
    d_ptr[i]   = &nullprep_design(*nprep[i]);
    off_ptr[i] = &nullprep_offset(*nprep[i]);
    bi_ptr[i]  = &nullprep_beta_init(*nprep[i]);
  }

  auto T_glmm = std::chrono::steady_clock::now();
  std::vector<FitNullResult> fits =
      multi_glmm_solver(paths, prep[0].cfg, d_ptr, off_ptr, bi_ptr, is_binary);
  {
    auto t = std::chrono::steady_clock::now();
    const double s = std::chrono::duration<double>(t - T_glmm).count();
    printf("[TIMER-MAIN] %-40s %8.2fs (%.2fs/trait)\n",
           "GLMM null model fitting (lockstep)", s, s / (double)P);
  }

  const bool have_plink =
      !paths[0].bed.empty() && !paths[0].bim.empty() && !paths[0].fam.empty();

  for (size_t i = 0; i < P; ++i) {
    fits[i] = engines[i].export_result(*nprep[i], std::move(fits[i]));

    if (prep[i].cfg.num_markers_for_vr > 0 && have_plink) {
      auto T_vr = std::chrono::steady_clock::now();
      VarianceRatioEngine vre(paths[i], prep[i].cfg, prep[i].chr);
      fits[i] = vre.run(fits[i], prep[i].design);
      auto t = std::chrono::steady_clock::now();
      printf("[TIMER-MAIN] %-40s %8.2fs\n", "Variance ratio estimation",
             std::chrono::duration<double>(t - T_vr).count());
    }
  }

  return fits;
}

} // namespace saige
