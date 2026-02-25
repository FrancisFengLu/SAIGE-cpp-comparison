// saige_null.cpp
// ------------------------------------------------------------------
// High-level orchestration for fitting the null model and (optionally)
// estimating variance ratios. This implements saige::fit_null().
// ------------------------------------------------------------------

#include "saige_null.hpp"
#include "preprocess_engine.hpp"
#include "null_model_engine.hpp"
#include "variance_ratio_engine.hpp"
#include "variance_ratio_compute.hpp"
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
  FitNullResult fit = nme.run(prep.design);

  // debug
  std::cout << "NullModelEngine Called" << std::endl; 
  //

  // --- Variance ratio (optional) ---
  // Heuristic: run VR when a positive number of markers is requested AND
  // PLINK files are present; your VR hook can relax/override this if desired.
  const bool have_plink =
      !paths.bed.empty() && !paths.bim.empty() && !paths.fam.empty();
  const bool do_vr =
      (prep.cfg.num_markers_for_vr > 0) && have_plink;

  if (do_vr) {
    VarianceRatioEngine vre(paths, prep.cfg, prep.chr);
      std::cout << "VarianceRatioEngine Called" << std::endl; 
    fit = vre.run(fit, prep.design);
  }

  return fit;
}

} // namespace saige
