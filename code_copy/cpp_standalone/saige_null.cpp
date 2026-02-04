// saige_null.cpp
// ------------------------------------------------------------------
// High-level orchestration for fitting the null model and (optionally)
// estimating variance ratios. This implements saige::fit_null().
// ------------------------------------------------------------------

#include "saige_null.hpp"
#include "preprocess_engine.hpp"
#include "null_model_engine.hpp"
#include "variance_ratio_engine.hpp"
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

static void neutral_vr_runner(const Paths& paths,
                              const FitNullConfig& cfg,
                              const LocoRanges& chr,
                              const FitNullResult& fit,
                              const Design& design,
                              std::string& out_vr_path,
                              std::string& out_marker_results_path)
{
  (void)chr; (void)fit; (void)design;

  // Decide output path
  std::string vr_out = paths.out_prefix_vr.empty()
                         ? (paths.out_prefix + ".variance_ratio.tsv")
                         : (paths.out_prefix_vr + ".variance_ratio.tsv");

  std::ofstream ofs(vr_out);
  // Header expected by VR validation
  ofs << "MAC_Lower\tMAC_Upper\tVarianceRatio\n";

  // Make a couple of bins with VR=1.0
  int lo  = cfg.vr_min_mac > 0 ? cfg.vr_min_mac : 1;
  int mid = std::max(lo, lo * 4 - 1);
  int hi  = std::max(mid + 1, (cfg.vr_max_mac > 0 ? cfg.vr_max_mac : mid + 10));

  ofs << lo     << '\t' << mid     << '\t' << 1.0 << '\n';
  ofs << mid+1  << '\t' << hi      << '\t' << 1.0 << '\n';
  ofs.close();

  out_vr_path = vr_out;
  out_marker_results_path.clear();
}

static inline void register_default_vr() {
  saige::register_vr_runner(&neutral_vr_runner);
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

  register_default_vr();
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
