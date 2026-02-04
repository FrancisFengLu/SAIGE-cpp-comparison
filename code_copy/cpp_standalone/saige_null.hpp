#pragma once
#include <string>
#include <vector>
#include <optional>

namespace saige {

// -------- Core data structures --------

struct LocoRanges {
  // 0-based inclusive marker ranges per chromosome (chr 1..22 order).
  // If a chromosome is absent, store -1 for both start/end at that slot.
  std::vector<int> start;
  std::vector<int> end;
  bool enabled{false};
};

struct FitNullConfig {
  // Trait
  std::string trait{"binary"};            // "binary" | "quantitative" | "survival"

  // Feature flags
  bool loco{true};
  bool lowmem_loco{false};
  bool use_sparse_grm_to_fit{false};
  bool use_sparse_grm_for_vr{false};
  bool covariate_qr{true};
  bool covariate_offset{false};
  bool inv_normalize{false};
  bool include_nonauto_for_vr{false};
  bool isDiagofKinSetAsOne{true};
  bool make_sparse_grm_only{false};  // NEW

  // Convergence / runtime
  double tol{0.02};
  int    maxiter{20};
  double tolPCG{1e-5};
  int    maxiterPCG{500};
  int    nrun{30};
  int    nthreads{1};

  // Diagnostics / CV thresholds
  double traceCVcutoff{0.0025};
  double ratio_cv_cutoff{0.001};

  // GRM marker QC
  double min_maf_grm{0.01};
  double max_miss_grm{0.15};
  double relatedness_cutoff{0.05};

  // VR controls
  int num_markers_for_vr{30};

  // VR min MAC
  double memory_chunk_gb;
  int vr_min_mac;
  int vr_max_mac;

  // Survival
  std::optional<int> event_time_bin_size; // if set, bin event times by this size
  bool pcg_for_uhat_surv = true;

  bool        female_only{false};
  bool        male_only{false};
  std::string sex_col;              // e.g. "Sex"
  std::string female_code{"1"};     // matches your R defaults
  std::string male_code{"0"};

  bool        overwrite_vr{false}; 

};

struct Paths {
  std::string bed;
  std::string bim;
  std::string fam;

  std::string sparse_grm;
  std::string sparse_grm_ids;

  std::string out_prefix;
  std::string out_prefix_vr; // <prefix>.varianceRatio.txt will use this prefix
};

// Numeric design matrix + vectors, row-aligned across all fields.
struct Design {
  // X stored row-major (n*p entries). Access with Eigen::Map for compute.
  std::vector<double> X;
  int n{0};
  int p{0};

  std::vector<double> y;          // length n
  std::vector<double> offset;     // optional, length n or empty
  std::vector<double> event_time; // optional, length n or empty

  std::vector<std::string> iid;   // length n, order must match rows
};

// Minimal score-null payload for SPA / downstream score tests
struct ScoreNullPack {
  int n{0}, p{0};
  std::vector<double> V;        // length n  (mu*(1-mu) for binary; 1/tau0 for quant)
  std::vector<double> S_a;      // length p  (colSums(X ⊙ residuals))
  std::vector<double> XVX;      // p x p row-major
  std::vector<double> XVX_inv;  // p x p row-major
};

struct FitNullResult {
  // Model parameters / summaries
  std::vector<double> theta;   // variance components
  std::vector<double> alpha;   // fixed effects (original design scale)
  std::vector<double> offset;  // final GLMM offset (length n)

  // Flags
  bool loco{false};
  bool lowmem_loco{false};

  // Artifacts (paths); may be empty if not written
  std::string model_rda_path;     // main model artifact path (json/rds/etc.)
  std::string vr_path;            // <prefix>.varianceRatio.txt
  std::string markers_out_path;   // <prefix>_<N>markers.SAIGE.results.txt

  // ... your existing fields (alpha, theta, offset, loco flags, paths, etc.)

  // Baseline score-null (R's obj.noK) in pure C++
  ScoreNullPack obj_noK;

  // Optional: per-chromosome LOCO score-nulls (only if you populate them)
  // Convention: index 0..21 => chr 1..22; entries may be empty (n=0) when not computed.
  std::vector<ScoreNullPack> loco_obj_noK;
};



// Optional library-level orchestrator (implemented in saige_null.cpp, if you use it)
FitNullResult fit_null(const FitNullConfig& cfg,
                       const Paths& paths,
                       const Design& design);

} // namespace saige
