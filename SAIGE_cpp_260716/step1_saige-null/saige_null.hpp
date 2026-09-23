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
  // Matches R SAIGE's default (LOCO=TRUE in R/SAIGE_fitGLMM_fast.R:1199 and
  // R/SAIGE_Test_main.R). LOCO auto-disables, with a log message, when it
  // cannot be conducted: fewer than 2 autosomes in the BIM (R/Util.R:53-57),
  // no BIM at all, or a sparse GRM used to fit the null model (main.cpp).
  // nullmodel.json's "loco" field reports what ACTUALLY happened, never the
  // config flag.
  bool loco{true};
  bool lowmem_loco{false};
  bool use_sparse_grm_to_fit{false};
  bool use_sparse_grm_for_vr{false};
  bool covariate_qr{true};
  bool covariate_offset{true};   // Docker step1_fitNULLGLMM.R CLI default is TRUE (function default is FALSE; CLI overrides)
  bool inv_normalize{false};
  // Stage-0 profiler for the sparse direct solve (see spsolve_prof). Off by
  // default; when off nothing is timed and no arithmetic changes.
  bool profile_spsolve{false};
  bool include_nonauto_for_vr{false};
  bool isDiagofKinSetAsOne{false};  // R default is FALSE; TRUE forces GRM diagonal to 1.0
  bool make_sparse_grm_only{false};  // NEW
  bool use_pcg_with_sparse_grm{false};  // false = direct solve (R default), true = PCG

  // Selective genotype load for the sparse-GRM fitting path.
  //
  // With use_sparse_grm_to_fit the GLMM fit never touches the genotype matrix
  // (psi*u goes through the sparse Psi, Sigma^-1 through the sparse solve), so
  // the only numerical consumer left is the variance-ratio marker loop, which
  // reads ~1000 drawn markers out of every QC-passing marker in the BED. When
  // this is on, step 1 decodes ONLY those drawn markers.
  //
  // OFF by default because it changes one output: .grm_diag.txt needs every
  // marker and is therefore NOT written under it (step 2 never reads that file;
  // it is a diagnostic). Everything else — the variance ratio, the VR marker
  // list, tau, the model artifacts — is bit-for-bit what the full load gives,
  // because the VR draw happens before any BED byte is read and a marker's VR
  // eligibility depends on that marker alone. Ignored (with a printed reason)
  // for any configuration that does need the matrix; see main.cpp.
  bool selective_geno_load{false};

  // Tier-2 multi-phenotype: advance the P AI-REML loops in lockstep so the
  // fixed-effect PCG solves of all still-active traits go through one batched
  // multi-Sigma solve. OFF by default because it is not bit-identical to
  // fitting each trait alone: psi*B in the batched kernel reduces in a
  // different order than psi*b, so tau lands a few ULP away (the same class of
  // difference the trace and variance-ratio block paths already carry). Tier-1
  // (this flag false) keeps the property that a P>1 run reproduces each
  // trait's solo run byte for byte.
  bool multi_lockstep{false};

  // Scheme C (optimization/missing_mt/SCHEME_C_DESIGN.md): load the UNION of a
  // group's sample sets once and mask each phenotype's missing rows in the
  // multiplication, instead of one full genotype load per distinct sample set.
  // OFF by default for the same reason as lockstep: the union has extra zero
  // rows, so the fp32 block sums land a few ulp away from the solo run's (§4).
  // Only pays when the sample sets overlap heavily — mask_min_coverage is the
  // whole cost model (§6): a phenotype joins a mask group only while every
  // member still covers at least this fraction of the group's union.
  bool   mask_missing{false};
  double mask_min_coverage{0.8};
  // Acceptance C2 only (SCHEME_C_DESIGN.md §5): "" | "freq" | "qc" | "corr"
  // removes exactly one of §1's steps from the per-trait rebuild, so the gate
  // can be shown to catch each of them. Anything but "" prints a loud line and
  // produces deliberately wrong numbers.
  std::string scheme_c_break{};

  // Convergence / runtime
  double tol{0.02};
  int    maxiter{20};
  double tolPCG{1e-5};
  int    maxiterPCG{500};
  int    nrun{30};
  // Trace-estimator RNG seed for AI-REML's stochastic tr(P·K) (Hutchinson probes).
  // -1 => use the builtin per-trait defaults (10 for binary GetTrace, 200 for
  // quantitative GetTrace_q), matching R SAIGE. Set >=0 to override and sweep
  // seeds (the trace estimate is unbiased, so different seeds give different tau
  // draws — useful for characterizing MC noise in the cpp↔R p-value comparison).
  int    trace_seed{-1};
  int    nthreads{1};

  // GPU acceleration of the K·u kernel (SAIGE_step1_fast.cpp::parallelCrossProd).
  // Enabled via fit.use_gpu: true in YAML or --gpu on CLI. Silent CPU fallback
  // when SAIGE was built without CUDA or no device is available at runtime.
  bool   use_gpu{false};

  // CPU-only blocked-GEMV K·u path (parallelCrossProd_blocked).
  // When true (and the sparse-Sigma path is not active), the SNP-by-SNP
  // dot+axpy loop is replaced by per-block GEMV calls. Off by default so
  // existing benchmarks remain bit-for-bit reproducible during A/B testing.
  bool   use_blocked_gemv{false};
  int    gemv_block_size{128};      // 64/128/256 are the suggested sweep points
  bool   gemv_verify{false};         // when true, run both paths and print rel/abs error

  // Diagnostics / CV thresholds
  double traceCVcutoff{0.0025};
  double ratio_cv_cutoff{0.001};

  // GRM marker QC
  double min_maf_grm{0.01};
  double max_miss_grm{0.15};
  double relatedness_cutoff{0.05};

  // VR controls
  int num_markers_for_vr{30};

  // step-2 knobs that step 1 writes into nullmodel.json. They used to be
  // hardcoded there, and nullmodel.json is step 2's ONLY source for them
  // (null_model_loader.cpp reads them off that file; there is no YAML override),
  // so a C++ user had no way to ask for the exact test or a different
  // imputation. Defaults are the previously hardcoded values, so an existing
  // config produces a byte-identical nullmodel.json.
  //
  // NOTE these two do NOT match the R CLI defaults: R defaults
  // --is_fastTest=FALSE and --impute_method=best_guess. Comparing C++ against
  // R without aligning them compares two different algorithms -- on mid (~1%
  // missing) the mean/best_guess difference alone moved log10 p by up to 1.1.
  bool        fast_test{true};              // R CLI: FALSE
  std::string impute_method{"mean"};        // R CLI: "best_guess"
  double      spa_cutoff{2.0};
  double      p_cutoff_for_firth{0.01};
  // Firth beta defaults to on for binary, off otherwise (the previous derived
  // behaviour). Set explicitly to override.
  int         firth_beta{-1};               // -1 = derive from trait, 0 = off, 1 = on

  // Categorical VR (R: --isCateVarianceRatio + --cateVarRatioMinMACVecExclude / --cateVarRatioMaxMACVecInclude)
  // When isCateVarianceRatio=true, markers are partitioned into MAC bins
  // (numCate = cateVarRatioMinMACVecExclude.size()). For bin i in [0..numCate-2]:
  //   MAC > cateVarRatioMinMACVecExclude[i]  AND  MAC <= cateVarRatioMaxMACVecInclude[i]
  // For the last bin (numCate-1): MAC > cateVarRatioMinMACVecExclude[numCate-1] (open upper)
  // Output is multi-line: "<vr>\tnull\t<k>", "<vr>\tnull_noXadj\t<k>", "<vr>\tsparse\t<k>" per bin k (1-based).
  bool   isCateVarianceRatio{false};
  std::vector<double> cateVarRatioMinMACVecExclude{10.0, 20.5};  // R default
  std::vector<double> cateVarRatioMaxMACVecInclude{20.5};        // R default (length = numCate - 1)
  std::vector<int>    cateVarRatioIndexVec{};                    // optional: 0 => default 1.0; 1 => estimate (per bin)

  // ---- Fused variance ratio (quantitative + sparse-GRM fit only) ----------
  // fit.fused_variance_ratio. Anchors the per-bin variance ratio on the exact
  // closed form tr(P Psi)/tr((I-H) Psi) -- every trace computed exactly by
  // connected components of the sparse Psi, no Hutchinson probes -- and spends
  // markers only on testing the per-bin correction delta = ratio / anchor
  // against 1.  See optimization/FUSED_VARIANCE_RATIO.md.  Off by default;
  // with it off nothing below is read and the output is byte-identical.
  bool   fused_variance_ratio{false};
  int    fused_vr_markers{10};      // markers per bin for the delta test (0 = pure closed form)
  int    fused_vr_max_markers{30};  // ceiling when the se rule keeps asking for more
  double fused_vr_delta_se{0.002};  // stop when se(delta)/delta <= this
  double fused_vr_delta_z{2.0};     // keep delta only when |delta-1| > z*se(delta)
  int    fused_vr_max_block{1024};  // refuse if a connected component is larger

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

  // Categorical covariate column names (must be subset of covar_cols)
  std::vector<std::string> q_covar_cols;

  // Skip model fitting (load pre-existing model file)
  bool skip_model_fitting{false};
  std::string model_file;   // path to pre-existing .rda/.json model (when skip_model_fitting=true)

  // Dry-run: validate inputs only, no genotype loading or solver
  bool dry_run{false};
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

  // Snapshot of the original full-covariate design (X_full has n*p_full entries,
  // row-major) BEFORE the covariate_offset path collapses X to an intercept-only
  // column. p_full=0 means "no snapshot taken / not in covariate_offset mode".
  // Step-2 reads X.arma / XVX.arma etc. to project g_tilde = g − X(X'X)⁻¹X'g; the
  // collapsed intercept-only X can only mean-center, which inflates Tstat at
  // markers in LD with the omitted covariates.
  std::vector<double> X_full;
  int p_full{0};

  // Initial-GLM fixed-effect coefficients on the FULL (uncollapsed) design,
  // length p_full. Only populated on the covariate_offset path. On that path
  // the GLMM re-estimates the INTERCEPT ONLY (all other covariate effects are
  // frozen inside the offset), so this is what lets us report a full-length
  // alpha in nullmodel.json instead of just the intercept. Reporting only —
  // nothing numerical is derived from it.
  std::vector<double> beta_full;

  std::vector<double> y;          // length n
  std::vector<double> offset;     // optional, length n or empty
  std::vector<double> event_time; // optional, length n or empty

  std::vector<std::string> iid;   // length n, order must match rows
};

// Score-null payload for Step 2 (matches R's obj.noK)
struct ScoreNullPack {
  int n{0}, p{0};
  std::string trait_type;         // "binary" | "quantitative" | "survival"

  // N-length vectors
  std::vector<double> V;          // mu*(1-mu) for binary; 1/tau0 for quant
  std::vector<double> mu;         // fitted values
  std::vector<double> res;        // residuals (y - mu)
  std::vector<double> y;          // phenotype

  // p-length vectors
  std::vector<double> S_a;        // colSums(X ⊙ residuals)

  // Matrices (row-major flat arrays)
  std::vector<double> XV;         // p x n
  std::vector<double> XVX;        // p x p
  std::vector<double> XVX_inv;    // p x p
  std::vector<double> XXVX_inv;   // n x p
  std::vector<double> XVX_inv_XV; // n x p
  std::vector<double> X;          // n x p (design matrix)
};

// Result of a LOCO batch run. `chroms` lists the autosomes (1-based) for which
// a chr<j>/ directory was actually written; `obj_noK` is parallel to it.
struct LocoBatchOut {
  std::vector<int>           chroms;
  std::vector<ScoreNullPack> obj_noK;
};

struct FitNullResult {
  // Model parameters / summaries
  std::vector<double> theta;   // variance components
  std::vector<double> alpha;   // fixed effects (original design scale)
  std::vector<double> offset;  // final GLMM offset (length n)

  // Converged linear predictor and fitted values (length n).
  // CRITICAL: these are the GLMM-converged eta/mu (they include the BLUP of the
  // random effect b̂); they must NOT be recomputed downstream as X*alpha + offset,
  // which would drop b̂ and corrupt step-2 .arma exports (mu, res, V, S_a).
  // Populated by binary_glmm_solver / quant_glmm_solver in glmm.cpp.
  std::vector<double> eta;     // length n; for binary == log-odds, for quant == fitted mean
  std::vector<double> mu;      // length n; for binary == 1/(1+exp(-eta)), for quant == eta

  // Flags
  bool loco{false};
  bool lowmem_loco{false};
  bool converged{false};       // Step 13: true if solver converged before maxiter
  int  iterations{0};          // number of outer iterations run

  // Artifacts (paths); may be empty if not written
  std::string model_rda_path;     // main model artifact path (json/rds/etc.)
  std::string vr_path;            // <prefix>.varianceRatio.txt
  std::string markers_out_path;   // <prefix>_<N>markers.SAIGE.results.txt

  // ... your existing fields (alpha, theta, offset, loco flags, paths, etc.)

  // Baseline score-null (R's obj.noK) in pure C++
  ScoreNullPack obj_noK;

  // Per-chromosome LOCO score-nulls, parallel to `loco_chroms`. Populated only
  // when the LOCO batch actually ran.
  std::vector<ScoreNullPack> loco_obj_noK;

  // Autosomes (1-based, ascending) for which chr<j>/ artifacts were written.
  // Empty when `loco` is false. Serialized to nullmodel.json as "loco_chroms".
  std::vector<int> loco_chroms;
};



// Optional library-level orchestrator (implemented in saige_null.cpp, if you use it)
// Tier-2 lockstep: fit P phenotypes together. Same per-trait outputs as P
// separate fit_null() calls, but the P AI-REML loops advance in lockstep so
// the fixed-effect PCG solves of all still-active traits are issued as ONE
// batched multi-Sigma solve (one pass of the packed 2-bit matrix per 8-column
// block instead of one per column). Enabled by fit.multi_lockstep; falls back
// to the per-trait path for P==1.
std::vector<FitNullResult> fit_null_multi(const FitNullConfig& cfg,
                                          const std::vector<Paths>& paths,
                                          const std::vector<Design>& designs);

FitNullResult fit_null(const FitNullConfig& cfg,
                       const Paths& paths,
                       const Design& design);

} // namespace saige
