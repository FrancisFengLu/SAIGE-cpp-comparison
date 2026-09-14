#pragma once
#include <memory>
#include <vector>
#include "saige_null.hpp"

namespace saige {

// Opaque per-trait state produced by NullModelEngine::prep() and consumed by
// solve() / export_result(). Defined in null_model_engine.cpp so this header
// stays free of Eigen (NullPrep holds the QR map).
struct NullPrep;

// ======= Solver / LOCO hook signatures =======

// Binary (incl. survival) null GLMM solver.
// Provide an implementation and register via register_binary_solver().
using BinarySolverFn = FitNullResult (*)(const Paths&,
                                         const FitNullConfig&,
                                         const Design&,
                                         const std::vector<double>& /*offset*/,
                                         const std::vector<double>& /*beta_init (optional)*/);

// Quantitative null GLMM solver.
// Provide an implementation and register via register_quant_solver().
using QuantSolverFn = FitNullResult (*)(const Paths&,
                                        const FitNullConfig&,
                                        const Design&,
                                        const std::vector<double>& /*offset*/,
                                        const std::vector<double>& /*beta_init (optional)*/);

// Batched LOCO runner across chromosomes.
// Provide an implementation and register via register_loco_batch();
// loco_engine.cpp supplies the default (saige::run_loco_batch).
//
// Returns true iff LOCO actually ran and the chr<j>/ artifacts were written.
// `alpha` is in the SAME space as `design.X` (i.e. before any QR back-transform)
// and `eta` is the converged full-genome linear predictor.
using LocoBatchFn = bool (*)(const Paths&,
                             const FitNullConfig&,
                             const LocoRanges&,
                             const Design&,
                             const std::vector<double>& /*theta*/,
                             const std::vector<double>& /*alpha (fit space)*/,
                             const std::vector<double>& /*offset*/,
                             const std::vector<double>& /*eta (converged)*/,
                             LocoBatchOut& /*out*/);

// Registration APIs (call these once during initialization).
void register_binary_solver(BinarySolverFn fn);
void register_quant_solver (QuantSolverFn  fn);
void register_loco_batch   (LocoBatchFn    fn);

// ======= NullModelEngine =======

class NullModelEngine {
public:
  NullModelEngine(const Paths& paths,
                  const FitNullConfig& cfg,
                  const LocoRanges& chr);

  // Fits the null model:
  //  - optional QR on covariates
  //  - baseline GLM (logistic/gaussian)
  //  - constructs final GLMM offset (offset-mode vs. explicit offset)
  //  - calls registered GLMM solver (binary or quantitative)
  //  - optional LOCO batch via registered hook
  //  - writes a compact model artifact path into FitNullResult
  FitNullResult run(const Design& design);

  // Three-phase form of run():
  //     auto p = prep(design);
  //     FitNullResult r = export_result(*p, solve(*p));
  // is exactly run(design). The split exists so a lockstep multi-phenotype
  // driver can prep P traits, solve all P together, then export each — see
  // solve_multi() below.
  std::shared_ptr<NullPrep> prep(const Design& design);
  FitNullResult             solve(NullPrep& prep);
  FitNullResult             export_result(NullPrep& prep, FitNullResult out);

private:
  Paths        paths_;
  FitNullConfig cfg_;
  LocoRanges   chr_;
};

// Read-only views into a NullPrep, so callers outside null_model_engine.cpp
// (the lockstep driver) can see what prep() produced without this header
// having to expose the Eigen-carrying NullPrep definition.
const Design&              nullprep_design   (const NullPrep&);
const std::vector<double>& nullprep_offset   (const NullPrep&);
const std::vector<double>& nullprep_beta_init(const NullPrep&);
bool                       nullprep_is_binary(const NullPrep&);

} // namespace saige
