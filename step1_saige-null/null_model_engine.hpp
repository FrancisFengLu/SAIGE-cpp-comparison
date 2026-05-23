#pragma once
#include "saige_null.hpp"

namespace saige {

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

// Batched LOCO runner across chromosomes (optional).
// Provide an implementation and register via register_loco_batch().
using LocoBatchFn = void (*)(const Paths&,
                             const FitNullConfig&,
                             const LocoRanges&,
                             const Design&,
                             const std::vector<double>& /*theta*/,
                             const std::vector<double>& /*alpha*/,
                             const std::vector<double>& /*offset*/);

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

private:
  Paths        paths_;
  FitNullConfig cfg_;
  LocoRanges   chr_;
};

} // namespace saige
