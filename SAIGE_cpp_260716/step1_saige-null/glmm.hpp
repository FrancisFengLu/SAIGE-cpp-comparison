#pragma once

#include <vector>
#include <string>
#include "saige_null.hpp"   // Paths, FitNullConfig, Design, FitNullResult

namespace saige {

// Register these adapters with NullModelEngine’s hooks.
void register_default_solvers();

// (Exported for unit tests, optional to use directly)
FitNullResult binary_glmm_solver(const Paths&,
                                 const FitNullConfig&,
                                 const Design&,
                                 const std::vector<double>& offset,
                                 const std::vector<double>& beta_init);

FitNullResult quant_glmm_solver (const Paths&,
                                 const FitNullConfig&,
                                 const Design&,
                                 const std::vector<double>& offset,
                                 const std::vector<double>& beta_init);

// If your null_model_engine.hpp doesn’t expose these, keep these forward decls:
using BinarySolverFn = FitNullResult (*)(const Paths&, const FitNullConfig&, const Design&,
                                         const std::vector<double>&, const std::vector<double>&);
using QuantSolverFn  = FitNullResult (*)(const Paths&, const FitNullConfig&, const Design&,
                                         const std::vector<double>&, const std::vector<double>&);

void register_binary_solver(BinarySolverFn);
void register_quant_solver (QuantSolverFn);

// ---- Tier-2: lockstep multi-phenotype AI-REML driver -----------------------
// Fits P traits that share one genotype object / GRM by advancing all of their
// AI-REML iterations together. The only thing actually shared is the psi*B
// product inside the fixed-effect (coefficient) solve: the P*(1+q) columns
// [Y_p | X_p] go through ONE getPCGofSigmaAndMatrix_multiSigma call instead of
// P separate ones, so the packed 2-bit matrix is streamed once per 8-column
// block rather than once per column. Everything else (IRLS build, AI score,
// trace probes, tau update, convergence test) stays per trait and reproduces
// the scalar driver's arithmetic exactly.
//
// All P designs must share n (the caller guarantees one sample set); traits may
// have different p. Trait p is FROZEN the moment its own convergence test
// fires: it leaves the batch, is finalized on the scalar path (the final
// Get_Coef re-solve), and its result no longer moves while the others iterate.
std::vector<FitNullResult>
multi_glmm_solver(const std::vector<Paths>& paths,
                  const FitNullConfig& cfg,
                  const std::vector<const Design*>& designs,
                  const std::vector<const std::vector<double>*>& offsets,
                  const std::vector<const std::vector<double>*>& beta_inits,
                  bool is_binary);

} // namespace saige
