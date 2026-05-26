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

} // namespace saige
