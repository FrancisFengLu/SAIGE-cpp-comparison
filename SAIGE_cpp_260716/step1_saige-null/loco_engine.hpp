#pragma once
#include "saige_null.hpp"
#include <vector>

namespace saige {

// Runs leave-one-chromosome-out re-solves of the fixed effects with tau held
// fixed at `theta`, and writes <out_prefix>/chr<j>/ for each autosome present.
//
// Returns true iff at least one chr<j>/ directory was written.
//
// `design`  : the design the GLMM was fit on (X = fit-space, possibly QR
//             transformed and/or collapsed to intercept-only). design.X_full /
//             design.p_full, when set, carry the original full-covariate X used
//             for the obj.noK export (matches R's `Xorig` under
//             isCovariateOffset).
// `alpha_fit`: fixed-effect coefficients in the SAME space as design.X (i.e.
//             BEFORE any QR back-transform).
// `eta_full`: converged full-genome linear predictor from the null fit.
bool run_loco_batch(const Paths& paths,
                    const FitNullConfig& cfg,
                    const LocoRanges& chr,
                    const Design& design,
                    const std::vector<double>& theta,
                    const std::vector<double>& alpha_fit,
                    const std::vector<double>& offset,
                    const std::vector<double>& eta_full,
                    LocoBatchOut& out);

// Registers run_loco_batch() with NullModelEngine's g_loco_batch hook.
void register_default_loco_batch();

} // namespace saige
