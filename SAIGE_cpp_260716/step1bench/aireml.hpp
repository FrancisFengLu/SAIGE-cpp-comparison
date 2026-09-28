// aireml.hpp -- a small, legible AI-REML whose only job is to give the three
// solvers identical work.
//
// It is NOT a reimplementation of SAIGE's null-model fit. It computes the same
// quantities in the same order (Sigma^-1 y, Sigma^-1 X, tr(P Psi), tr(P), the
// AI matrix, the same covariate projection), with the same Sigma, and stops on
// the same relative-change rule -- but it drops SAIGE's LOCO handling, its
// variance-ratio stage, its trace-CV growth loop and the stale-W nesting that
// makes upstream's M a mixed operator. What it must do, and what is gated, is
// reach the SAME tau whichever solver is plugged in, and land close to
// production on the same data.
#ifndef STEP1BENCH_AIREML_HPP
#define STEP1BENCH_AIREML_HPP

#include <armadillo>
#include <string>
#include <vector>

#include "sigma_solver.hpp"

namespace step1bench {

struct FitConfig {
    bool   binary      = false;       // false = quantitative
    bool   exactTrace  = false;       // false = Hutchinson
    int    probes      = 30;
    int    maxiter     = 20;
    double tol         = 0.02;        // relative change on tau
    double tolCoef     = 1e-5;        // inner IRLS, alpha relative change
    unsigned long long seed = 1;
};

struct FitResult {
    arma::fvec tau;
    arma::fvec alpha;
    int  iters = 0;
    bool converged = false;
    long long innerIRLS = 0;
    double    fitSeconds = 0.0;
    std::vector<std::vector<double>> traj;   // it, tau0, tau1, score..., relchg
};

// psi is the sparse GRM itself (both triangles, fp64) -- the same matrix the
// solvers were built from. It is used for Psi*v and, in the exact-trace path,
// for Psi*Sigma_iX.
FitResult run_aireml(SigmaSolver& sol, const arma::sp_mat& psi,
                     const arma::fvec& y, const arma::fmat& X,
                     const FitConfig& cfg);

}  // namespace step1bench

#endif  // STEP1BENCH_AIREML_HPP
