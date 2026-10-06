// gpu_spa.hpp — the binary-trait saddlepoint approximation on the device, for
// the (marker, trait) pairs the step-2 gate flags. Config key gpuSpa (needs
// useGPU + gpuBinary); default off. No CUDA types cross this boundary; the
// no-CUDA build links gpu_spa_stub.cpp and spaCreate() returns nullptr.
//
// ---------------------------------------------------------------------------
// What it computes
// ---------------------------------------------------------------------------
// SAIGEClass::getMarkerPval (saige_test.cpp) does, for a binary pair whose
// |Tstat|/sqrt(var1) exceeds SPA_Cutoff:
//
//   g~  = g - XXVX_inv (XV g)                       getadjGFast, carriers only
//   m1  = mu . g~
//   q   = Tstat / sqrt(var1/var2) + m1,  qinv mirrored about m1
//   fast variant iff |{g == 0}| / N >= 0.5:  sums over carriers only, the
//        non-carriers folded into NAmu = m1 - sum_c mu g~, NAsigma = var2 -
//        sum_c mu(1-mu) g~^2  (SPA_fast); otherwise sums over all N (SPA)
//   two Newton root-finds of K1(t) = q (resp. qinv) from t = 0, with the
//        prevJump halving on a sign flip, |dt| < eps^(1/4), maxiter 1000
//   two saddlepoint tail probabilities; the wrapper's rule that a
//        non-saddle withdraws convergence; |p1| + |p2|, or add_logp in the
//        log domain when the score p-value had underflowed
//
// spa_binary.cpp / spa.cpp are the reference; the kernel carries that control
// flow over statement by statement (gpu_spa.cu), one CUDA block per pair with
// the block's 256 threads sharing the length-N sums, so each pair iterates as
// many times as SAIGE's own loop would. What differs is the order of the
// additions inside each sum (tree instead of sequential) and CUDA's exp / log
// / erfc against glibc's and boost's; GPU_WHOLE_PIPELINE.md section 2.2
// measured 0 of 128 roots with a different iteration count and a p-value
// relative difference of 5e-14. fp64 throughout: in fp32 the squared
// denominator of K2 underflows.
//
// What stays on the host, after spaRun(): the quantile step that turns the
// SPA p-value into seBeta and can itself withdraw convergence, the p == 0
// rule, the "%.6E" / "%.1fE%d" string, and the Firth decision -- a pair whose
// SPA p-value asks for Firth is handed to the CPU scalar path whole (main.cpp,
// mainMarkerMTGpu). Survival traits are not handled (the gate refuses them).
//
// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------
// The genotype column is read from the reducer's resident packed rows and
// per-slot dosage tables (gpu_step2.hpp devicePacked / deviceLut), so g is the
// same decoded value the GEMM saw and the host never ships a dense vector.
// Per trait: mu (N), XV (p x N as the model stores it: sample i's p values are
// contiguous) and XXVX_inv (N x p column-major), uploaded once at spaCreate().
// Per pair: slot, trait, the fast / logp flags and Tstat, var1, var2, pval_noadj
// from the batch kernel -- the GPU GEMM's numbers, not a CPU recompute.
#pragma once

#include <cstddef>
#include <cstdint>

namespace saige {
namespace gpu2 {

struct Reducer;
struct Spa;   // opaque; defined in gpu_spa.cu

struct SpaTraitArgs {
    const double* mu       = nullptr;   // N
    const double* XV       = nullptr;   // p x N column-major
    const double* XXVX_inv = nullptr;   // N x p column-major
    int p = 0;                          // covariate count incl. intercept, <= SPA_PMAX
};

constexpr int SPA_PMAX = 8;

struct SpaCreateArgs {
    int device = 0;
    int N = 0;
    int nTraits = 0;
    const SpaTraitArgs* traits = nullptr;   // nTraits entries, indexed by the pair's trait
    int maxPairs = 0;                       // per spaRun() call
    double tol = 0.0;                       // |dt| convergence, eps^(1/4) in SAIGE
    int maxiter = 1000;
    int blocks = 256;                       // resident blocks; scratch is blocks x N doubles
};

struct SpaPairIn {
    int    slot;     // reducer slot of the marker (last reduce())
    int    trait;    // index into SpaCreateArgs::traits
    int    fast;     // 1: SPA_fast variant (|{g==0}|/N >= 0.5)
    int    logp;     // 1: pno is log(p) and the result is log(p)
    double Tstat, var1, var2, pno;
};

struct SpaPairOut {
    double pval;     // SPA p-value (log when logp), or pno when not converged
    int    conv;     // the wrapper's Isconverge
    int    s1, s2;   // isSaddle of each tail (-1: roots did not converge)
    int    niter1, niter2;
    double root1, root2;
    double m1;
};

// nullptr on any failure; the caller then keeps SPA on the CPU.
Spa* spaCreate(const SpaCreateArgs& t_args);
void spaDestroy(Spa* t_s);

// Pinned tables of maxPairs entries, owned here.
SpaPairIn*  spaIn(Spa* t_s);
SpaPairOut* spaOut(Spa* t_s);

// Solve pairs [0, t_nPairs) of spaIn() against the reducer's last reduce();
// results in spaOut(). Synchronous. false on a CUDA error (results undefined).
bool spaRun(Spa* t_s, const Reducer* t_r, int t_nPairs, int t_devSet = -1);  // t_devSet: the reducer's device set (-1 = last reduce())

// Cumulative device seconds in the kernel, and pairs solved, since spaCreate().
void spaTimings(const Spa* t_s, double* t_kernel, long long* t_pairs);
std::size_t spaDeviceBytes(const Spa* t_s);

}  // namespace gpu2
}  // namespace saige
