// gpu_firth.hpp — Firth-corrected logistic regression on the device for the
// binary (marker, trait) pairs whose final p-value asks for it. Config key
// gpuFirth (needs gpuSpa); default off. No CUDA types cross this boundary; the
// no-CUDA build links gpu_firth_stub.cpp and firthCreate() returns nullptr.
//
// ---------------------------------------------------------------------------
// What it computes
// ---------------------------------------------------------------------------
// SAIGEClass::getMarkerPval (saige_test.cpp), for a binary pair whose p-value
// (SPA's when SPA ran, the normal approximation's otherwise) is <= pCutoffforFirth:
//
//   g~  = g - XXVX_inv (XV g)                        getadjGFast, carriers only
//   fast_logistf_fit_simple(x = [1, g~], y, offset, firth = true, init = 0,
//                           maxit 50, maxstep 15, xconv = gconv = 1e-5)
//   seBeta = |beta| / |qnorm(p / 2)|                 the host, as before
//
// That fit is Newton on the Firth-penalised likelihood with X = [1, g~]: every
// quantity of one step -- the 2x2 Fisher information, the hat-matrix diagonal
// h_i = w_i x_i' F^-1 x_i and the adjusted score U* = X'((y - pi) + h (1/2 - pi))
// -- is a sum over the N samples of a function of (g~_i, y_i, offset_i, alpha,
// beta), so one pass over N with nine accumulators replaces the N x 2 weighted
// matrix, the QR and the hat vector (FIRTH_GPU.md). The step cap
// (max|delta| <= maxstep), the stopping rule (max|delta| <= xconv and
// max|U*| <= gconv, with U* from before the step), the maxit exit that SAIGE
// also reports as "converged", and the singular-information exit are carried
// over statement by statement, so each block iterates exactly as many times as
// the CPU would for that pair. One CUDA block per pair, the block's 256
// threads sharing each length-N sum; what differs from the CPU is the order of
// the additions inside each sum (tree instead of sequential), which the
// prototype measured at |d beta| <= 4.4e-15 on 24,870 converged pairs and
// <= 2.3e-12 on the MAC = 1 pairs that oscillate to maxit. fp64, --fmad=false.
//
// What stays on the host, after firthRun(): the flip sign on beta, seBeta from
// the p-value, the Firth counters and route bits. ER pairs (MAC <= the ER
// cutoff) and pairs that need the fast-test recompute keep the CPU scalar path.
//
// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------
// The genotype column is read from the reducer's resident packed rows and
// per-slot dosage tables (gpu_step2.hpp devicePacked / deviceLut), so g is the
// same decoded value the GEMM and the device SPA saw. Per trait: y (N),
// offset (N), XV (p x N as the model stores it: sample i's p values are
// contiguous) and XXVX_inv (N x p column-major), uploaded once at firthCreate().
// Per pair: slot and trait only (16 B).
#pragma once

#include <cstddef>
#include <cstdint>

namespace saige {
namespace gpu2 {

struct Reducer;
struct Firth;   // opaque; defined in gpu_firth.cu

struct FirthTraitArgs {
    const double* y        = nullptr;   // N
    const double* offset   = nullptr;   // N
    const double* XV       = nullptr;   // p x N column-major
    const double* XXVX_inv = nullptr;   // N x p column-major
    int p = 0;                          // covariate count incl. intercept, <= FIRTH_PMAX
    // FirthCreateArgs::ownSamples only: the trait's samples among the N (1 bit
    // per sample, bit i of word i >> 6); nullptr = all. y / offset / XV /
    // XXVX_inv are embedded at N with zeros outside; a sample outside the
    // mask has dosage 0 and is left out of the fit's sums.
    const uint64_t* mask   = nullptr;
};

constexpr int FIRTH_PMAX = 8;

struct FirthCreateArgs {
    int device = 0;
    int N = 0;
    int nTraits = 0;
    const FirthTraitArgs* traits = nullptr;   // nTraits entries, indexed by the pair's trait
    int maxPairs = 0;                         // per firthRun() call
    int    maxit   = 50;                      // fast_logistf_fit_simple's arguments
    double maxstep = 15.0;
    double xconv   = 1e-5;
    double gconv   = 1e-5;
    int blocks = 256;                         // resident blocks; scratch is blocks x N doubles
    // Per-trait sample sets (default off): masks per trait, and pair k's own
    // 4-entry dosage table in firthPairLut() at 4k instead of the slot's.
    int ownSamples = 0;
};

struct FirthPairIn {
    int slot;     // reducer slot of the marker (last reduce())
    int trait;    // index into FirthCreateArgs::traits
};

struct FirthPairOut {
    double beta;      // the fit's beta for g~ (the flipped allele when the column was flipped)
    double alpha;     // intercept
    double se;        // sqrt of the fit's cov(beta, beta); the host overwrites seBeta from p
    int    conv;      // fast_logistf_fit_simple's isfirthconverge (maxit also counts)
    int    strict;    // 1: stopped on the tolerances, 0: stopped at maxit
    int    niter;
    int    singular;  // 1: the information matrix could not be inverted (conv = 0)
};

// nullptr on any failure; the caller then keeps Firth on the CPU.
Firth* firthCreate(const FirthCreateArgs& t_args);
void   firthDestroy(Firth* t_f);

// Pinned tables of maxPairs entries, owned here.
FirthPairIn*  firthIn(Firth* t_f);
FirthPairOut* firthOut(Firth* t_f);
// ownSamples: pinned maxPairs x 4 doubles (pair k's table at 4k); else nullptr.
double*       firthPairLut(Firth* t_f);

// Fit pairs [0, t_nPairs) of firthIn() against the reducer's last reduce();
// results in firthOut(). Synchronous. false on a CUDA error (results undefined).
bool firthRun(Firth* t_f, const Reducer* t_r, int t_nPairs, int t_devSet = -1);  // t_devSet: the reducer's device set (-1 = last reduce())

// Cumulative device seconds in the kernel, and pairs fitted, since firthCreate().
void firthTimings(const Firth* t_f, double* t_kernel, long long* t_pairs);
std::size_t firthDeviceBytes(const Firth* t_f);
const char* firthLastError();

}  // namespace gpu2
}  // namespace saige
