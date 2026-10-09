// gpu_er.hpp — SAIGE's ER (exact test for binary traits, er_binary.cpp) on the
// device for the (marker, trait) pairs of the GPU path that take it. Config key
// gpuER; default off. No CUDA types cross this boundary; the no-CUDA build links
// gpu_er_stub.cpp and erCreate() returns nullptr.
//
// ---------------------------------------------------------------------------
// What it computes (er_binary.cpp, SKATExactBin_Work with the production
// arguments NResampling = 2e6, ExactMax = 1e4, epsilon = 1e-6, an empty resout)
// ---------------------------------------------------------------------------
// For a marker with k carriers (g_i != 0; at most ER::kDeviceMaxCarriers = 20,
// so all 2^k case/control assignments of the carriers are enumerated and no
// random number is drawn), with p_i = mu_i at the carriers and n, ncase the
// trait's sample and case counts:
//
//  1. Null distribution of the number of cases among the carriers: carriers
//     grouped by p (10 bins of width 0.1, p >= 1 clamped to 0.999), each bin's
//     odds = mean(p)/(1-mean(p)) relative to the non-carriers' odds; a weighted
//     hypergeometric over the groups (HyperGeo::Recursive, log scale) gives
//     prob[j], j = 0..k.
//  2. Every assignment S of j cases to the carriers: weight prod_{i in S} odds_i
//     normalised within j and scaled by prob[j]; statistic
//     T(S) = (sum_i g_i (1{i in S} - p_i))^2.
//  3. p = P(T >= T_obs) - P(T == T_obs) / 2 under that distribution, ties being
//     |T - T_obs| <= 1e-6; T_obs from the carriers that are cases.
//
// One thread per pair, every sum in the CPU's order: the non-carriers' mean of
// mu with armadillo's two interleaved accumulators, log C(n-k, ncase-j) as the
// CPU's running sum of log(n - d + 1) - log(d) (the k+1 values are prefixes of
// one loop), the hypergeometric leaves in depth-first order, the assignments in
// the CPU's order (lexicographic subsets of the cases, or of the controls for
// j > k/2 + 1). The CPU's arrays of 2^k weights and statistics are not stored;
// the three passes over them (stratum sums, total, p) each regenerate the
// values with the same arithmetic. exp and log are device ports of the glibc
// functions the CPU calls (gpu_er_glibc.cuh); log of an integer comes from a
// table the host fills with the host's log. fp64, --fmad=false.
//
// What stays on the host: the score test that decides whether the pair takes
// ER at all (unchanged, getMarkerPval with a PerMarkerCtx::erDefer), and after
// erRun() the statements that format the p-value and back-calculate seBeta
// (SAIGE::erFinish), the Firth decision and the route bits.
//
// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------
// Per trait, uploaded once at erCreate(): mu (the trait's own n samples, as
// SAIGEClass::m_mu), n, ncase. A log table log(0..maxN) from the host. Per pair:
// trait, k and the carriers (position in the trait's vector, ascending; dosage;
// case flag), 16 B + 13 B per carrier.
//
// Precision (config key gpuPrecisionER, gpu_precision.hpp): ErCreateArgs::
// precision, FP64 (default, everything above, bit-identical to the CPU) or
// FP32 (gpu_er_fp32.cu: the same assignments, strata and two kernels in fp32,
// log domain -- log Fisher weights, compensated log-sum-exp per stratum, one
// pass instead of three, no O(N) loop per pair; header of that file). Only
// modes erSupports() accepts may be passed. CONTRACT for any
// non-fp64 variant: ErPairOut::pval is handed back as double (the fp32 value
// converted); the host's erFinish ([0, 1] guard, formatting, seBeta) reads it
// unchanged. main.cpp runs the bit-for-bit startup self-check (gpuErSelfCheck:
// erMathCheck + 24 exact tests) only when ER is fp64; with any other mode it
// logs that the check is skipped, so an fp32 variant does not have to be, and
// cannot be, bit-identical. erMathCheck itself stays the fp64 glibc port.
#pragma once

#include <cstddef>
#include <cstdint>

#include "gpu_precision.hpp"

namespace saige {
namespace gpu2 {

struct Er;   // opaque; defined in gpu_er.cu

struct ErTraitArgs {
    const double* mu = nullptr;   // n
    int n = 0;
    int ncase = 0;
};

struct ErCreateArgs {
    int device = 0;
    int nTraits = 0;
    const ErTraitArgs* traits = nullptr;
    const double* logTable = nullptr;   // log((double)i), i = 0..maxN (entry 0 unused)
    int maxN = 0;                       // >= every trait's n
    Prec precision = Prec::FP64;        // gpuPrecisionER; see the contract above
};

// True when erCreate() accepts ErCreateArgs::precision = t_p. main.cpp asks
// first and stops the run on false ("useGPU: ER precision <mode> is not
// implemented yet").
bool erSupports(Prec t_p);

struct ErPairIn {
    int trait;     // index into ErCreateArgs::traits
    int k;         // carriers, 1..kMaxCarriers
    long long off; // first carrier in the erRun() carrier arrays
};

struct ErPairOut {
    double pval;   // pval[0] - pval_same[0] / 2, before erFinish's [0, 1] guard
};

constexpr int kErMaxCarriers = 20;   // == ER::kDeviceMaxCarriers

// nullptr on any failure; the caller then keeps ER on the CPU.
Er*  erCreate(const ErCreateArgs& t_args);
void erDestroy(Er* t_e);

// Exact p-values of pairs [0, t_nPairs): carrier c of pair q is entry
// in[q].off + c of carIdx / carG / carCase. Synchronous. false on a CUDA error.
bool erRun(Er* t_e, const ErPairIn* t_in, int t_nPairs,
           const uint32_t* t_carIdx, const double* t_carG, const unsigned char* t_carCase,
           long long t_nCar, ErPairOut* t_out);

// Cumulative device seconds (kernel only) and pairs since erCreate().
void erTimings(const Er* t_e, double* t_kernel, long long* t_pairs);
const char* erLastError();

// Test hook (tools/er_gpu_check): the device exp and log of x[0..n).
bool erMathCheck(int t_device, const double* t_x, long long t_n, double* t_exp, double* t_log);

}  // namespace gpu2
}  // namespace saige
