// stats_selfcheck.hpp -- start-up checks for the per-pair statistics on the
// device (config key gpuDeviceStats; gpu/gpu_step2.hpp "Per-pair statistics").
//
// The device kernel reproduces, bit for bit, what the host tail
// (scoreTestBatchMTBinPre, saige_mt.cpp) computes from the same C1 / C2
// columns: S = (g'res - S_a'Z) / tau0 and var2 = Z'XVX Z + (g^2)'mu2 - 2 W'Z.
// On the host the three length-p contractions are not under this code's
// control: XVX Z is an OpenBLAS dgemm, S_a'Z an OpenBLAS dgemv, and the two
// Hadamard column sums are armadillo's two-accumulator loop compiled by the
// build's own compiler, which contracts a*b+c into fma where it likes (g++
// does so for C++ even under -std=c++17). The order and the fma pattern of
// each decide the last bit. So:
//
//   1. identifyHostPatterns runs the real scoreTestBatchMTBinPre on random
//      data shaped like this run's traits and matches its S / var2 against
//      the candidate orders dotPattern / hadSumPattern enumerate (this file
//      is compiled with -ffp-contract=off, so "no fma" here means no fma;
//      the fma candidates call std::fma explicitly).
//   2. deviceSelfTest uploads random C1 / C2 / VR, runs the device kernel
//      with those patterns and requires every S / var2 to match the host
//      function bit for bit, every gate bit to match the host's decision
//      from boost's p, and the device p to be within pRelTol of boost's.
//
// Any failure leaves the switch off for the run (main.cpp prints why), as
// the gpuER self-check does. The candidate set is small on purpose: a length
// p <= 8 contraction has few plausible shapes, and an unlisted one is a
// clean failure, not a wrong number.
#ifndef SAIGE_STATS_SELFCHECK_HPP
#define SAIGE_STATS_SELFCHECK_HPP

#include <string>
#include <vector>

#include "saige_mt.hpp"
#include "gpu_step2.hpp"

namespace SAIGE {
namespace devstats {

// sum_k a[k] b[k], k < p:
//   0  no fma, left to right
//   1  accumulator from +0, one fma per term (OpenBLAS micro-kernel shape)
//   2  fma(a0, b0, a1 b1), then a chain of fma (g++'s contraction of t0+t1+t2)
//   3  (a0 b0 + a1 b1) unfused, then a chain of fma
//   4  fma(a0, b0, a1 b1), then unfused adds
double dotPattern(int t_pat, const double* t_a, const double* t_b, int t_p);
// armadillo's sum(x % y, 0) over p rows: even rows into val1, odd into val2,
// val1 + val2. 0: the product rounded then added; 1: fma into the accumulator.
double hadSumPattern(int t_pat, const double* t_x, const double* t_y, int t_p);

struct Patterns { int xz = -1, saz = -1, zxz = -1, gwz = -1; };

// S and var2 of one pair from its columns, in the given patterns, exactly as
// the device kernel forms them. XVX row-major with stride STATS_PMAX.
void refPairStats(const Patterns& t_pat, int t_p, const double* t_XVX, const double* t_Sa,
                  double t_tau0, const double* t_z, const double* t_w, double t_gr, double t_g2,
                  double& t_S, double& t_var2);

// Step 1. Empty string = found (t_out filled; t_detail says which); otherwise
// the reason. t_traits: internal indices of the binary traits to score.
std::string identifyHostPatterns(const MTContext& t_ctx, const std::vector<int>& t_traits,
                                 Patterns& t_out, std::string& t_detail);

// The device's per-trait table. t_enabled[b] = 0 marks binary trait b (by
// binIdx) as not scored on the device. t_oW / t_oR: C1 column offsets of the
// WXstack block and of RES (mainMarkerMTGpu).
void fillStatsTraits(const MTContext& t_ctx, const std::vector<int>& t_binTraits,
                     const std::vector<char>& t_enabled, int t_oW, int t_oR,
                     std::vector<saige::gpu2::StatsTrait>& t_out);

// Step 2. The reducer must have statsSetup done and no reduce() yet. Empty
// string = pass; t_detail carries the counts either way.
std::string deviceSelfTest(saige::gpu2::Reducer* t_R, const MTContext& t_ctx,
                           const std::vector<int>& t_binTraits, const std::vector<char>& t_enabled,
                           int t_K1, int t_K2, int t_oW, int t_oR, int t_maxSlots, int t_Bblk,
                           double t_statCutoff, double t_pRelTol, std::string& t_detail);

}  // namespace devstats
}  // namespace SAIGE

#endif  // SAIGE_STATS_SELFCHECK_HPP
