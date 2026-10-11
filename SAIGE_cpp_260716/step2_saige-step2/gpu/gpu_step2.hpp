// gpu_step2.hpp — the ONLY header main.cpp includes for the step-2 GPU path.
//
// No CUDA types cross this boundary, so main.cpp compiles identically whether
// or not the build has CUDA. With USE_CUDA unset the implementation is
// gpu_step2_stub.cpp, whose available() always says "built without CUDA" and
// whose create() always returns nullptr — every caller then takes the CPU path
// it would have taken anyway.
//
// ---------------------------------------------------------------------------
// What this computes, and why exactly this
// ---------------------------------------------------------------------------
// SAIGE::scoreTestBatchMT (saige_mt.cpp) reduces a block of markers against the
// stacked per-trait constants. For traits whose sample list is the union's,
// the only things it touches in sample space are
//
//     Zall  = Astack^T G      (sumP x B)       every trait
//     GWqnt = Xstack^T G      (sumPqnt x B)    quantitative traits
//     GWbin = WXstack^T G     (sumPbin x B)    binary traits  (WXstack = mu2 % X)
//     GR    = G^T RES         (B x P)          every trait
//     Gsq   = colsum(G % G)   (B)              quantitative traits
//     G2Mu2 = (G % G)^T MU2bin (B x nBin)      binary traits
//
// and everything after that is O(p^2) and O(P) per marker. So two GEMMs
//
//     C1 = G^T B1,         B1 = [ Astack | WXstack | Xstack_qnt | RES ]   (N x K1)
//     C2 = (G % G)^T B2,   B2 = MU2bin                                   (N x K2)
//
// are the entire sample-space cost; C1's columns split back into the first
// four blocks above and C2 is G2Mu2. Gsq is NOT computed by a GEMM: the device
// holds the marker's 2-bit code counts, so sum_i g_i^2 = sum_c count[c]*fd[c]^2
// exactly, in double, for free. G2Mu2 cannot be had that way (it is weighted
// by mu2 per sample), so for binary traits the marker is decoded a second time
// through the SQUARED dosage table and multiplied against MU2bin. A run with no
// binary trait passes K2 = 0 and this second GEMM does not exist.
//
// ---------------------------------------------------------------------------
// The host reads bytes; the device does the per-marker preprocessing
// ---------------------------------------------------------------------------
// The caller hands over the marker's PLINK 2-bit codes verbatim (packed(), 4
// samples per byte, sample i in bits 2*(i&3) of byte i>>2) and a 1-byte
// "this slot carries a marker" flag (valid()). Everything a single-trait
// reader would derive from the column is then formed on the device
// (prepSetup, PrepSlot / PrepPair):
//
//   * code counts of the column over every distinct sample list of the run
//     (group masks: group 0 is the union, i.e. every analysis sample; the
//     others are the traits' own lists) and over every binary trait's cases
//     -- three popcounts per 64-bit word of (column AND mask);
//   * from the union's counts: the pre-imputation ALT frequency and missing
//     rate, the QC filters, the allele flip, the imputed value for missing
//     calls, the MAC-gated .clean() zeroing, the post-imputation AF / allele
//     count -- i.e. PlinkClass::fusedPreStatsFromCounts + finalizeFusedStats
//     on the four counts, expression for expression -- and so the marker's
//     4-entry code -> dosage table (deviceLut / lutOut), which the decode
//     kernel applies; sum g and sum g^2 from the same counts;
//   * per (slot, trait), the variance ratio of the pair's MAC category
//     (SAIGEClass::computeVarianceRatio), and -- when the traits have their
//     own sample lists -- the trait's own statistics from ITS counts: QC,
//     flip, imputed value, its table fd_t, AF / MAC, and the affine map
//     g_t = a g + b 1_t + d m from the union column to the trait's vector
//     (saige_mt.hpp MTBlockAdj).
//
// So the per-marker work on the host is one fread and one memcpy, whatever P
// and however many distinct sample lists there are; the device's cost is one
// popcount pass per mask per slot plus O(markers x traits) scalar work.
//
// ---------------------------------------------------------------------------
// The genotype matrix never becomes fp32 on the host
// ---------------------------------------------------------------------------
// The decode kernel applies the per-marker table, so allele flip, missing-
// genotype imputation and the MAC-gated .clean() are all already baked in and
// the GPU reproduces the CPU's dosages exactly, cell for cell, with no separate
// missing-value path. That is the whole trick: 1/16 of the PCIe traffic of
// shipping fp32 dosages, and imputation is a per-marker constant rather than a
// per-cell special case.
//
// ---------------------------------------------------------------------------
// Streaming, not resident
// ---------------------------------------------------------------------------
// Step 2 reads every marker exactly once, so there is nothing for a resident
// matrix to amortise: uploading 12.5 GB once (3.7 s measured, V100 pinned-less)
// and then computing costs strictly more than streaming the same bytes chunk by
// chunk behind the compute. This reducer therefore always streams, one
// superblock (maxSlots markers) at a time: the superblock's packed codes go up
// in one copy and stay resident until the next reduce() -- the SPA kernel
// (gpu_spa.hpp) reads flagged markers' columns straight from there -- while the
// decode+GEMM runs in slots-per-pass batches on two streams. The "does the
// packed matrix fit in device memory" question does not arise.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "gpu_precision.hpp"

namespace saige {
namespace gpu2 {

struct Reducer;   // opaque; defined in gpu_step2.cu

// True if a usable device exists. `why` (optional) gets a one-line reason when
// the answer is false, for the startup line main() prints.
bool available(int t_device, std::string* t_why);

// One-line device description ("Tesla V100-SXM2-16GB, 16384 MiB, sm_70"), or
// "" when unavailable. For the startup banner only.
std::string describe(int t_device);

struct CreateArgs {
    int device = 0;
    int N = 0;                     // samples
    // Right operand of C1 = G^T B1: N x K1 column-major DOUBLE (the caller's
    // stacked constants as they already are); narrowed here when fp64 is
    // false. Copied to the device, so it may be freed on return.
    int K1 = 0; const double* B1 = nullptr;
    // Right operand of C2 = (G % G)^T B2, N x K2. K2 = 0: no second GEMM.
    int K2 = 0; const double* B2 = nullptr;
    int maxSlots = 0;              // markers per reduce() call, at most
    // Arithmetic of the decode and the GEMMs (config key gpuPrecisionScan,
    // gpu_precision.hpp). FP64 (default): double throughout. On a V100 this
    // costs about 40% more GEMM time than FP32 and about one extra pass of
    // device bandwidth in the decode -- a few percent of a real step-2 run --
    // and it removes the precision question rather than bounding it. FP32:
    // float decode (g - rint(mean g) per marker) and SGEMMs over 4096-sample
    // chunks, the chunks and the exact shift / centring corrections summed in
    // fp64. INT8: the Ozaki-style split (int8 genotypes plus a 0/1 column
    // per mean-imputed code, B1 / B2 split into int8Slices int8 slices with
    // per-column power-of-two scaling, int32 GEMMs, recombined in fp64). Every
    // mode hands back double results (outCd / outC2d). Only the modes
    // scanSupports() accepts may be passed; create() returns nullptr for any
    // other. See PRECISION in gpu_step2.cu and gpu_scan_lowp.cuh.
    Prec precision = Prec::FP64;
    // INT8 only: slices of the trait-side operand, kInt8SlicesMin..Max.
    int int8Slices = kInt8SlicesDefault;
    // Code-count masks: nMask rows of maskWords() 64-bit words each (the
    // layout popcount_af.hpp's pcBuildMask writes), 11 in every selected
    // sample's 2-bit field. The first nGrp masks are the sample groups
    // (group 0 = every analysis sample, the union; then each distinct own
    // sample list); the rest are the binary traits' case masks, in binary
    // trait order. nGrp >= 1 is required: the union's counts drive the
    // marker's table.
    int nMask = 0; const uint64_t* masks = nullptr;
    int nGrp = 1;
    // Pinned staging sets (packed rows + valid flags), each maxSlots wide. 1
    // is the original single set. More than one lets the caller fill set s+1
    // on the host while reduce(set s) and the host work after it run (the
    // read / compute overlap of config key gpuPrefetch, S2_PIPELINE.md).
    int stagingSets = 1;
    // fp64 decode with one 16-byte store per thread, lanes contiguous, instead
    // of two 16-byte stores per thread 32 bytes apart (default false; only
    // takes effect with fp64 and N % 4 == 0). Same output; on a V100 it
    // removes the read-modify-write that half-written sectors cost on ECC'd
    // HBM2 (S2_KERNEL_ROOFLINE.md).
    bool decodeX2 = false;
    // Device sets: resident packed rows + tables on the device, and pinned
    // host result buffers (C1, C2, counts, prep records), each maxSlots wide.
    // 1 is the original single set. More than one lets the caller keep
    // superblock k's device rows (for the device SPA / Firth) and host results
    // (for the host tail) alive while reduce() of a later superblock runs into
    // another set (config key gpuOverlap, S2_OVERLAP.md). The GEMM operands,
    // the decode buffers, the launch shapes and the cuBLAS calls are the same
    // whichever set is used, so the numbers are too.
    int deviceSets = 1;
};

// True when create() accepts CreateArgs::precision = t_p. main.cpp asks
// before create() and stops the run on false ("useGPU: scan precision <mode>
// is not implemented yet"); it never falls back to another precision.
bool scanSupports(Prec t_p);

// Returns nullptr on ANY failure (no device, allocation refused, ...). The
// caller must fall back to the CPU path; nothing is printed here.
Reducer* create(const CreateArgs& t_args);
void     destroy(Reducer* t_r);

// Words per packed row, = popcount_af.hpp's pcWords(N) = (N + 31) / 32.
int maskWords(int t_N);

// ---------------------------------------------------------------------------
// Per-marker preprocessing on the device (prepSetup)
// ---------------------------------------------------------------------------
// Per trait: which group mask counts its samples, its sample count, and the
// variance-ratio table the pair's MAC category selects (the table the run's
// context picks -- sparse, null or null_noXadj -- as SAIGEClass::
// computeVarianceRatio / computeSingleVarianceRatio would), plus the sparse
// table for the sparse-GRM statistics kernel.
struct PrepTrait {
    int grp = 0;            // mask index of the trait's sample list (0 = the union)
    int nObj = 0;           // the trait's sample count
    int isBin = 0;          // binary trait
    int single = 0;         // the table has one entry: always vr[0]
    int vrOff = 0;          // consts: the active table (nCat entries)
    int vrSpOff = -1;       // consts: the sparse table (nCat entries), -1 = none
    int vrAdjOff = -1;      // consts: the null (covariate-adjusted, dense) table (nCat entries) for the
                            // fast-test recompute of an isnoadjCov trait (StatsTrait::fastRc), -1 = none
    int catMinOff = 0;      // consts: cateVarRatioMinMACVecExclude (nCat)
    int catMaxOff = 0;      // consts: cateVarRatioMaxMACVecInclude (nCat)
    int nCat = 0;
};
struct PrepArgs {
    int P = 0;
    const PrepTrait* traits = nullptr;
    const double* consts = nullptr;
    int nConsts = 0;
    // The traits have their own sample lists: fill PrepPair per (slot, trait).
    int perTrait = 0;
    // code -> pre-flip dosage for codes 0, 2, 3 (code 1 is the missing call),
    // PlinkClass::m_genoMaps_*; refFirst: AlleleOrder ref-first (the
    // pre-imputation AF is 1 - the ALT count rate, as fusedPreStatsFromCounts).
    double dmap[4] = {2.0, -1.0, 1.0, 0.0};
    int refFirst = 0;
    int imputeCase = 1;     // 1 best_guess, 2 mean, 3 minor
    double zerodCutoff = 0.0, zerodMacCutoff = 0.0;
    double missCut = 0.15, minMAF = 0.0, minMAC = 0.5, minINFO = 0.0;
    double macER = 4.0;     // MACCutoffforER: PrepSlot::hi
};
// Per slot: the union column's statistics (PlinkClass::FusedMarkerStats
// after finalizeFusedStats, on the union's counts).
struct PrepSlot {
    double AF = 0.0;        // altFreq_post
    double AC = 0.0;        // altCounts_post
    double MR = 0.0;        // missingRate
    double ss = 0.0;        // sum_c count[c] fd[c]^2
    double gs = 0.0;        // sum_c count[c] fd[c]
    unsigned char valid = 0;   // the caller staged a marker in this slot
    unsigned char qc = 0;      // passed the pre- and post-imputation filters
    unsigned char flip = 0;    // altFreq > 0.5: the dosages count the other allele
    unsigned char hi = 0;      // MAC > macER
    unsigned char pad[4] = {0, 0, 0, 0};
};
// Per (slot, trait), PrepArgs::perTrait only: the trait's own statistics on
// its sample list, and its map from the union column (MTBlockAdj: a = flip
// agreement, b = 2 when a = -1, d and q the missing-cell terms).
struct PrepPair {
    double fd[4] = {0.0, 0.0, 0.0, 0.0};   // the trait's code -> final dosage table
    double AF = 0.0, AC = 0.0;             // post-imputation, on the trait's samples
    double d = 0.0, q = 0.0;               // fd_t[miss] - (a fd_u[miss] + b);  fd_t[miss]^2 - (a fd_u[miss] + b)^2
    unsigned char qc = 0;      // passed QC on the trait's samples
    unsigned char flip = 0;
    unsigned char aff = 0;     // fd_t[c] == a fd_u[c] + b for every non-missing code
    unsigned char pad[5] = {0, 0, 0, 0, 0};
};
bool prepSetup(Reducer* t_r, const PrepArgs& t_args);
std::string prepLastError();

// Pinned staging buffers owned by the reducer. The caller writes marker data
// straight into them, so no copy stands between the decode and the H2D.
//   packed()  t_maxSlots * bytesPerSlot() bytes; slot j at offset
//             j*bytesPerSlot(). bytesPerSlot() is (N+3)/4 rounded UP to a
//             multiple of 8 so every row is a whole number of 64-bit words;
//             the caller writes (N+3)/4 bytes and leaves the rest, which the
//             reducer zeroed once and the kernels never see as a sample.
//   valid()   t_maxSlots bytes; 1 where the caller staged a marker. A slot
//             with 0 gets an all-zero table (decodes to zeros) and
//             PrepSlot::valid = 0.
//   t_set     which staging set (0 .. stagingSets()-1); the single set of a
//             reducer created with stagingSets = 1 is set 0.
unsigned char* packed(Reducer* t_r, int t_set = 0);
unsigned char* valid(Reducer* t_r, int t_set = 0);
std::size_t    bytesPerSlot(const Reducer* t_r);
int            stagingSets(const Reducer* t_r);

// Reduce slots [0, t_nSlots) of staging set t_set. Returns false on any CUDA
// error, in which case the results are undefined and the caller must fall
// back for this batch. Synchronous: on return the host result buffers are
// filled and staging set t_set is no longer read by the device.
// t_devSet: the device set (0 .. deviceSets-1) the rows are uploaded into and
// the results come back into; -1 = 0. Only one thread may call reduce() at a
// time; the accessors below may be read from another thread for a set no
// reduce() in flight is writing.
bool reduce(Reducer* t_r, int t_nSlots, int t_set = 0, int t_devSet = -1);

// Results of the last reduce(): C1 is t_maxSlots x K1 column-major, so element
// (slot, k) sits at index (std::size_t)k * ldC() + slot; C2 likewise with K2
// columns. Pinned, owned here. Exactly one of each pair of accessors is
// non-null, per isFp64(), which says whether the RESULT buffers are double:
// true in every precision mode (FP32 and INT8 recombine in fp64); the float
// accessors are kept for the interface and return nullptr.
bool          isFp64(const Reducer* t_r);
const float*  outCf(const Reducer* t_r, int t_devSet = -1);
const double* outCd(const Reducer* t_r, int t_devSet = -1);
const float*  outC2f(const Reducer* t_r, int t_devSet = -1);
const double* outC2d(const Reducer* t_r, int t_devSet = -1);
std::size_t   ldC(const Reducer* t_r);         // == maxSlots
// Code counts of the last reduce(): for slot s and mask m, the four values at
// outCounts()[((std::size_t)s * nMask + m) * 4 + c] are the number of masked
// samples whose 2-bit code is c (c indexes PLINK codes 0..3, so c = 1 is the
// missing call). Mask 0 is the union, i.e. the column's own counts.
const uint32_t* outCounts(const Reducer* t_r, int t_devSet = -1);
// The device-derived per-marker tables and statistics of the last reduce()
// into device set t_devSet (prepSetup): the code -> dosage table (4 per
// slot, as the decode used it), the per-slot record, the per-pair records
// (index slot * P + t; nullptr without PrepArgs::perTrait) and the pair's
// variance ratio (slot * P + t).
const double*   lutOut(const Reducer* t_r, int t_devSet = -1);
const PrepSlot* prepSlots(const Reducer* t_r, int t_devSet = -1);
const PrepPair* prepPairs(const Reducer* t_r, int t_devSet = -1);
const double*   prepVr(const Reducer* t_r, int t_devSet = -1);

// Cumulative device-side timings in seconds since create(), for the end-to-end
// breakdown. Measured with events on the compute stream, so they are the GPU's
// own numbers and not wall clock around the call. t_popc covers the code
// counts and the preprocessing kernels.
void timings(const Reducer* t_r, double* t_h2d, double* t_decode, double* t_gemm,
             double* t_d2h, double* t_popc = nullptr);

// Peak device bytes allocated by this reducer.
std::size_t deviceBytes(const Reducer* t_r);

// For gpu_spa.hpp: the device address of the resident packed rows and tables
// of the last reduce() (row stride bytesPerSlot()), and the stream they were
// written on. Opaque pointers; only the SPA module dereferences them.
// t_devSet -1 = the set of the last reduce().
const void* devicePacked(const Reducer* t_r, int t_devSet = -1);
const void* deviceLut(const Reducer* t_r, int t_devSet = -1);
void*       deviceStream(const Reducer* t_r);
int         deviceSets(const Reducer* t_r);

// ---------------------------------------------------------------------------
// Per-pair statistics on the device (config key gpuDeviceStats, S2_DEVSTATS.md)
// ---------------------------------------------------------------------------
// After the GEMMs, one thread per (slot, trait) forms what the host tail
// (saige_mt.cpp scoreTestBatchMT*Pre) forms from the same C1 / C2 columns.
// Binary, with the covariate projection
//     S = (g'res - S_a'Z) / tau0,   var2 = Z'XVX Z + (g^2)'mu2 - 2 W'Z
// or, for a trait with isnoadjCov (R's scoreTestFast_noadjCov),
//     S = (g'res - 2 AF sum res) / tau0,
//     var2 = tau0 ((g^2)'mu2 - 4 AF g'mu2 + 4 AF^2 sum mu2)
// Quantitative:
//     S = (g'res - S_a'Z) / tau0,   var2 = tau0 Z'XVX Z + g'g - 2 W'Z
// or with isnoadjCov
//     S = (g'res - 2 AF sum res) / tau0,   var2 = g'g - 4 AF sum g + 4 AF^2 n
// (g'g and sum g from the code counts and the table, exact). Then var1 =
// var2 * VR, stat = S^2 / var1, StdStat = |S| / sqrt(var1), the chi-square(1)
// upper tail p = erfc(sqrt(stat / 2)), and the gate bits the host tail
// decides from them: SPA (StdStat > SPA_Cutoff, binary), Firth (is_Firth_beta
// and p <= pCutoffforFirth, binary), fast test (p < the cutoff). Plain fp64
// arithmetic; the host's own tail is the same algebra in another summation
// order, and the two agree to rounding.
//
// A p that underflows to 0 in double is reported in the log domain (STATS_LOGP:
// the value is log p, by the asymptotic erfc expansion the host's
// log_chisq1_uppertail uses in that range). A degenerate pair (var1 <=
// DBL_MIN, a negative / NaN / infinite stat) gets p = 1 and STATS_DEGEN;
// STATS_STAT0 says the host's format_score_result would have reset stat to 0
// (a NaN / infinite stat), which its seBeta reads.
//
// Traits with their own sample list (PrepArgs::perTrait, StatsTrait::own):
// the GEMMs ran on the union column g, and the trait's own vector is the exact
// affine image g_t = a g + b 1_t + d m (m = the column's missing-call
// indicator; saige_mt.hpp MTBlockAdj; a, d, q from PrepPair). The kernel
// applies the map term for term, as the host's scoreTestBatchMT*PreAdj do:
//     L(g_t) = a L(g) + b L(1_t) + d L(m)          L = A', W', res', mu2'
//     Q(g_t) = Q(g) + 2ab mu2'g + b^2 sum mu2 + q mu2'm       (binary)
//     g_t'g_t = sum_c count_t[c] fd_t[c]^2                    (quantitative)
// where L(1_t) are the per-trait constants (sumA / sumW / sumR / sumM) and
// L(m) -- the stack columns summed over the column's missing calls -- come
// from a second device pass: for every slot some pair of which has d != 0
// or q != 0, one block walks the slot's packed codes, and for every missing
// call adds that sample's row of B1 (kept row-major on the device for this)
// into a K1-wide sum (miss_sums, same column indexing as C1).
//
// The fast-test recompute (StatsTrait::fastRc, StatsArgs::fastRecompute):
// with isFastTest and isnoadjCov, the scalar path rescored a pair whose
// first-pass p is below pval_cutoff_for_fastTest (binary: MAC above the ER
// cutoff) with the covariate-adjusted statistic above and the NULL
// variance-ratio table (R's mainMarkerInCPP: isnoadjCov_cur = false, and
// flagSparseGRM_cur = false when the MAC is above the sparse category -- or
// there is no sparse GRM). The kernel forms that statistic for exactly the
// pairs it flagged STATS_FAST, from the same C1 / C2 columns, with var1 from
// the null table of the pair's MAC category (PrepTrait::vrAdjOff), and its
// p and SPA / Firth gate bits (no fast bit: the recompute is not retested),
// into the statsRc* buffers; a pair without one carries STATS_HOST there.
// The host then routes it (SPA / Firth on the device) as the scalar path
// routes its recompute; a pair whose recompute context is the sparse
// variance takes the sparse kernel's result instead (below).
//
// Sparse-GRM traits (StatsSparse, statsSparseSetup): a second kernel, run
// after the cross-term kernel (gpu_sparse.hpp), forms the exact sparse
// statistic of SAIGEClass::scoreTest with z = XV g (s2_gpu_sparse.hpp):
//     S    = (g'res - z'(Y'res)) / tau0
//     var2 = z'(Y'BY)z + g'Bg - 2 z'(BY)'g,   g'Bg = (g^2)'diag B + cross terms
// with the trait's own-list map applied as above (sumXV / sumBY / tr B are
// the L(1_t) here), and var1 from the SPARSE variance-ratio table of the
// pair's MAC category. For a trait whose first pass is the sparse variance
// (isFastTest false) this is the pair's statistic and gates; for a fast-test
// trait it is the recompute the scalar path would run on the pairs whose
// dense p crosses the cutoff.
enum StatsFlag : unsigned {
    STATS_HOST  = 1u,     // the host recomputes this pair (reason in bits 4..7)
    STATS_SPA   = 2u,     // StdStat > SPA_Cutoff
    STATS_FIRTH = 4u,     // is_Firth_beta and p <= pCutoffforFirth
    STATS_FAST  = 8u,     // isFastTest and p < pval_cutoff_for_fastTest (the host still applies its context test)
    STATS_REASON_SHIFT = 4,
    STATS_LOGP  = 16u,    // p underflowed: the p value is log p (format "%.1fE%d")
    STATS_DEGEN = 32u,    // var1 <= DBL_MIN or a bad stat: p = 1
    STATS_STAT0 = 64u     // with STATS_DEGEN: stat was NaN / infinite (format_score_result sets stat = 0)
};
enum StatsReason : unsigned {
    STATS_R_NONE = 0,
    STATS_R_TAIL = 1,     // (no longer used: the log-domain p is formed on the device)
    STATS_R_DEGEN = 2,    // (no longer used: the degenerate pair is formed on the device)
    STATS_R_OFF = 3,      // the trait is not scored by this kernel
    STATS_R_COUNT = 4
};
struct StatsTrait {
    int trait = 0;          // internal trait index (PrepTrait / PrepPair / prepVr index)
    int isBin = 1;          // binary (1) or quantitative (0)
    int p = 0;              // covariates incl. the intercept (any count)
    int rowZ = 0;           // first C1 column (row of C1^T) of this trait's A'g block
    int rowW = 0;           // first C1 column of (mu2 % X)'g (binary) / X'g (quantitative)
    int rowGR = 0;          // C1 column of g'res
    int rowGM = -1;         // C1 column of g'mu2 (binary isnoadjCov / own traits), -1 when absent
    int colG2 = -1;         // C2 column of (g^2)'mu2 (binary), -1 for quantitative
    int enabled = 1;        // 0: every pair of this trait is flagged STATS_HOST / STATS_R_OFF
    int isFirth = 0;        // is_Firth_beta
    int isFast = 0;         // isFastTest
    int noadj = 0;          // isnoadjCov: the centred score, no covariate block
    int fastRc = 0;         // isnoadjCov and isFastTest: the covariate-adjusted recompute of the pairs
                            // flagged STATS_FAST (needs rowZ / rowW / xvxOff / saOff and PrepTrait::vrAdjOff)
    double tau0 = 1.0, spaCut = 2.0, firthCut = 0.0, fastCut = 0.0;
    double sumR = 0.0;      // sum res (isnoadjCov, own)
    double sumM = 0.0;      // sum mu2 (binary) / n_t (quantitative) (isnoadjCov, own)
    // Offsets into StatsArgs::consts: XVX (row-major p x p) and S_a (p).
    int xvxOff = 0, saOff = 0;
    // The trait's sample list is not the union's: apply the pair's affine map
    // (PrepPair) and the missing-cell sums. Needs rowGM >= 0 (binary) and
    // p <= 32. sumAOff / sumWOff: offsets into consts of the trait's A / W
    // columns summed over its samples (p each).
    int own = 0;
    int sumAOff = 0, sumWOff = 0;
};
struct StatsArgs {
    int nTraits = 0;                       // = P (every trait, internal order)
    const StatsTrait* traits = nullptr;
    const double* consts = nullptr;        // the traits' XVX / S_a / sumA / sumW, nConsts doubles
    int nConsts = 0;
    // Some trait has own = 1: allocate the row-major copy of B1 (N x K1
    // doubles on the device) and the missing-cell sums, and run miss_sums
    // before the stats kernel. Needs the FP64 scan and PrepArgs::perTrait.
    int ownSets = 0;
    // Some trait has fastRc = 1: allocate the recompute buffers (statsRc*).
    int fastRecompute = 0;
};
// Allocate the buffers and upload the constants; from then on every reduce()
// also runs the stats kernel. false (with statsLastError()) on failure.
// Requires prepSetup.
bool statsSetup(Reducer* t_r, const StatsArgs& t_args);
// Stop running the stats kernel in reduce().
void statsDisable(Reducer* t_r);
std::string statsLastError();
// Results of the last reduce() into device set t_devSet, pair-major: index
// slot * nTraits + t. nullptr until statsSetup.
const double*        statsS(const Reducer* t_r, int t_devSet = -1);
const double*        statsVar2(const Reducer* t_r, int t_devSet = -1);
const double*        statsP(const Reducer* t_r, int t_devSet = -1);
const unsigned char* statsFlags(const Reducer* t_r, int t_devSet = -1);
// The fast-test recompute of the last reduce() (StatsArgs::fastRecompute),
// same indexing: S, var2, the null-table variance ratio, p and the flags
// (STATS_HOST where the kernel formed no recompute for the pair). nullptr
// without fastRecompute.
const double*        statsRcS(const Reducer* t_r, int t_devSet = -1);
const double*        statsRcVar2(const Reducer* t_r, int t_devSet = -1);
const double*        statsRcVr(const Reducer* t_r, int t_devSet = -1);
const double*        statsRcP(const Reducer* t_r, int t_devSet = -1);
const unsigned char* statsRcFlags(const Reducer* t_r, int t_devSet = -1);
// Cumulative kernel + D2H seconds of the stats stage, for the breakdown line.
double statsSeconds(const Reducer* t_r);

// Sparse-GRM statistics (see above). One record per sparse trait, in sparse
// trait order s (the cross-term kernel's trait order).
struct StatsSparse {
    int trait = 0;          // internal trait index
    int p = 0;
    int rowZ = 0;           // first C1 column of XV g
    int rowGW = 0;          // first C1 column of (BY)'g
    int rowGR = 0;          // C1 column of g'res
    int rowGB = -1;         // C1 column of g'diag B (own traits), -1 when absent
    int colQ = 0;           // C2 column of (g^2)'diag B
    int isBin = 1;
    int isFirth = 0;
    double tau0 = 1.0, spaCut = 2.0, firthCut = 0.0, trB = 0.0;
    double sumR = 0.0;
    int ybyOff = 0, yresOff = 0;     // consts: Y'BY (row-major p x p), Y'res (p)
    int own = 0;
    int sumXVOff = 0, sumBYOff = 0;  // consts: XV 1_t, (BY)'1_t (p each)
};
struct StatsSparseArgs {
    int nSparse = 0;
    const StatsSparse* traits = nullptr;
    const double* consts = nullptr;
    int nConsts = 0;
};
bool statsSparseSetup(Reducer* t_r, const StatsSparseArgs& t_args);
// Run the sparse kernel on slots [0, t_nSlots) of device set t_devSet, with
// t_cross the device address of the cross terms (slot * nSparse + s; the
// cross-term module's spqDeviceOut) -- after the reduce() into that set and
// the cross-term call, before the next reduce(). Synchronous. The device-side
// tables of the own-list traits for the cross-term kernel (slot x nSparse x 4,
// the trait's own fd or the union's, zero where the pair failed QC) are
// formed by sparseTables() from the prep records.
bool statsSparseRun(Reducer* t_r, int t_nSlots, int t_devSet, const void* t_cross);
const void* sparseTables(Reducer* t_r, int t_nSlots, int t_devSet);   // device address, or nullptr on failure
// Results, index slot * nSparse + s: S, var2, p, flags as the dense kernel's,
// and the sparse-table variance ratio the pair's MAC selects.
const double*        statsSpS(const Reducer* t_r, int t_devSet = -1);
const double*        statsSpVar2(const Reducer* t_r, int t_devSet = -1);
const double*        statsSpP(const Reducer* t_r, int t_devSet = -1);
const unsigned char* statsSpFlags(const Reducer* t_r, int t_devSet = -1);
const double*        statsSpVr(const Reducer* t_r, int t_devSet = -1);

// cudaSetDevice for a host thread other than the one that called create()
// (the gpuOverlap workers). false without CUDA or on error.
bool bindDevice(int t_device);
// Ask the driver to block, not spin, a host thread waiting on the device
// (cudaDeviceScheduleBlockingSync), so the gpuOverlap workers' waits leave
// their cores to the host tail. Returns the CUDA status as text ("" = ok).
std::string setBlockingSync(int t_device);

}  // namespace gpu2
}  // namespace saige
