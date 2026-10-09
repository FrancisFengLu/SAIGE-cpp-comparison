// gpu_step2.cu — CUDA implementation of the step-2 sample-space reduction.
// Contract, and why it is shaped this way: gpu_step2.hpp.
//
// Provenance: the decode kernel's structure (one block per marker, one byte =
// four samples per thread, one vector store, coalesced at both ends) is taken
// from the step-1 packed kernels in
// step1_saige-null/tools/gpu_matvec/gemv2bit.cu, and from the standalone
// prototype optimization/torchgwas2/scripts/step2core.cu. What is new here is
// the PER-MARKER code->dosage table: step 1 can hard-code 2-popcount(field)
// because its packer guarantees no missing code ever reaches the device,
// whereas step 2 must honour the marker's own allele flip, its imputed value
// for missing calls, and the MAC-gated .clean() zeroing. All three are
// functions of the 2-bit code alone, so they collapse into four floats per
// marker and the kernel stays branch-free -- no separate missing-value path,
// and imputation costs nothing. The binary-trait second GEMM reuses the same
// kernel with the squared table.
//
// Determinism: for a given (N, K1, K2, maxSlots, slotsPerPass, precision) every
// launch shape is fixed, so two runs of the same configuration produce
// bit-identical C. cuBLAS is left in PEDANTIC math mode so no TF32 or
// reduced-precision reduction can be substituted on a device newer than the
// V100 this was written on.
//
// PRECISION
// ---------
// FP64: decode and GEMM run in double. Measured on this V100 (14.28 TFLOP/s
// fp32, 7.19 fp64): fp64 costs
// about 40% more GEMM time and one extra pass of device bandwidth in the
// decode. In a real step-2 run the GEMM is single-digit seconds against
// a hundred-plus seconds of reading and decoding the .bed, so that is a few
// percent of wall clock and it buys agreement with the CPU path near 1e-15
// instead of near 1e-6. fp64 is therefore the default. fp32 stays available
// for devices where fp64 runs at 1/32 or 1/64 rate (L4, A10, consumer parts)
// and for measuring the fp32 error itself.
//
// The mode comes from CreateArgs::precision (config key gpuPrecisionScan,
// gpu_precision.hpp). FP64 is the code above, untouched. FP32 and INT8 are in
// gpu_scan_lowp.cuh: FP32 decodes g - rint(mean g) per marker in float, runs
// one SGEMM per kChunkF32 samples and sums the chunks (and the exact shift /
// centring corrections) in fp64; INT8 is the Ozaki-style split, int8 x int8
// -> int32 GEMMs per slice of B, recombined in fp64. Both hand back double
// results (outD), so the sparse-GRM variance (outCd) works in every mode.
//
// What the old fp32 (one SGEMM over all N samples, float results) got wrong,
// measured on bingpu_test bt full (N = 50,000, P = 8) against fp64 at full
// precision (outputFormat sgs):
//   * its 1.9e-2 max relative p.value difference is not a large score error:
//     it sits on MAC-5 markers of the 50%-prevalence traits whose score is at
//     the edge of its support (SPA p 9e-18 against a normal-approximation p
//     of 2e-2), where the saddle point runs off and p moves 2% for an 8e-8
//     relative move of the score; printed Tstat / var agree to 6 digits there.
//   * the real errors: Tstat 9e-6 of its standard deviation (one fp32 dot
//     product of 50,000 nearly cancelling terms), var 2.7e-6 and SE 2.5e-4
//     relative (the mu2 and W X columns are all-positive, so their sums carry
//     a large common part).
// The fp32 mode now (gpu_scan_lowp.cuh): 4096-sample chunks summed in fp64
// (Tstat 9e-6 -> 1.1e-6 of its sd), integer shift of g and centring of the
// columns with |mean| > sd, both undone exactly in fp64 (var -> 3e-8, SE ->
// 8e-7). What is left of the Tstat error is fp32 rounding of a sum whose size
// is far above its spread (A columns with a nonzero mean), magnified by the
// cancellation in S = g'res - S_a'(A'g); neither smaller chunks nor a hi/lo
// split of B move it, and centring those columns would cost a MAC-1 marker
// the relative precision of its one entry. INT8 is the mode for fp64-level
// agreement.
//
// The dosage TABLE is double in every mode (FP32 narrows g - r inside the
// kernel). Three of its four entries are exact small integers, but the fourth is
// the imputed mean 2*altFreq -- narrowing that on the host would put a 6e-8
// relative error into every missing cell before the reduction starts, which is
// a difference from the CPU's INPUT rather than from its arithmetic. It costs
// 32 bytes per marker to carry.

#include "gpu_step2.hpp"
#include "gpu_scan_lowp.cuh"

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace saige {
namespace gpu2 {

namespace {

// Never abort: a GPU failure must leave the caller free to run on the CPU.
#define CKR(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    lastErr = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)
#define CBR(x) do { cublasStatus_t s_ = (x); if (s_ != CUBLAS_STATUS_SUCCESS) { \
    lastErr = std::string(#x) + ": cuBLAS status " + std::to_string((int)s_); return false; } } while (0)

thread_local std::string lastErr;

// The marker's four dosages, held in registers.
template <typename T> struct Lut4 { T v0, v1, v2, v3; };

template <typename T>
__device__ __forceinline__ Lut4<T> loadLut(const double* __restrict__ lut, int m)
{
    const double2* d = reinterpret_cast<const double2*>(lut) + 2 * m;
    const double2 a = d[0], b = d[1];
    return Lut4<T>{(T)a.x, (T)a.y, (T)b.x, (T)b.y};
}

// code -> dosage. Two predicated selects; no local array (which would land in
// local memory) and no shared-memory lookup (which would bank-conflict across
// codes).
template <typename T>
__device__ __forceinline__ T pick(const Lut4<T>& L, unsigned c)
{
    const T a = (c & 1u) ? L.v1 : L.v0;
    const T b = (c & 1u) ? L.v3 : L.v2;
    return (c & 2u) ? b : a;
}

// double -> two 16-byte stores, aligned because N % 4 == 0 on this path, so
// every marker's column starts on a 16-byte boundary.
__device__ __forceinline__ void storeQuad(double* p, double a, double b, double c, double d)
{
    reinterpret_cast<double2*>(p)[0] = make_double2(a, b);
    reinterpret_cast<double2*>(p)[1] = make_double2(c, d);
}

// One block per marker. Thread b takes byte b (samples 4b..4b+3), so both the
// byte load and the value stores are fully coalesced.
template <typename T>
__global__ void __launch_bounds__(256)
decode_lut_x4(const uint8_t* __restrict__ packed, std::size_t bpv, int N,
              const double* __restrict__ lut, T* __restrict__ dG)
{
    const int m = blockIdx.x;
    const uint8_t* __restrict__ row = packed + (std::size_t)m * bpv;
    const Lut4<T> L = loadLut<T>(lut, m);
    T* __restrict__ o = dG + (std::size_t)m * N;
    const int nq = N >> 2;
    for (int b = threadIdx.x; b < nq; b += blockDim.x) {
        const unsigned p = row[b];
        storeQuad(o + 4 * b,
                  pick(L, p & 3u),        pick(L, (p >> 2) & 3u),
                  pick(L, (p >> 4) & 3u), pick(L, (p >> 6) & 3u));
    }
}

// fp64 variant with lane-contiguous stores (CreateArgs::decodeX2). The x4
// kernel's thread writes its 32 bytes as two 16-byte stores, so across a warp
// each store instruction touches 32 sectors by half; the ECC'd HBM2 of a V100
// turns a half-written sector into a read-modify-write, and ncu shows the
// kernel reading 0.55 byte and writing 1.7 bytes from DRAM per payload byte
// (S2_KERNEL_ROOFLINE.md). Here thread i takes samples 2i, 2i+1 -- one
// 16-byte store, lanes contiguous -- so a warp's store covers 16 whole sectors.
// Same table, same values, same positions: the output is identical.
__global__ void __launch_bounds__(256)
decode_lut_x2(const uint8_t* __restrict__ packed, std::size_t bpv, int N,
              const double* __restrict__ lut, double* __restrict__ dG)
{
    const int m = blockIdx.x;
    const uint8_t* __restrict__ row = packed + (std::size_t)m * bpv;
    const Lut4<double> L = loadLut<double>(lut, m);
    double* __restrict__ o = dG + (std::size_t)m * N;
    const int nh = N >> 1;
    for (int h = threadIdx.x; h < nh; h += blockDim.x) {
        const unsigned p = (unsigned)row[h >> 1] >> ((h & 1) * 4);
        *reinterpret_cast<double2*>(o + 2 * h) = make_double2(pick(L, p & 3u), pick(L, (p >> 2) & 3u));
    }
}

// General N. One sample per thread-iteration; used only when N % 4 != 0.
template <typename T>
__global__ void __launch_bounds__(256)
decode_lut_any(const uint8_t* __restrict__ packed, std::size_t bpv, int N,
               const double* __restrict__ lut, T* __restrict__ dG)
{
    const int m = blockIdx.x;
    const uint8_t* __restrict__ row = packed + (std::size_t)m * bpv;
    const Lut4<T> L = loadLut<T>(lut, m);
    T* __restrict__ o = dG + (std::size_t)m * N;
    for (int i = threadIdx.x; i < N; i += blockDim.x)
        o[i] = pick(L, (unsigned)((row[i >> 2] >> ((i & 3) * 2)) & 3u));
}

// Code counts: one block per slot; for each mask, every thread popcounts its
// share of the words of (column AND mask) and the block reduces the three
// counts. The column (one row, 12.5 KB at N = 50,000) stays in L1 across the
// masks; the masks (nMask * words * 8 B) are shared by every block and sit in
// L2. Same arithmetic as popcount_af.hpp's pcCountCodes, including the
// derivation of the fourth count from the mask population.
__device__ __forceinline__ unsigned warpSumU(unsigned v)
{
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

__global__ void __launch_bounds__(256)
count_codes(const uint8_t* __restrict__ packed, std::size_t bpv, int words,
            const uint64_t* __restrict__ masks, const uint32_t* __restrict__ maskPop,
            int nMask, uint32_t* __restrict__ out)
{
    __shared__ unsigned sh[3][8];
    const int s = blockIdx.x;
    const uint64_t* __restrict__ col = reinterpret_cast<const uint64_t*>(packed + (std::size_t)s * bpv);
    const uint64_t M55 = 0x5555555555555555ULL;
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int m = 0; m < nMask; ++m) {
        const uint64_t* __restrict__ mk = masks + (std::size_t)m * words;
        unsigned n11 = 0, n10 = 0, n01 = 0;
        for (int w = threadIdx.x; w < words; w += blockDim.x) {
            const uint64_t x  = col[w] & mk[w];
            const uint64_t lo = x & M55;
            const uint64_t hi = (x >> 1) & M55;
            n11 += (unsigned)__popcll(hi & lo);
            n10 += (unsigned)__popcll(hi & ~lo);
            n01 += (unsigned)__popcll(lo & ~hi);
        }
        n11 = warpSumU(n11); n10 = warpSumU(n10); n01 = warpSumU(n01);
        if (lane == 0) { sh[0][wid] = n11; sh[1][wid] = n10; sh[2][wid] = n01; }
        __syncthreads();
        if (threadIdx.x == 0) {
            unsigned a = 0, b = 0, c = 0;
            for (int k = 0; k < (int)(blockDim.x >> 5); ++k) { a += sh[0][k]; b += sh[1][k]; c += sh[2][k]; }
            uint32_t* o = out + ((std::size_t)s * nMask + m) * 4;
            // PLINK codes: 00 HOM_ALT, 01 MISSING, 10 HET, 11 HOM_REF
            o[3] = a; o[2] = b; o[1] = c; o[0] = maskPop[m] - a - b - c;
        }
        __syncthreads();
    }
}

// ---------------------------------------------------------------------------
// Per-pair statistics (gpu_step2.hpp "Per-pair statistics on the device").
// Every operation below goes through the IEEE intrinsics or fma(), so the
// result does not depend on nvcc's contraction setting; the orders are the
// ones stats_selfcheck.hpp enumerates and main.cpp identified on the host.
// ---------------------------------------------------------------------------
__device__ __forceinline__ double sdAdd(double a, double b) { return __dadd_rn(a, b); }
__device__ __forceinline__ double sdSub(double a, double b) { return __dsub_rn(a, b); }
__device__ __forceinline__ double sdMul(double a, double b) { return __dmul_rn(a, b); }
__device__ __forceinline__ double sdDiv(double a, double b) { return __ddiv_rn(a, b); }

// sum_k a[k] * b[k], k < p, in order `pat` (stats_selfcheck.hpp: dotPattern)
__device__ __forceinline__ double sdDot(int pat, const double* a, const double* b, int p)
{
    if (p == 1) return (pat == 1) ? fma(a[0], b[0], 0.0) : sdMul(a[0], b[0]);
    double acc;
    int k = 2;
    switch (pat) {
    case 1:   // accumulator from +0, one fma per term (a BLAS micro-kernel)
        acc = fma(a[0], b[0], 0.0);
        acc = fma(a[1], b[1], acc);
        for (; k < p; ++k) acc = fma(a[k], b[k], acc);
        break;
    case 2:   // first product fused into the second: fma(a0, b0, a1 b1), then a chain of fma
        acc = fma(a[0], b[0], sdMul(a[1], b[1]));
        for (; k < p; ++k) acc = fma(a[k], b[k], acc);
        break;
    case 3:   // (a0 b0 + a1 b1) unfused, then a chain of fma
        acc = sdAdd(sdMul(a[0], b[0]), sdMul(a[1], b[1]));
        for (; k < p; ++k) acc = fma(a[k], b[k], acc);
        break;
    case 4:   // fma(a0, b0, a1 b1), then unfused adds
        acc = fma(a[0], b[0], sdMul(a[1], b[1]));
        for (; k < p; ++k) acc = sdAdd(acc, sdMul(a[k], b[k]));
        break;
    default:  // 0: no fma anywhere, left to right
        acc = sdAdd(sdMul(a[0], b[0]), sdMul(a[1], b[1]));
        for (; k < p; ++k) acc = sdAdd(acc, sdMul(a[k], b[k]));
        break;
    }
    return acc;
}

// armadillo's column sum of the Hadamard product x % y over p rows: even rows
// into val1, odd rows into val2 (both from +0), then val1 + val2. pat 1: the
// term is fused into the accumulator (fma), pat 0: rounded first.
__device__ __forceinline__ double sdHadSum(int pat, const double* x, const double* y, int p)
{
    double v1 = 0.0, v2 = 0.0;
    int i = 0;
    if (pat == 1) {
        for (; i + 1 < p; i += 2) { v1 = fma(x[i], y[i], v1); v2 = fma(x[i + 1], y[i + 1], v2); }
        if (i < p) v1 = fma(x[i], y[i], v1);
    } else {
        for (; i + 1 < p; i += 2) { v1 = sdAdd(v1, sdMul(x[i], y[i])); v2 = sdAdd(v2, sdMul(x[i + 1], y[i + 1])); }
        if (i < p) v1 = sdAdd(v1, sdMul(x[i], y[i]));
    }
    return sdAdd(v1, v2);
}

// chi-square(1) upper tail as the host's score_vec.hpp forms it
__device__ __forceinline__ double sdTailP(double stat) { return erfc(sqrt(sdMul(stat, 0.5))); }

__constant__ double c_p10[20] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
                                 1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19};

__global__ void __launch_bounds__(256)
pair_stats(const double* __restrict__ C1, const double* __restrict__ C2, std::size_t ld, int nSlots,
           const StatsTrait* __restrict__ T, int nT, const double* __restrict__ vr,
           double statCutoff, double tol, int patXZ, int patSaz, int patZxz, int patGwz,
           double* __restrict__ oS, double* __restrict__ oV, double* __restrict__ oP,
           unsigned char* __restrict__ oF)
{
    const long long idx = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= (long long)nSlots * nT) return;
    const int b = (int)(idx / nSlots);
    const int slot = (int)(idx - (long long)b * nSlots);
    const std::size_t o = (std::size_t)slot * nT + b;
    const StatsTrait& t = T[b];
    if (!t.enabled) {
        oS[o] = 0.0; oV[o] = 0.0; oP[o] = 1.0;
        oF[o] = (unsigned char)(STATS_HOST | (STATS_R_OFF << STATS_REASON_SHIFT));
        return;
    }
    const int p = t.p;
    double z[STATS_PMAX], w[STATS_PMAX], xz[STATS_PMAX];
    for (int i = 0; i < p; ++i) {
        z[i] = C1[(std::size_t)(t.rowZ + i) * ld + slot];
        w[i] = C1[(std::size_t)(t.rowW + i) * ld + slot];
    }
    const double gr = C1[(std::size_t)t.rowGR * ld + slot];
    const double g2 = C2[(std::size_t)t.colG2 * ld + slot];
    for (int i = 0; i < p; ++i) xz[i] = sdDot(patXZ, &t.XVX[i * STATS_PMAX], z, p);
    const double zxz = sdHadSum(patZxz, z, xz, p);
    const double saz = sdDot(patSaz, t.Sa, z, p);
    const double gwz = sdHadSum(patGwz, w, z, p);
    // S = (GR - saz) / tau0;  var2 = zxz + G2Mu2 - 2.0 * gwz  (the host's expressions)
    const double S    = sdDiv(sdSub(gr, saz), t.tau0);
    const double var2 = sdSub(sdAdd(zxz, g2), sdMul(2.0, gwz));
    oS[o] = S; oV[o] = var2;

    const double var1 = sdMul(var2, vr[o]);
    const double stat = sdDiv(sdMul(S, S), var1);
    unsigned f = 0;
    double pv = 1.0;
    const double tiny = 2.2250738585072014e-308;   // DBL_MIN, format_score_result's test
    auto host = [&](unsigned reason) { if (!(f & STATS_HOST)) f |= STATS_HOST | (reason << STATS_REASON_SHIFT); };
    if (!(var1 > tiny) || !(stat >= 0.0) || !isfinite(stat)) {
        host(STATS_R_DEGEN);
    } else if (!(stat < statCutoff)) {
        host(STATS_R_TAIL);
    } else {
        pv = sdTailP(stat);
        const double plo = sdMul(pv, 1.0 - tol), phi = sdMul(pv, 1.0 + tol);
        const double sd = sdDiv(fabs(S), sqrt(var1));            // StdStat
        if (sd > t.spaCut) f |= STATS_SPA;
        if (t.isFirth) {
            if (phi >= t.firthCut && plo <= t.firthCut) host(STATS_R_FIRTH);
            else if (pv <= t.firthCut) f |= STATS_FIRTH;
        }
        // The 7 significant digits "%.6E" prints: m = p * 10^(6-e) in
        // [1e6, 1e7), rounded to an integer. Flag when p*(1 -+ tol) could
        // round differently (or cross a power of ten).
        int e = (int)floor(log10(pv));
        if (e > 0) e = 0;
        if (e < -13) e = -13;
        double m = sdMul(pv, c_p10[6 - e]);
        if (m < 1e6 && e > -13) { --e; m = sdMul(pv, c_p10[6 - e]); }
        else if (m >= 1e7 && e < 0) { ++e; m = sdMul(pv, c_p10[6 - e]); }
        const double mlo = m * (1.0 - tol) - 2e-8, mhi = m * (1.0 + tol) + 2e-8;
        if (mlo < 1e6 || mhi >= 1e7 || floor(mlo + 0.5) != floor(mhi + 0.5)) {
            host(STATS_R_PRINT);
        } else if (t.isFast) {
            // std::stod of the printed string against pval_cutoff_for_fastTest
            const double v7 = floor(m + 0.5) / c_p10[6 - e];
            if (fabs(v7 - t.fastCut) <= t.fastCut * 4e-15) host(STATS_R_FASTC);
            else if (v7 < t.fastCut) f |= STATS_FAST;
        }
    }
    oP[o] = pv;
    oF[o] = (unsigned char)f;
}

__global__ void erfc_tail(const double* __restrict__ stat, int n, double* __restrict__ p)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) p[i] = sdTailP(stat[i]);
}

}  // namespace

// ---------------------------------------------------------------------------

struct Reducer {
    int  N = 0, K1 = 0, K2 = 0, maxSlots = 0, nMask = 0, words = 0;
    Prec prec = Prec::FP64;        // CreateArgs::precision
    int  int8Slices = 0;           // INT8 only
    bool fp64 = true;              // decode + GEMM in double (prec == FP64)
    bool outD = true;              // result buffers hC1 / hC2 are double (FP64, INT8); isFp64()
    bool decodeX2 = false;
    std::size_t bpv = 0;           // padded row stride, bytes (multiple of 8)
    std::size_t esz = 0;           // sizeof(T)
    int slotsPerPass = 0;          // markers per device pass; two are resident

    // pinned host staging: nSets sets of (packed, lut); hLut2 is one buffer,
    // filled from the set being reduced
    int            nSets   = 1;
    std::vector<unsigned char*> hPacked;
    std::vector<double*>        hLut;
    double*        hLut2   = nullptr;  // squared table, filled in reduce()
    // pinned host results, one per device set
    int            nDev    = 1;
    int            lastDev = 0;       // the device set of the last reduce()
    std::vector<void*>     hC1;
    std::vector<void*>     hC2;
    std::vector<uint32_t*> hCnt;

    // device
    void*          dB1  = nullptr;      // N x K1
    void*          dB2  = nullptr;      // N x K2
    void*          dC1  = nullptr;      // maxSlots x K1
    void*          dC2  = nullptr;      // maxSlots x K2
    std::vector<unsigned char*> dPk;    // per device set: maxSlots x bpv, resident for the superblock
    std::vector<double*>        dLut;   // per device set: maxSlots x 4
    double*        dLut2 = nullptr;     // maxSlots x 4
    uint64_t*      dMask = nullptr;     // nMask x words
    uint32_t*      dMaskPop = nullptr;  // nMask
    uint32_t*      dCnt = nullptr;      // maxSlots x nMask x 4
    void*          dG[2]   = {nullptr, nullptr};

    cudaStream_t   st[2] = {nullptr, nullptr};
    cudaEvent_t    evUp = nullptr;      // packed + tables landed
    cudaEvent_t    ev[2][4] = {{nullptr, nullptr, nullptr, nullptr},
                               {nullptr, nullptr, nullptr, nullptr}};
    bool           evLive[2] = {false, false};
    cublasHandle_t cub = nullptr;

    // ---- FP32 / INT8 (gpu_scan_lowp.cuh); unused in FP64 ----
    std::size_t cesz = 8;               // bytes per result element (dC*, hC*)
    double*  dCs1 = nullptr;            // FP32: colsum(B1), fp64, K1
    double*  dCs2 = nullptr;            // FP32: colsum(B2), fp64, K2
    double*  dShift1 = nullptr;         // FP32: per-slot integer shift of g, maxSlots
    double*  dShift2 = nullptr;         // FP32: per-slot integer shift of g^2, maxSlots
    void*    dPart[2] = {nullptr, nullptr};  // per stream: fp32 chunk / int32 slice product
    int      Np = 0;                    // INT8: padded sample count (multiple of 16)
    int      capX = 0;                  // INT8: indicator columns per pass, at most
    int8_t*  dBs1 = nullptr;            // INT8: slices of B1, int8Slices x Np x K1
    int8_t*  dBs2 = nullptr;            // INT8: slices of B2
    int*     dE1 = nullptr;             // INT8: per-column exponent of B1, K1
    int*     dE2 = nullptr;             // INT8: per-column exponent of B2, K2
    // INT8 indicator columns of one reduce(): per slot (global) its first
    // indicator (pass-relative) and count; per indicator (global, in slot
    // order) its slot (pass-relative), code and the two table deltas.
    int*     hXf = nullptr; int* hXn = nullptr; int* hXs = nullptr; int* hXc = nullptr;
    double*  hXd1 = nullptr; double* hXd2 = nullptr;
    int*     dXf = nullptr; int* dXn = nullptr; int* dXs = nullptr; int* dXc = nullptr;
    double*  dXd1 = nullptr; double* dXd2 = nullptr;

    double tH2D = 0, tDec = 0, tGemm = 0, tD2H = 0, tPopc = 0;
    std::size_t devBytes = 0;

    // ---- per-pair statistics (statsSetup; gpu_step2.hpp "Per-pair statistics") ----
    bool       stats = false;
    int        nStat = 0;                 // binary traits (= nMask)
    StatsTrait* dStatT = nullptr;         // device: nStat
    double     statCutoff = 0.0, pRelTol = 2e-14;
    int        patXZ = 0, patSaz = 0, patZxz = 0, patGwz = 0;
    std::vector<double*> hVr;             // pinned, per staging set: maxSlots x nStat
    double*    dVr = nullptr;             // device: maxSlots x nStat
    double*    dStS = nullptr;            // device: maxSlots x nStat, pair-major
    double*    dStV = nullptr;
    double*    dStP = nullptr;
    unsigned char* dStF = nullptr;
    std::vector<double*>        hStS, hStV, hStP;   // pinned, per device set
    std::vector<unsigned char*> hStF;
    double     tStats = 0;
};

namespace {

// Fold buffer `b`'s recorded stage times into the totals. Only called where the
// stream has already been synchronised, so it never adds a stall of its own.
// ev[b][0..3]: start, after decode(G), after GEMM1, after decode(G^2)+GEMM2.
void harvest(Reducer* r, int b)
{
    if (!r->evLive[b]) return;
    float ms = 0;
    if (cudaEventElapsedTime(&ms, r->ev[b][0], r->ev[b][1]) == cudaSuccess) r->tDec  += ms * 1e-3;
    if (cudaEventElapsedTime(&ms, r->ev[b][1], r->ev[b][2]) == cudaSuccess) r->tGemm += ms * 1e-3;
    if (r->K2 > 0 && cudaEventElapsedTime(&ms, r->ev[b][2], r->ev[b][3]) == cudaSuccess) r->tGemm += ms * 1e-3;
    r->evLive[b] = false;
}

}  // namespace

bool available(int t_device, std::string* t_why)
{
    int n = 0;
    cudaError_t e = cudaGetDeviceCount(&n);
    if (e != cudaSuccess) {
        if (t_why) *t_why = std::string("cudaGetDeviceCount: ") + cudaGetErrorString(e);
        return false;
    }
    if (n <= 0) { if (t_why) *t_why = "no CUDA device"; return false; }
    if (t_device < 0 || t_device >= n) {
        if (t_why) *t_why = "gpuDevice " + std::to_string(t_device) + " out of range (" +
                            std::to_string(n) + " device(s) present)";
        return false;
    }
    cudaDeviceProp p;
    if (cudaGetDeviceProperties(&p, t_device) != cudaSuccess) {
        if (t_why) *t_why = "cudaGetDeviceProperties failed";
        return false;
    }
    return true;
}

std::string describe(int t_device)
{
    cudaDeviceProp p;
    if (cudaGetDeviceProperties(&p, t_device) != cudaSuccess) return "";
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s, %.0f MiB, sm_%d%d", p.name,
                  (double)p.totalGlobalMem / (1024.0 * 1024.0), p.major, p.minor);
    return std::string(buf);
}

int maskWords(int t_N) { return (t_N + 31) / 32; }

bool scanSupports(Prec t_p)
{
    switch (t_p) {
        case Prec::FP64: return true;
        case Prec::FP32: return true;
        // The sparse-GRM cross terms (gpu_sparse.cu, spqSupports) are part of
        // this stage too and have their own switch.
        case Prec::INT8: return true;
    }
    return false;
}

Reducer* create(const CreateArgs& a)
{
    if (a.N <= 0 || a.K1 <= 0 || a.maxSlots <= 0 || a.B1 == nullptr) return nullptr;
    if (a.K2 > 0 && a.B2 == nullptr) return nullptr;
    if (a.nMask > 0 && a.masks == nullptr) return nullptr;
    // ---- precision dispatch (the one place the mode is decided) ----
    if (!scanSupports(a.precision)) {
        lastErr = std::string("scan precision ") + precName(a.precision) + " is not implemented yet";
        return nullptr;
    }
    if (a.precision == Prec::INT8 &&
        (a.int8Slices < kInt8SlicesMin || a.int8Slices > kInt8SlicesMax)) {
        lastErr = "int8Slices out of range";
        return nullptr;
    }
    if (cudaSetDevice(a.device) != cudaSuccess) return nullptr;

    Reducer* r = new Reducer();
    r->prec = a.precision;
    r->int8Slices = (a.precision == Prec::INT8) ? a.int8Slices : 0;
    const bool fp64 = (a.precision == Prec::FP64);
    const bool i8   = (a.precision == Prec::INT8);
    // int8: every partial sum is bounded by 4 * 64 * N (gpu_scan_lowp.cuh)
    if (i8 && (long long)a.N * lowp::kI8MaxA * 64 >= (1LL << 31)) {
        lastErr = "int8 scan: N too large for the int32 accumulation bound";
        return nullptr;
    }
    r->N = a.N; r->K1 = a.K1; r->K2 = a.K2; r->maxSlots = a.maxSlots; r->fp64 = fp64;
    r->outD = true;                     // every mode hands back double results
    r->cesz = sizeof(double);
    r->decodeX2 = a.decodeX2 && fp64 && ((a.N & 3) == 0);
    r->nMask = a.nMask;
    r->words = maskWords(a.N);
    r->bpv = (std::size_t)r->words * 8;           // >= (N+3)/4, whole 64-bit words
    r->esz = fp64 ? sizeof(double) : sizeof(float);
    r->Np = i8 ? ((a.N + 15) / 16) * 16 : a.N;
    const int Kmax = a.K1 > a.K2 ? a.K1 : a.K2;

    // dG is the big one: slotsPerPass * N elements, twice over. Keep each
    // buffer near 512 MB and never wider than the caller's batch. FP32 / INT8
    // also keep the per-stream partial product (slots x K, 4 bytes) near 256 MB;
    // INT8's dG is int8, Np x (2 * slots + indicator columns), capX = slots.
    const std::size_t perSlot = i8 ? (std::size_t)r->Np * 3 : (std::size_t)a.N * r->esz;
    long long sp = (long long)((512ull << 20) / perSlot);
    if (!fp64) {
        const long long spp = (long long)((256ull << 20) / ((std::size_t)Kmax * 4 * (i8 ? 2 : 1)));
        if (sp > spp) sp = spp;
    }
    if (sp < 1) sp = 1;
    if (sp > a.maxSlots) sp = a.maxSlots;
    if (sp > 4096) sp = 4096;
    r->slotsPerPass = (int)sp;
    r->capX = i8 ? std::max((int)sp, 4) : 0;

    auto fail = [&]() -> Reducer* { destroy(r); return nullptr; };

    const std::size_t nSlots = (std::size_t)a.maxSlots;
    r->nSets = a.stagingSets < 1 ? 1 : a.stagingSets;
    r->hPacked.assign((std::size_t)r->nSets, nullptr);
    r->hLut.assign((std::size_t)r->nSets, nullptr);
    for (int s = 0; s < r->nSets; ++s) {
        if (cudaHostAlloc((void**)&r->hPacked[s], nSlots * r->bpv, cudaHostAllocDefault) != cudaSuccess) return fail();
        std::memset(r->hPacked[s], 0, nSlots * r->bpv);   // the padding bytes stay zero for good
        if (cudaHostAlloc((void**)&r->hLut[s], nSlots * 4 * sizeof(double), cudaHostAllocDefault) != cudaSuccess) return fail();
    }
    if (cudaHostAlloc((void**)&r->hLut2, nSlots * 4 * sizeof(double), cudaHostAllocDefault) != cudaSuccess) return fail();
    r->nDev = a.deviceSets < 1 ? 1 : a.deviceSets;
    r->hC1.assign((std::size_t)r->nDev, nullptr);
    r->hC2.assign((std::size_t)r->nDev, nullptr);
    r->hCnt.assign((std::size_t)r->nDev, nullptr);
    r->dPk.assign((std::size_t)r->nDev, nullptr);
    r->dLut.assign((std::size_t)r->nDev, nullptr);
    for (int d = 0; d < r->nDev; ++d) {
        if (cudaHostAlloc(&r->hC1[d], nSlots * a.K1 * r->cesz, cudaHostAllocDefault) != cudaSuccess) return fail();
        if (a.K2 > 0 && cudaHostAlloc(&r->hC2[d], nSlots * a.K2 * r->cesz, cudaHostAllocDefault) != cudaSuccess) return fail();
        if (a.nMask > 0 && cudaHostAlloc((void**)&r->hCnt[d], nSlots * a.nMask * 4 * sizeof(uint32_t), cudaHostAllocDefault) != cudaSuccess) return fail();
    }
    if (i8) {
        auto ha = [&](void** p, std::size_t n) { return cudaHostAlloc(p, n, cudaHostAllocDefault) == cudaSuccess; };
        if (!ha((void**)&r->hXf, nSlots * sizeof(int)) || !ha((void**)&r->hXn, nSlots * sizeof(int)) ||
            !ha((void**)&r->hXs, nSlots * 4 * sizeof(int)) || !ha((void**)&r->hXc, nSlots * 4 * sizeof(int)) ||
            !ha((void**)&r->hXd1, nSlots * 4 * sizeof(double)) || !ha((void**)&r->hXd2, nSlots * 4 * sizeof(double)))
            return fail();
    }

    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    const int nSl = i8 ? a.int8Slices : 1;
    const std::size_t bEsz = i8 ? (std::size_t)nSl : r->esz;   // bytes per B element (all slices)
    if (!dev(&r->dB1, (std::size_t)r->Np * a.K1 * bEsz)) return fail();
    if (!dev(&r->dC1, nSlots * a.K1 * r->cesz)) return fail();
    if (a.K2 > 0) {
        if (!dev(&r->dB2, (std::size_t)r->Np * a.K2 * bEsz)) return fail();
        if (!dev(&r->dC2, nSlots * a.K2 * r->cesz)) return fail();
    }
    for (int d = 0; d < r->nDev; ++d) {
        if (!dev((void**)&r->dPk[d],  nSlots * r->bpv)) return fail();
        if (!dev((void**)&r->dLut[d], nSlots * 4 * sizeof(double))) return fail();
    }
    if (!dev((void**)&r->dLut2, nSlots * 4 * sizeof(double))) return fail();
    const std::size_t gBytes = i8 ? (std::size_t)r->Np * (2 * (std::size_t)r->slotsPerPass + r->capX)
                                  : (std::size_t)r->slotsPerPass * perSlot;
    for (int b = 0; b < 2; ++b)
        if (!dev(&r->dG[b], gBytes)) return fail();
    if (!fp64) {
        const std::size_t rows = (std::size_t)r->slotsPerPass + (std::size_t)r->capX;
        for (int b = 0; b < 2; ++b)
            if (!dev(&r->dPart[b], rows * Kmax * 4)) return fail();
    }
    if (a.precision == Prec::FP32) {
        if (!dev((void**)&r->dCs1, (std::size_t)2 * a.K1 * sizeof(double))) return fail();
        if (!dev((void**)&r->dShift1, nSlots * 2 * sizeof(double))) return fail();
        if (a.K2 > 0) {
            if (!dev((void**)&r->dCs2, (std::size_t)2 * a.K2 * sizeof(double))) return fail();
            if (!dev((void**)&r->dShift2, nSlots * 2 * sizeof(double))) return fail();
        }
    }
    if (i8) {
        if (!dev((void**)&r->dE1, (std::size_t)a.K1 * sizeof(int))) return fail();
        if (a.K2 > 0 && !dev((void**)&r->dE2, (std::size_t)a.K2 * sizeof(int))) return fail();
        if (!dev((void**)&r->dXf, nSlots * sizeof(int)) || !dev((void**)&r->dXn, nSlots * sizeof(int)) ||
            !dev((void**)&r->dXs, nSlots * 4 * sizeof(int)) || !dev((void**)&r->dXc, nSlots * 4 * sizeof(int)) ||
            !dev((void**)&r->dXd1, nSlots * 4 * sizeof(double)) || !dev((void**)&r->dXd2, nSlots * 4 * sizeof(double)))
            return fail();
    }
    if (a.nMask > 0) {
        if (!dev((void**)&r->dMask, (std::size_t)a.nMask * r->words * sizeof(uint64_t))) return fail();
        if (!dev((void**)&r->dMaskPop, (std::size_t)a.nMask * sizeof(uint32_t))) return fail();
        if (!dev((void**)&r->dCnt, nSlots * a.nMask * 4 * sizeof(uint32_t))) return fail();
        if (cudaMemcpy(r->dMask, a.masks, (std::size_t)a.nMask * r->words * sizeof(uint64_t),
                       cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        std::vector<uint32_t> pop(a.nMask, 0);
        for (int m = 0; m < a.nMask; ++m) {
            unsigned long long c = 0;
            for (int w = 0; w < r->words; ++w)
                c += (unsigned long long)__builtin_popcountll(a.masks[(std::size_t)m * r->words + w] & 0x5555555555555555ULL);
            pop[m] = (uint32_t)c;
        }
        if (cudaMemcpy(r->dMaskPop, pop.data(), pop.size() * sizeof(uint32_t),
                       cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    }
    r->devBytes = db;

    // Right operands. FP64: as given. FP32: narrowed, plus colsum in fp64 for
    // the shift. INT8: the slices and per-column exponents.
    auto upload = [&](void* dst, const double* src, std::size_t n) -> bool {
        if (fp64) return cudaMemcpy(dst, src, n * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
        std::vector<float> bf(n);
        for (std::size_t i = 0; i < n; ++i) bf[i] = (float)src[i];
        return cudaMemcpy(dst, bf.data(), n * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess;
    };
    // FP32: column k goes up as float(b_ik - mean_k) and dCs gets
    // (colsum_k, mean_k) in fp64, 2 x K.
    auto uploadCentered = [&](void* dst, double* dCs, const double* src, int K) -> bool {
        std::vector<double> cs((std::size_t)2 * K, 0.0);
        std::vector<float> bf((std::size_t)a.N * K);
        for (int k = 0; k < K; ++k) {
            const double* c = src + (std::size_t)k * a.N;
            double t = 0.0;
            for (int i = 0; i < a.N; ++i) t += c[i];
            // Centre only a column whose entries cluster around a nonzero
            // mean (|mean| > sd: mu2, W X's intercept). Centring a column of
            // small entries around a larger mean would cost those entries
            // their relative precision (a MAC-1 marker reads one of them).
            double mu = t / (double)a.N, ss = 0.0;
            for (int i = 0; i < a.N; ++i) ss += (c[i] - mu) * (c[i] - mu);
            if (!(mu * mu * (double)a.N > ss)) mu = 0.0;
            cs[(std::size_t)k] = t; cs[(std::size_t)K + k] = mu;
            for (int i = 0; i < a.N; ++i) bf[(std::size_t)k * a.N + i] = (float)(c[i] - mu);
        }
        return cudaMemcpy(dst, bf.data(), bf.size() * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess &&
               cudaMemcpy(dCs, cs.data(), cs.size() * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
    };
    // B[:,k] = 2^(E_k-6) sum_s d_s 2^(-7s) + O(2^(E_k-6-7 nSl)), |d_s| <= 64;
    // every step exact in fp64 (power-of-two scalings, x - rint(x)).
    auto slices = [&](int8_t* dst, int* dE, const double* src, int K) -> bool {
        const int Np = r->Np;
        std::vector<int8_t> sl((std::size_t)nSl * Np * K, 0);
        std::vector<int> E((std::size_t)K, 0);
        for (int k = 0; k < K; ++k) {
            const double* c = src + (std::size_t)k * a.N;
            double mx = 0.0;
            for (int i = 0; i < a.N; ++i) mx = std::max(mx, std::fabs(c[i]));
            if (!(mx > 0.0) || !std::isfinite(mx)) { if (!std::isfinite(mx)) return false; continue; }
            int e = 0; std::frexp(mx, &e);              // mx < 2^e
            E[(std::size_t)k] = e;
            for (int i = 0; i < a.N; ++i) {
                double x = std::ldexp(c[i], 6 - e);     // |x| < 64
                for (int q = 0; q < nSl; ++q) {
                    const double d = std::rint(x);      // |d| <= 64
                    sl[((std::size_t)q * K + k) * Np + i] = (int8_t)d;
                    x = (x - d) * 128.0;                // |x - d| <= 1/2, exact
                }
            }
        }
        return cudaMemcpy(dst, sl.data(), sl.size(), cudaMemcpyHostToDevice) == cudaSuccess &&
               cudaMemcpy(dE, E.data(), E.size() * sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess;
    };
    if (i8) {
        if (!slices((int8_t*)r->dB1, r->dE1, a.B1, a.K1)) return fail();
        if (a.K2 > 0 && !slices((int8_t*)r->dB2, r->dE2, a.B2, a.K2)) return fail();
        r->dBs1 = (int8_t*)r->dB1; r->dBs2 = (int8_t*)r->dB2;
    } else {
        if (fp64) {
            if (!upload(r->dB1, a.B1, (std::size_t)a.N * a.K1)) return fail();
            if (a.K2 > 0 && !upload(r->dB2, a.B2, (std::size_t)a.N * a.K2)) return fail();
        } else {
            if (!uploadCentered(r->dB1, r->dCs1, a.B1, a.K1)) return fail();
            if (a.K2 > 0 && !uploadCentered(r->dB2, r->dCs2, a.B2, a.K2)) return fail();
        }
    }

    for (int b = 0; b < 2; ++b) {
        if (cudaStreamCreate(&r->st[b]) != cudaSuccess) return fail();
        for (int k = 0; k < 4; ++k)
            if (cudaEventCreate(&r->ev[b][k]) != cudaSuccess) return fail();
    }
    if (cudaEventCreate(&r->evUp) != cudaSuccess) return fail();
    if (cublasCreate(&r->cub) != CUBLAS_STATUS_SUCCESS) return fail();
    // No TF32, no split-k reduced-precision accumulation: the tolerance this
    // path is validated at assumes plain IEEE multiply-add in T.
    cublasSetMathMode(r->cub, CUBLAS_PEDANTIC_MATH);
    return r;
}

void destroy(Reducer* r)
{
    if (!r) return;
    if (r->cub) cublasDestroy(r->cub);
    for (int b = 0; b < 2; ++b) {
        for (int k = 0; k < 4; ++k) if (r->ev[b][k]) cudaEventDestroy(r->ev[b][k]);
        if (r->st[b]) cudaStreamDestroy(r->st[b]);
        if (r->dG[b]) cudaFree(r->dG[b]);
    }
    if (r->evUp) cudaEventDestroy(r->evUp);
    if (r->dB1) cudaFree(r->dB1);
    if (r->dB2) cudaFree(r->dB2);
    if (r->dC1) cudaFree(r->dC1);
    if (r->dC2) cudaFree(r->dC2);
    for (unsigned char* p : r->dPk) if (p) cudaFree(p);
    for (double* p : r->dLut) if (p) cudaFree(p);
    if (r->dLut2) cudaFree(r->dLut2);
    if (r->dMask) cudaFree(r->dMask);
    if (r->dMaskPop) cudaFree(r->dMaskPop);
    if (r->dCnt) cudaFree(r->dCnt);
    for (unsigned char* p : r->hPacked) if (p) cudaFreeHost(p);
    for (double* p : r->hLut) if (p) cudaFreeHost(p);
    if (r->hLut2)   cudaFreeHost(r->hLut2);
    for (void* p : r->hC1) if (p) cudaFreeHost(p);
    for (void* p : r->hC2) if (p) cudaFreeHost(p);
    for (uint32_t* p : r->hCnt) if (p) cudaFreeHost(p);
    for (int b = 0; b < 2; ++b) if (r->dPart[b]) cudaFree(r->dPart[b]);
    for (void* p : {(void*)r->dCs1, (void*)r->dCs2, (void*)r->dShift1, (void*)r->dShift2,
                    (void*)r->dE1, (void*)r->dE2, (void*)r->dXf, (void*)r->dXn, (void*)r->dXs,
                    (void*)r->dXc, (void*)r->dXd1, (void*)r->dXd2})
        if (p) cudaFree(p);
    for (void* p : {(void*)r->hXf, (void*)r->hXn, (void*)r->hXs, (void*)r->hXc,
                    (void*)r->hXd1, (void*)r->hXd2})
        if (p) cudaFreeHost(p);
    for (void* p : {(void*)r->dStatT, (void*)r->dVr, (void*)r->dStS, (void*)r->dStV, (void*)r->dStP, (void*)r->dStF})
        if (p) cudaFree(p);
    for (double* p : r->hVr) if (p) cudaFreeHost(p);
    for (double* p : r->hStS) if (p) cudaFreeHost(p);
    for (double* p : r->hStV) if (p) cudaFreeHost(p);
    for (double* p : r->hStP) if (p) cudaFreeHost(p);
    for (unsigned char* p : r->hStF) if (p) cudaFreeHost(p);
    delete r;
}

unsigned char*  packed(Reducer* r, int s)      { return (r && s >= 0 && s < r->nSets) ? r->hPacked[s] : nullptr; }
double*         lut(Reducer* r, int s)         { return (r && s >= 0 && s < r->nSets) ? r->hLut[s] : nullptr; }
std::size_t     bytesPerSlot(const Reducer* r) { return r ? r->bpv : 0; }
int             stagingSets(const Reducer* r)  { return r ? r->nSets : 0; }
bool            isFp64(const Reducer* r)       { return r ? r->outD : false; }
namespace {
// -1 = the set of the last reduce(); out of range -> -1 (callers get nullptr)
inline int devSetOf(const Reducer* r, int d) {
    if (d < 0) return r->lastDev;
    return d < r->nDev ? d : -1;
}
}  // namespace
const float*    outCf(const Reducer* r, int d)  { const int k = r ? devSetOf(r, d) : -1; return (k >= 0 && !r->outD) ? (const float*)r->hC1[k] : nullptr; }
const double*   outCd(const Reducer* r, int d)  { const int k = r ? devSetOf(r, d) : -1; return (k >= 0 &&  r->outD) ? (const double*)r->hC1[k] : nullptr; }
const float*    outC2f(const Reducer* r, int d) { const int k = r ? devSetOf(r, d) : -1; return (k >= 0 && !r->outD && r->K2 > 0) ? (const float*)r->hC2[k] : nullptr; }
const double*   outC2d(const Reducer* r, int d) { const int k = r ? devSetOf(r, d) : -1; return (k >= 0 &&  r->outD && r->K2 > 0) ? (const double*)r->hC2[k] : nullptr; }
std::size_t     ldC(const Reducer* r)          { return r ? (std::size_t)r->maxSlots : 0; }
const uint32_t* outCounts(const Reducer* r, int d) { const int k = r ? devSetOf(r, d) : -1; return (k >= 0 && r->nMask > 0) ? r->hCnt[k] : nullptr; }
std::size_t     deviceBytes(const Reducer* r)  { return r ? r->devBytes : 0; }
const void*     devicePacked(const Reducer* r, int d) { const int k = r ? devSetOf(r, d) : -1; return k >= 0 ? (const void*)r->dPk[k] : nullptr; }
const void*     deviceLut(const Reducer* r, int d)    { const int k = r ? devSetOf(r, d) : -1; return k >= 0 ? (const void*)r->dLut[k] : nullptr; }
void*           deviceStream(const Reducer* r) { return r ? (void*)r->st[0] : nullptr; }
int             deviceSets(const Reducer* r)   { return r ? r->nDev : 0; }

// ---------------------------------------------------------------------------
// Per-pair statistics
// ---------------------------------------------------------------------------
namespace {
thread_local std::string statsErr;

// The kernel on nS slots of the device result buffers into device set d's
// pinned stats buffers, on stream s. Asynchronous; the caller synchronises.
bool statsLaunch(Reducer* r, int nS, int d, cudaStream_t s)
{
    const long long nPair = (long long)nS * r->nStat;
    if (nPair <= 0) return true;
    const int blocks = (int)((nPair + 255) / 256);
    pair_stats<<<blocks, 256, 0, s>>>((const double*)r->dC1, (const double*)r->dC2, (std::size_t)r->maxSlots, nS,
                                      r->dStatT, r->nStat, r->dVr, r->statCutoff, r->pRelTol,
                                      r->patXZ, r->patSaz, r->patZxz, r->patGwz,
                                      r->dStS, r->dStV, r->dStP, r->dStF);
    CKR(cudaGetLastError());
    const std::size_t n = (std::size_t)nPair;
    CKR(cudaMemcpyAsync(r->hStS[d], r->dStS, n * sizeof(double), cudaMemcpyDeviceToHost, s));
    CKR(cudaMemcpyAsync(r->hStV[d], r->dStV, n * sizeof(double), cudaMemcpyDeviceToHost, s));
    CKR(cudaMemcpyAsync(r->hStP[d], r->dStP, n * sizeof(double), cudaMemcpyDeviceToHost, s));
    CKR(cudaMemcpyAsync(r->hStF[d], r->dStF, n, cudaMemcpyDeviceToHost, s));
    return true;
}
}  // namespace

std::string statsLastError() { return statsErr; }

bool statsSetup(Reducer* r, const StatsArgs& a)
{
    if (!r) { statsErr = "no reducer"; return false; }
    if (r->stats) { statsErr = "statsSetup called twice"; return false; }
    if (a.nTraits <= 0 || a.traits == nullptr) { statsErr = "no traits"; return false; }
    if (a.nTraits != r->nMask) { statsErr = "nTraits must equal the reducer's mask count"; return false; }
    if (!r->outD || r->K2 <= 0) { statsErr = "the reducer has no fp64 C2 results"; return false; }
    for (int b = 0; b < a.nTraits; ++b) {
        const StatsTrait& t = a.traits[b];
        if (t.p < 1 || t.p > STATS_PMAX) { statsErr = "trait p out of range"; return false; }
        if (t.rowZ < 0 || t.rowZ + t.p > r->K1 || t.rowW < 0 || t.rowW + t.p > r->K1 ||
            t.rowGR < 0 || t.rowGR >= r->K1 || t.colG2 < 0 || t.colG2 >= r->K2) {
            statsErr = "trait column indices outside C1 / C2"; return false;
        }
    }
    if (a.patXZ < 0 || a.patXZ > 4 || a.patSaz < 0 || a.patSaz > 4 ||
        a.patZxz < 0 || a.patZxz > 1 || a.patGwz < 0 || a.patGwz > 1) { statsErr = "pattern id out of range"; return false; }
    const std::size_t nPair = (std::size_t)r->maxSlots * (std::size_t)a.nTraits;
    r->nStat = a.nTraits;
    r->statCutoff = a.statCutoff; r->pRelTol = a.pRelTol;
    r->patXZ = a.patXZ; r->patSaz = a.patSaz; r->patZxz = a.patZxz; r->patGwz = a.patGwz;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; statsErr = "cudaMalloc failed"; return false; }
        r->devBytes += n; return true;
    };
    auto pin = [&](void** p, std::size_t n) {
        if (cudaHostAlloc(p, n, cudaHostAllocDefault) != cudaSuccess) { *p = nullptr; statsErr = "cudaHostAlloc failed"; return false; }
        return true;
    };
    if (!dev((void**)&r->dStatT, (std::size_t)a.nTraits * sizeof(StatsTrait))) return false;
    if (cudaMemcpy(r->dStatT, a.traits, (std::size_t)a.nTraits * sizeof(StatsTrait), cudaMemcpyHostToDevice) != cudaSuccess) {
        statsErr = "trait upload failed"; return false;
    }
    if (!dev((void**)&r->dVr, nPair * sizeof(double))) return false;
    if (!dev((void**)&r->dStS, nPair * sizeof(double))) return false;
    if (!dev((void**)&r->dStV, nPair * sizeof(double))) return false;
    if (!dev((void**)&r->dStP, nPair * sizeof(double))) return false;
    if (!dev((void**)&r->dStF, nPair)) return false;
    r->hVr.assign((std::size_t)r->nSets, nullptr);
    for (int s = 0; s < r->nSets; ++s) {
        if (!pin((void**)&r->hVr[s], nPair * sizeof(double))) return false;
        for (std::size_t i = 0; i < nPair; ++i) r->hVr[s][i] = 1.0;
    }
    r->hStS.assign((std::size_t)r->nDev, nullptr); r->hStV.assign((std::size_t)r->nDev, nullptr);
    r->hStP.assign((std::size_t)r->nDev, nullptr); r->hStF.assign((std::size_t)r->nDev, nullptr);
    for (int d = 0; d < r->nDev; ++d) {
        if (!pin((void**)&r->hStS[d], nPair * sizeof(double))) return false;
        if (!pin((void**)&r->hStV[d], nPair * sizeof(double))) return false;
        if (!pin((void**)&r->hStP[d], nPair * sizeof(double))) return false;
        if (!pin((void**)&r->hStF[d], nPair)) return false;
    }
    r->stats = true;
    return true;
}

void statsDisable(Reducer* r) { if (r) r->stats = false; }

double* statsVr(Reducer* r, int s)
{
    return (r && r->stats && s >= 0 && s < r->nSets) ? r->hVr[(std::size_t)s] : nullptr;
}
const double*        statsS(const Reducer* r, int d)     { const int k = (r && r->stats) ? devSetOf(r, d) : -1; return k >= 0 ? r->hStS[(std::size_t)k] : nullptr; }
const double*        statsVar2(const Reducer* r, int d)  { const int k = (r && r->stats) ? devSetOf(r, d) : -1; return k >= 0 ? r->hStV[(std::size_t)k] : nullptr; }
const double*        statsP(const Reducer* r, int d)     { const int k = (r && r->stats) ? devSetOf(r, d) : -1; return k >= 0 ? r->hStP[(std::size_t)k] : nullptr; }
const unsigned char* statsFlags(const Reducer* r, int d) { const int k = (r && r->stats) ? devSetOf(r, d) : -1; return k >= 0 ? r->hStF[(std::size_t)k] : nullptr; }
double               statsSeconds(const Reducer* r)      { return r ? r->tStats : 0.0; }

bool statsSelfTest(Reducer* r, const double* C1, const double* C2, const double* vr, int nS)
{
    if (!r || !r->stats) { statsErr = "stats not set up"; return false; }
    if (nS <= 0 || nS > r->maxSlots) { statsErr = "nSlots out of range"; return false; }
    const std::size_t ldb = (std::size_t)r->maxSlots * sizeof(double);
    if (cudaMemcpy2D(r->dC1, ldb, C1, ldb, (std::size_t)nS * sizeof(double), (std::size_t)r->K1, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy2D(r->dC2, ldb, C2, ldb, (std::size_t)nS * sizeof(double), (std::size_t)r->K2, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(r->dVr, vr, (std::size_t)nS * r->nStat * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) {
        statsErr = "self-test upload failed"; return false;
    }
    r->lastDev = 0;
    if (!statsLaunch(r, nS, 0, r->st[0])) { statsErr = lastErr; return false; }
    if (cudaStreamSynchronize(r->st[0]) != cudaSuccess) { statsErr = "self-test kernel failed"; return false; }
    return true;
}

bool statsErfcTest(Reducer* r, const double* stat, int n, double* p)
{
    if (!r || n <= 0) return false;
    double* ds = nullptr; double* dp = nullptr;
    if (cudaMalloc((void**)&ds, (std::size_t)n * sizeof(double)) != cudaSuccess) return false;
    if (cudaMalloc((void**)&dp, (std::size_t)n * sizeof(double)) != cudaSuccess) { cudaFree(ds); return false; }
    bool ok = cudaMemcpy(ds, stat, (std::size_t)n * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
    if (ok) {
        erfc_tail<<<(n + 255) / 256, 256>>>(ds, n, dp);
        ok = cudaGetLastError() == cudaSuccess && cudaDeviceSynchronize() == cudaSuccess;
    }
    if (ok) ok = cudaMemcpy(p, dp, (std::size_t)n * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
    cudaFree(ds); cudaFree(dp);
    return ok;
}


bool bindDevice(int t_device) { return cudaSetDevice(t_device) == cudaSuccess; }

std::string setBlockingSync(int t_device)
{
    cudaError_t e = cudaSetDevice(t_device);
    if (e == cudaSuccess) e = cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
    return e == cudaSuccess ? std::string() : std::string(cudaGetErrorString(e));
}

void timings(const Reducer* r, double* h2d, double* dec, double* gemm, double* d2h, double* popc)
{
    if (!r) return;
    if (h2d)  *h2d  = r->tH2D;
    if (dec)  *dec  = r->tDec;
    if (gemm) *gemm = r->tGemm;
    if (d2h)  *d2h  = r->tD2H;
    if (popc) *popc = r->tPopc;
}

namespace {

// Launch shape of acc_f32 / acc_i8: x over the sc slots, y over the K columns.
inline dim3 accGrid(int sc, int K)
{
    const int gx = (sc + 255) / 256;
    return dim3((unsigned)(gx < 64 ? gx : 64), (unsigned)(K < 65535 ? K : 65535));
}

// FP32 pass (gpu_scan_lowp.cuh): shifted float decode, one SGEMM per
// kChunkF32 samples into dPart[b], each chunk added into the fp64 result.
// Records ev[b][1..3] as the fp64 pass does.
bool passF32(Reducer* r, int b, cudaStream_t s, const unsigned char* pk, const double* lu,
             const double* lu2, int s0, int sc)
{
    const float onef = 1.f, zerof = 0.f;
    float* g = (float*)r->dG[b];
    float* part = (float*)r->dPart[b];
    auto gemm = [&](const float* B, int K, double* C, const double* shift, const double* cs, int chunk) -> bool {
        // floor(N / chunk) chunks of `chunk` samples, the last one also taking
        // the remainder (a narrow GEMM of its own costs more than the
        // remainder's share: cuBLAS picks small tiles for it).
        const int nc = std::max(1, r->N / chunk);
        for (int c = 0; c < nc; ++c) {
            const int c0 = c * chunk;
            const int kk = (c == nc - 1) ? r->N - c0 : chunk;
            CBR(cublasSgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, K, kk,
                            &onef, g + c0, r->N, B + c0, r->N, &zerof, part, sc));
            lowp::acc_f32<<<accGrid(sc, K), 256, 0, s>>>(part, sc, K, C, (std::size_t)r->maxSlots,
                                                                   c0 == 0 ? 1 : 0, shift, cs);
            CKR(cudaGetLastError());
        }
        return true;
    };
    lowp::decode_f32_shift<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g, r->dShift1 + 2 * (std::size_t)s0);
    CKR(cudaGetLastError());
    CKR(cudaEventRecord(r->ev[b][1], s));
    if (!gemm((const float*)r->dB1, r->K1, (double*)r->dC1 + s0, r->dShift1 + 2 * (std::size_t)s0, r->dCs1,
              lowp::kChunkF32)) return false;
    CKR(cudaEventRecord(r->ev[b][2], s));
    if (r->K2 > 0) {
        lowp::decode_f32_shift<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g, r->dShift2 + 2 * (std::size_t)s0);
        CKR(cudaGetLastError());
        if (!gemm((const float*)r->dB2, r->K2, (double*)r->dC2 + s0, r->dShift2 + 2 * (std::size_t)s0, r->dCs2,
                  lowp::kChunkF32Sq)) return false;
        CKR(cudaEventRecord(r->ev[b][3], s));
    }
    return true;
}

// INT8 pass (gpu_scan_lowp.cuh): A = [g_int (sc) | indicators (nx) | g2_int
// (sc, binary only)], one int8 GEMM per slice of B into dPart[b], each slice
// recombined into the fp64 result, smallest slice first.
bool passI8(Reducer* r, int b, cudaStream_t s, const unsigned char* pk, const double* lu,
            const double* lu2, int s0, int sc, int x0, int nx)
{
    const int onei = 1, zeroi = 0;
    const int Np = r->Np, nSl = r->int8Slices;
    int8_t* A = (int8_t*)r->dG[b];
    int* part = (int*)r->dPart[b];
    const int sc2 = (r->K2 > 0) ? sc : 0;
    lowp::decode_i8<<<sc + nx + sc2, 256, 0, s>>>(pk, r->bpv, r->N, Np, lu, lu2, sc, nx,
                                                 r->dXs + x0, r->dXc + x0, A);
    CKR(cudaGetLastError());
    CKR(cudaEventRecord(r->ev[b][1], s));
    auto gemm = [&](const int8_t* Aop, int rows, int rowOff, int xOff, const int8_t* Bs, int K,
                    const int* E, double* C, const double* delta) -> bool {
        for (int q = nSl - 1; q >= 0; --q) {
            CBR(cublasGemmEx(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, rows, K, Np,
                             &onei, Aop, CUDA_R_8I, Np,
                             Bs + (std::size_t)q * Np * K, CUDA_R_8I, Np,
                             &zeroi, part, CUDA_R_32I, rows,
                             CUBLAS_COMPUTE_32I, CUBLAS_GEMM_DEFAULT));
            lowp::acc_i8<<<accGrid(sc, K), 256, 0, s>>>(part, rows, rowOff, xOff, sc, K, E, q,
                                                                  C, (std::size_t)r->maxSlots,
                                                                  q == nSl - 1 ? 1 : 0,
                                                                  r->dXf + s0, r->dXn + s0, delta);
            CKR(cudaGetLastError());
        }
        return true;
    };
    if (!gemm(A, sc + nx, 0, sc, r->dBs1, r->K1, r->dE1, (double*)r->dC1 + s0, r->dXd1 + x0)) return false;
    CKR(cudaEventRecord(r->ev[b][2], s));
    if (r->K2 > 0) {
        if (!gemm(A + (std::size_t)sc * Np, nx + sc, nx, 0, r->dBs2, r->K2, r->dE2,
                  (double*)r->dC2 + s0, r->dXd2 + x0)) return false;
        CKR(cudaEventRecord(r->ev[b][3], s));
    }
    return true;
}

}  // namespace

bool reduce(Reducer* r, int t_nSlots, int t_set, int t_devSet)
{
    if (!r) return false;
    if (t_nSlots <= 0) return true;
    if (t_nSlots > r->maxSlots) { lastErr = "nSlots > maxSlots"; return false; }
    if (t_set < 0 || t_set >= r->nSets) { lastErr = "staging set out of range"; return false; }
    if (t_devSet < 0) t_devSet = 0;
    if (t_devSet >= r->nDev) { lastErr = "device set out of range"; return false; }
    const unsigned char* hPk = r->hPacked[(std::size_t)t_set];
    const double*        hLu = r->hLut[(std::size_t)t_set];
    unsigned char* const dPk  = r->dPk[(std::size_t)t_devSet];
    double* const        dLut = r->dLut[(std::size_t)t_devSet];
    void* const          hC1  = r->hC1[(std::size_t)t_devSet];
    void* const          hC2  = r->hC2[(std::size_t)t_devSet];
    uint32_t* const      hCnt = r->hCnt[(std::size_t)t_devSet];
    r->lastDev = t_devSet;

    const double oned = 1.0, zerod = 0.0;
    const bool x4 = ((r->N & 3) == 0);
    const std::size_t nS = (std::size_t)t_nSlots;

    // Squared table for the second GEMM: fd[c]*fd[c] in double, the product the
    // CPU kernel's Gb % Gb forms per cell.
    if (r->K2 > 0)
        for (std::size_t i = 0; i < nS * 4; ++i) r->hLut2[i] = hLu[i] * hLu[i];

    // Passes: slotsPerPass markers each; INT8 also caps the indicator columns
    // of a pass at capX and builds them here (gpu_scan_lowp.cuh).
    struct Pass { int s0, sc, x0, nx; };
    std::vector<Pass> passes;
    int nxTot = 0;
    if (r->prec == Prec::INT8) {
        Pass cur{0, 0, 0, 0};
        for (int j = 0; j < t_nSlots; ++j) {
            int nxj = 0;
            double d1[4], d2[4];
            for (int c = 0; c < 4; ++c) {
                const double v = hLu[(std::size_t)j * 4 + c], n1 = std::rint(v);
                const double v2 = (r->K2 > 0) ? r->hLut2[(std::size_t)j * 4 + c] : 0.0, n2 = std::rint(v2);
                if (!(std::fabs(n1) <= lowp::kI8MaxA) || !(std::fabs(n2) <= lowp::kI8MaxA)) {
                    lastErr = "int8 scan: dosage table entry outside [-4, 4]"; return false;
                }
                d1[c] = v - n1; d2[c] = v2 - n2;
                if (d1[c] != 0.0 || d2[c] != 0.0) ++nxj;
            }
            if (cur.sc > 0 && (cur.sc + 1 > r->slotsPerPass || cur.nx + nxj > r->capX)) {
                passes.push_back(cur);
                cur = Pass{j, 0, nxTot, 0};
            }
            r->hXf[j] = cur.nx; r->hXn[j] = nxj;
            for (int c = 0; c < 4; ++c) {
                if (d1[c] == 0.0 && d2[c] == 0.0) continue;
                r->hXs[nxTot] = cur.sc; r->hXc[nxTot] = c;
                r->hXd1[nxTot] = d1[c]; r->hXd2[nxTot] = d2[c];
                ++nxTot; ++cur.nx;
            }
            ++cur.sc;
        }
        passes.push_back(cur);
    } else {
        for (int s0 = 0; s0 < t_nSlots; s0 += r->slotsPerPass)
            passes.push_back(Pass{s0, std::min(r->slotsPerPass, t_nSlots - s0), 0, 0});
    }

    // Everything from the previous reduce() has been harvested already; both
    // streams are idle at the top of a call.
    cudaEvent_t a0 = nullptr, a1 = nullptr;
    const bool timed = (cudaEventCreate(&a0) == cudaSuccess) && (cudaEventCreate(&a1) == cudaSuccess);
    if (timed) cudaEventRecord(a0, r->st[0]);
    CKR(cudaMemcpyAsync(dPk, hPk, nS * r->bpv, cudaMemcpyHostToDevice, r->st[0]));
    CKR(cudaMemcpyAsync(dLut, hLu, nS * 4 * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
    if (r->K2 > 0)
        CKR(cudaMemcpyAsync(r->dLut2, r->hLut2, nS * 4 * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
    if (r->prec == Prec::INT8) {
        CKR(cudaMemcpyAsync(r->dXf, r->hXf, nS * sizeof(int), cudaMemcpyHostToDevice, r->st[0]));
        CKR(cudaMemcpyAsync(r->dXn, r->hXn, nS * sizeof(int), cudaMemcpyHostToDevice, r->st[0]));
        if (nxTot > 0) {
            const std::size_t nx = (std::size_t)nxTot;
            CKR(cudaMemcpyAsync(r->dXs, r->hXs, nx * sizeof(int), cudaMemcpyHostToDevice, r->st[0]));
            CKR(cudaMemcpyAsync(r->dXc, r->hXc, nx * sizeof(int), cudaMemcpyHostToDevice, r->st[0]));
            CKR(cudaMemcpyAsync(r->dXd1, r->hXd1, nx * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
            CKR(cudaMemcpyAsync(r->dXd2, r->hXd2, nx * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
        }
    }
    if (r->stats)
        CKR(cudaMemcpyAsync(r->dVr, r->hVr[(std::size_t)t_set], nS * r->nStat * sizeof(double),
                            cudaMemcpyHostToDevice, r->st[0]));
    CKR(cudaEventRecord(r->evUp, r->st[0]));
    if (timed) {
        cudaEventRecord(a1, r->st[0]);
    }
    // Code counts on the second stream, behind the upload.
    cudaEvent_t p0 = nullptr, p1 = nullptr;
    if (r->nMask > 0) {
        CKR(cudaStreamWaitEvent(r->st[1], r->evUp, 0));
        if (timed) { cudaEventCreate(&p0); cudaEventCreate(&p1); cudaEventRecord(p0, r->st[1]); }
        count_codes<<<t_nSlots, 256, 0, r->st[1]>>>(dPk, r->bpv, r->words, r->dMask,
                                                    r->dMaskPop, r->nMask, r->dCnt);
        CKR(cudaGetLastError());
        if (timed) cudaEventRecord(p1, r->st[1]);
        CKR(cudaMemcpyAsync(hCnt, r->dCnt, nS * r->nMask * 4 * sizeof(uint32_t),
                            cudaMemcpyDeviceToHost, r->st[1]));
    }

    for (int pass = 0; pass < (int)passes.size(); ++pass) {
        const int s0 = passes[(std::size_t)pass].s0;
        const int sc = passes[(std::size_t)pass].sc;
        const int b  = pass & 1;
        cudaStream_t s = r->st[b];

        // Buffer b is still in flight from two passes ago; wait, then bank its
        // times (free: the wait had to happen anyway).
        CKR(cudaStreamSynchronize(s));
        harvest(r, b);
        CKR(cudaStreamWaitEvent(s, r->evUp, 0));

        const unsigned char* pk = dPk + (std::size_t)s0 * r->bpv;
        const double* lu  = dLut  + (std::size_t)s0 * 4;
        const double* lu2 = r->dLut2 + (std::size_t)s0 * 4;

        CKR(cudaEventRecord(r->ev[b][0], s));
        if (r->prec != Prec::FP64) {
            CBR(cublasSetStream(r->cub, s));
            if (r->prec == Prec::FP32) { if (!passF32(r, b, s, pk, lu, lu2, s0, sc)) return false; }
            else if (!passI8(r, b, s, pk, lu, lu2, s0, sc, passes[(std::size_t)pass].x0,
                             passes[(std::size_t)pass].nx)) return false;
            r->evLive[b] = true;
            continue;
        }
        {
            double* g = (double*)r->dG[b];
            if (r->decodeX2) decode_lut_x2<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
            else if (x4) decode_lut_x4<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
            else    decode_lut_any<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
        }
        CKR(cudaGetLastError());
        CKR(cudaEventRecord(r->ev[b][1], s));

        CBR(cublasSetStream(r->cub, s));
        // C1(sc x K1) = dG^T (sc x N) * dB1 (N x K1), into rows [s0, s0+sc) of
        // the maxSlots x K1 result.
        CBR(cublasDgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K1, r->N,
                        &oned, (const double*)r->dG[b], r->N,
                        (const double*)r->dB1, r->N,
                        &zerod, (double*)r->dC1 + s0, r->maxSlots));
        CKR(cudaEventRecord(r->ev[b][2], s));

        if (r->K2 > 0) {
            // Same buffer, squared table, second right operand.
            {
                double* g = (double*)r->dG[b];
                if (r->decodeX2) decode_lut_x2<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                else if (x4) decode_lut_x4<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                else    decode_lut_any<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                CKR(cudaGetLastError());
                CBR(cublasDgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K2, r->N,
                                &oned, (const double*)r->dG[b], r->N,
                                (const double*)r->dB2, r->N,
                                &zerod, (double*)r->dC2 + s0, r->maxSlots));
            }
            CKR(cudaEventRecord(r->ev[b][3], s));
        }
        r->evLive[b] = true;
    }

    CKR(cudaStreamSynchronize(r->st[0]));
    CKR(cudaStreamSynchronize(r->st[1]));
    harvest(r, 0);
    harvest(r, 1);
    if (timed) {
        float ms = 0;
        if (cudaEventElapsedTime(&ms, a0, a1) == cudaSuccess) r->tH2D += ms * 1e-3;
        if (p0 && p1 && cudaEventElapsedTime(&ms, p0, p1) == cudaSuccess) r->tPopc += ms * 1e-3;
    }

    cudaEvent_t z0 = nullptr, z1 = nullptr;
    const bool timed2 = timed && (cudaEventCreate(&z0) == cudaSuccess) && (cudaEventCreate(&z1) == cudaSuccess);
    if (timed2) cudaEventRecord(z0, r->st[0]);
    // Only the first t_nSlots rows of each of the K columns are live.
    CKR(cudaMemcpy2DAsync(hC1, (std::size_t)r->maxSlots * r->cesz,
                          r->dC1, (std::size_t)r->maxSlots * r->cesz,
                          nS * r->cesz, (std::size_t)r->K1,
                          cudaMemcpyDeviceToHost, r->st[0]));
    if (r->K2 > 0)
        CKR(cudaMemcpy2DAsync(hC2, (std::size_t)r->maxSlots * r->cesz,
                              r->dC2, (std::size_t)r->maxSlots * r->cesz,
                              nS * r->cesz, (std::size_t)r->K2,
                              cudaMemcpyDeviceToHost, r->st[0]));
    if (timed2) cudaEventRecord(z1, r->st[0]);
    // Per-pair statistics (statsSetup): both streams are idle here, so every
    // C1 / C2 column is final; the kernel and its D2H queue behind the result
    // copies on the same stream.
    cudaEvent_t y0 = nullptr, y1 = nullptr;
    const bool timed3 = r->stats && timed && (cudaEventCreate(&y0) == cudaSuccess) && (cudaEventCreate(&y1) == cudaSuccess);
    if (r->stats) {
        if (timed3) cudaEventRecord(y0, r->st[0]);
        if (!statsLaunch(r, t_nSlots, t_devSet, r->st[0])) return false;
        if (timed3) cudaEventRecord(y1, r->st[0]);
    }
    CKR(cudaStreamSynchronize(r->st[0]));
    if (timed2) {
        float ms = 0;
        if (cudaEventElapsedTime(&ms, z0, z1) == cudaSuccess) r->tD2H += ms * 1e-3;
    }
    if (timed3) {
        float ms = 0;
        if (cudaEventElapsedTime(&ms, y0, y1) == cudaSuccess) r->tStats += ms * 1e-3;
    }
    if (y0) cudaEventDestroy(y0);
    if (y1) cudaEventDestroy(y1);
    if (a0) cudaEventDestroy(a0);
    if (a1) cudaEventDestroy(a1);
    if (p0) cudaEventDestroy(p0);
    if (p1) cudaEventDestroy(p1);
    if (z0) cudaEventDestroy(z0);
    if (z1) cudaEventDestroy(z1);
    return true;
}

}  // namespace gpu2
}  // namespace saige
