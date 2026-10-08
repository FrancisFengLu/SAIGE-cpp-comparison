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
// Decode and GEMM both run in one type T -- float or double -- chosen at
// create(). Measured on this V100 (14.28 TFLOP/s fp32, 7.19 fp64): fp64 costs
// about 40% more GEMM time and one extra pass of device bandwidth in the
// decode. In a real step-2 run the GEMM is single-digit seconds against
// a hundred-plus seconds of reading and decoding the .bed, so that is a few
// percent of wall clock and it buys agreement with the CPU path near 1e-15
// instead of near 1e-6. fp64 is therefore the default. fp32 stays available
// for devices where fp64 runs at 1/32 or 1/64 rate (L4, A10, consumer parts)
// and for measuring the fp32 error itself.
//
// The mode comes from CreateArgs::precision (config key gpuPrecisionScan,
// gpu_precision.hpp): FP64 and FP32 are the two instantiations of T above.
// INT8 (the Ozaki-style split) is not implemented yet: scanSupports() says no
// and create() refuses it. Its plug-in points are marked TODO(precision:scan)
// -- the operand split / upload in create() and the GEMM block in reduce().
//
// The dosage TABLE is double in both modes and narrowed inside the kernel for
// fp32. Three of its four entries are exact small integers, but the fourth is
// the imputed mean 2*altFreq -- narrowing that on the host would put a 6e-8
// relative error into every missing cell before the reduction starts, which is
// a difference from the CPU's INPUT rather than from its arithmetic. It costs
// 32 bytes per marker to carry.

#include "gpu_step2.hpp"

#include <cuda_runtime.h>
#include <cublas_v2.h>

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

// float -> one 16-byte store; double -> two. Both aligned because N % 4 == 0
// on this path, so every marker's column starts on a 16-byte boundary.
__device__ __forceinline__ void storeQuad(float* p, float a, float b, float c, float d)
{
    *reinterpret_cast<float4*>(p) = make_float4(a, b, c, d);
}
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

    double tH2D = 0, tDec = 0, tGemm = 0, tD2H = 0, tPopc = 0;
    std::size_t devBytes = 0;
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
        // TODO(precision:scan): return true once the int8 split below is in.
        // The sparse-GRM cross terms (gpu_sparse.cu, spqSupports) are part of
        // this stage too and have their own switch.
        case Prec::INT8: return false;
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
    // TODO(precision:scan): INT8 -- the decode can stay in double or go
    // straight to int8 (genotypes are exact small integers; the imputed mean
    // of a missing call is not, see the LUT note above), B1 / B2 are split
    // into int8Slices slices at upload, and the results are double (outD is
    // already true for INT8, so isFp64() says so and main.cpp reads outCd /
    // outC2d; hC1 / hC2 / dC1 / dC2 must then be sized in double, not esz).
    const bool fp64 = (a.precision == Prec::FP64);
    r->N = a.N; r->K1 = a.K1; r->K2 = a.K2; r->maxSlots = a.maxSlots; r->fp64 = fp64;
    r->outD = (a.precision != Prec::FP32);
    r->decodeX2 = a.decodeX2 && fp64 && ((a.N & 3) == 0);
    r->nMask = a.nMask;
    r->words = maskWords(a.N);
    r->bpv = (std::size_t)r->words * 8;           // >= (N+3)/4, whole 64-bit words
    r->esz = fp64 ? sizeof(double) : sizeof(float);

    // dG is the big one: slotsPerPass * N elements, twice over. Keep each
    // buffer near 512 MB and never wider than the caller's batch.
    const std::size_t perSlot = (std::size_t)a.N * r->esz;
    long long sp = (long long)((512ull << 20) / perSlot);
    if (sp < 1) sp = 1;
    if (sp > a.maxSlots) sp = a.maxSlots;
    if (sp > 4096) sp = 4096;
    r->slotsPerPass = (int)sp;

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
        if (cudaHostAlloc(&r->hC1[d], nSlots * a.K1 * r->esz, cudaHostAllocDefault) != cudaSuccess) return fail();
        if (a.K2 > 0 && cudaHostAlloc(&r->hC2[d], nSlots * a.K2 * r->esz, cudaHostAllocDefault) != cudaSuccess) return fail();
        if (a.nMask > 0 && cudaHostAlloc((void**)&r->hCnt[d], nSlots * a.nMask * 4 * sizeof(uint32_t), cudaHostAllocDefault) != cudaSuccess) return fail();
    }

    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    if (!dev(&r->dB1, (std::size_t)a.N * a.K1 * r->esz)) return fail();
    if (!dev(&r->dC1, nSlots * a.K1 * r->esz)) return fail();
    if (a.K2 > 0) {
        if (!dev(&r->dB2, (std::size_t)a.N * a.K2 * r->esz)) return fail();
        if (!dev(&r->dC2, nSlots * a.K2 * r->esz)) return fail();
    }
    for (int d = 0; d < r->nDev; ++d) {
        if (!dev((void**)&r->dPk[d],  nSlots * r->bpv)) return fail();
        if (!dev((void**)&r->dLut[d], nSlots * 4 * sizeof(double))) return fail();
    }
    if (!dev((void**)&r->dLut2, nSlots * 4 * sizeof(double))) return fail();
    for (int b = 0; b < 2; ++b)
        if (!dev(&r->dG[b], (std::size_t)r->slotsPerPass * a.N * r->esz)) return fail();
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

    auto upload = [&](void* dst, const double* src, std::size_t n) -> bool {
        if (fp64) return cudaMemcpy(dst, src, n * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess;
        std::vector<float> bf(n);
        for (std::size_t i = 0; i < n; ++i) bf[i] = (float)src[i];
        return cudaMemcpy(dst, bf.data(), n * sizeof(float), cudaMemcpyHostToDevice) == cudaSuccess;
    };
    if (!upload(r->dB1, a.B1, (std::size_t)a.N * a.K1)) return fail();
    if (a.K2 > 0 && !upload(r->dB2, a.B2, (std::size_t)a.N * a.K2)) return fail();

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

    const float  onef = 1.f, zerof = 0.f;
    const double oned = 1.0, zerod = 0.0;
    const bool x4 = ((r->N & 3) == 0);
    const std::size_t nS = (std::size_t)t_nSlots;

    // Squared table for the second GEMM: fd[c]*fd[c] in double, the product the
    // CPU kernel's Gb % Gb forms per cell.
    if (r->K2 > 0)
        for (std::size_t i = 0; i < nS * 4; ++i) r->hLut2[i] = hLu[i] * hLu[i];

    // Everything from the previous reduce() has been harvested already; both
    // streams are idle at the top of a call.
    cudaEvent_t a0 = nullptr, a1 = nullptr;
    const bool timed = (cudaEventCreate(&a0) == cudaSuccess) && (cudaEventCreate(&a1) == cudaSuccess);
    if (timed) cudaEventRecord(a0, r->st[0]);
    CKR(cudaMemcpyAsync(dPk, hPk, nS * r->bpv, cudaMemcpyHostToDevice, r->st[0]));
    CKR(cudaMemcpyAsync(dLut, hLu, nS * 4 * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
    if (r->K2 > 0)
        CKR(cudaMemcpyAsync(r->dLut2, r->hLut2, nS * 4 * sizeof(double), cudaMemcpyHostToDevice, r->st[0]));
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

    int pass = 0;
    for (int s0 = 0; s0 < t_nSlots; s0 += r->slotsPerPass, ++pass) {
        const int sc = (t_nSlots - s0 < r->slotsPerPass) ? (t_nSlots - s0) : r->slotsPerPass;
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
        if (r->fp64) {
            double* g = (double*)r->dG[b];
            if (r->decodeX2) decode_lut_x2<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
            else if (x4) decode_lut_x4<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
            else    decode_lut_any<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
        } else {
            float* g = (float*)r->dG[b];
            if (x4) decode_lut_x4<float><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
            else    decode_lut_any<float><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu, g);
        }
        CKR(cudaGetLastError());
        CKR(cudaEventRecord(r->ev[b][1], s));

        CBR(cublasSetStream(r->cub, s));
        // TODO(precision:scan): the INT8 GEMMs (slices x cublasGemmEx int8 ->
        // int32, recombined in fp64 into dC1 / dC2) go here, as a third branch
        // beside the fp64 / fp32 ones below and for the second GEMM.
        // C1(sc x K1) = dG^T (sc x N) * dB1 (N x K1), into rows [s0, s0+sc) of
        // the maxSlots x K1 result.
        if (r->fp64) {
            CBR(cublasDgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K1, r->N,
                            &oned, (const double*)r->dG[b], r->N,
                            (const double*)r->dB1, r->N,
                            &zerod, (double*)r->dC1 + s0, r->maxSlots));
        } else {
            CBR(cublasSgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K1, r->N,
                            &onef, (const float*)r->dG[b], r->N,
                            (const float*)r->dB1, r->N,
                            &zerof, (float*)r->dC1 + s0, r->maxSlots));
        }
        CKR(cudaEventRecord(r->ev[b][2], s));

        if (r->K2 > 0) {
            // Same buffer, squared table, second right operand.
            if (r->fp64) {
                double* g = (double*)r->dG[b];
                if (r->decodeX2) decode_lut_x2<<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                else if (x4) decode_lut_x4<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                else    decode_lut_any<double><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                CKR(cudaGetLastError());
                CBR(cublasDgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K2, r->N,
                                &oned, (const double*)r->dG[b], r->N,
                                (const double*)r->dB2, r->N,
                                &zerod, (double*)r->dC2 + s0, r->maxSlots));
            } else {
                float* g = (float*)r->dG[b];
                if (x4) decode_lut_x4<float><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                else    decode_lut_any<float><<<sc, 256, 0, s>>>(pk, r->bpv, r->N, lu2, g);
                CKR(cudaGetLastError());
                CBR(cublasSgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K2, r->N,
                                &onef, (const float*)r->dG[b], r->N,
                                (const float*)r->dB2, r->N,
                                &zerof, (float*)r->dC2 + s0, r->maxSlots));
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
    CKR(cudaMemcpy2DAsync(hC1, (std::size_t)r->maxSlots * r->esz,
                          r->dC1, (std::size_t)r->maxSlots * r->esz,
                          nS * r->esz, (std::size_t)r->K1,
                          cudaMemcpyDeviceToHost, r->st[0]));
    if (r->K2 > 0)
        CKR(cudaMemcpy2DAsync(hC2, (std::size_t)r->maxSlots * r->esz,
                              r->dC2, (std::size_t)r->maxSlots * r->esz,
                              nS * r->esz, (std::size_t)r->K2,
                              cudaMemcpyDeviceToHost, r->st[0]));
    if (timed2) cudaEventRecord(z1, r->st[0]);
    CKR(cudaStreamSynchronize(r->st[0]));
    if (timed2) {
        float ms = 0;
        if (cudaEventElapsedTime(&ms, z0, z1) == cudaSuccess) r->tD2H += ms * 1e-3;
    }
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
