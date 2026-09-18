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
// and imputation costs nothing.
//
// Determinism: for a given (N, K, maxSlots, slotsPerPass, precision) every
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

}  // namespace

// ---------------------------------------------------------------------------

struct Reducer {
    int  N = 0, K = 0, maxSlots = 0;
    bool fp64 = true;
    std::size_t bpv = 0;
    std::size_t esz = 0;           // sizeof(T)
    int slotsPerPass = 0;          // markers per device pass; two are resident

    // pinned host staging
    unsigned char* hPacked = nullptr;
    double*        hLut    = nullptr;
    void*          hC      = nullptr;

    // device
    void*          dB   = nullptr;      // N x K
    void*          dC   = nullptr;      // maxSlots x K
    unsigned char* dPk[2]  = {nullptr, nullptr};
    double*        dLut[2] = {nullptr, nullptr};
    void*          dG[2]   = {nullptr, nullptr};

    cudaStream_t   st[2] = {nullptr, nullptr};
    cudaEvent_t    ev[2][4] = {{nullptr, nullptr, nullptr, nullptr},
                               {nullptr, nullptr, nullptr, nullptr}};
    bool           evLive[2] = {false, false};
    cublasHandle_t cub = nullptr;

    double tH2D = 0, tDec = 0, tGemm = 0, tD2H = 0;
    std::size_t devBytes = 0;
};

namespace {

// Fold buffer `b`'s recorded stage times into the totals. Only called where the
// stream has already been synchronised, so it never adds a stall of its own.
void harvest(Reducer* r, int b)
{
    if (!r->evLive[b]) return;
    float ms = 0;
    if (cudaEventElapsedTime(&ms, r->ev[b][0], r->ev[b][1]) == cudaSuccess) r->tH2D  += ms * 1e-3;
    if (cudaEventElapsedTime(&ms, r->ev[b][1], r->ev[b][2]) == cudaSuccess) r->tDec  += ms * 1e-3;
    if (cudaEventElapsedTime(&ms, r->ev[b][2], r->ev[b][3]) == cudaSuccess) r->tGemm += ms * 1e-3;
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

Reducer* create(int t_device, int t_N, int t_K, const double* t_B, int t_maxSlots,
                bool t_fp64)
{
    if (t_N <= 0 || t_K <= 0 || t_maxSlots <= 0 || t_B == nullptr) return nullptr;
    if (cudaSetDevice(t_device) != cudaSuccess) return nullptr;

    Reducer* r = new Reducer();
    r->N = t_N; r->K = t_K; r->maxSlots = t_maxSlots; r->fp64 = t_fp64;
    r->bpv = (std::size_t)((t_N + 3) / 4);
    r->esz = t_fp64 ? sizeof(double) : sizeof(float);

    // dG is the big one: slotsPerPass * N elements, twice over. Keep each
    // buffer near 512 MB and never wider than the caller's batch.
    const std::size_t perSlot = (std::size_t)t_N * r->esz;
    long long sp = (long long)((512ull << 20) / perSlot);
    if (sp < 1) sp = 1;
    if (sp > t_maxSlots) sp = t_maxSlots;
    if (sp > 4096) sp = 4096;
    r->slotsPerPass = (int)sp;

    auto fail = [&]() -> Reducer* { destroy(r); return nullptr; };

    if (cudaHostAlloc((void**)&r->hPacked, (std::size_t)t_maxSlots * r->bpv,
                      cudaHostAllocDefault) != cudaSuccess) return fail();
    if (cudaHostAlloc((void**)&r->hLut, (std::size_t)t_maxSlots * 4 * sizeof(double),
                      cudaHostAllocDefault) != cudaSuccess) return fail();
    if (cudaHostAlloc(&r->hC, (std::size_t)t_maxSlots * t_K * r->esz,
                      cudaHostAllocDefault) != cudaSuccess) return fail();

    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    if (!dev(&r->dB, (std::size_t)t_N * t_K * r->esz)) return fail();
    if (!dev(&r->dC, (std::size_t)t_maxSlots * t_K * r->esz)) return fail();
    for (int b = 0; b < 2; ++b) {
        if (!dev((void**)&r->dPk[b],  (std::size_t)r->slotsPerPass * r->bpv)) return fail();
        if (!dev((void**)&r->dLut[b], (std::size_t)r->slotsPerPass * 4 * sizeof(double))) return fail();
        if (!dev(&r->dG[b],           (std::size_t)r->slotsPerPass * t_N * r->esz)) return fail();
    }
    r->devBytes = db;

    const std::size_t nB = (std::size_t)t_N * t_K;
    if (t_fp64) {
        if (cudaMemcpy(r->dB, t_B, nB * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
            return fail();
    } else {
        std::vector<float> bf(nB);
        for (std::size_t i = 0; i < nB; ++i) bf[i] = (float)t_B[i];
        if (cudaMemcpy(r->dB, bf.data(), nB * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess)
            return fail();
    }

    for (int b = 0; b < 2; ++b) {
        if (cudaStreamCreate(&r->st[b]) != cudaSuccess) return fail();
        for (int k = 0; k < 4; ++k)
            if (cudaEventCreate(&r->ev[b][k]) != cudaSuccess) return fail();
    }
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
        if (r->dPk[b])  cudaFree(r->dPk[b]);
        if (r->dLut[b]) cudaFree(r->dLut[b]);
        if (r->dG[b])   cudaFree(r->dG[b]);
    }
    if (r->dB) cudaFree(r->dB);
    if (r->dC) cudaFree(r->dC);
    if (r->hPacked) cudaFreeHost(r->hPacked);
    if (r->hLut)    cudaFreeHost(r->hLut);
    if (r->hC)      cudaFreeHost(r->hC);
    delete r;
}

unsigned char* packed(Reducer* r)             { return r ? r->hPacked : nullptr; }
double*        lut(Reducer* r)                { return r ? r->hLut : nullptr; }
std::size_t    bytesPerSlot(const Reducer* r) { return r ? r->bpv : 0; }
bool           isFp64(const Reducer* r)       { return r ? r->fp64 : false; }
const float*   outCf(const Reducer* r)        { return (r && !r->fp64) ? (const float*)r->hC : nullptr; }
const double*  outCd(const Reducer* r)        { return (r &&  r->fp64) ? (const double*)r->hC : nullptr; }
std::size_t    ldC(const Reducer* r)          { return r ? (std::size_t)r->maxSlots : 0; }
std::size_t    deviceBytes(const Reducer* r)  { return r ? r->devBytes : 0; }

void timings(const Reducer* r, double* h2d, double* dec, double* gemm, double* d2h)
{
    if (!r) return;
    if (h2d)  *h2d  = r->tH2D;
    if (dec)  *dec  = r->tDec;
    if (gemm) *gemm = r->tGemm;
    if (d2h)  *d2h  = r->tD2H;
}

bool reduce(Reducer* r, int t_nSlots)
{
    if (!r) return false;
    if (t_nSlots <= 0) return true;
    if (t_nSlots > r->maxSlots) { lastErr = "nSlots > maxSlots"; return false; }

    const float  onef = 1.f, zerof = 0.f;
    const double oned = 1.0, zerod = 0.0;
    const bool x4 = ((r->N & 3) == 0);

    int pass = 0;
    for (int s0 = 0; s0 < t_nSlots; s0 += r->slotsPerPass, ++pass) {
        const int sc = (t_nSlots - s0 < r->slotsPerPass) ? (t_nSlots - s0) : r->slotsPerPass;
        const int b  = pass & 1;
        cudaStream_t s = r->st[b];

        // Buffer b is still in flight from two passes ago; wait, then bank its
        // times (free: the wait had to happen anyway).
        CKR(cudaStreamSynchronize(s));
        harvest(r, b);

        CKR(cudaEventRecord(r->ev[b][0], s));
        CKR(cudaMemcpyAsync(r->dPk[b], r->hPacked + (std::size_t)s0 * r->bpv,
                            (std::size_t)sc * r->bpv, cudaMemcpyHostToDevice, s));
        CKR(cudaMemcpyAsync(r->dLut[b], r->hLut + (std::size_t)s0 * 4,
                            (std::size_t)sc * 4 * sizeof(double),
                            cudaMemcpyHostToDevice, s));
        CKR(cudaEventRecord(r->ev[b][1], s));

        if (r->fp64) {
            double* g = (double*)r->dG[b];
            if (x4) decode_lut_x4<double><<<sc, 256, 0, s>>>(r->dPk[b], r->bpv, r->N, r->dLut[b], g);
            else    decode_lut_any<double><<<sc, 256, 0, s>>>(r->dPk[b], r->bpv, r->N, r->dLut[b], g);
        } else {
            float* g = (float*)r->dG[b];
            if (x4) decode_lut_x4<float><<<sc, 256, 0, s>>>(r->dPk[b], r->bpv, r->N, r->dLut[b], g);
            else    decode_lut_any<float><<<sc, 256, 0, s>>>(r->dPk[b], r->bpv, r->N, r->dLut[b], g);
        }
        CKR(cudaGetLastError());
        CKR(cudaEventRecord(r->ev[b][2], s));

        CBR(cublasSetStream(r->cub, s));
        // C(sc x K) = dG^T (sc x N) * dB (N x K), into rows [s0, s0+sc) of the
        // maxSlots x K result.
        if (r->fp64) {
            CBR(cublasDgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K, r->N,
                            &oned, (const double*)r->dG[b], r->N,
                            (const double*)r->dB, r->N,
                            &zerod, (double*)r->dC + s0, r->maxSlots));
        } else {
            CBR(cublasSgemm(r->cub, CUBLAS_OP_T, CUBLAS_OP_N, sc, r->K, r->N,
                            &onef, (const float*)r->dG[b], r->N,
                            (const float*)r->dB, r->N,
                            &zerof, (float*)r->dC + s0, r->maxSlots));
        }
        CKR(cudaEventRecord(r->ev[b][3], s));
        r->evLive[b] = true;
    }

    CKR(cudaStreamSynchronize(r->st[0]));
    CKR(cudaStreamSynchronize(r->st[1]));
    harvest(r, 0);
    harvest(r, 1);

    cudaEvent_t a = nullptr, z = nullptr;
    const bool timed = (cudaEventCreate(&a) == cudaSuccess) &&
                       (cudaEventCreate(&z) == cudaSuccess);
    if (timed) cudaEventRecord(a, r->st[0]);
    // Only the first t_nSlots rows of each of the K columns are live.
    CKR(cudaMemcpy2DAsync(r->hC, (std::size_t)r->maxSlots * r->esz,
                          r->dC,  (std::size_t)r->maxSlots * r->esz,
                          (std::size_t)t_nSlots * r->esz, (std::size_t)r->K,
                          cudaMemcpyDeviceToHost, r->st[0]));
    if (timed) cudaEventRecord(z, r->st[0]);
    CKR(cudaStreamSynchronize(r->st[0]));
    if (timed) {
        float ms = 0;
        if (cudaEventElapsedTime(&ms, a, z) == cudaSuccess) r->tD2H += ms * 1e-3;
    }
    if (a) cudaEventDestroy(a);
    if (z) cudaEventDestroy(z);
    return true;
}

}  // namespace gpu2
}  // namespace saige
