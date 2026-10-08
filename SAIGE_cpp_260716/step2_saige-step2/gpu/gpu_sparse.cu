// gpu_sparse.cu — see gpu_sparse.hpp.
#include "gpu_sparse.hpp"
#include "gpu_step2.hpp"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstring>
#include <vector>
#include <string>

namespace saige {
namespace gpu2 {

namespace {

constexpr int NT = 256;          // threads per block (one block per slot)
constexpr int TR = 8;            // traits accumulated per pass over the pairs
constexpr std::size_t SHMAX = 40960;   // stage the slot's column in shared memory up to this

thread_local std::string lastErrSpq;
#define CKQ(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    lastErrSpq = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)

__device__ __forceinline__ double warpSum(double v)
{
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

__device__ __forceinline__ double decode(const unsigned char* col, const double* L, int i)
{
    return L[(col[i >> 2] >> ((i & 3) * 2)) & 3u];
}

// One block per slot. Each thread owns pairs k = tid, tid + NT, ...; for every
// pass of TR traits it accumulates w[t][k] * g_i g_j in registers, then the
// block reduces each trait (warp shuffle, then the warps' partials in order).
__global__ void __launch_bounds__(NT)
spq_kernel(const unsigned char* __restrict__ pk, std::size_t bpv,
           const double* __restrict__ lut,
           long long nPairs, const int* __restrict__ pi, const int* __restrict__ pj,
           const double* __restrict__ w, int nTr, int useShared,
           double* __restrict__ out)
{
    extern __shared__ unsigned char shCol[];
    __shared__ double L[4];
    __shared__ double part[TR][NT / 32];
    const int s = blockIdx.x;
    const unsigned char* gcol = pk + (std::size_t)s * bpv;
    if (threadIdx.x < 4) L[threadIdx.x] = lut[(std::size_t)s * 4 + threadIdx.x];
    const unsigned char* col = gcol;
    if (useShared) {
        for (std::size_t b = threadIdx.x; b < bpv; b += NT) shCol[b] = gcol[b];
        col = shCol;
    }
    __syncthreads();
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int t0 = 0; t0 < nTr; t0 += TR) {
        const int nt = (nTr - t0 < TR) ? (nTr - t0) : TR;
        double acc[TR];
        #pragma unroll
        for (int k = 0; k < TR; k++) acc[k] = 0.0;
        for (long long p = threadIdx.x; p < nPairs; p += NT) {
            const double gi = decode(col, L, pi[p]);
            if (gi == 0.0) continue;
            const double gg = gi * decode(col, L, pj[p]);
            if (gg == 0.0) continue;
            const double* wp = w + (std::size_t)t0 * nPairs + p;
            #pragma unroll
            for (int k = 0; k < TR; k++)
                if (k < nt) acc[k] += gg * wp[(std::size_t)k * nPairs];
        }
        #pragma unroll
        for (int k = 0; k < TR; k++) {
            const double v = warpSum(acc[k]);
            if (lane == 0) part[k][wid] = v;
        }
        __syncthreads();
        if (threadIdx.x < nt) {
            double v = 0.0;
            for (int q = 0; q < NT / 32; q++) v += part[threadIdx.x][q];
            out[(std::size_t)s * nTr + t0 + threadIdx.x] = v;
        }
        __syncthreads();
    }
}

// gpuOwnSampleSets: as spq_kernel, but every trait decodes the slot's codes
// through its own 4-entry table (tl, nSlots x nTr x 4). The codes of a pair
// are read once; a pair is skipped when no trait of the pass gives it a
// nonzero product. Per trait the sum is (L_t[ci] * L_t[cj]) * w, in the same
// (thread, pair) order and tree as spq_kernel.
__global__ void __launch_bounds__(NT)
spq_kernel_own(const unsigned char* __restrict__ pk, std::size_t bpv,
               const double* __restrict__ tl,
               long long nPairs, const int* __restrict__ pi, const int* __restrict__ pj,
               const double* __restrict__ w, int nTr, int useShared,
               double* __restrict__ out)
{
    extern __shared__ unsigned char shCol[];
    __shared__ double L[TR][4];
    __shared__ unsigned char nz[16];
    __shared__ double part[TR][NT / 32];
    const int s = blockIdx.x;
    const unsigned char* gcol = pk + (std::size_t)s * bpv;
    const unsigned char* col = gcol;
    if (useShared) {
        for (std::size_t b = threadIdx.x; b < bpv; b += NT) shCol[b] = gcol[b];
        col = shCol;
    }
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int t0 = 0; t0 < nTr; t0 += TR) {
        const int nt = (nTr - t0 < TR) ? (nTr - t0) : TR;
        if (threadIdx.x < TR * 4) {
            const int k = threadIdx.x >> 2, c = threadIdx.x & 3;
            L[k][c] = (k < nt) ? tl[((std::size_t)s * nTr + t0 + k) * 4 + c] : 0.0;
        }
        __syncthreads();
        if (threadIdx.x < 16) {
            const int ci = threadIdx.x >> 2, cj = threadIdx.x & 3;
            unsigned char any = 0;
            for (int k = 0; k < nt; k++) if (L[k][ci] * L[k][cj] != 0.0) any = 1;
            nz[threadIdx.x] = any;
        }
        __syncthreads();
        double acc[TR];
        #pragma unroll
        for (int k = 0; k < TR; k++) acc[k] = 0.0;
        for (long long p = threadIdx.x; p < nPairs; p += NT) {
            const int i = pi[p], j = pj[p];
            const unsigned ci = (col[i >> 2] >> ((i & 3) * 2)) & 3u;
            const unsigned cj = (col[j >> 2] >> ((j & 3) * 2)) & 3u;
            if (!nz[ci * 4 + cj]) continue;
            const double* wp = w + (std::size_t)t0 * nPairs + p;
            #pragma unroll
            for (int k = 0; k < TR; k++)
                if (k < nt) acc[k] += (L[k][ci] * L[k][cj]) * wp[(std::size_t)k * nPairs];
        }
        #pragma unroll
        for (int k = 0; k < TR; k++) {
            const double v = warpSum(acc[k]);
            if (lane == 0) part[k][wid] = v;
        }
        __syncthreads();
        if (threadIdx.x < nt) {
            double v = 0.0;
            for (int q = 0; q < NT / 32; q++) v += part[threadIdx.x][q];
            out[(std::size_t)s * nTr + t0 + threadIdx.x] = v;
        }
        __syncthreads();
    }
}

}  // namespace

struct SpQuad {
    int N = 0, nTr = 0, maxSlots = 0;
    long long nPairs = 0;
    int* dI = nullptr; int* dJ = nullptr; double* dW = nullptr;
    double* dOut = nullptr;
    std::vector<double*> hOut;             // one pinned result set per outSets
    int lastOut = 0;
    double* dTl = nullptr;                 // spqRunOwn: nSlots x nTr x 4, allocated on first use
    cudaStream_t st = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    double tKernel = 0.0;
    long long nSlotsDone = 0;
    std::size_t devBytes = 0;
};

bool spqSupports(Prec t_p)
{
    switch (t_p) {
        case Prec::FP64: return true;
        // TODO(precision:scan): the sparse-GRM cross terms in the scan's fp32 /
        // int8 mode. Plug the variant in at spq_kernel / spq_kernel_own's
        // launch in spqRun / spqRunOwn (results stay double) and return true.
        case Prec::FP32: return false;
        case Prec::INT8: return false;
    }
    return false;
}

SpQuad* spqCreate(const SpQuadCreateArgs& a)
{
    if (a.N <= 0 || a.nTraits <= 0 || a.maxSlots <= 0 || a.nPairs < 0) return nullptr;
    if (a.nPairs > 0 && (!a.pi || !a.pj || !a.w)) return nullptr;
    // ---- precision dispatch (the one place the mode is decided) ----
    if (!spqSupports(a.precision)) {
        lastErrSpq = std::string("scan precision ") + precName(a.precision) +
                     " is not implemented yet (sparse-GRM cross terms)";
        return nullptr;
    }
    if (cudaSetDevice(a.device) != cudaSuccess) return nullptr;
    SpQuad* q = new SpQuad();
    q->N = a.N; q->nTr = a.nTraits; q->maxSlots = a.maxSlots; q->nPairs = a.nPairs;
    auto fail = [&]() -> SpQuad* { spqDestroy(q); return nullptr; };
    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (n == 0) n = 8;
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    const std::size_t np = (std::size_t)a.nPairs;
    if (!dev((void**)&q->dI, np * sizeof(int))) return fail();
    if (!dev((void**)&q->dJ, np * sizeof(int))) return fail();
    if (!dev((void**)&q->dW, np * a.nTraits * sizeof(double))) return fail();
    if (!dev((void**)&q->dOut, (std::size_t)a.maxSlots * a.nTraits * sizeof(double))) return fail();
    q->hOut.assign((std::size_t)(a.outSets < 1 ? 1 : a.outSets), nullptr);
    for (double*& h : q->hOut)
        if (cudaHostAlloc((void**)&h, (std::size_t)a.maxSlots * a.nTraits * sizeof(double),
                          cudaHostAllocDefault) != cudaSuccess) return fail();
    q->devBytes = db;
    if (np > 0) {
        if (cudaMemcpy(q->dI, a.pi, np * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(q->dJ, a.pj, np * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(q->dW, a.w, np * a.nTraits * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    }
    if (cudaStreamCreate(&q->st) != cudaSuccess) return fail();
    if (cudaEventCreate(&q->e0) != cudaSuccess || cudaEventCreate(&q->e1) != cudaSuccess) return fail();
    return q;
}

void spqDestroy(SpQuad* q)
{
    if (!q) return;
    if (q->e0) cudaEventDestroy(q->e0);
    if (q->e1) cudaEventDestroy(q->e1);
    if (q->st) cudaStreamDestroy(q->st);
    if (q->dI) cudaFree(q->dI);
    if (q->dJ) cudaFree(q->dJ);
    if (q->dW) cudaFree(q->dW);
    if (q->dOut) cudaFree(q->dOut);
    if (q->dTl) cudaFree(q->dTl);
    for (double* h : q->hOut) if (h) cudaFreeHost(h);
    delete q;
}

const double* spqOut(const SpQuad* q, int t_set)
{
    if (!q) return nullptr;
    const int k = t_set < 0 ? q->lastOut : t_set;
    return k < (int)q->hOut.size() ? q->hOut[(std::size_t)k] : nullptr;
}
std::size_t spqDeviceBytes(const SpQuad* q) { return q ? q->devBytes : 0; }
const char* spqLastError() { return lastErrSpq.c_str(); }
void spqTimings(const SpQuad* q, double* t_kernel, long long* t_slots)
{
    if (!q) return;
    if (t_kernel) *t_kernel = q->tKernel;
    if (t_slots)  *t_slots = q->nSlotsDone;
}

bool spqRun(SpQuad* q, const Reducer* r, int nSlots, int t_set)
{
    if (!q || !r) return false;
    if (t_set >= (int)q->hOut.size()) { lastErrSpq = "result set out of range"; return false; }
    double* const hOut = q->hOut[(std::size_t)(t_set < 0 ? 0 : t_set)];
    q->lastOut = t_set < 0 ? 0 : t_set;
    if (nSlots <= 0) return true;
    if (nSlots > q->maxSlots) { lastErrSpq = "nSlots > maxSlots"; return false; }
    const std::size_t nOut = (std::size_t)nSlots * q->nTr;
    if (q->nPairs == 0) {
        std::memset(hOut, 0, nOut * sizeof(double));
        return true;
    }
    const unsigned char* dPk = (const unsigned char*)devicePacked(r, t_set);
    const double* dLut = (const double*)deviceLut(r, t_set);
    if (!dPk || !dLut) { lastErrSpq = "reducer has no resident packed rows"; return false; }
    const std::size_t bpv = bytesPerSlot(r);
    const int useShared = (bpv <= SHMAX) ? 1 : 0;
    // The reducer's last reduce() is complete (it synchronises before
    // returning), so its resident rows are safe to read on our own stream.
    CKQ(cudaEventRecord(q->e0, q->st));
    spq_kernel<<<nSlots, NT, useShared ? bpv : 0, q->st>>>(dPk, bpv, dLut, q->nPairs, q->dI, q->dJ, q->dW,
                                                          q->nTr, useShared, q->dOut);
    CKQ(cudaGetLastError());
    CKQ(cudaEventRecord(q->e1, q->st));
    CKQ(cudaMemcpyAsync(hOut, q->dOut, nOut * sizeof(double), cudaMemcpyDeviceToHost, q->st));
    CKQ(cudaStreamSynchronize(q->st));
    float ms = 0;
    if (cudaEventElapsedTime(&ms, q->e0, q->e1) == cudaSuccess) q->tKernel += ms * 1e-3;
    q->nSlotsDone += nSlots;
    return true;
}

bool spqRunOwn(SpQuad* q, const Reducer* r, int nSlots, const double* tlut, int t_set)
{
    if (!q || !r || !tlut) return false;
    if (t_set >= (int)q->hOut.size()) { lastErrSpq = "result set out of range"; return false; }
    double* const hOut = q->hOut[(std::size_t)(t_set < 0 ? 0 : t_set)];
    q->lastOut = t_set < 0 ? 0 : t_set;
    if (nSlots <= 0) return true;
    if (nSlots > q->maxSlots) { lastErrSpq = "nSlots > maxSlots"; return false; }
    const std::size_t nOut = (std::size_t)nSlots * q->nTr;
    if (q->nPairs == 0) {
        std::memset(hOut, 0, nOut * sizeof(double));
        return true;
    }
    const unsigned char* dPk = (const unsigned char*)devicePacked(r, t_set);
    if (!dPk) { lastErrSpq = "reducer has no resident packed rows"; return false; }
    if (!q->dTl) {
        const std::size_t nb = (std::size_t)q->maxSlots * q->nTr * 4 * sizeof(double);
        CKQ(cudaMalloc((void**)&q->dTl, nb));
        q->devBytes += nb;
    }
    const std::size_t bpv = bytesPerSlot(r);
    const int useShared = (bpv <= SHMAX) ? 1 : 0;
    CKQ(cudaMemcpyAsync(q->dTl, tlut, nOut * 4 * sizeof(double), cudaMemcpyHostToDevice, q->st));
    CKQ(cudaEventRecord(q->e0, q->st));
    spq_kernel_own<<<nSlots, NT, useShared ? bpv : 0, q->st>>>(dPk, bpv, q->dTl, q->nPairs, q->dI, q->dJ,
                                                              q->dW, q->nTr, useShared, q->dOut);
    CKQ(cudaGetLastError());
    CKQ(cudaEventRecord(q->e1, q->st));
    CKQ(cudaMemcpyAsync(hOut, q->dOut, nOut * sizeof(double), cudaMemcpyDeviceToHost, q->st));
    CKQ(cudaStreamSynchronize(q->st));
    float ms = 0;
    if (cudaEventElapsedTime(&ms, q->e0, q->e1) == cudaSuccess) q->tKernel += ms * 1e-3;
    q->nSlotsDone += nSlots;
    return true;
}

}  // namespace gpu2
}  // namespace saige
