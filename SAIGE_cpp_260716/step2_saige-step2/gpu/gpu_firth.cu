// gpu_firth.cu — one CUDA block per Firth-flagged (marker, trait) pair; the
// block's threads share every length-N sum, so fast_logistf_fit_simple's
// Newton loop runs unchanged per pair. Contract and provenance: gpu_firth.hpp.
// Reference: saige_test.cpp (SAIGEClass::fast_logistf_fit_simple, getadjGFast)
// and the validated prototype gpu/firth_proto/firth_gpu.py on branch firth-gpu
// (FIRTH_GPU.md). Compiled with --fmad=false so a*b+c rounds twice, as the CPU
// build (-std=c++17, no contraction) does.

#include "gpu_firth.hpp"
#include "gpu_step2.hpp"

#include <cuda_runtime.h>
#include <math_constants.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace saige {
namespace gpu2 {

namespace {

#define NT 256
#define NACC 9

thread_local std::string lastErrFirth;
#define CKF(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    lastErrFirth = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)

__device__ __forceinline__ double warpSum(double v)
{
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}
// Block-wide sum; every thread gets the result. sh needs NT/32 + 1 doubles.
__device__ double blockSum(double v, double* sh)
{
    v = warpSum(v);
    if ((threadIdx.x & 31) == 0) sh[threadIdx.x >> 5] = v;
    __syncthreads();
    double r = 0.0;
    if (threadIdx.x < 32) {
        r = (threadIdx.x < NT / 32) ? sh[threadIdx.x] : 0.0;
        r = warpSum(r);
    }
    if (threadIdx.x == 0) sh[NT / 32] = r;
    __syncthreads();
    r = sh[NT / 32];
    __syncthreads();
    return r;
}
// The nine accumulators of one Newton step summed across the block in a fixed
// shape (warp shuffle tree, then thread 0 adds the 8 warp partials in order),
// so a pair's fit is bit-reproducible run to run. sh needs NACC*(NT/32)+NACC.
__device__ void blockSum9(double* acc, double* sh)
{
    const int lane = threadIdx.x & 31, wid = threadIdx.x >> 5;
    for (int q = 0; q < NACC; q++) {
        const double v = warpSum(acc[q]);
        if (lane == 0) sh[q * (NT / 32) + wid] = v;
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        for (int q = 0; q < NACC; q++) {
            double s = 0.0;
            for (int w = 0; w < NT / 32; w++) s += sh[q * (NT / 32) + w];
            sh[NACC * (NT / 32) + q] = s;
        }
    }
    __syncthreads();
    for (int q = 0; q < NACC; q++) acc[q] = sh[NACC * (NT / 32) + q];
    __syncthreads();
}

__device__ __forceinline__ double dose(const unsigned char* col, const double4& L, int i)
{
    const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
    return (c == 0) ? L.x : (c == 1 ? L.y : (c == 2 ? L.z : L.w));
}

// armadillo's inv_sympd on the 2x2 Fisher matrix [[s0, s1], [s1, s2]]: the
// tiny-matrix closed form when it applies, else LAPACK's Cholesky, which fails
// only on a non-positive pivot. Returns false when the CPU's inv_sympd would.
__device__ __forceinline__ bool inv2x2(double s0, double s1, double s2, double& A, double& B, double& C)
{
    const double eps = 2.220446049250313e-16;
    const double det = s0 * s2 - s1 * s1;
    const bool closed = (s0 > 0.0) && (s2 > 0.0) && (det >= eps) && (det <= 1.0 / eps);
    if (!closed) {
        // dpotrf on the 2x2: l11 = sqrt(s0) needs s0 > 0; the second pivot
        // s2 - (s1/l11)^2 must be > 0. NaN anywhere fails both tests.
        if (!(s0 > 0.0)) return false;
        const double l21 = s1 / sqrt(s0);
        if (!((s2 - l21 * l21) > 0.0)) return false;
        if (!(det > 0.0) || !isfinite(det)) return false;
    }
    A = s2 / det; B = -s1 / det; C = s0 / det;
    return true;
}

// OWN (FirthCreateArgs::ownSamples): pair k's dosage table at plut + 4k, and
// a sample outside the trait's mask has dosage 0 and is left out of every
// Newton-step sum (its y / offset are embedded as zeros but never read there).
template <bool OWN>
__global__ void __launch_bounds__(NT)
firth_pairs(const unsigned char* __restrict__ packed, std::size_t bpv, const double* __restrict__ lut, int N,
            const double* __restrict__ plut, const uint64_t* __restrict__ masks, int maskWords,
            const double* __restrict__ Y, const double* __restrict__ OFF,
            const double* __restrict__ XV, const double* __restrict__ XX,
            const int* __restrict__ pOfTrait, std::size_t traitStride,
            const FirthPairIn* __restrict__ in, int nPairs, double* __restrict__ scratch,
            int maxit, double maxstep, double xconv, double gconv,
            FirthPairOut* __restrict__ out)
{
    __shared__ double sh[NACC * (NT / 32) + NACC];
    __shared__ double bsh[FIRTH_PMAX];
    double* gt = scratch + (std::size_t)blockIdx.x * N;
    for (int k = blockIdx.x; k < nPairs; k += gridDim.x) {
        const FirthPairIn pin = in[k];
        const unsigned char* col = packed + (std::size_t)pin.slot * bpv;
        double4 L;
        {
            const double2* d = reinterpret_cast<const double2*>(OWN ? plut : lut) +
                               2 * (std::size_t)(OWN ? k : pin.slot);
            const double2 a = d[0], b = d[1];
            L = make_double4(a.x, a.y, b.x, b.y);
        }
        const int t = pin.trait;
        const uint64_t* mk = OWN ? masks + (std::size_t)t * maskWords : nullptr;
        auto in = [&](int i) -> bool { return !OWN || ((mk[i >> 6] >> (i & 63)) & 1ull); };
        const int p = pOfTrait[t];
        const double* y   = Y   + (std::size_t)t * N;
        const double* off = OFF + (std::size_t)t * N;
        const double* xv  = XV  + (std::size_t)t * traitStride;   // sample i: xv[i*p + j]
        const double* xx  = XX  + (std::size_t)t * traitStride;   // xx[j*N + i]

        // pass A: b = XV g over carriers (getadjGFast's loop over iIndex)
        double a[FIRTH_PMAX];
        for (int j = 0; j < FIRTH_PMAX; ++j) a[j] = 0.0;
        for (int i = threadIdx.x; i < N; i += NT) {
            const double g = in(i) ? dose(col, L, i) : 0.0;
            if (g == 0.0) continue;
            const double* x = xv + (std::size_t)i * p;
            for (int j = 0; j < p; ++j) a[j] += x[j] * g;
        }
        for (int j = 0; j < p; ++j) {
            const double s = blockSum(a[j], sh);
            if (threadIdx.x == 0) bsh[j] = s;
        }
        __syncthreads();
        // pass B: g~ = g - XXVX_inv b
        for (int i = threadIdx.x; i < N; i += NT) {
            const double g = in(i) ? dose(col, L, i) : 0.0;
            double proj = 0.0;
            for (int j = 0; j < p; ++j) proj += xx[(std::size_t)j * N + i] * bsh[j];
            gt[i] = g - proj;
        }
        __syncthreads();

        // ---- fast_logistf_fit_simple(x = [1, g~], y, offset, firth, init = 0) ----
        double alpha = 0.0, beta = 0.0, cov11 = 0.0;   // XX_covs starts as zeros
        int iter = 0, flag = 0, strict = 0, singular = 0;
        while (iter <= maxit) {
            double acc[NACC];
            for (int q = 0; q < NACC; q++) acc[q] = 0.0;
            for (int i = threadIdx.x; i < N; i += NT) {
                if (!in(i)) continue;
                const double g = gt[i];
                // legacy: pi = 1 / (exp(-x*beta - offset) + 1)
                const double pi = 1.0 / (exp(-(alpha + beta * g) - off[i]) + 1.0);
                const double w = pi * (1.0 - pi);
                const double W2 = sqrt(w), c = g * W2;        // XW2 columns
                const double r = y[i] - pi;
                const double ad = w * (0.5 - pi);              // Firth term weight
                acc[0] += W2 * W2;  acc[1] += W2 * c;  acc[2] += c * c;   // Fisher = XW2' XW2
                acc[3] += r;        acc[4] += g * r;                      // X' (y - pi)
                acc[5] += ad;       acc[6] += ad * g;  acc[7] += ad * g * g;  acc[8] += ad * g * g * g;
            }
            blockSum9(acc, sh);
            double A, Bc, C;
            if (!inv2x2(acc[0], acc[1], acc[2], A, Bc, C)) { singular = 1; break; }
            // U* = X' ((y - pi) + h (0.5 - pi)),  h_i = w_i (A + 2 Bc g_i + C g_i^2)
            const double U0 = acc[3] + (A * acc[5] + 2.0 * Bc * acc[6] + C * acc[7]);
            const double U1 = acc[4] + (A * acc[6] + 2.0 * Bc * acc[7] + C * acc[8]);
            double d0 = A * U0 + Bc * U1, d1 = Bc * U0 + C * U1;      // delta = XX_covs * U*
            double mxd = fmax(fabs(d0), fabs(d1));
            const double mx = mxd / maxstep;
            if (mx > 1.0) { d0 /= mx; d1 /= mx; mxd = fmax(fabs(d0), fabs(d1)); }
            cov11 = C;
            iter++;
            alpha += d0; beta += d1;
            const bool small = (mxd <= xconv) && (fabs(U0) <= gconv) && (fabs(U1) <= gconv);
            if (iter == maxit || small) { flag = 1; strict = small ? 1 : 0; break; }
        }
        if (threadIdx.x == 0) {
            FirthPairOut o;
            o.beta = beta; o.alpha = alpha; o.se = sqrt(cov11);
            o.conv = flag; o.strict = strict; o.niter = iter; o.singular = singular;
            out[k] = o;
        }
        __syncthreads();
    }
}

}  // namespace

struct Firth {
    int N = 0, nTraits = 0, maxPairs = 0, blocks = 256, maxit = 50;
    double maxstep = 15.0, xconv = 1e-5, gconv = 1e-5;
    std::size_t traitStride = 0;
    FirthPairIn*  hIn  = nullptr;
    FirthPairOut* hOut = nullptr;
    FirthPairIn*  dIn  = nullptr;
    FirthPairOut* dOut = nullptr;
    double* dY = nullptr; double* dOff = nullptr; double* dXV = nullptr; double* dXX = nullptr;
    int*    dP  = nullptr;
    double* dScratch = nullptr;
    int own = 0, maskWords = 0;
    double* hPLut = nullptr; double* dPLut = nullptr; uint64_t* dMask = nullptr;
    cudaStream_t st = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    double tKernel = 0.0;
    long long nPairsDone = 0;
    std::size_t devBytes = 0;
};

Firth* firthCreate(const FirthCreateArgs& a)
{
    if (a.N <= 0 || a.nTraits <= 0 || a.traits == nullptr || a.maxPairs <= 0) return nullptr;
    if (a.maxit <= 0 || !(a.maxstep > 0.0)) return nullptr;
    if (cudaSetDevice(a.device) != cudaSuccess) return nullptr;
    int pMax = 0;
    for (int t = 0; t < a.nTraits; ++t) {
        if (a.traits[t].p <= 0 || a.traits[t].p > FIRTH_PMAX) return nullptr;
        if (!a.traits[t].y || !a.traits[t].offset || !a.traits[t].XV || !a.traits[t].XXVX_inv) return nullptr;
        if (a.traits[t].p > pMax) pMax = a.traits[t].p;
    }
    Firth* s = new Firth();
    s->N = a.N; s->nTraits = a.nTraits; s->maxPairs = a.maxPairs;
    s->blocks = a.blocks > 0 ? a.blocks : 256;
    s->maxit = a.maxit; s->maxstep = a.maxstep; s->xconv = a.xconv; s->gconv = a.gconv;
    s->traitStride = (std::size_t)a.N * pMax;
    auto fail = [&]() -> Firth* { firthDestroy(s); return nullptr; };
    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    if (cudaHostAlloc((void**)&s->hIn,  (std::size_t)a.maxPairs * sizeof(FirthPairIn),  cudaHostAllocDefault) != cudaSuccess) return fail();
    if (cudaHostAlloc((void**)&s->hOut, (std::size_t)a.maxPairs * sizeof(FirthPairOut), cudaHostAllocDefault) != cudaSuccess) return fail();
    if (!dev((void**)&s->dIn,  (std::size_t)a.maxPairs * sizeof(FirthPairIn)))  return fail();
    if (!dev((void**)&s->dOut, (std::size_t)a.maxPairs * sizeof(FirthPairOut))) return fail();
    if (!dev((void**)&s->dY,   (std::size_t)a.nTraits * a.N * sizeof(double))) return fail();
    if (!dev((void**)&s->dOff, (std::size_t)a.nTraits * a.N * sizeof(double))) return fail();
    if (!dev((void**)&s->dXV, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail();
    if (!dev((void**)&s->dXX, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail();
    if (!dev((void**)&s->dP,  (std::size_t)a.nTraits * sizeof(int))) return fail();
    if (!dev((void**)&s->dScratch, (std::size_t)s->blocks * a.N * sizeof(double))) return fail();
    if (a.ownSamples) {
        s->own = 1;
        s->maskWords = (a.N + 63) / 64;
        if (cudaHostAlloc((void**)&s->hPLut, (std::size_t)a.maxPairs * 4 * sizeof(double), cudaHostAllocDefault) != cudaSuccess) return fail();
        if (!dev((void**)&s->dPLut, (std::size_t)a.maxPairs * 4 * sizeof(double))) return fail();
        if (!dev((void**)&s->dMask, (std::size_t)a.nTraits * s->maskWords * sizeof(uint64_t))) return fail();
        std::vector<uint64_t> all((std::size_t)s->maskWords, ~0ull);
        for (int t = 0; t < a.nTraits; ++t) {
            const uint64_t* m = a.traits[t].mask ? a.traits[t].mask : all.data();
            if (cudaMemcpy(s->dMask + (std::size_t)t * s->maskWords, m, (std::size_t)s->maskWords * sizeof(uint64_t),
                           cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        }
    }
    s->devBytes = db;
    std::vector<int> pv(a.nTraits);
    for (int t = 0; t < a.nTraits; ++t) {
        const FirthTraitArgs& T = a.traits[t];
        pv[t] = T.p;
        if (cudaMemcpy(s->dY   + (std::size_t)t * a.N, T.y,      (std::size_t)a.N * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(s->dOff + (std::size_t)t * a.N, T.offset, (std::size_t)a.N * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(s->dXV + (std::size_t)t * s->traitStride, T.XV, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(s->dXX + (std::size_t)t * s->traitStride, T.XXVX_inv, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    }
    if (cudaMemcpy(s->dP, pv.data(), pv.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    if (cudaStreamCreate(&s->st) != cudaSuccess) return fail();
    if (cudaEventCreate(&s->e0) != cudaSuccess || cudaEventCreate(&s->e1) != cudaSuccess) return fail();
    return s;
}

void firthDestroy(Firth* s)
{
    if (!s) return;
    if (s->e0) cudaEventDestroy(s->e0);
    if (s->e1) cudaEventDestroy(s->e1);
    if (s->st) cudaStreamDestroy(s->st);
    if (s->dIn) cudaFree(s->dIn);
    if (s->dOut) cudaFree(s->dOut);
    if (s->dY) cudaFree(s->dY);
    if (s->dOff) cudaFree(s->dOff);
    if (s->dXV) cudaFree(s->dXV);
    if (s->dXX) cudaFree(s->dXX);
    if (s->dP) cudaFree(s->dP);
    if (s->dScratch) cudaFree(s->dScratch);
    if (s->dPLut) cudaFree(s->dPLut);
    if (s->dMask) cudaFree(s->dMask);
    if (s->hPLut) cudaFreeHost(s->hPLut);
    if (s->hIn) cudaFreeHost(s->hIn);
    if (s->hOut) cudaFreeHost(s->hOut);
    delete s;
}

FirthPairIn*  firthIn(Firth* s)  { return s ? s->hIn : nullptr; }
FirthPairOut* firthOut(Firth* s) { return s ? s->hOut : nullptr; }
double*       firthPairLut(Firth* s) { return (s && s->own) ? s->hPLut : nullptr; }
std::size_t firthDeviceBytes(const Firth* s) { return s ? s->devBytes : 0; }
const char* firthLastError() { return lastErrFirth.c_str(); }
void firthTimings(const Firth* s, double* t_kernel, long long* t_pairs)
{
    if (!s) return;
    if (t_kernel) *t_kernel = s->tKernel;
    if (t_pairs)  *t_pairs = s->nPairsDone;
}

bool firthRun(Firth* s, const Reducer* r, int nPairs)
{
    if (!s || !r) return false;
    if (nPairs <= 0) return true;
    if (nPairs > s->maxPairs) { lastErrFirth = "nPairs > maxPairs"; return false; }
    const unsigned char* dPk = (const unsigned char*)devicePacked(r);
    const double* dLut = (const double*)deviceLut(r);
    if (!dPk || !dLut) { lastErrFirth = "reducer has no resident packed rows"; return false; }
    // The reducer's last reduce() is complete (it synchronises before
    // returning), so its resident rows are safe to read on our own stream.
    CKF(cudaMemcpyAsync(s->dIn, s->hIn, (std::size_t)nPairs * sizeof(FirthPairIn), cudaMemcpyHostToDevice, s->st));
    if (s->own) CKF(cudaMemcpyAsync(s->dPLut, s->hPLut, (std::size_t)nPairs * 4 * sizeof(double), cudaMemcpyHostToDevice, s->st));
    const int grid = nPairs < s->blocks ? nPairs : s->blocks;
    CKF(cudaEventRecord(s->e0, s->st));
    if (s->own)
        firth_pairs<true><<<grid, NT, 0, s->st>>>(dPk, bytesPerSlot(r), dLut, s->N, s->dPLut, s->dMask, s->maskWords,
                                            s->dY, s->dOff, s->dXV, s->dXX, s->dP,
                                            s->traitStride, s->dIn, nPairs, s->dScratch,
                                            s->maxit, s->maxstep, s->xconv, s->gconv, s->dOut);
    else
        firth_pairs<false><<<grid, NT, 0, s->st>>>(dPk, bytesPerSlot(r), dLut, s->N, nullptr, nullptr, 0,
                                             s->dY, s->dOff, s->dXV, s->dXX, s->dP,
                                             s->traitStride, s->dIn, nPairs, s->dScratch,
                                             s->maxit, s->maxstep, s->xconv, s->gconv, s->dOut);
    CKF(cudaGetLastError());
    CKF(cudaEventRecord(s->e1, s->st));
    CKF(cudaMemcpyAsync(s->hOut, s->dOut, (std::size_t)nPairs * sizeof(FirthPairOut), cudaMemcpyDeviceToHost, s->st));
    CKF(cudaStreamSynchronize(s->st));
    float ms = 0;
    if (cudaEventElapsedTime(&ms, s->e0, s->e1) == cudaSuccess) s->tKernel += ms * 1e-3;
    s->nPairsDone += nPairs;
    return true;
}

}  // namespace gpu2
}  // namespace saige
