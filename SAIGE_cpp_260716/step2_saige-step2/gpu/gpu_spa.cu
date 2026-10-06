// gpu_spa.cu — one CUDA block per flagged (marker, trait) pair; the block's
// threads share every length-N sum, so SAIGE's Newton loop runs unchanged per
// pair. Contract and provenance: gpu_spa.hpp. Reference: spa_binary.cpp
// (getroot_K1_Binom / getroot_K1_fast_Binom / Get_Saddle_Prob_*_Binom) and
// spa.cpp (SPA / SPA_fast). Compiled with --fmad=false so a*b+c rounds twice,
// as the CPU build (-std=c++17, no contraction) does.

#include "gpu_spa.hpp"
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

thread_local std::string lastErrSpa;
#define CKS(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    lastErrSpa = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)

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

__device__ __forceinline__ int sgn(double x) { return (x > 0) - (x < 0); }   // arma::sign

__device__ __forceinline__ double dose(const unsigned char* col, const double4& L, int i)
{
    const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
    return (c == 0) ? L.x : (c == 1 ? L.y : (c == 2 ? L.z : L.w));
}

struct PairCtx {
    const unsigned char* col;
    double4 L;
    const double* mu;
    const double* gt;      // g~ for this block
    int N;
    int fast;
    double NAmu, NAsigma;
};

// mode 0: K1 sum  sum mu g~ / ((1-mu) e^{-g~ t} + mu)           (+ NAmu + NAsigma t)
// mode 1: K2 sum  sum (1-mu) mu g~^2 e^{-g~ t} / (...)^2        (+ NAsigma); non-finite terms skipped (sum_arma1)
// mode 2: Korg    sum log(1 - mu + mu e^{g~ t})                 (+ NAmu t + NAsigma t^2 / 2)
// fast: carriers only (g != 0), exactly the index set iIndex the CPU gathers.
__device__ double kpass(int mode, double t, const PairCtx& P, double* sh)
{
    double acc = 0.0;
    for (int i = threadIdx.x; i < P.N; i += NT) {
        if (P.fast && dose(P.col, P.L, i) == 0.0) continue;
        const double g = P.gt[i], m = P.mu[i];
        if (mode == 0) {
            acc += m * g / ((1.0 - m) * exp(-g * t) + m);
        } else if (mode == 1) {
            const double e = exp(-g * t);
            const double d = (1.0 - m) * e + m;
            const double term = (1.0 - m) * m * (g * g * e) / (d * d);
            if (isfinite(term)) acc += term;
        } else {
            acc += log(1.0 - m + m * exp(g * t));
        }
    }
    double s = blockSum(acc, sh);
    if (P.fast) {
        if (mode == 0)      s += P.NAmu + P.NAsigma * t;
        else if (mode == 1) s += P.NAsigma;
        else                s += P.NAmu * t + 0.5 * P.NAsigma * t * t;
    }
    return s;
}

// getroot_K1_Binom (fast == 0) / getroot_K1_fast_Binom (fast == 1), from init 0.
__device__ void getroot(const PairCtx& P, double q, double gpos, double gneg, double tol, int maxiter,
                        double* sh, double* root, int* niter, int* conv)
{
    if (q >= gpos || q <= gneg) { *root = CUDART_INF; *niter = 0; *conv = 1; return; }
    double t = 0.0;
    double K1 = kpass(0, t, P, sh) - q;
    double prevJump = CUDART_INF;
    int rep = 1, c = 1;
    double tnew = t, newK1 = K1;
    while (rep <= maxiter) {
        const double K2 = kpass(1, t, P, sh);
        if (!P.fast && (!isfinite(K2) || fabs(K2) < 1e-15)) { c = 0; break; }   // non-fast only
        tnew = t - K1 / K2;
        if (P.fast ? isnan(tnew) : !isfinite(tnew)) { c = 0; break; }
        if (fabs(tnew - t) < tol) { c = 1; break; }
        if (rep == maxiter) { c = 0; break; }
        newK1 = kpass(0, tnew, P, sh) - q;
        const bool flipped = P.fast ? ((K1 * newK1) < 0) : (sgn(K1) != sgn(newK1));
        if (flipped) {
            if (fabs(tnew - t) > (prevJump - tol)) {
                tnew = t + (double)sgn(newK1 - K1) * prevJump / 2.0;
                newK1 = kpass(0, tnew, P, sh) - q;
                prevJump = prevJump / 2.0;
            } else {
                prevJump = fabs(tnew - t);
            }
        }
        rep = rep + 1;
        t = tnew;
        K1 = newK1;
    }
    *root = t; *niter = rep; *conv = c;
}

__device__ __forceinline__ double PhiUpper(double z) { return 0.5 * erfc(z / sqrt(2.0)); }
__device__ __forceinline__ double PhiLower(double z) { return 0.5 * erfc(-z / sqrt(2.0)); }

// Get_Saddle_Prob_Binom / Get_Saddle_Prob_fast_Binom.
__device__ double saddle(const PairCtx& P, double zeta, double q, int logp, double* sh, int* isSaddle)
{
    const double k1 = kpass(2, zeta, P, sh);
    const double k2 = kpass(1, zeta, P, sh);
    const double temp1 = zeta * q - k1;
    *isSaddle = 0;
    bool flagrun = false;
    double w = 0.0, v = 0.0;
    if (isfinite(k1) && isfinite(k2) && temp1 >= 0 && k2 >= 0) {
        w = (double)sgn(zeta) * sqrt(2.0 * temp1);
        v = zeta * sqrt(k2);
        if (w != 0) flagrun = true;
    }
    if (flagrun) {
        const double Ztest = w + (1.0 / w) * log(v / w);
        *isSaddle = 1;
        if (Ztest > 0) return logp ? log(PhiUpper(Ztest)) : PhiUpper(Ztest);
        else           return logp ? -log(PhiLower(Ztest)) : -PhiLower(Ztest);
    }
    return logp ? -CUDART_INF : 0.0;
}

__device__ __forceinline__ double addLogp(double p1, double p2)   // UTIL.cpp add_logp
{
    p1 = -fabs(p1); p2 = -fabs(p2);
    const double mx = fmax(p1, p2), mn = fmin(p1, p2);
    return mx + log(1.0 + exp(mn - mx));
}

__global__ void __launch_bounds__(NT)
spa_pairs(const unsigned char* __restrict__ packed, std::size_t bpv, const double* __restrict__ lut, int N,
          const double* __restrict__ MU, const double* __restrict__ XV, const double* __restrict__ XX,
          const int* __restrict__ pOfTrait, std::size_t traitStride,
          const SpaPairIn* __restrict__ in, int nPairs, double* __restrict__ scratch,
          double tol, int maxiter, SpaPairOut* __restrict__ out)
{
    __shared__ double sh[NT / 32 + 1];
    __shared__ double bsh[SPA_PMAX];
    double* gt = scratch + (std::size_t)blockIdx.x * N;
    for (int k = blockIdx.x; k < nPairs; k += gridDim.x) {
        const SpaPairIn pin = in[k];
        PairCtx P;
        P.col = packed + (std::size_t)pin.slot * bpv;
        {
            const double2* d = reinterpret_cast<const double2*>(lut) + 2 * (std::size_t)pin.slot;
            const double2 a = d[0], b = d[1];
            P.L = make_double4(a.x, a.y, b.x, b.y);
        }
        const int t = pin.trait;
        const int p = pOfTrait[t];
        P.mu = MU + (std::size_t)t * N;
        const double* xv = XV + (std::size_t)t * traitStride;   // sample i: xv[i*p + j]
        const double* xx = XX + (std::size_t)t * traitStride;   // xx[j*N + i]
        P.gt = gt; P.N = N; P.fast = pin.fast;

        // pass A: b = XV g over carriers (getadjGFast's loop over iIndex)
        double a[SPA_PMAX];
        for (int j = 0; j < SPA_PMAX; ++j) a[j] = 0.0;
        for (int i = threadIdx.x; i < N; i += NT) {
            const double g = dose(P.col, P.L, i);
            if (g == 0.0) continue;
            const double* x = xv + (std::size_t)i * p;
            for (int j = 0; j < p; ++j) a[j] += x[j] * g;
        }
        for (int j = 0; j < p; ++j) {
            const double s = blockSum(a[j], sh);
            if (threadIdx.x == 0) bsh[j] = s;
        }
        __syncthreads();

        // pass B: g~ = g - XXVX_inv b; gpos / gneg; m1; the fast variant's
        // carrier sums.
        double ap = 0.0, an = 0.0, m1 = 0.0, nbmu = 0.0, nbsig = 0.0;
        for (int i = threadIdx.x; i < N; i += NT) {
            const double g = dose(P.col, P.L, i);
            double proj = 0.0;
            for (int j = 0; j < p; ++j) proj += xx[(std::size_t)j * N + i] * bsh[j];
            const double v = g - proj;
            gt[i] = v;
            if (v > 0) ap += v; else if (v < 0) an += v;
            const double m = P.mu[i];
            m1 += m * v;
            if (g != 0.0) { nbmu += v * m; nbsig += m * (1.0 - m) * v * v; }
        }
        const double gpos = blockSum(ap, sh);
        const double gneg = blockSum(an, sh);
        m1    = blockSum(m1, sh);
        nbmu  = blockSum(nbmu, sh);
        nbsig = blockSum(nbsig, sh);
        P.NAmu = m1 - nbmu;
        P.NAsigma = pin.var2 - nbsig;

        // q, qinv as getMarkerPval forms them for a binary trait
        const double q = pin.Tstat / sqrt(pin.var1 / pin.var2) + m1;
        double qinv;
        if ((q - m1) > 0)       qinv = -1.0 * fabs(q - m1) + m1;
        else if ((q - m1) == 0) qinv = m1;
        else                    qinv = fabs(q - m1) + m1;

        double r1, r2; int n1, n2, c1, c2, s1 = -1, s2 = -1, conv = 0; double pv;
        getroot(P, q,    gpos, gneg, tol, maxiter, sh, &r1, &n1, &c1);
        getroot(P, qinv, gpos, gneg, tol, maxiter, sh, &r2, &n2, &c2);
        if (c1 && c2) {
            // spa.cpp SPA / SPA_fast: a non-saddle tail withdraws convergence
            double p1 = saddle(P, r1, q,    pin.logp, sh, &s1);
            double p2 = saddle(P, r2, qinv, pin.logp, sh, &s2);
            conv = 1;
            if (!s1) { conv = 0; p1 = pin.logp ? pin.pno - log(2.0) : pin.pno / 2; }
            if (!s2) { conv = 0; p2 = pin.logp ? pin.pno - log(2.0) : pin.pno / 2; }
            pv = pin.logp ? addLogp(p1, p2) : fabs(p1) + fabs(p2);
        } else {
            pv = pin.pno; conv = 0;
        }
        if (threadIdx.x == 0) {
            SpaPairOut o;
            o.pval = pv; o.conv = conv; o.s1 = s1; o.s2 = s2;
            o.niter1 = n1; o.niter2 = n2; o.root1 = r1; o.root2 = r2; o.m1 = m1;
            out[k] = o;
        }
        __syncthreads();
    }
}

}  // namespace

struct Spa {
    int N = 0, nTraits = 0, maxPairs = 0, blocks = 256, maxiter = 1000;
    double tol = 0.0;
    std::size_t traitStride = 0;
    SpaPairIn*  hIn  = nullptr;
    SpaPairOut* hOut = nullptr;
    SpaPairIn*  dIn  = nullptr;
    SpaPairOut* dOut = nullptr;
    double* dMu = nullptr; double* dXV = nullptr; double* dXX = nullptr;
    int*    dP  = nullptr;
    double* dScratch = nullptr;
    cudaStream_t st = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    double tKernel = 0.0;
    long long nPairsDone = 0;
    std::size_t devBytes = 0;
};

Spa* spaCreate(const SpaCreateArgs& a)
{
    if (a.N <= 0 || a.nTraits <= 0 || a.traits == nullptr || a.maxPairs <= 0) return nullptr;
    if (cudaSetDevice(a.device) != cudaSuccess) return nullptr;
    int pMax = 0;
    for (int t = 0; t < a.nTraits; ++t) {
        if (a.traits[t].p <= 0 || a.traits[t].p > SPA_PMAX) return nullptr;
        if (!a.traits[t].mu || !a.traits[t].XV || !a.traits[t].XXVX_inv) return nullptr;
        if (a.traits[t].p > pMax) pMax = a.traits[t].p;
    }
    Spa* s = new Spa();
    s->N = a.N; s->nTraits = a.nTraits; s->maxPairs = a.maxPairs;
    s->blocks = a.blocks > 0 ? a.blocks : 256; s->maxiter = a.maxiter; s->tol = a.tol;
    s->traitStride = (std::size_t)a.N * pMax;
    auto fail = [&]() -> Spa* { spaDestroy(s); return nullptr; };
    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        if (cudaMalloc(p, n) != cudaSuccess) { *p = nullptr; return false; }
        db += n; return true;
    };
    if (cudaHostAlloc((void**)&s->hIn,  (std::size_t)a.maxPairs * sizeof(SpaPairIn),  cudaHostAllocDefault) != cudaSuccess) return fail();
    if (cudaHostAlloc((void**)&s->hOut, (std::size_t)a.maxPairs * sizeof(SpaPairOut), cudaHostAllocDefault) != cudaSuccess) return fail();
    if (!dev((void**)&s->dIn,  (std::size_t)a.maxPairs * sizeof(SpaPairIn)))  return fail();
    if (!dev((void**)&s->dOut, (std::size_t)a.maxPairs * sizeof(SpaPairOut))) return fail();
    if (!dev((void**)&s->dMu, (std::size_t)a.nTraits * a.N * sizeof(double))) return fail();
    if (!dev((void**)&s->dXV, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail();
    if (!dev((void**)&s->dXX, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail();
    if (!dev((void**)&s->dP,  (std::size_t)a.nTraits * sizeof(int))) return fail();
    if (!dev((void**)&s->dScratch, (std::size_t)s->blocks * a.N * sizeof(double))) return fail();
    s->devBytes = db;
    std::vector<int> pv(a.nTraits);
    for (int t = 0; t < a.nTraits; ++t) {
        const SpaTraitArgs& T = a.traits[t];
        pv[t] = T.p;
        if (cudaMemcpy(s->dMu + (std::size_t)t * a.N, T.mu, (std::size_t)a.N * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(s->dXV + (std::size_t)t * s->traitStride, T.XV, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        if (cudaMemcpy(s->dXX + (std::size_t)t * s->traitStride, T.XXVX_inv, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    }
    if (cudaMemcpy(s->dP, pv.data(), pv.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
    if (cudaStreamCreate(&s->st) != cudaSuccess) return fail();
    if (cudaEventCreate(&s->e0) != cudaSuccess || cudaEventCreate(&s->e1) != cudaSuccess) return fail();
    return s;
}

void spaDestroy(Spa* s)
{
    if (!s) return;
    if (s->e0) cudaEventDestroy(s->e0);
    if (s->e1) cudaEventDestroy(s->e1);
    if (s->st) cudaStreamDestroy(s->st);
    if (s->dIn) cudaFree(s->dIn);
    if (s->dOut) cudaFree(s->dOut);
    if (s->dMu) cudaFree(s->dMu);
    if (s->dXV) cudaFree(s->dXV);
    if (s->dXX) cudaFree(s->dXX);
    if (s->dP) cudaFree(s->dP);
    if (s->dScratch) cudaFree(s->dScratch);
    if (s->hIn) cudaFreeHost(s->hIn);
    if (s->hOut) cudaFreeHost(s->hOut);
    delete s;
}

SpaPairIn*  spaIn(Spa* s)  { return s ? s->hIn : nullptr; }
SpaPairOut* spaOut(Spa* s) { return s ? s->hOut : nullptr; }
std::size_t spaDeviceBytes(const Spa* s) { return s ? s->devBytes : 0; }
void spaTimings(const Spa* s, double* t_kernel, long long* t_pairs)
{
    if (!s) return;
    if (t_kernel) *t_kernel = s->tKernel;
    if (t_pairs)  *t_pairs = s->nPairsDone;
}

bool spaRun(Spa* s, const Reducer* r, int nPairs, int t_devSet)
{
    if (!s || !r) return false;
    if (nPairs <= 0) return true;
    if (nPairs > s->maxPairs) { lastErrSpa = "nPairs > maxPairs"; return false; }
    const unsigned char* dPk = (const unsigned char*)devicePacked(r, t_devSet);
    const double* dLut = (const double*)deviceLut(r, t_devSet);
    if (!dPk || !dLut) { lastErrSpa = "reducer has no resident packed rows"; return false; }
    // The reducer's last reduce() is complete (it synchronises before
    // returning), so its resident rows are safe to read on our own stream.
    CKS(cudaMemcpyAsync(s->dIn, s->hIn, (std::size_t)nPairs * sizeof(SpaPairIn), cudaMemcpyHostToDevice, s->st));
    const int grid = nPairs < s->blocks ? nPairs : s->blocks;
    CKS(cudaEventRecord(s->e0, s->st));
    spa_pairs<<<grid, NT, 0, s->st>>>(dPk, bytesPerSlot(r), dLut, s->N, s->dMu, s->dXV, s->dXX, s->dP,
                                      s->traitStride, s->dIn, nPairs, s->dScratch, s->tol, s->maxiter, s->dOut);
    CKS(cudaGetLastError());
    CKS(cudaEventRecord(s->e1, s->st));
    CKS(cudaMemcpyAsync(s->hOut, s->dOut, (std::size_t)nPairs * sizeof(SpaPairOut), cudaMemcpyDeviceToHost, s->st));
    CKS(cudaStreamSynchronize(s->st));
    float ms = 0;
    if (cudaEventElapsedTime(&ms, s->e0, s->e1) == cudaSuccess) s->tKernel += ms * 1e-3;
    s->nPairsDone += nPairs;
    return true;
}

}  // namespace gpu2
}  // namespace saige
