// gpu_spa.cu — one CUDA block per flagged (marker, trait) pair; the block's
// threads share every length-N sum, so SAIGE's Newton loop runs unchanged per
// pair. Contract and provenance: gpu_spa.hpp. Reference: spa_binary.cpp
// (getroot_K1_Binom / getroot_K1_fast_Binom / Get_Saddle_Prob_*_Binom) and
// spa.cpp (SPA / SPA_fast). Compiled with --fmad=false so a*b+c rounds twice,
// as the CPU build (-std=c++17, no contraction) does.

#include "gpu_spa.hpp"
#include "gpu_step2.hpp"
#include "spa_gpu/spa_fp32_terms.cuh"

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
// Covariate columns of b = XV g formed per pass over the samples (pass A):
// SPA_PTILE per-thread accumulators, p / SPA_PTILE passes. One pass for
// p <= SPA_PTILE, which is the loop as it was before p > 8 was allowed; for
// any p, column j is summed over the same samples in the same order, so b
// does not depend on the tiling.
constexpr int SPA_PTILE = 8;

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
    // fp32 variant only: (s g~, s min(mu, 1 - mu)) per sample, s = -1 where
    // mu > 1/2, and the m1 of those stored floats, for centring (see kpass32).
    const float2* gt32;
    double m1;
    float4 L32;            // the dosage table as floats
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

// ---------------------------------------------------------------------------
// fp32 variant (gpuPrecisionSPA: fp32, gpuSpaImpl: own). Same passes, same
// Newton loop; the per-sample work of kpass runs in float and the sums are
// centred so the large cancellations happen in double (the same rewrite as
// spa_gpu/spa_gpu.cu's fp32 variant, where it is derived):
//   mode 0: sum g~ (p - mu)            (+ NAsigma t)          = K1 - m1
//   mode 1: sum g~^2 p (1 - p)         (+ NAsigma)            = K2
//   mode 2: sum log(1-mu+mu e^x) - mu x (+ NAsigma t^2 / 2)   = Korg - t m1
// with x = g~ t, p = mu e^x / (1 - mu + mu e^x); only exp(-|x|) <= 1 is
// formed, so nothing overflows float. Stored mu is min(mu, 1 - mu), negative
// when flipped. fp64's own overflow of Korg (exp(g~ t) = inf for
// g~ t > 709.78, which makes the tail "not a saddle") is reproduced: mode 2
// returns +inf then. Per-thread sums are Neumaier-compensated floats, merged
// in float-float across the block (blockSumF) and read out in double.
// ---------------------------------------------------------------------------
using spa_fp32::NSum;

__device__ __forceinline__ float dose32(const unsigned char* col, const float4& L, int i)
{
    const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
    return (c == 0) ? L.x : (c == 1 ? L.y : (c == 2 ? L.z : L.w));
}
// Block sum of one compensated float sum, float throughout; every thread gets
// it. sh: the fp64 path's NT/32 + 1 doubles, used as 2 NT/32 floats.
__device__ __forceinline__ NSum blockSumF(NSum v, double* sh)
{
    NSum r[1] = {v};
    spa_fp32::blockReduceF<1, NT / 32>(r, reinterpret_cast<float*>(sh));
    return r[0];
}

__device__ double kpass32(int mode, double t, const PairCtx& P, double* sh)
{
    NSum acc;
    bool inf = false;
    const float th = (float)t, tl = (float)(t - (double)th);   // t = th + tl: g~ t to float precision
    for (int i = threadIdx.x; i < P.N; i += NT) {
        if (P.fast && dose32(P.col, P.L32, i) == 0.0f) continue;
        const float2 e2 = P.gt32[i];
        const float g = e2.x, m = fabsf(e2.y), a = 1.0f - m;   // m <= 1/2
        const float x = g * th + g * tl;
        if (mode == 2 && (signbit(e2.y) ? -x : x) > 709.782712893384f) inf = true;   // fp64 exp overflow
        float p, q, pm, h = 0.0f;
        if (x >= 0.0f) {
            const float e = expf(-x), em = expm1f(-x);
            const float den = a * e + m;
            p = m / den; q = (a * e) / den; pm = -(m * a * em) / den;
            if (mode == 2) h = spa_fp32::hTerm(m, a, x, em);
        } else {
            const float e = expf(x), em = expm1f(x);
            const float den = a + m * e;
            p = (m * e) / den; q = a / den; pm = (m * a * em) / den;
            if (mode == 2) h = spa_fp32::hTerm(m, a, x, 0.0f);
        }
        if (mode == 0) {
            acc.add(g * pm);
        } else if (mode == 1) {
            const float term = (g * g) * (p * q);
            if (isfinite(term)) acc.add(term);
        } else {
            acc.add(h);
        }
    }
    inf = __syncthreads_or(mode == 2 && inf);
    const double sv = blockSumF(acc, sh).val();
    double s = inf ? CUDART_INF : sv;
    if (P.fast) {
        if (mode == 0)      s += P.NAsigma * t;
        else if (mode == 1) s += P.NAsigma;
        else                s += 0.5 * P.NAsigma * t * t;
    }
    return s;
}

// The pass of the chosen arithmetic. Modes 0 and 2 return K1 - m1 and
// Korg - t m1 in fp32 (centred), K1 and Korg in fp64.
template <bool F32>
__device__ __forceinline__ double kp(int mode, double t, const PairCtx& P, double* sh)
{
    if constexpr (F32) return kpass32(mode, t, P, sh); else return kpass(mode, t, P, sh);
}

// getroot_K1_Binom (fast == 0) / getroot_K1_fast_Binom (fast == 1), from init 0.
// F32: K1 is formed as (K1 - m1) - (q - m1).
template <bool F32>
__device__ void getroot(const PairCtx& P, double q, double gpos, double gneg, double tol, int maxiter,
                        double* sh, double* root, int* niter, int* conv)
{
    if (q >= gpos || q <= gneg) { *root = CUDART_INF; *niter = 0; *conv = 1; return; }
    const double qc = F32 ? q - P.m1 : q;   // what the pass's sum is compared with
    double t = 0.0;
    double K1 = kp<F32>(0, t, P, sh) - qc;
    double prevJump = CUDART_INF;
    int rep = 1, c = 1;
    double tnew = t, newK1 = K1;
    while (rep <= maxiter) {
        const double K2 = kp<F32>(1, t, P, sh);
        if (!P.fast && (!isfinite(K2) || fabs(K2) < 1e-15)) { c = 0; break; }   // non-fast only
        tnew = t - K1 / K2;
        if (P.fast ? isnan(tnew) : !isfinite(tnew)) { c = 0; break; }
        // fp32: max(tol, 2^-17 max(|t|, |tnew|)) -- see gpu_spa.hpp
        const double tolE = F32 ? fmax(tol, (1.0 / 131072) * fmax(fabs(t), fabs(tnew))) : tol;
        if (fabs(tnew - t) < tolE) { c = 1; break; }
        if (rep == maxiter) { c = 0; break; }
        newK1 = kp<F32>(0, tnew, P, sh) - qc;
        const bool flipped = P.fast ? ((K1 * newK1) < 0) : (sgn(K1) != sgn(newK1));
        if (flipped) {
            if (fabs(tnew - t) > (prevJump - tolE)) {
                tnew = t + (double)sgn(newK1 - K1) * prevJump / 2.0;
                newK1 = kp<F32>(0, tnew, P, sh) - qc;
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

// Get_Saddle_Prob_Binom / Get_Saddle_Prob_fast_Binom. F32: k1 = Korg - zeta m1,
// so temp1 = zeta (q - m1) - k1, in double.
template <bool F32>
__device__ double saddle(const PairCtx& P, double zeta, double q, int logp, double* sh, int* isSaddle)
{
    const double k1 = kp<F32>(2, zeta, P, sh);
    const double k2 = kp<F32>(1, zeta, P, sh);
    const double temp1 = F32 ? zeta * (q - P.m1) - k1 : zeta * q - k1;
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

template <bool F32>
__global__ void __launch_bounds__(NT)
spa_pairs(const unsigned char* __restrict__ packed, std::size_t bpv, const double* __restrict__ lut, int N,
          const double* __restrict__ MU, const double* __restrict__ XV, const double* __restrict__ XX,
          const int* __restrict__ pOfTrait, std::size_t traitStride,
          const float* __restrict__ MU32, const float* __restrict__ MU32lo,
          const float* __restrict__ XV32, const float* __restrict__ XX32,
          const float* __restrict__ XV32lo, const float* __restrict__ XX32lo,
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

        double gpos, gneg, m1;
        if constexpr (F32) {
            // gpuPrecisionSPA fp32: passes A and B in float (float copies of XV,
            // XXVX_inv and mu; compensated sums; float block sums), same outputs.
            const float* mu32 = MU32 + (std::size_t)t * N;
            const float* xv32 = XV32 + (std::size_t)t * traitStride;
            const float* xx32 = XX32 + (std::size_t)t * traitStride;
            const float* xv32lo = XV32lo + (std::size_t)t * traitStride;
            const float* xx32lo = XX32lo + (std::size_t)t * traitStride;
            P.L32 = make_float4((float)P.L.x, (float)P.L.y, (float)P.L.z, (float)P.L.w);
            __shared__ float bsh32[SPA_PMAX], bsh32lo[SPA_PMAX];   // b = hi + lo
            // pass A: b = XV g over carriers, SPA_PTILE columns per pass
            for (int j0 = 0; j0 < p; j0 += SPA_PTILE) {
                const int pt = (p - j0 < SPA_PTILE) ? p - j0 : SPA_PTILE;
                NSum a[SPA_PTILE];
                for (int i = threadIdx.x; i < N; i += NT) {
                    const float g = dose32(P.col, P.L32, i);
                    if (g == 0.0f) continue;
                    const float* x = xv32 + (std::size_t)i * p + j0;
                    const float* xl = xv32lo + (std::size_t)i * p + j0;
                    for (int j = 0; j < pt; ++j) { a[j].addProd(x[j], g); a[j].c += xl[j] * g; }
                }
                for (int j = 0; j < pt; ++j) {
                    const NSum s = blockSumF(a[j], sh);
                    if (threadIdx.x == 0) { bsh32[j0 + j] = s.s + s.c; bsh32lo[j0 + j] = s.c - (bsh32[j0 + j] - s.s); }
                }
            }
            __syncthreads();
            // pass B: g~, gpos / gneg, m1 = sum mu g~ from exactly the stored
            // floats (error-free products; also the centring of kpass32), the
            // fast variant's carrier sums; store (s g~, s min(mu, 1 - mu)).
            // b is read from the block's shared copy (any p).
            NSum sp, sn, sm, snbmu, snbsig;
            for (int i = threadIdx.x; i < N; i += NT) {
                const float g = dose32(P.col, P.L32, i);
                float v, vl = 0.0f;
                if (g != 0.0f) {
                    spa_fp32::carrierGtP(g, xx32 + i, xx32lo + i, (std::size_t)N, bsh32, bsh32lo, p, &v, &vl);
                } else {
                    float proj = 0.0f;
                    for (int j = 0; j < p; ++j) proj += xx32[(std::size_t)j * N + i] * bsh32[j];
                    v = g - proj;
                }
                const float ms = mu32[i], m = fabsf(ms);
                const bool flip = signbit(ms);   // stored mu is 1 - mu
                reinterpret_cast<float2*>(gt)[i] = make_float2(flip ? -v : v, ms);
                if (v > 0) sp.add(v); else if (v < 0) sn.add(v);
                if (flip) { sm.add(v); sm.addProd(-m, v); } else sm.addProd(m, v);
                if (g != 0.0f) {
                    if (flip) { snbmu.add(v); snbmu.addProd(-m, v); } else snbmu.addProd(m, v);
                    spa_fp32::addVarTerm(snbsig, m, MU32lo[(std::size_t)t * N + i], v, vl);
                }
            }
            gpos = blockSumF(sp, sh).val();
            gneg = blockSumF(sn, sh).val();
            m1   = blockSumF(sm, sh).val();
            const double nbmu = blockSumF(snbmu, sh).val();
            const double nbsig = blockSumF(snbsig, sh).val();
            P.NAmu = m1 - nbmu;
            P.NAsigma = pin.var2 - nbsig;
            P.gt32 = reinterpret_cast<const float2*>(gt);
            P.m1 = m1;
        } else {
            // pass A: b = XV g over carriers (getadjGFast's loop over iIndex),
            // SPA_PTILE columns per pass
            for (int j0 = 0; j0 < p; j0 += SPA_PTILE) {
                const int pt = (p - j0 < SPA_PTILE) ? p - j0 : SPA_PTILE;
                double a[SPA_PTILE];
                for (int j = 0; j < SPA_PTILE; ++j) a[j] = 0.0;
                for (int i = threadIdx.x; i < N; i += NT) {
                    const double g = dose(P.col, P.L, i);
                    if (g == 0.0) continue;
                    const double* x = xv + (std::size_t)i * p + j0;
                    for (int j = 0; j < pt; ++j) a[j] += x[j] * g;
                }
                for (int j = 0; j < pt; ++j) {
                    const double s = blockSum(a[j], sh);
                    if (threadIdx.x == 0) bsh[j0 + j] = s;
                }
            }
            __syncthreads();

            // pass B: g~ = g - XXVX_inv b; gpos / gneg; m1; the fast variant's
            // carrier sums.
            double ap = 0.0, an = 0.0, nbmu = 0.0, nbsig = 0.0;
            m1 = 0.0;
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
            gpos = blockSum(ap, sh);
            gneg = blockSum(an, sh);
            m1    = blockSum(m1, sh);
            nbmu  = blockSum(nbmu, sh);
            nbsig = blockSum(nbsig, sh);
            P.NAmu = m1 - nbmu;
            P.NAsigma = pin.var2 - nbsig;
        }

        // q, qinv as getMarkerPval forms them for a binary trait
        const double q = pin.Tstat / sqrt(pin.var1 / pin.var2) + m1;
        double qinv;
        if ((q - m1) > 0)       qinv = -1.0 * fabs(q - m1) + m1;
        else if ((q - m1) == 0) qinv = m1;
        else                    qinv = fabs(q - m1) + m1;

        double r1, r2; int n1, n2, c1, c2, s1 = -1, s2 = -1, conv = 0; double pv;
        getroot<F32>(P, q,    gpos, gneg, tol, maxiter, sh, &r1, &n1, &c1);
        getroot<F32>(P, qinv, gpos, gneg, tol, maxiter, sh, &r2, &n2, &c2);
        if (c1 && c2) {
            // spa.cpp SPA / SPA_fast: a non-saddle tail withdraws convergence
            double p1 = saddle<F32>(P, r1, q,    pin.logp, sh, &s1);
            double p2 = saddle<F32>(P, r2, qinv, pin.logp, sh, &s2);
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
    Prec prec = Prec::FP64;   // SpaCreateArgs::precision
    int N = 0, nTraits = 0, maxPairs = 0, blocks = 256, maxiter = 1000;
    double tol = 0.0;
    std::size_t traitStride = 0;
    SpaPairIn*  hIn  = nullptr;
    SpaPairOut* hOut = nullptr;
    SpaPairIn*  dIn  = nullptr;
    SpaPairOut* dOut = nullptr;
    double* dMu = nullptr; double* dXV = nullptr; double* dXX = nullptr;
    float*  dMu32 = nullptr; float* dXV32 = nullptr; float* dXX32 = nullptr;   // FP32 only
    float*  dMu32lo = nullptr; float* dXV32lo = nullptr; float* dXX32lo = nullptr;
    int*    dP  = nullptr;
    double* dScratch = nullptr;
    cudaStream_t st = nullptr;
    cudaEvent_t e0 = nullptr, e1 = nullptr;
    double tKernel = 0.0;
    long long nPairsDone = 0;
    std::size_t devBytes = 0;
};

bool spaSupports(Prec t_p)
{
    switch (t_p) {
        case Prec::FP64: return true;
        case Prec::FP32: return true;    // spa_pairs<true>
        case Prec::INT8: return false;   // not a mode of this stage
    }
    return false;
}

Spa* spaCreate(const SpaCreateArgs& a)
{
    lastErrSpa.clear();
    if (a.N <= 0 || a.nTraits <= 0 || a.traits == nullptr || a.maxPairs <= 0) { lastErrSpa = "bad arguments"; return nullptr; }
    // ---- precision dispatch ----
    if (!spaSupports(a.precision)) {
        lastErrSpa = std::string("SPA precision ") + precName(a.precision) + " is not implemented yet";
        return nullptr;
    }
    if (cudaSetDevice(a.device) != cudaSuccess) return nullptr;
    int pMax = 0;
    for (int t = 0; t < a.nTraits; ++t) {
        if (a.traits[t].p <= 0 || a.traits[t].p > SPA_PMAX) {
            lastErrSpa = "trait " + std::to_string(t) + " has p = " + std::to_string(a.traits[t].p) +
                         " covariate columns; the device SPA takes 1.." + std::to_string(SPA_PMAX);
            return nullptr;
        }
        if (!a.traits[t].mu || !a.traits[t].XV || !a.traits[t].XXVX_inv) return nullptr;
        if (a.traits[t].p > pMax) pMax = a.traits[t].p;
    }
    Spa* s = new Spa();
    s->prec = a.precision;
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
    if (a.precision == Prec::FP32) {
        // float copies for passes A and B; mu as min(mu, 1 - mu), negated
        // where mu > 1/2 (1 - mu formed in double)
        auto dev32 = [&](float** p, std::size_t n) {
            if (cudaMalloc((void**)p, n * sizeof(float)) != cudaSuccess) { *p = nullptr; return false; }
            s->devBytes += n * sizeof(float); return true;
        };
        if (!dev32(&s->dMu32, (std::size_t)a.nTraits * a.N) || !dev32(&s->dMu32lo, (std::size_t)a.nTraits * a.N) ||
            !dev32(&s->dXV32lo, (std::size_t)a.nTraits * s->traitStride) || !dev32(&s->dXX32lo, (std::size_t)a.nTraits * s->traitStride) || !dev32(&s->dXV32, (std::size_t)a.nTraits * s->traitStride) ||
            !dev32(&s->dXX32, (std::size_t)a.nTraits * s->traitStride)) return fail();
        std::vector<float> f;
        for (int t = 0; t < a.nTraits; ++t) {
            const SpaTraitArgs& T = a.traits[t];
            f.resize((std::size_t)a.N);
            for (int i = 0; i < a.N; ++i) f[i] = T.mu[i] > 0.5 ? -(float)(1.0 - T.mu[i]) : (float)T.mu[i];
            if (cudaMemcpy(s->dMu32 + (std::size_t)t * a.N, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
            for (int i = 0; i < a.N; ++i) { const double mm = T.mu[i] > 0.5 ? 1.0 - T.mu[i] : T.mu[i]; f[i] = (float)(mm - (double)(float)mm); }
            if (cudaMemcpy(s->dMu32lo + (std::size_t)t * a.N, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
            f.resize((std::size_t)a.N * T.p);
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)T.XV[i];
            if (cudaMemcpy(s->dXV32 + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)(T.XV[i] - (double)(float)T.XV[i]);
            if (cudaMemcpy(s->dXV32lo + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)T.XXVX_inv[i];
            if (cudaMemcpy(s->dXX32 + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)(T.XXVX_inv[i] - (double)(float)T.XXVX_inv[i]);
            if (cudaMemcpy(s->dXX32lo + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail();
        }
    }
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
    if (s->dMu32) cudaFree(s->dMu32);
    if (s->dMu32lo) cudaFree(s->dMu32lo);
    if (s->dXV32lo) cudaFree(s->dXV32lo);
    if (s->dXX32lo) cudaFree(s->dXX32lo);
    if (s->dXV32) cudaFree(s->dXV32);
    if (s->dXX32) cudaFree(s->dXX32);
    if (s->dScratch) cudaFree(s->dScratch);
    if (s->hIn) cudaFreeHost(s->hIn);
    if (s->hOut) cudaFreeHost(s->hOut);
    delete s;
}

SpaPairIn*  spaIn(Spa* s)  { return s ? s->hIn : nullptr; }
SpaPairOut* spaOut(Spa* s) { return s ? s->hOut : nullptr; }
std::size_t spaDeviceBytes(const Spa* s) { return s ? s->devBytes : 0; }
const char* spaLastError() { return lastErrSpa.c_str(); }
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
    // ---- precision dispatch: the kernel variant ----
    if (s->prec == Prec::FP64) {
        spa_pairs<false><<<grid, NT, 0, s->st>>>(dPk, bytesPerSlot(r), dLut, s->N, s->dMu, s->dXV, s->dXX, s->dP,
                                          s->traitStride, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, s->dIn, nPairs, s->dScratch, s->tol, s->maxiter, s->dOut);
    } else if (s->prec == Prec::FP32) {
        // Newton tolerance: the caller's, floored at 1e-5 (gpu_spa.hpp).
        const double tol32 = s->tol < 1e-5 ? 1e-5 : s->tol;
        spa_pairs<true><<<grid, NT, 0, s->st>>>(dPk, bytesPerSlot(r), dLut, s->N, s->dMu, s->dXV, s->dXX, s->dP,
                                                s->traitStride, s->dMu32, s->dMu32lo, s->dXV32, s->dXX32, s->dXV32lo, s->dXX32lo, s->dIn, nPairs, s->dScratch, tol32, s->maxiter, s->dOut);
    } else {
        lastErrSpa = std::string("SPA precision ") + precName(s->prec) + " has no kernel";
        return false;
    }
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
