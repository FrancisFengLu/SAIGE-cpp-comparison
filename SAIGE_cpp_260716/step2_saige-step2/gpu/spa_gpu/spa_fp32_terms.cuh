// spa_fp32_terms.cuh -- per-sample float arithmetic shared by the fp32
// variants of both device SPA implementations (gpuPrecisionSPA: fp32;
// spa_gpu.cu for gpuSpaImpl: lib, ../gpu_spa.cu for gpuSpaImpl: own).
// Derivation of the centred sums they feed: spa_gpu.hpp, "FP32".
#pragma once

#include <cuda_runtime.h>

#include <cstddef>

namespace saige {
namespace spa_fp32 {

// Neumaier-compensated float sum; read out in double.
struct NSum {
    float s = 0.0f, c = 0.0f;
    __device__ __forceinline__ void add(float v)
    {
        const float t = s + v;
        if (fabsf(s) >= fabsf(v)) c += (s - t) + v; else c += (v - t) + s;
        s = t;
    }
    __device__ __forceinline__ double val() const { return (double)s + (double)c; }
    // Error-free product m * v (TwoProd via an explicit fma, unaffected by
    // --fmad=false) added with its rounding error: the sum is then exact to
    // the compensated-sum level, not to one rounding per product.
    __device__ __forceinline__ void addProd(float m, float v)
    {
        const float pr = m * v;
        add(pr);
        c += __fmaf_rn(m, v, -pr);
    }
    // this += (s2, c2), another compensated sum (TwoSum on the leading parts).
    __device__ __forceinline__ void merge(float s2, float c2)
    {
        const float t = s + s2, bp = t - s;
        const float e = (s - (t - bp)) + (s2 - bp);
        s = t;
        c += c2 + e;
    }
};

// Block reduction of K compensated float sums, float throughout (no fp64
// instruction): warp shuffle-down tree, then the NWARP warp totals merged in a
// fixed order, so the result is deterministic and every thread gets it.
// shf: 2 * K * NWARP floats of shared memory.
template <int K, int NWARP>
__device__ __forceinline__ void blockReduceF(NSum (&v)[K], float* shf)
{
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    #pragma unroll
    for (int k = 0; k < K; ++k) {
        #pragma unroll
        for (int o = 16; o > 0; o >>= 1) {
            const float s2 = __shfl_down_sync(0xffffffffu, v[k].s, o);
            const float c2 = __shfl_down_sync(0xffffffffu, v[k].c, o);
            v[k].merge(s2, c2);
        }
    }
    if (lane == 0) {
        #pragma unroll
        for (int k = 0; k < K; ++k) { shf[(2 * k) * NWARP + warp] = v[k].s; shf[(2 * k + 1) * NWARP + warp] = v[k].c; }
    }
    __syncthreads();
    #pragma unroll
    for (int k = 0; k < K; ++k) {
        NSum r; r.s = shf[(2 * k) * NWARP]; r.c = shf[(2 * k + 1) * NWARP];
        #pragma unroll
        for (int w = 1; w < NWARP; ++w) r.merge(shf[(2 * k) * NWARP + w], shf[(2 * k + 1) * NWARP + w]);
        v[k] = r;
    }
    __syncthreads();
}

// Carrier g~ = g - XXVX_inv b to about float^2 accuracy, all float
// instructions: XXVX_inv as hi + lo floats (xh[j*N], xl[j*N]), b as hi + lo
// (bh, bl), products error-free. Returns g~ = vh + vl, |vl| <= ulp(vh)/2.
// Only carriers need it: their g~ feed NAsigma = var2 - sum_c mu(1-mu) g~^2,
// which cancels by up to ~1e3 for rare markers (see addVarTerm), so a
// float-accurate g~ there moves the root.
template <int PMAXT>
__device__ __forceinline__ void carrierGt(float g, const float* xh, const float* xl, std::size_t ld,
                                          const float (&bh)[PMAXT], const float (&bl)[PMAXT], int p,
                                          float* vh, float* vl)
{
    NSum pr;
    #pragma unroll
    for (int j = 0; j < PMAXT; ++j) if (j < p) {
        const float h = xh[(std::size_t)j * ld], l = xl[(std::size_t)j * ld];
        pr.addProd(h, bh[j]);
        pr.c += h * bl[j] + l * bh[j];
    }
    const float ph = pr.s + pr.c, pl = pr.c - (ph - pr.s);
    // (g - ph) exactly, then - pl, renormalised
    const float d = g - ph, bp = d - g;
    const float e = ((g - (d - bp)) + (-ph - bp)) - pl;
    const float r = d + e;
    *vh = r;
    *vl = e - (r - d);
}

// sum += mu (1 - mu) v^2 to about float^2 accuracy, all float instructions:
// mu = mh + ml (|mh| <= 1/2, the stored hi / lo pair), v = vh + vl (from
// carrierGt). Used for NAsigma = var2 - sum_c mu (1-mu) g~^2, where the
// subtraction cancels by up to ~1e3 for very rare markers, so a float-accurate
// sum would move the root (and can move it across fp64's Korg overflow
// threshold).
__device__ __forceinline__ void addVarTerm(NSum& sum, float mh, float ml, float vh, float vl)
{
    // a = 1 - mu = ah + al
    const float ah = 1.0f - mh, bq = ah - 1.0f;
    const float al = ((1.0f - (ah - bq)) + (-mh - bq)) - ml;
    // w = mu (1 - mu) = wh + wl
    const float wh = mh * ah;
    const float wl = __fmaf_rn(mh, ah, -wh) + (mh * al + ml * ah);
    // s = v^2 = sh + sl
    const float sh = vh * vh;
    const float sl = __fmaf_rn(vh, vh, -sh) + 2.0f * vh * vl;
    sum.addProd(wh, sh);
    sum.c += wh * sl + wl * sh;
}

// Warp sum of an int (carrier counts), shuffle-down; lane 0 holds the total.
__device__ __forceinline__ int warpSumI(int v)
{
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

// log1p(u) - u. |u| < 1/4: -u^2/(2+u) + 2 s^3 (1/3 + s^2/5 + ...) with
// s = u/(2+u) (log1p(u) = 2 atanh(s)); truncation below 1e-9 of the result.
__device__ __forceinline__ float log1pmx(float u)
{
    if (fabsf(u) < 0.25f) {
        const float s = u / (2.0f + u), s2 = s * s;
        const float tail = 2.0f * s * s2 * (1.0f / 3 + s2 * (1.0f / 5 + s2 * (1.0f / 7 + s2 * (1.0f / 9 + s2 * (1.0f / 11)))));
        return -(u * u) / (2.0f + u) + tail;
    }
    return log1pf(u) - u;
}

// expm1(x) - x. |x| < 1/2: Taylor to x^8 (next term below 1e-7 of the result).
__device__ __forceinline__ float expm1mx(float x)
{
    if (fabsf(x) < 0.5f)
        return (x * x) * (0.5f + x * (1.0f / 6 + x * (1.0f / 24 + x * (1.0f / 120 + x * (1.0f / 720 +
               x * (1.0f / 5040 + x * (1.0f / 40320)))))));
    return expm1f(x) - x;
}

// h = log(1 - mu + mu e^x) - mu x, for mu = m <= 1/2, a = 1 - m, and
// em = expm1(-x) when x >= 0 (unused otherwise).
//   |x| < 0.1      Bernoulli cumulants, m a x^2 [1/2 + (1-2m) x/6
//                  + (1-6ma) x^2/24 + (1-2m)(1-12ma) x^3/120]; next term < 1e-7
//   0.1 <= |x| < 2 [log1p(u) - u] + m [expm1(x) - x], u = m expm1(x): the
//                  brackets are about -m^2 x^2/2 and m x^2/2, so the sum loses
//                  at most a factor 1/(1-m) <= 2, where log1p(u) - m x loses
//                  ~ 2/((1-m)|x|) and a x + log1p(a em) ~ 2/(m|x|)
//   other x < 80   log1p(m expm1(x)) - m x: loses at most a small factor there
//                  (the split form would subtract two numbers of size e^x)
//   x >= 80        a x + log1p(a em): m expm1(x) may overflow float, and
//                  nothing cancels
__device__ __forceinline__ float hTerm(float m, float a, float x, float em)
{
    const float ax = fabsf(x);
    if (ax < 0.1f) {
        const float ma = m * a, b = 1.0f - 2.0f * m;
        const float c = 0.5f + x * (b * (1.0f / 6) + x * ((1.0f - 6.0f * ma) * (1.0f / 24) +
                        x * (b * (1.0f - 12.0f * ma) * (1.0f / 120))));
        return ma * (x * x) * c;
    }
    if (ax < 2.0f) return log1pmx(m * expm1f(x)) + m * expm1mx(x);
    if (x < 80.0f) return log1pf(m * expm1f(x)) - m * x;
    return a * x + log1pf(a * em);
}

}  // namespace spa_fp32
}  // namespace saige
