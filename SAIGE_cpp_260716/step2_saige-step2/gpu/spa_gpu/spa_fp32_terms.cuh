// spa_fp32_terms.cuh -- per-sample float arithmetic shared by the fp32
// variants of both device SPA implementations (gpuPrecisionSPA: fp32;
// spa_gpu.cu for gpuSpaImpl: lib, ../gpu_spa.cu for gpuSpaImpl: own).
// Derivation of the centred sums they feed: spa_gpu.hpp, "FP32".
#pragma once

#include <cuda_runtime.h>

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
};

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
