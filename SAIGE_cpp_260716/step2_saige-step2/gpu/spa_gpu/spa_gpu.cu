// spa_gpu.cu -- one CUDA block per flagged (marker, trait) pair; the block's
// 256 threads share every sum, so SAIGE's scalar Newton loop runs unchanged
// per pair. Contract and provenance: spa_gpu.hpp. Reference: spa_binary.cpp
// (getroot_K1_Binom / getroot_K1_fast_Binom / Get_Saddle_Prob_*_Binom),
// spa.cpp (SPA / SPA_fast), saige_test.cpp (getMarkerPval: g~, m1, q, qinv,
// the variant rule, NAmu / NAsigma), UTIL.cpp (add_logp, sum_arma1).
//
// Build with --fmad=false: a*b+c must round twice, as the CPU's -std=c++17
// build does, or the per-term values drift by an ulp before the sums start.
//
// Layout of a pair's work
//   pass A  b = XV g over carriers, and the carrier count          (reads packed, XV)
//   pass B  g~ = g - XXVX_inv b; gpos, gneg, m1 and the fast
//           variant's carrier sums; stores (g~, mu) interleaved    (reads packed, XX, mu)
//           into this block's scratch -- all N for the full
//           variant, carriers only (stream-compacted, index
//           order kept) for the fast one
//   roots   each Newton step is ONE pass over the stored pairs
//           producing K1(t)-q and K2(t) together; SAIGE evaluates
//           K2(t) at the top of an iteration and K1(tnew) at the
//           bottom, which is the same two numbers one step apart
//   tails   one pass per tail producing Korg(zeta) and K2(zeta)
// Passes A and B give each warp a contiguous eighth of the samples, so the
// fast variant's compaction needs no block-wide scan: a warp's write base is
// the carrier count of the warps before it, which pass A's reduction already
// has. Every other pass is thread-strided over the contiguous store.
//
// Determinism: for a given N and block count the launch shape and every
// reduction order are fixed, so two runs give bit-identical results.

#include "spa_gpu.hpp"

#include <cuda_runtime.h>
#include <math_constants.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "erfc_boost53.cuh"
#include "spa_fp32_terms.cuh"

namespace saige {
namespace spa_gpu {

namespace {

constexpr int NT    = 256;
constexpr int NWARP = NT / 32;
constexpr int NACC  = PMAX + 1;   // widest reduction: p partial dot products + the carrier count
// gpuPrecisionSPA fp32: floor of the Newton |dt| tolerance (spa_gpu.hpp).
constexpr double kTolFp32 = 1e-5;

thread_local std::string g_lastErr;
#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    g_lastErr = std::string(#x) + ": " + cudaGetErrorString(e_); return false; } } while (0)

// Boost's erfc, ported: erfc_boost53.cuh. mode 1: scaled port (default),
// 2: literal port, anything else: CUDA libm.
__device__ __forceinline__ double erfcSel(double z, int mode)
{
    if (mode == 1) return erfc53::erfImp53s(z, true);
    if (mode == 2) return erfc53::erfImp53(z, true);
    return erfc(z);
}

// boost::math::cdf(complement(normal(0,1), x)) and cdf(normal(0,1), x):
// erfc(+-(x - 0) / (1 * root_two)) / 2, with the +-infinity shortcuts.
constexpr double ROOT_TWO = 1.4142135623730951;   // constants::root_two<double>()
__device__ __forceinline__ double phiUpper(double x, int mode)
{
    if (isinf(x)) return x < 0 ? 1.0 : 0.0;
    const double diff = (x - 0.0) / (1.0 * ROOT_TWO);
    return erfcSel(diff, mode) / 2;
}
__device__ __forceinline__ double phiLower(double x, int mode)
{
    if (isinf(x)) return x < 0 ? 0.0 : 1.0;
    const double diff = (x - 0.0) / (1.0 * ROOT_TWO);
    return erfcSel(-diff, mode) / 2;
}

// ---------------------------------------------------------------------------
// Reductions. Every thread receives the block total; the 8 warp partials are
// added in a fixed order, so the result is deterministic.
// ---------------------------------------------------------------------------
__device__ __forceinline__ double warpSum(double v)
{
    #pragma unroll
    for (int o = 16; o > 0; o >>= 1) v += __shfl_down_sync(0xffffffffu, v, o);
    return v;
}

template <int K>
__device__ __forceinline__ void blockReduce(double (&v)[K], double (*sh)[NWARP])
{
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5;
    #pragma unroll
    for (int k = 0; k < K; ++k) v[k] = warpSum(v[k]);
    if (lane == 0) {
        #pragma unroll
        for (int k = 0; k < K; ++k) sh[k][warp] = v[k];
    }
    __syncthreads();
    #pragma unroll
    for (int k = 0; k < K; ++k) {
        double s = 0.0;
        #pragma unroll
        for (int w = 0; w < NWARP; ++w) s += sh[k][w];
        v[k] = s;
    }
    __syncthreads();
}

__device__ __forceinline__ int sgn(double x) { return (x > 0) - (x < 0); }   // arma::sign

// std::max / std::min as written in add_logp (first argument wins on NaN).
__device__ __forceinline__ double stdMax(double a, double b) { return (a < b) ? b : a; }
__device__ __forceinline__ double stdMin(double a, double b) { return (b < a) ? b : a; }
constexpr double LOG2 = 0.6931471805599453;   // std::log(2)

__device__ __forceinline__ double addLogp(double p1, double p2)   // UTIL.cpp add_logp
{
    p1 = -fabs(p1);
    p2 = -fabs(p2);
    const double maxp = stdMax(p1, p2);
    const double minp = stdMin(p1, p2);
    return maxp + log(1 + exp(minp - maxp));
}


// ---------------------------------------------------------------------------
// Kernel
// ---------------------------------------------------------------------------
struct KParams {
    // genotype source
    const unsigned char* packed; std::size_t bpv; const double* lut;
    const double* dense; std::size_t ld;
    // traits
    int N;
    const double* MU;            // nTraits x N
    // gpuPrecisionSPA fp32 only: float copies made at create(). MU32 holds
    // min(mu, 1 - mu), negated where mu > 1/2 (the stored form of terms32).
    const float* MU32; const float* XV32; const float* XX32;
    const float* MU32lo;         // min(mu, 1 - mu) - |MU32|, for the NAsigma sum
    const float* XV32lo; const float* XX32lo;   // XV - XV32, XXVX_inv - XX32 (carriers only)
    const double* XV;            // per trait: traitStride doubles, sample i at i*p
    const double* XX;            // per trait: traitStride doubles, column j at j*N
    const int*    pOfTrait;
    std::size_t   traitStride;
    // pairs
    const PairIn* in; int nPairs;
    double* scratch;             // blocks x 2N doubles
    double tol; int maxiter; int erfcMode;
    PairOut* out;
    int* counter;                // dynamicPairs: next pair index (zeroed per run)
    // ownSamples only (spaBody<.., OWN = true>): pair k's dosage table at
    // plut + 4k instead of the slot's, and trait t's sample mask (1 bit per
    // sample, bit i of word i >> 6) at masks + t * maskWords; a sample outside
    // the mask has dosage 0.
    const double*   plut;
    const uint64_t* masks;
    int             maskWords;
};

struct GenoCol {
    const unsigned char* col;    // packed row, or nullptr
    double4 L;                   // its dosage table
    const double* dcol;          // dense column, or nullptr
    const uint64_t* mk;          // OWN: the trait's sample mask
    template <bool OWN>
    __device__ __forceinline__ double dose(int i) const
    {
        if (OWN && !((mk[i >> 6] >> (i & 63)) & 1ull)) return 0.0;
        if (dcol) return dcol[i];
        const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
        return (c & 2u) ? ((c & 1u) ? L.w : L.z) : ((c & 1u) ? L.y : L.x);
    }
    float4 L32;                  // F32: the dosage table as floats
    template <bool OWN>
    __device__ __forceinline__ float dose32(int i) const
    {
        if (OWN && !((mk[i >> 6] >> (i & 63)) & 1ull)) return 0.0f;
        if (dcol) return (float)dcol[i];   // dense input: tests only
        const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
        return (c & 2u) ? ((c & 1u) ? L32.w : L32.z) : ((c & 1u) ? L32.y : L32.x);
    }
};

// Per-pair state the root / tail passes need.
struct Ctx {
    const double2* buf;   // (g~, mu), nEff entries
    int nEff;
    bool fast;
    double NAmu, NAsigma;
    // fp32 variant only (gpuPrecisionSPA: fp32): the same entries as float2
    // (s g~, s min(mu, 1 - mu)) with s = -1 where mu > 1/2, and the m1 used
    // for centring (that of the stored floats; see spaBody).
    const float2* buf32;
    double m1;
};

// ---------------------------------------------------------------------------
// fp32 variant (gpuPrecisionSPA: fp32). Same control flow, same block
// structure; the per-sample work of the root and tail passes (exp / log, the
// K', K'' and K sums) runs in float, and so do passes A and B (float copies
// of XV, XXVX_inv, mu; see spa_gpu.hpp "FP32") and the block reductions. Only
// the per-pair scalars stay double: the Newton step (t, K1, K2, prevJump), q,
// NAsigma, and the saddlepoint tail (w, v, Ztest, the erfc port, the log
// domain), which reads the float sums as doubles. Mathematically identical rewrites, chosen so that
// float cannot overflow and the big cancellations happen in double:
//
//   p(x)  = mu e^x / (1 - mu + mu e^x),  x = g~ t  (the tilted mean)
//   K'(t) - q    = sum g~ (p - mu) + [fast: NAsigma t] - (q - m1)
//   K''(t)       = sum g~^2 p (1 - p) + [fast: NAsigma]
//   K(t) - t m1  = sum h,  h = log(1 - mu + mu e^x) - mu x  (+ [fast: NAsigma t^2 / 2])
//   zeta q - K(zeta) = zeta (q - m1) - (K(zeta) - zeta m1)
//
// m1 = sum mu g~ is the double from pass B, so neither K' - q nor zeta q - K
// is formed as a difference of two large float numbers. Per term only
// exp(-|x|) <= 1 is formed (p, 1 - p, p - mu and h are written with it and
// expm1 / log1p), so e^{x} never overflows float (it would from x ~ 88). The
// stored mu is min(mu, 1 - mu) with g~ negated when mu > 1/2 (every term is
// invariant under mu -> 1 - mu, g~ -> -g~), so 1 - mu is formed without
// cancellation and a mu within 1e-8 of 1 is not lost; the flip is kept in the
// sign bit of the stored mu. Per-thread sums are Neumaier-compensated floats,
// merged across the block in float-float (blockReduceF) and read out in double.
//
// One fp64 overflow is part of SAIGE's control flow and is reproduced: in
// Korg's log(1 - mu + mu exp(g~ t)) the exp overflows to +inf for
// g~ t > log(DBL_MAX) = 709.78, Korg is then +inf and the tail is "not a
// saddle" (p = pno / 2, convergence withdrawn, Is.SPA false). That happens on
// real data (very rare markers, roots of several hundred). The fp32 tail pass
// returns H = +inf in exactly that case, so the route matches fp64; every
// other fp64 overflow (exp(-g~ t) in K', K'') only produces the limit value,
// which the rewritten terms give directly.
// ---------------------------------------------------------------------------
using spa_fp32::NSum;
using spa_fp32::blockReduceF;

// exp(x) overflows double for x above this (CUDA and glibc agree).
constexpr float kLogDblMax = 709.782712893384f;

// Per-term K' (centred), K'' and, when H, h, for one stored entry
// (s g~, sign(s) min(mu, 1 - mu)). With H, *ovf is set when the fp64 Korg term
// of the unflipped sample overflows (see above).
// t as a float pair th + tl (th = float(t), tl = float(t - th)), so g~ t is
// formed to float precision even when t is large: float(t) alone would put an
// error of up to |t| 2^-24 into every x.
struct T2 { float h, l; };
__device__ __forceinline__ T2 split(double t) { const float h = (float)t; return T2{h, (float)(t - (double)h)}; }

template <bool H>
__device__ __forceinline__ void terms32(float g, float ms, T2 t, float* k1c, float* k2, float* h, bool* ovf)
{
    const float m = fabsf(ms);
    const float a = 1.0f - m;     // m <= 1/2: no cancellation
    const float x = g * t.h + g * t.l;
    if (H) *ovf = (signbit(ms) ? -x : x) > kLogDblMax;
    float p, q, pm;
    if (x >= 0.0f) {
        const float e = expf(-x), em = expm1f(-x);
        const float den = a * e + m;                 // 1 - mu + mu e^x, times e^-x
        p  = m / den;
        q  = (a * e) / den;
        pm = -(m * a * em) / den;                    // p - mu
        if (H) *h = spa_fp32::hTerm(m, a, x, em);   // forms without cancellation: spa_fp32_terms.cuh
    } else {
        const float e = expf(x), em = expm1f(x);
        const float den = a + m * e;                 // 1 - mu + mu e^x
        p  = (m * e) / den;
        q  = a / den;
        pm = (m * a * em) / den;
        if (H) *h = spa_fp32::hTerm(m, a, x, 0.0f);
    }
    *k1c = g * pm;
    *k2  = (g * g) * (p * q);
}

// One pass: K1(t) - q and K2(t), fp32 per sample (the fp32 passK1K2).
__device__ void passK1K2_32(const Ctx& C, double t, double q, double (*sh)[NWARP], double* K1, double* K2)
{
    NSum s0, s1;
    const T2 tf = split(t);
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const float2 e = C.buf32[i];
        float k1c, k2, h;
        terms32<false>(e.x, e.y, tf, &k1c, &k2, &h, nullptr);
        s0.add(k1c);
        if (isfinite(k2)) s1.add(k2);
    }
    NSum r[2] = {s0, s1};
    blockReduceF<2, NWARP>(r, reinterpret_cast<float*>(&sh[0][0]));
    const double a[2] = {r[0].val(), r[1].val()};
    const double dq = q - C.m1;
    if (C.fast) { *K1 = a[0] + C.NAsigma * t - dq; *K2 = a[1] + C.NAsigma; }
    else        { *K1 = a[0] - dq;                 *K2 = a[1]; }
}

// One pass: H = K(zeta) - zeta m1 and K2(zeta), fp32 per sample.
__device__ void passK0K2_32(const Ctx& C, double t, double (*sh)[NWARP], double* H, double* K2)
{
    NSum s0, s1;
    bool inf = false;
    const T2 tf = split(t);
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const float2 e = C.buf32[i];
        float k1c, k2, h; bool ovf;
        terms32<true>(e.x, e.y, tf, &k1c, &k2, &h, &ovf);
        s0.add(h); inf |= ovf;
        if (isfinite(k2)) s1.add(k2);
    }
    inf = __syncthreads_or(inf);
    NSum r[2] = {s0, s1};
    blockReduceF<2, NWARP>(r, reinterpret_cast<float*>(&sh[0][0]));
    const double a[2] = {inf ? CUDART_INF : r[0].val(), r[1].val()};
    if (C.fast) { *H = a[0] + 0.5 * C.NAsigma * (t * t); *K2 = a[1] + C.NAsigma; }
    else        { *H = a[0];                             *K2 = a[1]; }
}

__device__ void passK1K2x2_32(const Ctx& C, bool e1, double t1, double q1, bool e2, double t2, double q2,
                              double (*sh)[NWARP], double* K1a, double* K2a, double* K1b, double* K2b)
{
    NSum s[4];
    const T2 tf1 = split(t1), tf2 = split(t2);
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const float2 e = C.buf32[i];
        float k1c, k2, h;
        if (e1) { terms32<false>(e.x, e.y, tf1, &k1c, &k2, &h, nullptr); s[0].add(k1c); if (isfinite(k2)) s[1].add(k2); }
        if (e2) { terms32<false>(e.x, e.y, tf2, &k1c, &k2, &h, nullptr); s[2].add(k1c); if (isfinite(k2)) s[3].add(k2); }
    }
    blockReduceF<4, NWARP>(s, reinterpret_cast<float*>(&sh[0][0]));
    const double a[4] = {s[0].val(), s[1].val(), s[2].val(), s[3].val()};
    if (e1) {
        const double dq = q1 - C.m1;
        if (C.fast) { *K1a = a[0] + C.NAsigma * t1 - dq; *K2a = a[1] + C.NAsigma; }
        else        { *K1a = a[0] - dq; *K2a = a[1]; }
    }
    if (e2) {
        const double dq = q2 - C.m1;
        if (C.fast) { *K1b = a[2] + C.NAsigma * t2 - dq; *K2b = a[3] + C.NAsigma; }
        else        { *K1b = a[2] - dq; *K2b = a[3]; }
    }
}

__device__ void passK0K2x2_32(const Ctx& C, double t1, double t2, double (*sh)[NWARP],
                              double* Ha, double* K2a, double* Hb, double* K2b)
{
    NSum s[4];
    bool inf1 = false, inf2 = false;
    const T2 tf1 = split(t1), tf2 = split(t2);
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const float2 e = C.buf32[i];
        float k1c, k2, h; bool ovf;
        terms32<true>(e.x, e.y, tf1, &k1c, &k2, &h, &ovf); s[0].add(h); inf1 |= ovf; if (isfinite(k2)) s[1].add(k2);
        terms32<true>(e.x, e.y, tf2, &k1c, &k2, &h, &ovf); s[2].add(h); inf2 |= ovf; if (isfinite(k2)) s[3].add(k2);
    }
    inf1 = __syncthreads_or(inf1);
    inf2 = __syncthreads_or(inf2);
    blockReduceF<4, NWARP>(s, reinterpret_cast<float*>(&sh[0][0]));
    const double a[4] = {inf1 ? CUDART_INF : s[0].val(), s[1].val(), inf2 ? CUDART_INF : s[2].val(), s[3].val()};
    if (C.fast) {
        *Ha = a[0] + 0.5 * C.NAsigma * (t1 * t1);  *K2a = a[1] + C.NAsigma;
        *Hb = a[2] + 0.5 * C.NAsigma * (t2 * t2);  *K2b = a[3] + C.NAsigma;
    } else {
        *Ha = a[0]; *K2a = a[1];
        *Hb = a[2]; *K2b = a[3];
    }
}

// One pass: K1(t) - q and K2(t). Per-term arithmetic associated as the
// Armadillo expressions in spa_binary.cpp evaluate it:
//   K1: (mu % g) / ((1 - mu) % exp(-g*t) + mu)
//   K2: ((1 - mu) % mu % (pow(g,2) % exp(-g*t))) / pow((1 - mu) % exp(-g*t) + mu, 2),
//       non-finite terms skipped (sum_arma1)
__device__ void passK1K2(const Ctx& C, double t, double q, double (*sh)[NWARP], double* K1, double* K2)
{
    double a[2] = {0.0, 0.0};
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const double2 e = C.buf[i];
        const double g = e.x, m = e.y;
        const double ex = exp(-g * t);
        const double d  = (1 - m) * ex + m;
        a[0] += (m * g) / d;
        const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
        if (isfinite(term)) a[1] += term;
    }
    blockReduce<2>(a, sh);
    if (C.fast) {
        const double temp3 = C.NAmu + C.NAsigma * t;
        *K1 = a[0] + temp3 - q;        // K1_adj_fast_Binom: sum + temp3 - q
        *K2 = a[1] + C.NAsigma;        // K2_fast_Binom
    } else {
        *K1 = a[0] - q;                // K1_adj_Binom
        *K2 = a[1];                    // K2_Binom
    }
}

// One pass: Korg(zeta) and K2(zeta), for Get_Saddle_Prob_*_Binom.
//   Korg: log(1 - mu + mu % exp(g*t))   (+ NAmu*t + 0.5*NAsigma*pow(t,2))
__device__ void passK0K2(const Ctx& C, double t, double (*sh)[NWARP], double* K0, double* K2)
{
    double a[2] = {0.0, 0.0};
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const double2 e = C.buf[i];
        const double g = e.x, m = e.y;
        a[0] += log(1 - m + m * exp(g * t));
        const double ex = exp(-g * t);
        const double d  = (1 - m) * ex + m;
        const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
        if (isfinite(term)) a[1] += term;
    }
    blockReduce<2>(a, sh);
    if (C.fast) {
        *K0 = a[0] + C.NAmu * t + 0.5 * C.NAsigma * (t * t);
        *K2 = a[1] + C.NAsigma;
    } else {
        *K0 = a[0];
        *K2 = a[1];
    }
}

// Newton tolerance at the step t -> tnew. fp64: tol. fp32: max(tol, kRelTolFp32
// max(|t|, |tnew|)) -- the float sums resolve K' to about 1e-7 of its terms,
// so a step is only known to about 1e-7 |t|; with a purely absolute tol a root
// at |t| ~ 1e4 (very rare markers, full-N variant) could never converge.
constexpr double kRelTolFp32 = 1.0 / 131072;   // 2^-17, ~64 float ulps
template <bool F32>
__device__ __forceinline__ double tolAt(double tol, double t, double tnew)
{
    if constexpr (F32) return fmax(tol, kRelTolFp32 * fmax(fabs(t), fabs(tnew)));
    else return tol;
}

// The pass of the chosen arithmetic (F32: gpuPrecisionSPA fp32). In fp32 the
// "K0" a tail pass returns is H = K(zeta) - zeta m1 (see passK0K2_32).
template <bool F32>
__device__ __forceinline__ void PK1K2(const Ctx& C, double t, double q, double (*sh)[NWARP], double* K1, double* K2)
{
    if constexpr (F32) passK1K2_32(C, t, q, sh, K1, K2); else passK1K2(C, t, q, sh, K1, K2);
}

// getroot_K1_Binom (full) / getroot_K1_fast_Binom (fast), init 0, with the
// two variants' own tests. Every thread runs the same scalar code on the
// same block-wide sums, so there is no divergence and no broadcast.
template <bool F32>
__device__ void getroot(const Ctx& C, double q, double gpos, double gneg, double tol, int maxiter,
                        double (*sh)[NWARP], double* root, int* niter, int* conv, int* reason)
{
    *reason = RS_OK;
    if (q >= gpos || q <= gneg) { *root = CUDART_INF; *niter = 0; *conv = 1; return; }
    double t = 0.0;
    double K1, K2;
    PK1K2<F32>(C, t, q, sh, &K1, &K2);         // K1_eval at init; K2 at the same t for the loop top
    double prevJump = CUDART_INF;
    int rep = 1;
    int c = 1;
    while (rep <= maxiter) {
        // K2_eval = K2(t): already in K2 (evaluated with the current K1 at this t)
        if (!C.fast && (!isfinite(K2) || fabs(K2) < 1e-15)) { c = 0; *reason = RS_K2_GUARD; break; }
        double tnew = t - K1 / K2;
        if (C.fast ? isnan(tnew) : !isfinite(tnew)) { c = 0; *reason = RS_TNEW; break; }
        const double tolE = tolAt<F32>(tol, t, tnew);
        if (fabs(tnew - t) < tolE) { c = 1; break; }
        if (rep == maxiter) { c = 0; *reason = RS_MAXITER; break; }
        double newK1, newK2;
        PK1K2<F32>(C, tnew, q, sh, &newK1, &newK2);
        const bool flipped = C.fast ? ((K1 * newK1) < 0) : (sgn(K1) != sgn(newK1));
        if (flipped) {
            if (fabs(tnew - t) > (prevJump - tolE)) {
                tnew = t + (double)sgn(newK1 - K1) * prevJump / 2;
                PK1K2<F32>(C, tnew, q, sh, &newK1, &newK2);
                prevJump = prevJump / 2;
            } else {
                prevJump = fabs(tnew - t);
            }
        }
        rep = rep + 1;
        t = tnew;
        K1 = newK1;
        K2 = newK2;
    }
    *root = t; *niter = rep; *conv = c;
}

// ---------------------------------------------------------------------------
// Fused passes (CreateArgs::fusedRoots). One pass over the stored pairs
// serves both Newton solves (or both tails) of a pair. Per root the per-thread
// accumulation is the same sequence of additions as passK1K2 / passK0K2, and
// blockReduce adds the warp partials in the same order, so each root's sums
// are bit-identical to the unfused pass; e1 / e2 are block-uniform, so there is
// no divergence, and a root that is not wanted this pass costs no arithmetic.
// ---------------------------------------------------------------------------
__device__ void passK1K2x2(const Ctx& C, bool e1, double t1, double q1, bool e2, double t2, double q2,
                           double (*sh)[NWARP], double* K1a, double* K2a, double* K1b, double* K2b)
{
    double a[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const double2 e = C.buf[i];
        const double g = e.x, m = e.y;
        if (e1) {
            const double ex = exp(-g * t1);
            const double d  = (1 - m) * ex + m;
            a[0] += (m * g) / d;
            const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
            if (isfinite(term)) a[1] += term;
        }
        if (e2) {
            const double ex = exp(-g * t2);
            const double d  = (1 - m) * ex + m;
            a[2] += (m * g) / d;
            const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
            if (isfinite(term)) a[3] += term;
        }
    }
    blockReduce<4>(a, sh);
    if (e1) {
        if (C.fast) { const double temp3 = C.NAmu + C.NAsigma * t1; *K1a = a[0] + temp3 - q1; *K2a = a[1] + C.NAsigma; }
        else        { *K1a = a[0] - q1; *K2a = a[1]; }
    }
    if (e2) {
        if (C.fast) { const double temp3 = C.NAmu + C.NAsigma * t2; *K1b = a[2] + temp3 - q2; *K2b = a[3] + C.NAsigma; }
        else        { *K1b = a[2] - q2; *K2b = a[3]; }
    }
}

__device__ void passK0K2x2(const Ctx& C, double t1, double t2, double (*sh)[NWARP],
                           double* K0a, double* K2a, double* K0b, double* K2b)
{
    double a[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = threadIdx.x; i < C.nEff; i += NT) {
        const double2 e = C.buf[i];
        const double g = e.x, m = e.y;
        {
            a[0] += log(1 - m + m * exp(g * t1));
            const double ex = exp(-g * t1);
            const double d  = (1 - m) * ex + m;
            const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
            if (isfinite(term)) a[1] += term;
        }
        {
            a[2] += log(1 - m + m * exp(g * t2));
            const double ex = exp(-g * t2);
            const double d  = (1 - m) * ex + m;
            const double term = ((1 - m) * m * (g * g * ex)) / (d * d);
            if (isfinite(term)) a[3] += term;
        }
    }
    blockReduce<4>(a, sh);
    if (C.fast) {
        *K0a = a[0] + C.NAmu * t1 + 0.5 * C.NAsigma * (t1 * t1);  *K2a = a[1] + C.NAsigma;
        *K0b = a[2] + C.NAmu * t2 + 0.5 * C.NAsigma * (t2 * t2);  *K2b = a[3] + C.NAsigma;
    } else {
        *K0a = a[0]; *K2a = a[1];
        *K0b = a[2]; *K2b = a[3];
    }
}

// Both Newton solves of a pair in lockstep: getroot()'s control flow kept per
// root (its own t, K1, K2, prevJump, rep, failure reason), every pass shared.
// A root finishes when getroot() would have; the other continues alone, its
// passes then carrying one root's arithmetic.
template <bool F32>
__device__ __forceinline__ void PK1K2x2(const Ctx& C, bool e1, double t1, double q1, bool e2, double t2, double q2,
                                        double (*sh)[NWARP], double* K1a, double* K2a, double* K1b, double* K2b)
{
    if constexpr (F32) passK1K2x2_32(C, e1, t1, q1, e2, t2, q2, sh, K1a, K2a, K1b, K2b);
    else               passK1K2x2(C, e1, t1, q1, e2, t2, q2, sh, K1a, K2a, K1b, K2b);
}

struct RootSt {
    double q, t, K1, K2, prevJump, tn, nK1, nK2;
    int rep, c, reason;
    bool live, eval, inf;
};

template <bool F32>
__device__ void getroot2(const Ctx& C, double q1, double q2, double gpos, double gneg, double tol, int maxiter,
                         double (*sh)[NWARP],
                         double* root1, int* niter1, int* conv1, int* reason1,
                         double* root2, int* niter2, int* conv2, int* reason2)
{
    RootSt R[2];
    R[0].q = q1; R[1].q = q2;
    #pragma unroll
    for (int r = 0; r < 2; ++r) {
        RootSt& S = R[r];
        S.reason = RS_OK; S.t = 0.0; S.K1 = 0.0; S.K2 = 0.0; S.prevJump = CUDART_INF;
        S.tn = 0.0; S.nK1 = 0.0; S.nK2 = 0.0; S.rep = 1; S.c = 1; S.eval = false;
        S.inf  = (S.q >= gpos || S.q <= gneg);
        S.live = !S.inf;
    }
    if (R[0].live || R[1].live)
        PK1K2x2<F32>(C, R[0].live, 0.0, q1, R[1].live, 0.0, q2, sh, &R[0].K1, &R[0].K2, &R[1].K1, &R[1].K2);
    while (R[0].live || R[1].live) {
        #pragma unroll
        for (int r = 0; r < 2; ++r) {
            RootSt& S = R[r];
            S.eval = false;
            if (!S.live) continue;
            if (!C.fast && (!isfinite(S.K2) || fabs(S.K2) < 1e-15)) { S.c = 0; S.reason = RS_K2_GUARD; S.live = false; continue; }
            const double tnew = S.t - S.K1 / S.K2;
            if (C.fast ? isnan(tnew) : !isfinite(tnew)) { S.c = 0; S.reason = RS_TNEW; S.live = false; continue; }
            if (fabs(tnew - S.t) < tolAt<F32>(tol, S.t, tnew)) { S.c = 1; S.live = false; continue; }
            if (S.rep == maxiter) { S.c = 0; S.reason = RS_MAXITER; S.live = false; continue; }
            S.tn = tnew; S.eval = true;
        }
        if (!(R[0].eval || R[1].eval)) continue;
        PK1K2x2<F32>(C, R[0].eval, R[0].tn, q1, R[1].eval, R[1].tn, q2, sh, &R[0].nK1, &R[0].nK2, &R[1].nK1, &R[1].nK2);
        bool re0 = false, re1 = false;
        #pragma unroll
        for (int r = 0; r < 2; ++r) {
            RootSt& S = R[r];
            if (!S.eval) continue;
            const bool flipped = C.fast ? ((S.K1 * S.nK1) < 0) : (sgn(S.K1) != sgn(S.nK1));
            if (flipped) {
                if (fabs(S.tn - S.t) > (S.prevJump - tolAt<F32>(tol, S.t, S.tn))) {
                    S.tn = S.t + (double)sgn(S.nK1 - S.K1) * S.prevJump / 2;
                    if (r == 0) re0 = true; else re1 = true;
                } else {
                    S.prevJump = fabs(S.tn - S.t);
                }
            }
        }
        if (re0 || re1) {
            PK1K2x2<F32>(C, re0, R[0].tn, q1, re1, R[1].tn, q2, sh, &R[0].nK1, &R[0].nK2, &R[1].nK1, &R[1].nK2);
            if (re0) R[0].prevJump = R[0].prevJump / 2;
            if (re1) R[1].prevJump = R[1].prevJump / 2;
        }
        #pragma unroll
        for (int r = 0; r < 2; ++r) {
            RootSt& S = R[r];
            if (!S.eval) continue;
            S.rep = S.rep + 1;
            S.t = S.tn; S.K1 = S.nK1; S.K2 = S.nK2;
        }
    }
    *reason1 = R[0].reason; *reason2 = R[1].reason;
    if (R[0].inf) { *root1 = CUDART_INF; *niter1 = 0; *conv1 = 1; } else { *root1 = R[0].t; *niter1 = R[0].rep; *conv1 = R[0].c; }
    if (R[1].inf) { *root2 = CUDART_INF; *niter2 = 0; *conv2 = 1; } else { *root2 = R[1].t; *niter2 = R[1].rep; *conv2 = R[1].c; }
}

// The scalar tail of Get_Saddle_Prob_*_Binom once Korg(zeta) and K2(zeta) are
// known; shared by the unfused and the fused path.
// temp1 = zeta q - Korg(zeta); k1fin = isfinite(Korg(zeta)).
__device__ double saddleTailT(double zeta, double temp1, bool k1fin, double k2, int logp, int erfcMode, int* isSaddle)
{
    *isSaddle = 0;
    bool flagrun = false;
    double w = 0.0, v = 0.0;
    if (k1fin && isfinite(k2) && temp1 >= 0 && k2 >= 0) {
        w = (double)sgn(zeta) * sqrt(2 * temp1);
        v = zeta * sqrt(k2);
        if (w != 0) flagrun = true;
    }
    if (!flagrun) return logp ? -CUDART_INF : 0.0;
    const double Ztest = w + (1 / w) * log(v / w);
    *isSaddle = 1;
    double pval0;
    if (Ztest > 0) {
        pval0 = phiUpper(Ztest, erfcMode);
        if (logp) pval0 = log(pval0);
        return pval0;
    } else {
        pval0 = phiLower(Ztest, erfcMode);
        if (logp) pval0 = log(pval0);
        return -pval0;
    }
}
__device__ __forceinline__ double saddleTail(double zeta, double q, double k1, double k2, int logp, int erfcMode, int* isSaddle)
{
    return saddleTailT(zeta, zeta * q - k1, isfinite(k1), k2, logp, erfcMode, isSaddle);
}
// fp32: H = K(zeta) - zeta m1, so zeta q - K(zeta) = zeta (q - m1) - H in double.
__device__ __forceinline__ double saddleTail32(const Ctx& C, double zeta, double q, double H, double k2, int logp,
                                               int erfcMode, int* isSaddle)
{
    return saddleTailT(zeta, zeta * (q - C.m1) - H, isfinite(H), k2, logp, erfcMode, isSaddle);
}

// Get_Saddle_Prob_Binom / Get_Saddle_Prob_fast_Binom.
template <bool F32>
__device__ double saddle(const Ctx& C, double zeta, double q, int logp, int erfcMode,
                         double (*sh)[NWARP], int* isSaddle)
{
    double k1, k2;
    if constexpr (F32) {
        passK0K2_32(C, zeta, sh, &k1, &k2);
        return saddleTail32(C, zeta, q, k1, k2, logp, erfcMode, isSaddle);
    } else {
        passK0K2(C, zeta, sh, &k1, &k2);
        return saddleTail(zeta, q, k1, k2, logp, erfcMode, isSaddle);
    }
}

// Both tails of a pair from one pass (fusedRoots).
template <bool F32>
__device__ void saddle2(const Ctx& C, double z1, double q1, double z2, double q2, int logp, int erfcMode,
                        double (*sh)[NWARP], double* p1, int* s1, double* p2, int* s2)
{
    double k1a, k2a, k1b, k2b;
    if constexpr (F32) {
        passK0K2x2_32(C, z1, z2, sh, &k1a, &k2a, &k1b, &k2b);
        *p1 = saddleTail32(C, z1, q1, k1a, k2a, logp, erfcMode, s1);
        *p2 = saddleTail32(C, z2, q2, k1b, k2b, logp, erfcMode, s2);
    } else {
        passK0K2x2(C, z1, z2, sh, &k1a, &k2a, &k1b, &k2b);
        *p1 = saddleTail(z1, q1, k1a, k2a, logp, erfcMode, s1);
        *p2 = saddleTail(z2, q2, k1b, k2b, logp, erfcMode, s2);
    }
}

template <bool FUSED, bool DYN, bool OWN, bool F32>
__device__ __forceinline__ void spaBody(const KParams& P)
{
    __shared__ double sh[NACC][NWARP];
    __shared__ int shK;
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const int N = P.N;
    double2* buf = reinterpret_cast<double2*>(P.scratch) + (std::size_t)blockIdx.x * N;
    float2* buf32 = reinterpret_cast<float2*>(buf);   // F32: the same region, 8 bytes per entry
    // this warp's contiguous share of the samples (passes A and B)
    const int segLo = (int)(((long long)N * warp) / NWARP);
    const int segHi = (int)(((long long)N * (warp + 1)) / NWARP);

    for (int k = blockIdx.x; ; k += gridDim.x) {
        if (DYN) {
            // next pair from the shared counter; every thread reads shK before
            // the end-of-pair barrier, so the next write cannot race it
            if (tid == 0) shK = atomicAdd(P.counter, 1);
            __syncthreads();
            k = shK;
        }
        if (k >= P.nPairs) break;
        const PairIn pin = P.in[k];
        GenoCol G;
        if (P.dense) {
            G.dcol = P.dense + (std::size_t)pin.slot * P.ld; G.col = nullptr;
            G.L = make_double4(0, 0, 0, 0);
        } else {
            G.dcol = nullptr;
            G.col = P.packed + (std::size_t)pin.slot * P.bpv;
            const double2* d = reinterpret_cast<const double2*>(OWN ? P.plut : P.lut) +
                               2 * (std::size_t)(OWN ? k : pin.slot);
            const double2 a = d[0], b = d[1];
            G.L = make_double4(a.x, a.y, b.x, b.y);
        }
        if constexpr (F32) G.L32 = make_float4((float)G.L.x, (float)G.L.y, (float)G.L.z, (float)G.L.w);
        const int t = pin.trait;
        G.mk = OWN ? P.masks + (std::size_t)t * P.maskWords : nullptr;
        const int p = P.pOfTrait[t];
        const double* mu = P.MU + (std::size_t)t * N;
        const double* xv = P.XV + (std::size_t)t * P.traitStride;
        const double* xx = P.XX + (std::size_t)t * P.traitStride;

        int nnz; bool fast;
        double gpos, gneg, m1, NAmu, NAsigma, cen;   // cen: centring of the fp32 passes
        if constexpr (F32) {
            // ---- gpuPrecisionSPA fp32: passes A and B in float (float copies of
            //      XV, XXVX_inv and mu; compensated sums; float block reductions).
            //      Same carrier set, same compaction, same outputs as below.
            const float* mu32 = P.MU32 + (std::size_t)t * N;
            const float* mu32lo = P.MU32lo + (std::size_t)t * N;
            const float* xv32 = P.XV32 + (std::size_t)t * P.traitStride;
            const float* xx32 = P.XX32 + (std::size_t)t * P.traitStride;
            const float* xv32lo = P.XV32lo + (std::size_t)t * P.traitStride;
            const float* xx32lo = P.XX32lo + (std::size_t)t * P.traitStride;
            float* shf = reinterpret_cast<float*>(&sh[0][0]);
            __shared__ int shC[NWARP];
            // pass A: b = XV g over carriers, carrier count
            NSum accA[PMAX];
            int cnt = 0;
            for (int i = segLo + lane; i < segHi; i += 32) {
                const float g = G.template dose32<OWN>(i);
                if (g != 0.0f) {
                    const float* x = xv32 + (std::size_t)i * p;
                    const float* xl = xv32lo + (std::size_t)i * p;
                    #pragma unroll
                    for (int j = 0; j < PMAX; ++j) if (j < p) { accA[j].addProd(x[j], g); accA[j].c += xl[j] * g; }
                    ++cnt;
                }
            }
            cnt = spa_fp32::warpSumI(cnt);
            if (lane == 0) shC[warp] = cnt;
            blockReduceF<PMAX, NWARP>(accA, shf);   // its barriers publish shC
            int wbase = 0;
            nnz = 0;
            for (int w = 0; w < NWARP; ++w) { if (w < warp) wbase += shC[w]; nnz += shC[w]; }
            // the CPU's rule |{g==0}| / N >= 0.5, in integers
            fast = (pin.fast > 0) || (pin.fast < 0 && 2LL * (N - nnz) >= (long long)N);
            float b[PMAX], bl[PMAX];   // b = b + bl (float-float), for the carriers' g~
            #pragma unroll
            for (int j = 0; j < PMAX; ++j) { b[j] = accA[j].s + accA[j].c; bl[j] = accA[j].c - (b[j] - accA[j].s); }

            // pass B: g~ = g - XXVX_inv b; gpos, gneg, m1 = sum mu g~, the fast
            // variant's carrier sums; store (s g~, s min(mu, 1 - mu)). m1 is summed
            // from exactly the stored floats (error-free products), so it is also
            // the centring of the root and tail passes.
            NSum sB[5];   // gpos, gneg, m1, sum_c mu g~, sum_c mu(1-mu) g~^2
            int base = wbase;
            for (int i0 = segLo; i0 < segHi; i0 += 32) {
                const int i = i0 + lane;
                const bool valid = i < segHi;
                float g = 0.0f, v = 0.0f, ms = 0.0f;
                if (valid) {
                    g = G.template dose32<OWN>(i);
                    float vl = 0.0f;
                    if (g != 0.0f) {
                        spa_fp32::carrierGt<PMAX>(g, xx32 + i, xx32lo + i, (std::size_t)N, b, bl, p, &v, &vl);
                    } else {
                        float proj = 0.0f;
                        #pragma unroll
                        for (int j = 0; j < PMAX; ++j) if (j < p) proj += xx32[(std::size_t)j * N + i] * b[j];
                        v = g - proj;
                    }
                    ms = mu32[i];
                    const float m = fabsf(ms);
                    if (v > 0) sB[0].add(v); else if (v < 0) sB[1].add(v);
                    // mu v with mu = m, or 1 - m where flipped (sign bit of ms)
                    const bool flip = signbit(ms);
                    if (flip) { sB[2].add(v); sB[2].addProd(-m, v); } else sB[2].addProd(m, v);
                    if (g != 0.0f) {
                        if (flip) { sB[3].add(v); sB[3].addProd(-m, v); } else sB[3].addProd(m, v);
                        spa_fp32::addVarTerm(sB[4], m, mu32lo[i], v, vl);
                    }
                }
                const float2 e32 = make_float2(signbit(ms) ? -v : v, ms);
                if (!fast) {
                    if (valid) buf32[i] = e32;
                } else {
                    const bool carrier = valid && (g != 0.0f);
                    const unsigned mask = __ballot_sync(0xffffffffu, carrier);
                    if (carrier) buf32[base + __popc(mask & ((1u << lane) - 1u))] = e32;
                    base += __popc(mask);
                }
            }
            blockReduceF<5, NWARP>(sB, shf);   // its barriers also complete buf32
            gpos = sB[0].val(); gneg = sB[1].val(); m1 = sB[2].val();
            NAmu    = fast ? (m1 - sB[3].val()) : 0.0;
            NAsigma = fast ? (pin.var2 - sB[4].val()) : 0.0;
            cen = m1;
        } else {
            // ---- pass A: b = XV g over carriers (getadjGFast's loop), carrier count
            double acc[NACC];
            #pragma unroll
            for (int j = 0; j < NACC; ++j) acc[j] = 0.0;
            for (int i = segLo + lane; i < segHi; i += 32) {
                const double g = G.template dose<OWN>(i);
                if (g != 0.0) {
                    const double* x = xv + (std::size_t)i * p;
                    #pragma unroll
                    for (int j = 0; j < PMAX; ++j) if (j < p) acc[j] += x[j] * g;
                    acc[PMAX] += 1.0;
                }
            }
            // reduce; keep this warp's exclusive carrier prefix for the compaction
            {
                #pragma unroll
                for (int j = 0; j < NACC; ++j) acc[j] = warpSum(acc[j]);
                if (lane == 0) {
                    #pragma unroll
                    for (int j = 0; j < NACC; ++j) sh[j][warp] = acc[j];
                }
            }
            __syncthreads();
            int wbase = 0;
            {
                #pragma unroll
                for (int j = 0; j < NACC; ++j) {
                    double s = 0.0;
                    #pragma unroll
                    for (int w = 0; w < NWARP; ++w) s += sh[j][w];
                    acc[j] = s;
                }
                for (int w = 0; w < warp; ++w) wbase += (int)sh[PMAX][w];
            }
            __syncthreads();
            nnz = (int)acc[PMAX];
            // the CPU's rule: p_iIndexComVecSize = double(|{g==0}|) / m_n >= 0.5
            fast = (pin.fast > 0) || (pin.fast < 0 && ((double)(N - nnz) / (double)N) >= 0.5);
            double b[PMAX];
            #pragma unroll
            for (int j = 0; j < PMAX; ++j) b[j] = acc[j];

            // ---- pass B: g~, its positive / negative sums, m1, the fast variant's
            //      carrier sums; store (g~, mu)
            double s[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // gpos, gneg, m1, sum_c g~ mu, sum_c mu(1-mu) g~^2
            int base = wbase;
            for (int i0 = segLo; i0 < segHi; i0 += 32) {
                const int i = i0 + lane;
                const bool valid = i < segHi;
                double g = 0.0, v = 0.0, m = 0.0;
                if (valid) {
                    g = G.template dose<OWN>(i);
                    double proj = 0.0;
                    #pragma unroll
                    for (int j = 0; j < PMAX; ++j) if (j < p) proj += xx[(std::size_t)j * N + i] * b[j];
                    v = g - proj;
                    m = mu[i];
                    if (v > 0) s[0] += v; else if (v < 0) s[1] += v;
                    s[2] += m * v;
                    if (g != 0.0) { s[3] += v * m; s[4] += m * (1 - m) * (v * v); }
                }
                if (!fast) {
                    if (valid) buf[i] = make_double2(v, m);
                } else {
                    const bool carrier = valid && (g != 0.0);
                    const unsigned mask = __ballot_sync(0xffffffffu, carrier);
                    if (carrier) buf[base + __popc(mask & ((1u << lane) - 1u))] = make_double2(v, m);
                    base += __popc(mask);
                }
            }
            blockReduce<5>(s, sh);
            __syncthreads();   // buf complete before the passes read it
            gpos = s[0]; gneg = s[1]; m1 = s[2];
            NAmu    = fast ? (m1 - s[3]) : 0.0;           // NAmu = m1 - dot(gNB, muNB)
            NAsigma = fast ? (pin.var2 - s[4]) : 0.0;     // NAsigma = var2 - sum(muNB % (1-muNB) % pow(gNB,2))
            cen = m1;
        }

        Ctx C;
        C.buf = buf; C.nEff = fast ? nnz : N; C.fast = fast;
        C.NAmu = NAmu; C.NAsigma = NAsigma;
        C.buf32 = buf32; C.m1 = cen;

        // q, qinv as getMarkerPval forms them for a binary trait
        const double q = pin.Tstat / sqrt(pin.var1 / pin.var2) + m1;
        double qinv;
        if ((q - m1) > 0)       qinv = -1 * fabs(q - m1) + m1;
        else if ((q - m1) == 0) qinv = m1;
        else                    qinv = fabs(q - m1) + m1;

        unsigned status = (fast ? ST_FAST : 0u) | (pin.logp ? ST_LOGP : 0u);
        double r1, r2; int n1, n2, c1, c2, rs1, rs2;
        if (FUSED) {
            getroot2<F32>(C, q, qinv, gpos, gneg, P.tol, P.maxiter, sh, &r1, &n1, &c1, &rs1, &r2, &n2, &c2, &rs2);
        } else {
            getroot<F32>(C, q,    gpos, gneg, P.tol, P.maxiter, sh, &r1, &n1, &c1, &rs1);
            getroot<F32>(C, qinv, gpos, gneg, P.tol, P.maxiter, sh, &r2, &n2, &c2, &rs2);
        }
        if (n1 == 0 && isinf(r1)) status |= ST_ROOT1_INF;
        if (n2 == 0 && isinf(r2)) status |= ST_ROOT2_INF;
        if (!c1) status |= ST_ROOT1_FAIL;
        if (!c2) status |= ST_ROOT2_FAIL;

        double pv, p1 = 0.0, p2 = 0.0; int conv, s1 = -1, s2 = -1;
        if (c1 && c2) {
            // spa.cpp SPA / SPA_fast: a tail that is not a saddle withdraws convergence
            if (FUSED) {
                saddle2<F32>(C, r1, q, r2, qinv, pin.logp, P.erfcMode, sh, &p1, &s1, &p2, &s2);
            } else {
                p1 = saddle<F32>(C, r1, q,    pin.logp, P.erfcMode, sh, &s1);
                p2 = saddle<F32>(C, r2, qinv, pin.logp, P.erfcMode, sh, &s2);
            }
            conv = 1;
            if (!s1) { conv = 0; status |= ST_SADDLE1_FAIL; p1 = pin.logp ? pin.pno - LOG2 : pin.pno / 2; }
            if (!s2) { conv = 0; status |= ST_SADDLE2_FAIL; p2 = pin.logp ? pin.pno - LOG2 : pin.pno / 2; }
            pv = pin.logp ? addLogp(p1, p2) : fabs(p1) + fabs(p2);
        } else {
            pv = pin.pno; conv = 0; status |= ST_PNO;
        }
        if (tid == 0) {
            PairOut o;
            o.pval = pv; o.conv = conv; o.status = status;
            o.reason1 = rs1; o.reason2 = rs2;
            o.niter1 = n1; o.niter2 = n2; o.root1 = r1; o.root2 = r2;
            o.p1 = p1; o.p2 = p2; o.m1 = m1; o.q = q; o.qinv = qinv;
            o.nnz = nnz; o.s1 = s1; o.s2 = s2;
            P.out[k] = o;
        }
        __syncthreads();   // buf and sh are reused by the next pair
    }
}

// minBlocksPerSM 0 (default): the kernel as it was before the hint existed,
// __launch_bounds__(NT) only -- ptxas picks 128 registers for the unfused
// kernel (2 blocks of 8 warps per SM on a V100) and 194 for the fused one.
// 1..4: __launch_bounds__(NT, MINB). An explicit 1 is NOT the same as no hint:
// ptxas then takes 162 registers for the unfused kernel (1 block per SM).
// 2 / 3 / 4 cap the registers at 128 / 80 / 64 and spill the rest.
template <bool FUSED, bool DYN, bool OWN, bool F32>
__global__ void __launch_bounds__(NT) spaKernel0(const KParams P) { spaBody<FUSED, DYN, OWN, F32>(P); }
template <bool FUSED, bool DYN, int MINB, bool OWN, bool F32>
__global__ void __launch_bounds__(NT, MINB) spaKernel(const KParams P) { spaBody<FUSED, DYN, OWN, F32>(P); }

template <bool F, bool D, bool O, bool F32>
void launchSpaM(int minb, int grid, cudaStream_t st, const KParams& P)
{
    switch (minb) {
        case 1:  spaKernel<F, D, 1, O, F32><<<grid, NT, 0, st>>>(P); break;
        case 2:  spaKernel<F, D, 2, O, F32><<<grid, NT, 0, st>>>(P); break;
        case 3:  spaKernel<F, D, 3, O, F32><<<grid, NT, 0, st>>>(P); break;
        case 4:  spaKernel<F, D, 4, O, F32><<<grid, NT, 0, st>>>(P); break;
        default: spaKernel0<F, D, O, F32><<<grid, NT, 0, st>>>(P); break;
    }
}
template <bool O, bool F32>
void launchSpaO(int fused, int dyn, int minb, int grid, cudaStream_t st, const KParams& P)
{
    if (fused) { if (dyn) launchSpaM<true, true, O, F32>(minb, grid, st, P);  else launchSpaM<true, false, O, F32>(minb, grid, st, P); }
    else       { if (dyn) launchSpaM<false, true, O, F32>(minb, grid, st, P); else launchSpaM<false, false, O, F32>(minb, grid, st, P); }
}
template <bool F32>
void launchSpa(int fused, int dyn, int minb, int own, int grid, cudaStream_t st, const KParams& P)
{
    if (own) launchSpaO<true, F32>(fused, dyn, minb, grid, st, P);
    else     launchSpaO<false, F32>(fused, dyn, minb, grid, st, P);
}

__global__ void erfcDebugKernel(const double* z, int n, double* a, double* b, double* c)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    a[i] = erfc53::erfImp53s(z[i], true);
    b[i] = erfc53::erfImp53(z[i], true);
    c[i] = erfc(z[i]);
}

}  // namespace

// ---------------------------------------------------------------------------

struct Spa {
    saige::gpu2::Prec prec = saige::gpu2::Prec::FP64;   // CreateArgs::precision
    int N = 0, nTraits = 0, maxPairs = 0, blocks = 256, maxiter = 1000, erfcMode = 1;
    int fused = 0, dyn = 0, minb = 0;
    int* dCounter = nullptr;
    // ownSamples: per-pair dosage tables (pinned + device) and per-trait masks
    int own = 0, maskWords = 0;
    double*   hPLut = nullptr;
    double*   dPLut = nullptr;
    uint64_t* dMask = nullptr;
    double tol = 0.0;
    std::size_t traitStride = 0;
    PairIn*  hIn  = nullptr;
    PairOut* hOut = nullptr;
    PairIn*  dIn  = nullptr;
    PairOut* dOut = nullptr;
    double* dMu = nullptr; double* dXV = nullptr; double* dXX = nullptr;
    float*  dMu32 = nullptr; float* dXV32 = nullptr; float* dXX32 = nullptr;   // FP32 only (KParams)
    float*  dMu32lo = nullptr; float* dXV32lo = nullptr; float* dXX32lo = nullptr;
    int*    dP  = nullptr;
    double* dScratch = nullptr;
    // library-owned genotype buffers (uploadPacked / uploadDense)
    unsigned char* dGPk = nullptr; double* dGLut = nullptr; double* dGDense = nullptr;
    std::size_t capPk = 0, capLut = 0, capDense = 0;
    cudaStream_t st = nullptr;
    cudaEvent_t ev[4] = {nullptr, nullptr, nullptr, nullptr};
    double tKernel = 0.0, tH2D = 0.0, tD2H = 0.0;
    long long nPairsDone = 0;
    std::size_t devBytes = 0;
};

bool supports(saige::gpu2::Prec t_p)
{
    switch (t_p) {
        case saige::gpu2::Prec::FP64: return true;
        case saige::gpu2::Prec::FP32: return true;    // spaBody<.., F32 = true>
        case saige::gpu2::Prec::INT8: return false;   // not a mode of this stage
    }
    return false;
}

Spa* create(const CreateArgs& a)
{
    g_lastErr.clear();
    if (a.N <= 0 || a.nTraits <= 0 || a.traits == nullptr || a.maxPairs <= 0) {
        g_lastErr = "create: bad arguments"; return nullptr;
    }
    // ---- precision dispatch ----
    if (!supports(a.precision)) {
        g_lastErr = std::string("SPA precision ") + saige::gpu2::precName(a.precision) + " is not implemented yet";
        return nullptr;
    }
    if (cudaSetDevice(a.device) != cudaSuccess) { g_lastErr = "cudaSetDevice failed"; return nullptr; }
    int pMax = 0;
    for (int t = 0; t < a.nTraits; ++t) {
        const TraitArgs& T = a.traits[t];
        if (T.p <= 0 || T.p > PMAX || !T.mu || !T.XV || !T.XXVX_inv) {
            g_lastErr = "create: trait " + std::to_string(t) + " malformed (p, mu, XV, XXVX_inv)";
            return nullptr;
        }
        if (T.p > pMax) pMax = T.p;
    }
    Spa* s = new Spa();
    s->prec = a.precision;
    s->N = a.N; s->nTraits = a.nTraits; s->maxPairs = a.maxPairs;
    s->blocks = a.blocks > 0 ? a.blocks : 256;
    s->maxiter = a.maxiter; s->tol = a.tol; s->erfcMode = a.erfcMode;
    s->fused = a.fusedRoots ? 1 : 0; s->dyn = a.dynamicPairs ? 1 : 0;
    s->minb = (a.minBlocksPerSM >= 1 && a.minBlocksPerSM <= 4) ? a.minBlocksPerSM : 0;
    s->traitStride = (std::size_t)a.N * pMax;

    auto fail = [&](const char* what) -> Spa* {
        if (g_lastErr.empty()) g_lastErr = what;
        destroy(s); return nullptr;
    };
    std::size_t db = 0;
    auto dev = [&](void** p, std::size_t n) {
        cudaError_t e = cudaMalloc(p, n);
        if (e != cudaSuccess) { *p = nullptr; g_lastErr = std::string("cudaMalloc: ") + cudaGetErrorString(e); return false; }
        db += n; return true;
    };
    if (cudaHostAlloc((void**)&s->hIn,  (std::size_t)a.maxPairs * sizeof(PairIn),  cudaHostAllocDefault) != cudaSuccess) return fail("cudaHostAlloc in");
    if (cudaHostAlloc((void**)&s->hOut, (std::size_t)a.maxPairs * sizeof(PairOut), cudaHostAllocDefault) != cudaSuccess) return fail("cudaHostAlloc out");
    if (!dev((void**)&s->dIn,  (std::size_t)a.maxPairs * sizeof(PairIn)))  return fail("");
    if (!dev((void**)&s->dOut, (std::size_t)a.maxPairs * sizeof(PairOut))) return fail("");
    if (!dev((void**)&s->dMu, (std::size_t)a.nTraits * a.N * sizeof(double))) return fail("");
    if (!dev((void**)&s->dXV, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail("");
    if (!dev((void**)&s->dXX, (std::size_t)a.nTraits * s->traitStride * sizeof(double))) return fail("");
    if (!dev((void**)&s->dP,  (std::size_t)a.nTraits * sizeof(int))) return fail("");
    if (!dev((void**)&s->dScratch, (std::size_t)s->blocks * 2 * a.N * sizeof(double))) return fail("");
    if (!dev((void**)&s->dCounter, sizeof(int))) return fail("");
    if (a.ownSamples) {
        s->own = 1;
        s->maskWords = (a.N + 63) / 64;
        if (cudaHostAlloc((void**)&s->hPLut, (std::size_t)a.maxPairs * 4 * sizeof(double), cudaHostAllocDefault) != cudaSuccess) return fail("cudaHostAlloc pair luts");
        if (!dev((void**)&s->dPLut, (std::size_t)a.maxPairs * 4 * sizeof(double))) return fail("");
        if (!dev((void**)&s->dMask, (std::size_t)a.nTraits * s->maskWords * sizeof(uint64_t))) return fail("");
        std::vector<uint64_t> all((std::size_t)s->maskWords, ~0ull);
        for (int t = 0; t < a.nTraits; ++t) {
            const uint64_t* m = a.traits[t].mask ? a.traits[t].mask : all.data();
            if (cudaMemcpy(s->dMask + (std::size_t)t * s->maskWords, m, (std::size_t)s->maskWords * sizeof(uint64_t),
                           cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy mask");
        }
    }
    s->devBytes = db;
    std::vector<int> pv(a.nTraits);
    for (int t = 0; t < a.nTraits; ++t) {
        const TraitArgs& T = a.traits[t];
        pv[t] = T.p;
        if (cudaMemcpy(s->dMu + (std::size_t)t * a.N, T.mu, (std::size_t)a.N * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy mu");
        if (cudaMemcpy(s->dXV + (std::size_t)t * s->traitStride, T.XV, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XV");
        if (cudaMemcpy(s->dXX + (std::size_t)t * s->traitStride, T.XXVX_inv, (std::size_t)a.N * T.p * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XXVX_inv");
    }
    if (cudaMemcpy(s->dP, pv.data(), pv.size() * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy p");
    if (a.precision == saige::gpu2::Prec::FP32) {
        // float copies for passes A and B; mu in the stored form of terms32:
        // min(mu, 1 - mu), negated where mu > 1/2 (1 - mu formed in double)
        std::size_t db32 = 0;
        auto dev32 = [&](float** p, std::size_t n) {
            if (cudaMalloc((void**)p, n * sizeof(float)) != cudaSuccess) { *p = nullptr; return false; }
            db32 += n * sizeof(float); return true;
        };
        if (!dev32(&s->dMu32, (std::size_t)a.nTraits * a.N) || !dev32(&s->dMu32lo, (std::size_t)a.nTraits * a.N) ||
            !dev32(&s->dXV32lo, (std::size_t)a.nTraits * s->traitStride) || !dev32(&s->dXX32lo, (std::size_t)a.nTraits * s->traitStride) || !dev32(&s->dXV32, (std::size_t)a.nTraits * s->traitStride) ||
            !dev32(&s->dXX32, (std::size_t)a.nTraits * s->traitStride)) return fail("cudaMalloc fp32 copies");
        s->devBytes += db32;
        std::vector<float> f;
        for (int t = 0; t < a.nTraits; ++t) {
            const TraitArgs& T = a.traits[t];
            f.resize((std::size_t)a.N);
            for (int i = 0; i < a.N; ++i) f[i] = T.mu[i] > 0.5 ? -(float)(1.0 - T.mu[i]) : (float)T.mu[i];
            if (cudaMemcpy(s->dMu32 + (std::size_t)t * a.N, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy mu32");
            for (int i = 0; i < a.N; ++i) { const double mm = T.mu[i] > 0.5 ? 1.0 - T.mu[i] : T.mu[i]; f[i] = (float)(mm - (double)(float)mm); }
            if (cudaMemcpy(s->dMu32lo + (std::size_t)t * a.N, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy mu32lo");
            f.resize((std::size_t)a.N * T.p);
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)T.XV[i];
            if (cudaMemcpy(s->dXV32 + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XV32");
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)(T.XV[i] - (double)(float)T.XV[i]);
            if (cudaMemcpy(s->dXV32lo + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XV32lo");
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)T.XXVX_inv[i];
            if (cudaMemcpy(s->dXX32 + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XX32");
            for (std::size_t i = 0; i < f.size(); ++i) f[i] = (float)(T.XXVX_inv[i] - (double)(float)T.XXVX_inv[i]);
            if (cudaMemcpy(s->dXX32lo + (std::size_t)t * s->traitStride, f.data(), f.size() * sizeof(float), cudaMemcpyHostToDevice) != cudaSuccess) return fail("memcpy XX32lo");
        }
    }
    if (cudaStreamCreate(&s->st) != cudaSuccess) return fail("cudaStreamCreate");
    for (int i = 0; i < 4; ++i) if (cudaEventCreate(&s->ev[i]) != cudaSuccess) return fail("cudaEventCreate");
    return s;
}

void destroy(Spa* s)
{
    if (!s) return;
    for (int i = 0; i < 4; ++i) if (s->ev[i]) cudaEventDestroy(s->ev[i]);
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
    if (s->dCounter) cudaFree(s->dCounter);
    if (s->dPLut) cudaFree(s->dPLut);
    if (s->dMask) cudaFree(s->dMask);
    if (s->hPLut) cudaFreeHost(s->hPLut);
    if (s->dGPk) cudaFree(s->dGPk);
    if (s->dGLut) cudaFree(s->dGLut);
    if (s->dGDense) cudaFree(s->dGDense);
    if (s->hIn) cudaFreeHost(s->hIn);
    if (s->hOut) cudaFreeHost(s->hOut);
    delete s;
}

PairIn*  in(Spa* s)  { return s ? s->hIn : nullptr; }
PairOut* out(Spa* s) { return s ? s->hOut : nullptr; }
double*  pairLut(Spa* s) { return (s && s->own) ? s->hPLut : nullptr; }
std::size_t deviceBytes(const Spa* s) { return s ? s->devBytes : 0; }
const char* lastError() { return g_lastErr.c_str(); }

void timings(const Spa* s, double* k, double* h2d, double* d2h, long long* pairs)
{
    if (!s) return;
    if (k)     *k = s->tKernel;
    if (h2d)   *h2d = s->tH2D;
    if (d2h)   *d2h = s->tD2H;
    if (pairs) *pairs = s->nPairsDone;
}

namespace {
bool grow(void** p, std::size_t* cap, std::size_t need, std::size_t* devBytes)
{
    if (need <= *cap) return true;
    if (*p) { cudaFree(*p); *devBytes -= *cap; *p = nullptr; *cap = 0; }
    CK(cudaMalloc(p, need));
    *cap = need; *devBytes += need;
    return true;
}
}  // namespace

bool uploadPacked(Spa* s, const unsigned char* packed, std::size_t bpv, const double* lut, int nSlots, Geno* view)
{
    if (!s || !packed || !lut || nSlots <= 0 || !view) { g_lastErr = "uploadPacked: bad arguments"; return false; }
    if (!grow((void**)&s->dGPk, &s->capPk, (std::size_t)nSlots * bpv, &s->devBytes)) return false;
    if (!grow((void**)&s->dGLut, &s->capLut, (std::size_t)nSlots * 4 * sizeof(double), &s->devBytes)) return false;
    CK(cudaMemcpy(s->dGPk, packed, (std::size_t)nSlots * bpv, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(s->dGLut, lut, (std::size_t)nSlots * 4 * sizeof(double), cudaMemcpyHostToDevice));
    *view = Geno();
    view->packed = s->dGPk; view->bpv = bpv; view->lut = s->dGLut;
    return true;
}

bool uploadDense(Spa* s, const double* g, int nSlots, Geno* view)
{
    if (!s || !g || nSlots <= 0 || !view) { g_lastErr = "uploadDense: bad arguments"; return false; }
    const std::size_t n = (std::size_t)nSlots * s->N * sizeof(double);
    if (!grow((void**)&s->dGDense, &s->capDense, n, &s->devBytes)) return false;
    CK(cudaMemcpy(s->dGDense, g, n, cudaMemcpyHostToDevice));
    *view = Geno();
    view->dense = s->dGDense; view->ld = (std::size_t)s->N;
    return true;
}

bool run(Spa* s, const Geno& geno, int nPairs)
{
    if (!s) { g_lastErr = "run: null Spa"; return false; }
    if (nPairs <= 0) return true;
    if (nPairs > s->maxPairs) { g_lastErr = "run: nPairs > maxPairs"; return false; }
    const bool hasPk = geno.packed && geno.lut && geno.bpv > 0;
    const bool hasDn = geno.dense && geno.ld >= (std::size_t)s->N;
    if (hasPk == hasDn) { g_lastErr = "run: Geno must set exactly one of packed / dense"; return false; }

    KParams P;
    P.packed = hasPk ? geno.packed : nullptr; P.bpv = geno.bpv; P.lut = hasPk ? geno.lut : nullptr;
    P.dense = hasDn ? geno.dense : nullptr; P.ld = geno.ld;
    P.MU32lo = s->dMu32lo; P.XV32lo = s->dXV32lo; P.XX32lo = s->dXX32lo;
    P.MU32 = s->dMu32; P.XV32 = s->dXV32; P.XX32 = s->dXX32;
    P.N = s->N; P.MU = s->dMu; P.XV = s->dXV; P.XX = s->dXX; P.pOfTrait = s->dP; P.traitStride = s->traitStride;
    P.in = s->dIn; P.nPairs = nPairs; P.scratch = s->dScratch;
    P.tol = s->tol; P.maxiter = s->maxiter; P.erfcMode = s->erfcMode; P.out = s->dOut;
    P.counter = s->dCounter;
    P.plut = s->dPLut; P.masks = s->dMask; P.maskWords = s->maskWords;

    const int grid = nPairs < s->blocks ? nPairs : s->blocks;
    CK(cudaEventRecord(s->ev[0], s->st));
    CK(cudaMemcpyAsync(s->dIn, s->hIn, (std::size_t)nPairs * sizeof(PairIn), cudaMemcpyHostToDevice, s->st));
    if (s->dyn) CK(cudaMemsetAsync(s->dCounter, 0, sizeof(int), s->st));
    if (s->own) CK(cudaMemcpyAsync(s->dPLut, s->hPLut, (std::size_t)nPairs * 4 * sizeof(double), cudaMemcpyHostToDevice, s->st));
    CK(cudaEventRecord(s->ev[1], s->st));
    // ---- precision dispatch: the kernel variant ----
    switch (s->prec) {
        case saige::gpu2::Prec::FP64:
            launchSpa<false>(s->fused, s->dyn, s->minb, s->own, grid, s->st, P);
            break;
        case saige::gpu2::Prec::FP32:
            // Newton tolerance: the caller's, floored at kTolFp32 (spa_gpu.hpp).
            if (P.tol < kTolFp32) P.tol = kTolFp32;
            launchSpa<true>(s->fused, s->dyn, s->minb, s->own, grid, s->st, P);
            break;
        default:
            g_lastErr = std::string("run: SPA precision ") + saige::gpu2::precName(s->prec) + " has no kernel";
            return false;
    }
    CK(cudaGetLastError());
    CK(cudaEventRecord(s->ev[2], s->st));
    CK(cudaMemcpyAsync(s->hOut, s->dOut, (std::size_t)nPairs * sizeof(PairOut), cudaMemcpyDeviceToHost, s->st));
    CK(cudaEventRecord(s->ev[3], s->st));
    CK(cudaStreamSynchronize(s->st));
    float ms = 0;
    if (cudaEventElapsedTime(&ms, s->ev[0], s->ev[1]) == cudaSuccess) s->tH2D += ms * 1e-3;
    if (cudaEventElapsedTime(&ms, s->ev[1], s->ev[2]) == cudaSuccess) s->tKernel += ms * 1e-3;
    if (cudaEventElapsedTime(&ms, s->ev[2], s->ev[3]) == cudaSuccess) s->tD2H += ms * 1e-3;
    s->nPairsDone += nPairs;
    return true;
}

bool debugErfc(int device, const double* z, int n, double* a, double* b, double* c)
{
    if (n <= 0) return true;
    CK(cudaSetDevice(device));
    double *dz = nullptr, *da = nullptr, *db = nullptr, *dc = nullptr;
    CK(cudaMalloc(&dz, n * sizeof(double)));
    CK(cudaMalloc(&da, n * sizeof(double)));
    CK(cudaMalloc(&db, n * sizeof(double)));
    CK(cudaMalloc(&dc, n * sizeof(double)));
    CK(cudaMemcpy(dz, z, n * sizeof(double), cudaMemcpyHostToDevice));
    erfcDebugKernel<<<(n + 255) / 256, 256>>>(dz, n, da, db, dc);
    CK(cudaGetLastError());
    CK(cudaMemcpy(a, da, n * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(b, db, n * sizeof(double), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(c, dc, n * sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(dz); cudaFree(da); cudaFree(db); cudaFree(dc);
    return true;
}

}  // namespace spa_gpu
}  // namespace saige
