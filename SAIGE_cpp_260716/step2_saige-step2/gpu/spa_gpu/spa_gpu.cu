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

namespace saige {
namespace spa_gpu {

namespace {

constexpr int NT    = 256;
constexpr int NWARP = NT / 32;
constexpr int NACC  = PMAX + 1;   // widest reduction: p partial dot products + the carrier count

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
    const double* XV;            // per trait: traitStride doubles, sample i at i*p
    const double* XX;            // per trait: traitStride doubles, column j at j*N
    const int*    pOfTrait;
    std::size_t   traitStride;
    // pairs
    const PairIn* in; int nPairs;
    double* scratch;             // blocks x 2N doubles
    double tol; int maxiter; int erfcMode;
    PairOut* out;
};

struct GenoCol {
    const unsigned char* col;    // packed row, or nullptr
    double4 L;                   // its dosage table
    const double* dcol;          // dense column, or nullptr
    __device__ __forceinline__ double dose(int i) const
    {
        if (dcol) return dcol[i];
        const unsigned c = (col[i >> 2] >> (2 * (i & 3))) & 3u;
        return (c & 2u) ? ((c & 1u) ? L.w : L.z) : ((c & 1u) ? L.y : L.x);
    }
};

// Per-pair state the root / tail passes need.
struct Ctx {
    const double2* buf;   // (g~, mu), nEff entries
    int nEff;
    bool fast;
    double NAmu, NAsigma;
};

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

// getroot_K1_Binom (full) / getroot_K1_fast_Binom (fast), init 0, with the
// two variants' own tests. Every thread runs the same scalar code on the
// same block-wide sums, so there is no divergence and no broadcast.
__device__ void getroot(const Ctx& C, double q, double gpos, double gneg, double tol, int maxiter,
                        double (*sh)[NWARP], double* root, int* niter, int* conv, int* reason)
{
    *reason = RS_OK;
    if (q >= gpos || q <= gneg) { *root = CUDART_INF; *niter = 0; *conv = 1; return; }
    double t = 0.0;
    double K1, K2;
    passK1K2(C, t, q, sh, &K1, &K2);         // K1_eval at init; K2 at the same t for the loop top
    double prevJump = CUDART_INF;
    int rep = 1;
    int c = 1;
    while (rep <= maxiter) {
        // K2_eval = K2(t): already in K2 (evaluated with the current K1 at this t)
        if (!C.fast && (!isfinite(K2) || fabs(K2) < 1e-15)) { c = 0; *reason = RS_K2_GUARD; break; }
        double tnew = t - K1 / K2;
        if (C.fast ? isnan(tnew) : !isfinite(tnew)) { c = 0; *reason = RS_TNEW; break; }
        if (fabs(tnew - t) < tol) { c = 1; break; }
        if (rep == maxiter) { c = 0; *reason = RS_MAXITER; break; }
        double newK1, newK2;
        passK1K2(C, tnew, q, sh, &newK1, &newK2);
        const bool flipped = C.fast ? ((K1 * newK1) < 0) : (sgn(K1) != sgn(newK1));
        if (flipped) {
            if (fabs(tnew - t) > (prevJump - tol)) {
                tnew = t + (double)sgn(newK1 - K1) * prevJump / 2;
                passK1K2(C, tnew, q, sh, &newK1, &newK2);
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

// Get_Saddle_Prob_Binom / Get_Saddle_Prob_fast_Binom.
__device__ double saddle(const Ctx& C, double zeta, double q, int logp, int erfcMode,
                         double (*sh)[NWARP], int* isSaddle)
{
    double k1, k2;
    passK0K2(C, zeta, sh, &k1, &k2);
    const double temp1 = zeta * q - k1;
    *isSaddle = 0;
    bool flagrun = false;
    double w = 0.0, v = 0.0;
    if (isfinite(k1) && isfinite(k2) && temp1 >= 0 && k2 >= 0) {
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

__global__ void __launch_bounds__(NT)
spaKernel(const KParams P)
{
    __shared__ double sh[NACC][NWARP];
    const int tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const int N = P.N;
    double2* buf = reinterpret_cast<double2*>(P.scratch) + (std::size_t)blockIdx.x * N;
    // this warp's contiguous share of the samples (passes A and B)
    const int segLo = (int)(((long long)N * warp) / NWARP);
    const int segHi = (int)(((long long)N * (warp + 1)) / NWARP);

    for (int k = blockIdx.x; k < P.nPairs; k += gridDim.x) {
        const PairIn pin = P.in[k];
        GenoCol G;
        if (P.dense) {
            G.dcol = P.dense + (std::size_t)pin.slot * P.ld; G.col = nullptr;
            G.L = make_double4(0, 0, 0, 0);
        } else {
            G.dcol = nullptr;
            G.col = P.packed + (std::size_t)pin.slot * P.bpv;
            const double2* d = reinterpret_cast<const double2*>(P.lut) + 2 * (std::size_t)pin.slot;
            const double2 a = d[0], b = d[1];
            G.L = make_double4(a.x, a.y, b.x, b.y);
        }
        const int t = pin.trait;
        const int p = P.pOfTrait[t];
        const double* mu = P.MU + (std::size_t)t * N;
        const double* xv = P.XV + (std::size_t)t * P.traitStride;
        const double* xx = P.XX + (std::size_t)t * P.traitStride;

        // ---- pass A: b = XV g over carriers (getadjGFast's loop), carrier count
        double acc[NACC];
        #pragma unroll
        for (int j = 0; j < NACC; ++j) acc[j] = 0.0;
        for (int i = segLo + lane; i < segHi; i += 32) {
            const double g = G.dose(i);
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
        const int nnz = (int)acc[PMAX];
        // the CPU's rule: p_iIndexComVecSize = double(|{g==0}|) / m_n >= 0.5
        const bool fast = (pin.fast > 0) || (pin.fast < 0 && ((double)(N - nnz) / (double)N) >= 0.5);
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
                g = G.dose(i);
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
        const double gpos = s[0], gneg = s[1], m1 = s[2];

        Ctx C;
        C.buf = buf; C.nEff = fast ? nnz : N; C.fast = fast;
        C.NAmu    = fast ? (m1 - s[3]) : 0.0;           // NAmu = m1 - dot(gNB, muNB)
        C.NAsigma = fast ? (pin.var2 - s[4]) : 0.0;     // NAsigma = var2 - sum(muNB % (1-muNB) % pow(gNB,2))

        // q, qinv as getMarkerPval forms them for a binary trait
        const double q = pin.Tstat / sqrt(pin.var1 / pin.var2) + m1;
        double qinv;
        if ((q - m1) > 0)       qinv = -1 * fabs(q - m1) + m1;
        else if ((q - m1) == 0) qinv = m1;
        else                    qinv = fabs(q - m1) + m1;

        unsigned status = (fast ? ST_FAST : 0u) | (pin.logp ? ST_LOGP : 0u);
        double r1, r2; int n1, n2, c1, c2, rs1, rs2;
        getroot(C, q,    gpos, gneg, P.tol, P.maxiter, sh, &r1, &n1, &c1, &rs1);
        getroot(C, qinv, gpos, gneg, P.tol, P.maxiter, sh, &r2, &n2, &c2, &rs2);
        if (n1 == 0 && isinf(r1)) status |= ST_ROOT1_INF;
        if (n2 == 0 && isinf(r2)) status |= ST_ROOT2_INF;
        if (!c1) status |= ST_ROOT1_FAIL;
        if (!c2) status |= ST_ROOT2_FAIL;

        double pv, p1 = 0.0, p2 = 0.0; int conv, s1 = -1, s2 = -1;
        if (c1 && c2) {
            // spa.cpp SPA / SPA_fast: a tail that is not a saddle withdraws convergence
            p1 = saddle(C, r1, q,    pin.logp, P.erfcMode, sh, &s1);
            p2 = saddle(C, r2, qinv, pin.logp, P.erfcMode, sh, &s2);
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
    int N = 0, nTraits = 0, maxPairs = 0, blocks = 256, maxiter = 1000, erfcMode = 1;
    double tol = 0.0;
    std::size_t traitStride = 0;
    PairIn*  hIn  = nullptr;
    PairOut* hOut = nullptr;
    PairIn*  dIn  = nullptr;
    PairOut* dOut = nullptr;
    double* dMu = nullptr; double* dXV = nullptr; double* dXX = nullptr;
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

Spa* create(const CreateArgs& a)
{
    g_lastErr.clear();
    if (a.N <= 0 || a.nTraits <= 0 || a.traits == nullptr || a.maxPairs <= 0) {
        g_lastErr = "create: bad arguments"; return nullptr;
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
    s->N = a.N; s->nTraits = a.nTraits; s->maxPairs = a.maxPairs;
    s->blocks = a.blocks > 0 ? a.blocks : 256;
    s->maxiter = a.maxiter; s->tol = a.tol; s->erfcMode = a.erfcMode;
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
    if (s->dScratch) cudaFree(s->dScratch);
    if (s->dGPk) cudaFree(s->dGPk);
    if (s->dGLut) cudaFree(s->dGLut);
    if (s->dGDense) cudaFree(s->dGDense);
    if (s->hIn) cudaFreeHost(s->hIn);
    if (s->hOut) cudaFreeHost(s->hOut);
    delete s;
}

PairIn*  in(Spa* s)  { return s ? s->hIn : nullptr; }
PairOut* out(Spa* s) { return s ? s->hOut : nullptr; }
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
    P.N = s->N; P.MU = s->dMu; P.XV = s->dXV; P.XX = s->dXX; P.pOfTrait = s->dP; P.traitStride = s->traitStride;
    P.in = s->dIn; P.nPairs = nPairs; P.scratch = s->dScratch;
    P.tol = s->tol; P.maxiter = s->maxiter; P.erfcMode = s->erfcMode; P.out = s->dOut;

    const int grid = nPairs < s->blocks ? nPairs : s->blocks;
    CK(cudaEventRecord(s->ev[0], s->st));
    CK(cudaMemcpyAsync(s->dIn, s->hIn, (std::size_t)nPairs * sizeof(PairIn), cudaMemcpyHostToDevice, s->st));
    CK(cudaEventRecord(s->ev[1], s->st));
    spaKernel<<<grid, NT, 0, s->st>>>(P);
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
