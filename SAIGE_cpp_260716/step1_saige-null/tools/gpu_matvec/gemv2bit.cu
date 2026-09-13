// gemv2bit.cu — tier-4 kernels: 2-bit-resident GRM matvec with the
// standardization lifted out by a rank-one decomposition.
//
// Ported from the SAIGE R package, gpu-opt branch, src/gpuSymMatMult.cu
// (class gpuSymMatMult: gm_unpack16 / gm_field2f / gm_sumx / gm_pass1 /
// gm_pass1_finish / gm_calcC / gm_pass2 / gm_pass2_finish, set_matrix_packed,
// sym_gemv_2bit_range). The math is copied unchanged; see gemv2bit.hpp for the
// identity and the encoding contract.
//
// Deliberate differences from the R original, all mechanical:
//   1. `int rank` (an MPI rank, used only for log prefixes and
//      cudaSetDevice(rank % count)) is gone. This binary is single-process, so
//      rank was always 0 and rank%count always selects device 0. Device choice
//      is now SAIGE_GPU_DEVICE (default 0).
//   2. The final scale by 1/M_pass moved INTO gm_pass2_finish. The R package
//      divides on the host in SAIGE_fitGLMM_fast.cpp; the C++ facade contract
//      (gpu_matvec.hpp) says matvec returns K·u = (A Aᵀ / M_pass)·u, and doing
//      it in the kernel avoids an extra length-N host pass.
//   3. Error paths return false/nullptr instead of std::cerr + EXIT_FAILURE —
//      the facade's contract is a silent CPU fallback, not an abort.
//   4. No cuBLAS handle. The R class created one in set_matrix_packed; the
//      2-bit path never used it.
//   5. Padding: the R version rounds W32 up to MM_MT/16 and M up to MM_MT so
//      the fp16/wmma batch kernels need no bounds checks. That path is NOT
//      ported, so W32 = ⌈N/16⌉ and the marker count is used as-is — which is
//      what shrinks the packed allocation from Mpad·W32·4 to M·W32·4.
//   6. Comments translated from Chinese; measured numbers kept verbatim.
//
// NOT ported (stay in the R package): sym_gemm_2bit (fp16 + wmma tensor-core
// batch path — mm_pack_X rounds the RHS to fp16 and moves tau by ~2e-4, and
// the C++ trunk has no multi-RHS entry point yet), sym_gemm_2bit_fp32 (fp32
// multi-column, only useful once parallelCrossProdMat is wired up),
// set_matrix / sym_sgemv / sym_sgemm / sym_hgemm (dense fp32 and fp16 cuBLAS —
// tier 1/2 in gpu_matvec.cu already cover that).
#include "gemv2bit.hpp"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace saige::gpu::g2b {

namespace {

#define G2B_CHECK_RET(expr, ret_on_fail)                                       \
  do {                                                                         \
    cudaError_t _s = (expr);                                                   \
    if (_s != cudaSuccess) {                                                   \
      std::fprintf(stderr, "[gemv2bit] CUDA error at %s:%d: %s\n",             \
                   __FILE__, __LINE__, cudaGetErrorString(_s));                \
      return ret_on_fail;                                                      \
    }                                                                          \
  } while (0)

#define G2B_CHECK_FALSE(expr) G2B_CHECK_RET(expr, false)

// ---------------------------------------------------------------- tuning
#define G_NTHREAD   256
#define G_NWARP     (G_NTHREAD/32)
// Pass-1 shared x tile. 128 uint32 = 2048 samples → 8704 B of shared, which on
// V100 (96 KB/SM) fits exactly 8 blocks = 2048 threads = 100 % occupancy.
// 256 would need 17408 B → only 5 blocks = 62 %, measured 10 % slower.
// x reuse is ACROSS markers (G_MTILE1 markers share one tile), independent of
// tile width, so halving it costs no reuse — only grid.y doubles and Ypart
// grows from 2 MB to 4 MB.
#define G_TILE_W    128
#define G_SPAD      17     // 17 floats per 16 samples in shared → bank-offset
#define G_MTILE1    32     // markers per pass-1 block (sharing one x tile)
#define G_MTILE2    256    // pass-2 marker slice

// g = 2 − popc(field), whole-word parallel. Each field's popc ≤ 2 so the
// subtraction never borrows across fields.
__device__ __forceinline__ unsigned int gm_unpack16(unsigned int p)
{
    return 0xAAAAAAAAu - (p & 0x55555555u) - ((p >> 1) & 0x55555555u);
}

// Extract the 2-bit field at bit s and turn it into a float WITHOUT an integer
// →float conversion instruction. 8388608.0f = 2^23 has bit pattern 0x4B000000,
// whose mantissa LSB has weight exactly 1, so reading 0x4B000000|f as a float
// gives 8388608+f (exact for f ≤ 2); subtracting back yields f exactly.
// This dodges I2F — Volta's conversion unit issues only 16/SM/cycle, a quarter
// of FP32. (a & 3) | magic is one LOP3 (immLut 0xEA = (a&b)|c) with magic in a
// register; as two literals the compiler would have to split it in two.
// So each genotype costs SHF + LOP3 (INT pipe) + FADD + FFMA (FP pipe).
// Bit-identical to the I2F version: both f and 8388608+f are exact in fp32.
// Measured (optimization/bench/10_gemv2bit.cu): both passes 2.38 → 2.23 ms.
__device__ __forceinline__ float gm_field2f(unsigned int p, int s, unsigned int magic)
{
    unsigned int f;
    asm("lop3.b32 %0, %1, 3, %2, 0xEA;" : "=r"(f) : "r"(p >> s), "r"(magic));
    return __int_as_float(f) - 8388608.0f;
}

// Pass 1: Ypart[tile][j] = Σ_{i in tile} g[i,j] · x[i]
// One warp per marker; the G_MTILE1 markers in a block share one shared-memory
// x tile — otherwise every marker re-reads all of x and L2 traffic (M·N·4)
// exceeds the matrix itself.
// The second __launch_bounds__ argument caps registers at 32 to max out
// occupancy; measured both passes 52 % → 58 % of roofline. This kernel is
// latency-bound, so more warps beat more registers (64/128/256 registers made
// no difference; squeezing to 32 was 10 % faster).
// j0/jn is the marker range (LOCO uses one chromosome's slice; full GRM passes
// j0=0, jn=M). Ypart is indexed by the LOCAL index j−j0 so finish/calcC/pass2
// can all walk the same contiguous span.
__global__ void __launch_bounds__(G_NTHREAD, 8)
gm_pass1(const unsigned int* __restrict__ Ap, int W32,
         const float* __restrict__ x, float* __restrict__ Ypart,
         int j0, int jn, int ypitch)
{
    __shared__ float sx[G_TILE_W*G_SPAD];

    int w0 = blockIdx.y * G_TILE_W;
    int nw = W32 - w0; if (nw > G_TILE_W) nw = G_TILE_W;
    if (nw <= 0) return;

    for (int t = threadIdx.x; t < nw*16; t += G_NTHREAD)
        sx[(t >> 4)*G_SPAD + (t & 15)] = x[w0*16 + t];
    __syncthreads();

    unsigned int magic = 0x4B000000u;
    int warp = threadIdx.x >> 5, lane = threadIdx.x & 31;
    int jb = blockIdx.x * G_MTILE1;
    for (int jj = warp; jj < G_MTILE1; jj += G_NWARP) {
        int jl = jb + jj;                       // local index within the range
        if (jl >= jn) break;
        const unsigned int *row = Ap + (size_t)(j0 + jl)*W32 + w0;
        // 4 accumulators: with one, the 16 FMAs form a single dependency chain
        // and the full latency is exposed.
        float a0 = 0.f, a1 = 0.f, a2 = 0.f, a3 = 0.f;
        for (int w = lane; w < nw; w += 32) {
            unsigned int g = gm_unpack16(row[w]);
            const float *xs = sx + w*G_SPAD;
            #pragma unroll
            for (int b = 0; b < 16; b += 4) {
                a0 += gm_field2f(g, 2*b    , magic) * xs[b  ];
                a1 += gm_field2f(g, 2*b + 2, magic) * xs[b+1];
                a2 += gm_field2f(g, 2*b + 4, magic) * xs[b+2];
                a3 += gm_field2f(g, 2*b + 6, magic) * xs[b+3];
            }
        }
        float acc = (a0 + a1) + (a2 + a3);
        #pragma unroll
        for (int off = 16; off; off >>= 1) acc += __shfl_down_sync(0xffffffffu, acc, off);
        if (lane == 0) Ypart[(size_t)blockIdx.y*ypitch + jl] = acc;
    }
}

// Pass 2: Zpart[tile][i] = Σ_{j in tile} g[i,j] · w[j]
// One uint32 (16 samples) per thread, accumulating along the marker axis:
// consecutive threads read consecutive uint32 → one coalesced load per marker,
// and w[j] is a broadcast.
__global__ void __launch_bounds__(G_NTHREAD, 8)
gm_pass2(const unsigned int* __restrict__ Ap, int W32,
         const float* __restrict__ w, float* __restrict__ Zpart,
         int j0, int jn, int Npad)
{
    int wi = blockIdx.x*blockDim.x + threadIdx.x;
    if (wi >= W32) return;
    int jl0 = blockIdx.y * G_MTILE2;
    int jl1 = jl0 + G_MTILE2; if (jl1 > jn) jl1 = jn;

    unsigned int magic = 0x4B000000u;
    float acc[16];
    #pragma unroll
    for (int b = 0; b < 16; ++b) acc[b] = 0.f;
    for (int jl = jl0; jl < jl1; ++jl) {
        unsigned int g = gm_unpack16(Ap[(size_t)(j0 + jl)*W32 + wi]);
        float wj = w[jl];
        #pragma unroll
        for (int b = 0; b < 16; ++b) acc[b] += gm_field2f(g, 2*b, magic) * wj;
    }
    float *out = Zpart + (size_t)blockIdx.y*Npad + wi*16;
    #pragma unroll
    for (int b = 0; b < 16; ++b) out[b] = acc[b];
}

// S = Σ_i x[i], single-block reduction (N elements, a few microseconds).
// The result stays in device memory — it never goes back to the host.
__global__ void gm_sumx(const float* __restrict__ x, int N, float* __restrict__ out)
{
    __shared__ float s[G_NTHREAD];
    float acc = 0.f;
    for (int i = threadIdx.x; i < N; i += G_NTHREAD) acc += x[i];
    s[threadIdx.x] = acc; __syncthreads();
    for (int k = G_NTHREAD/2; k; k >>= 1) {
        if (threadIdx.x < k) s[threadIdx.x] += s[threadIdx.x+k];
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = s[0];
}

// Pass-1 reduction + rank-one correction.
// M here is the range length jn; freq/invStd are offset to j0 by the caller.
// ypitch is Ypart's row stride and is always the FULL marker count (not jn):
// the full-GRM call and a LOCO range call share one buffer, so the stride has
// to agree or any future cross-call reuse (warm start, say) reads skewed data.
__global__ void gm_pass1_finish(const float* __restrict__ Ypart, int ntile, int M, int ypitch,
                                const float* __restrict__ freq, const float* __restrict__ invStd,
                                const float* __restrict__ Sp,
                                float* __restrict__ Y, float* __restrict__ w)
{
    int j = blockIdx.x*blockDim.x + threadIdx.x;
    if (j >= M) return;
    float s = 0.f;
    for (int t = 0; t < ntile; ++t) s += Ypart[(size_t)t*ypitch + j];
    float y = invStd[j] * (s - 2.f*freq[j]*(*Sp));
    Y[j] = y;
    w[j] = invStd[j] * y;
}

// C = Σ_j 2 freq[j] w[j]
__global__ void gm_calcC(const float* __restrict__ freq, const float* __restrict__ w,
                         int M, float* __restrict__ out)
{
    __shared__ float s[G_NTHREAD];
    float acc = 0.f;
    for (int j = threadIdx.x; j < M; j += G_NTHREAD) acc += 2.f*freq[j]*w[j];
    s[threadIdx.x] = acc; __syncthreads();
    for (int k = G_NTHREAD/2; k; k >>= 1) {
        if (threadIdx.x < k) s[threadIdx.x] += s[threadIdx.x+k];
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = s[0];
}

// Pass-2 reduction + rank-one correction, then the facade's 1/M_pass scale.
// (The R original stops at `Z[i] = s − C`; the division happens on the host
// there. Folding inv_M in here keeps it off the host critical path.)
__global__ void gm_pass2_finish(const float* __restrict__ Zpart, int ntile, int N, int Npad,
                                const float* __restrict__ Cp, float inv_M,
                                float* __restrict__ Z)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= N) return;
    float s = 0.f;
    for (int t = 0; t < ntile; ++t) s += Zpart[(size_t)t*Npad + i];
    Z[i] = (s - (*Cp)) * inv_M;
}

// ======================= fp32 small-batch path (ncol 2..8) ================
//
// Motive: the trunk's block-PCG (getCrossprodMatAndKinMat, used by GetTrace's
// nrun probes and by the batched VR solves) asks for ψ·B with k columns at
// once. Served by k separate single-column matvecs, the 2-bit matrix is
// re-streamed k times; these kernels read it once per chunk of ≤8 columns.
//
// The fp16/wmma batch path from the R package is deliberately NOT used here:
// mm_pack_X rounds the RHS to fp16 and moved tau by ~2e-4. This family is
// fp32 end to end.
//
// Microbenchmark (optimization/bench/12_gemv2bit_mc.cu, V100), speedup vs n
// separate single-column calls: ncol=2 1.5x / 4 1.9x / 8 2.1x (same on mid and
// big2). Against an fp64 reference: mc 2.863e-6 vs per-column 2.848e-6 — no
// extra error, only a different summation order. Two runs bit-identical.
// The "theoretical 4x" is unreachable: at ncol=8 the FFMA lower bound for both
// passes is 4.3 ms and the decode (SHF+LOP3+FADD) plus shared/global loads add
// ~65 % of issue slots on top, so the instruction-bound ceiling is about 2.4x.
//
// pass1_mc keeps the warp-per-marker shape of the single-column version, plus:
//   - shared holds NCOL per-column planes, SPAD=17 skew inside each plane (all
//     lanes of one instruction read the same column, so the bank analysis is
//     identical to the single-column case);
//   - staged loading: STAGE=256/NCOL uint32 per stage (shared stays 17408 B)
//     while a block covers GMC_CHUNK1=512 uint32, carrying accumulators in
//     registers across stages. Without staging the tile shrinks with NCOL and
//     Ypart grows to 1.4 GB at big2/ncol=8; staged it is a flat 90 MB.
//   - MPW=4 markers per warp: one shared read feeds 4 FMAs (FFMA:LDS = 4:1;
//     undiluted, the LDS issue slot becomes the new bottleneck at ncol=8).
//     Accumulators are 4·NCOL ≤ 32. Row pointers are kept as base+index —
//     storing 4 pointers (8 registers) measurably spilled the kernel.
// pass2_mc: the single-column version gives each thread one uint32 and 16
//   accumulators; 16·NCOL would be 128 registers at ncol=8. Instead T=NCOL/2
//   threads share a uint32 (2 words / 4 half-words / 8 quarter-words) so the
//   accumulator count is a constant (16/T)·NCOL = 32. T=NCOL/2 rather than
//   T=NCOL because a warp's coalesced load then still covers 256/NCOL = 32 B
//   at ncol=8, one whole sector; T=NCOL leaves 16 B and doubles DRAM read
//   amplification. w is stored marker-major (w[jl*NCOL+c]) so a row of NCOL
//   weights is one or two float4/float2 loads — pure issue-slot savings.
//   The marker slice mtile2 is chosen at RUNTIME to fill the GPU (see
//   gmc_pick_g2y), not as a compile-time constant: a fixed 4096 leaves only
//   130 blocks on mid, 33 k threads, which does not saturate a V100.
#define GMC_CHUNK1 512
#define GMC_MPW    4
#define GMC_NCMAX  8
#define GMC_G2YCAP 64      // Zpart tile cap (64·8·Npad·4 B; 410 MB on big2)

template<int NCOL>
__global__ void __launch_bounds__(G_NTHREAD, 4)
gm_pass1_mc(const unsigned int* __restrict__ Ap, int W32,
            const float* __restrict__ X, int Npad, float* __restrict__ Ypart,
            int j0, int jn, int ypitch)
{
    constexpr int STAGE = 256/NCOL;
    constexpr int PLANE = STAGE*G_SPAD;
    constexpr int MPW   = GMC_MPW;
    constexpr int MT    = G_NWARP*MPW;
    __shared__ float sx[NCOL*PLANE];

    int w0 = blockIdx.y * GMC_CHUNK1;
    int wend = W32 - w0; if (wend > GMC_CHUNK1) wend = GMC_CHUNK1;
    if (wend <= 0) return;

    unsigned int magic = 0x4B000000u;
    int warp = threadIdx.x >> 5, lane = threadIdx.x & 31;
    int jb = blockIdx.x * MT + warp*MPW;

    // Out-of-range markers are clamped to jn-1 (legal memory, result never
    // written). Clamping beats branching: the g[r] load is innermost and a
    // per-r predicate would break up the coalesced access.
    const unsigned int *row0 = Ap + (size_t)j0*W32 + w0;
    int jrow[MPW];
    #pragma unroll
    for (int r = 0; r < MPW; ++r) {
        int jl = jb + r; if (jl >= jn) jl = jn - 1;
        jrow[r] = jl;
    }

    float acc[MPW][NCOL];
    #pragma unroll
    for (int r = 0; r < MPW; ++r)
        #pragma unroll
        for (int c = 0; c < NCOL; ++c) acc[r][c] = 0.f;

    for (int s0 = 0; s0 < wend; s0 += STAGE) {
        int nw = wend - s0; if (nw > STAGE) nw = STAGE;
        #pragma unroll
        for (int c = 0; c < NCOL; ++c)
            for (int t = threadIdx.x; t < nw*16; t += G_NTHREAD)
                sx[c*PLANE + (t >> 4)*G_SPAD + (t & 15)] = X[(size_t)c*Npad + (w0 + s0)*16 + t];
        __syncthreads();

        for (int w = lane; w < nw; w += 32) {
            unsigned int g[MPW];
            #pragma unroll
            for (int r = 0; r < MPW; ++r) g[r] = gm_unpack16(row0[(size_t)jrow[r]*W32 + s0 + w]);
            const float *xs = sx + w*G_SPAD;
            #pragma unroll
            for (int b = 0; b < 16; ++b) {
                float f[MPW];
                #pragma unroll
                for (int r = 0; r < MPW; ++r) f[r] = gm_field2f(g[r], 2*b, magic);
                #pragma unroll
                for (int c = 0; c < NCOL; ++c) {
                    float xv = xs[c*PLANE + b];      // one LDS feeds MPW FMAs
                    #pragma unroll
                    for (int r = 0; r < MPW; ++r) acc[r][c] += f[r]*xv;
                }
            }
        }
        __syncthreads();
    }

    #pragma unroll
    for (int r = 0; r < MPW; ++r) {
        int jl = jb + r;
        if (jl >= jn) break;
        #pragma unroll
        for (int c = 0; c < NCOL; ++c) {
            float v = acc[r][c];
            #pragma unroll
            for (int off = 16; off; off >>= 1) v += __shfl_down_sync(0xffffffffu, v, off);
            if (lane == 0)
                Ypart[((size_t)blockIdx.y*NCOL + c)*ypitch + jl] = v;
        }
    }
}

template<int NCOL>
__global__ void __launch_bounds__(G_NTHREAD, 4)
gm_pass2_mc(const unsigned int* __restrict__ Ap, int W32,
            const float* __restrict__ w, float* __restrict__ Zpart,
            int j0, int jn, int Npad, int mtile2)
{
    constexpr int T   = NCOL/2 > 0 ? NCOL/2 : 1;
    constexpr int SPT = 16/T;
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int wi  = tid / T, sub = tid % T;
    if (wi >= W32) return;
    int jl0 = blockIdx.y * mtile2;
    int jl1 = jl0 + mtile2; if (jl1 > jn) jl1 = jn;

    unsigned int magic = 0x4B000000u;
    float acc[SPT][NCOL];
    #pragma unroll
    for (int s = 0; s < SPT; ++s)
        #pragma unroll
        for (int c = 0; c < NCOL; ++c) acc[s][c] = 0.f;

    for (int jl = jl0; jl < jl1; ++jl) {
        unsigned int g = gm_unpack16(Ap[(size_t)(j0 + jl)*W32 + wi]);
        float wv[NCOL];
        if (NCOL == 8) {
            float4 lo = ((const float4*)(w + (size_t)jl*NCOL))[0];
            float4 hi = ((const float4*)(w + (size_t)jl*NCOL))[1];
            wv[0]=lo.x; wv[1]=lo.y; wv[2]=lo.z; wv[3]=lo.w;
            wv[4]=hi.x; wv[5]=hi.y; wv[6]=hi.z; wv[7]=hi.w;
        } else if (NCOL == 4) {
            float4 lo = *(const float4*)(w + (size_t)jl*NCOL);
            wv[0]=lo.x; wv[1]=lo.y; wv[2]=lo.z; wv[3]=lo.w;
        } else {
            float2 lo = *(const float2*)(w + (size_t)jl*NCOL);
            wv[0]=lo.x; wv[1]=lo.y;
        }
        #pragma unroll
        for (int s = 0; s < SPT; ++s) {
            float f = gm_field2f(g, 2*(sub*SPT + s), magic);
            #pragma unroll
            for (int c = 0; c < NCOL; ++c) acc[s][c] += f*wv[c];
        }
    }
    #pragma unroll
    for (int c = 0; c < NCOL; ++c) {
        float *out = Zpart + ((size_t)blockIdx.y*NCOL + c)*Npad + wi*16 + sub*SPT;
        #pragma unroll
        for (int s = 0; s < SPT; ++s) out[s] = acc[s][c];
    }
}

// S[c] = Σ_i X[i,c] (X's tail is zero-filled, so summing to Npad is the same)
__global__ void gm_colsum_mc(const float* __restrict__ X, int Npad, float* __restrict__ S)
{
    __shared__ float s[G_NTHREAD];
    int c = blockIdx.x;
    float acc = 0.f;
    for (int i = threadIdx.x; i < Npad; i += G_NTHREAD) acc += X[(size_t)c*Npad + i];
    s[threadIdx.x] = acc; __syncthreads();
    for (int k = G_NTHREAD/2; k; k >>= 1) {
        if (threadIdx.x < k) s[threadIdx.x] += s[threadIdx.x+k];
        __syncthreads();
    }
    if (threadIdx.x == 0) S[c] = s[0];
}

// Pass-1 reduction + rank-one correction; emits w marker-major [jl*NCOL+c]
// (the layout pass 2's vectorized load wants).
template<int NCOL>
__global__ void gm_pass1_finish_mc(const float* __restrict__ Ypart, int ntile, int M, int ypitch,
                                   const float* __restrict__ freq, const float* __restrict__ invStd,
                                   const float* __restrict__ S, float* __restrict__ w)
{
    int t = blockIdx.x*blockDim.x + threadIdx.x;
    if (t >= M*NCOL) return;
    int j = t / NCOL, c = t % NCOL;
    float s = 0.f;
    for (int ti = 0; ti < ntile; ++ti) s += Ypart[((size_t)ti*NCOL + c)*ypitch + j];
    float y = invStd[j] * (s - 2.f*freq[j]*S[c]);
    w[t] = invStd[j] * y;
}

// C[c] = Σ_j 2 freq[j] w[j,c]
template<int NCOL>
__global__ void gm_calcC_mc(const float* __restrict__ freq, const float* __restrict__ w,
                            int M, float* __restrict__ C)
{
    __shared__ float s[G_NTHREAD];
    int c = blockIdx.x;
    float acc = 0.f;
    for (int j = threadIdx.x; j < M; j += G_NTHREAD) acc += 2.f*freq[j]*w[(size_t)j*NCOL + c];
    s[threadIdx.x] = acc; __syncthreads();
    for (int k = G_NTHREAD/2; k; k >>= 1) {
        if (threadIdx.x < k) s[threadIdx.x] += s[threadIdx.x+k];
        __syncthreads();
    }
    if (threadIdx.x == 0) C[c] = s[0];
}

// Pass-2 reduction + rank-one correction + the 1/M_pass scale. Z column-major,
// ld = N. (Same inv_M deviation from the R original as gm_pass2_finish.)
template<int NCOL>
__global__ void gm_pass2_finish_mc(const float* __restrict__ Zpart, int ntile, int N, int Npad,
                                   const float* __restrict__ C, float inv_M,
                                   float* __restrict__ Z)
{
    size_t t = blockIdx.x*(size_t)blockDim.x + threadIdx.x;
    if (t >= (size_t)N*NCOL) return;
    int c = t / N, i = t % N;
    float s = 0.f;
    for (int ti = 0; ti < ntile; ++ti) s += Zpart[((size_t)ti*NCOL + c)*Npad + i];
    Z[t] = (s - C[c]) * inv_M;
}

// Pass-2 marker slice count: g2x is fixed by W32·T, so pad g2y until there are
// ~1536 blocks to fill the GPU, capped by GMC_G2YCAP (the Zpart buffer size)
// and by a floor of 256 markers per slice.
int gmc_pick_g2y(int W32, int T, int jn)
{
    int g2x = (int)(((long long)W32*T + G_NTHREAD - 1)/G_NTHREAD);
    int g2y = (1536 + g2x - 1)/g2x;
    int maxy = (jn + 255)/256;
    if (g2y > maxy) g2y = maxy;
    if (g2y > GMC_G2YCAP) g2y = GMC_G2YCAP;
    if (g2y < 1) g2y = 1;
    return g2y;
}

// ---------------------------------------------------------------- geometry
struct Geom {
    int W32;    // uint32 per marker row = ⌈N/16⌉
    int Npad;   // W32*16 ≥ N; x's tail is permanently zero
    int g1y;    // pass-1 sample tiles  = ⌈W32 / G_TILE_W⌉
    int g2y;    // pass-2 marker tiles  = ⌈M / G_MTILE2⌉  (upper bound)
};

Geom geom_of(int N, int M)
{
    Geom g;
    g.W32  = (N + 15)/16;
    g.Npad = g.W32*16;
    g.g1y  = (g.W32 + G_TILE_W - 1)/G_TILE_W;
    g.g2y  = (M + G_MTILE2 - 1)/G_MTILE2;
    return g;
}

}  // namespace

// ---------------------------------------------------------------- Ctx
struct Ctx {
    int N = 0, M = 0;
    int W32 = 0, Npad = 0, g1y = 0, g2y = 0;
    std::size_t bytes = 0;

    unsigned int* Ap     = nullptr;   // M × W32 uint32, packed 2-bit
    float*        freq   = nullptr;   // M
    float*        invStd = nullptr;   // M
    float*        xv     = nullptr;   // Npad (tail permanently 0)
    float*        Zv     = nullptr;   // Npad
    float*        Ypart  = nullptr;   // g1y × M
    float*        Zpart  = nullptr;   // g2y × Npad
    float*        Yv     = nullptr;   // M
    float*        wv     = nullptr;   // M
    float*        Sd     = nullptr;   // scalar, stays resident
    float*        Cd     = nullptr;   // scalar, stays resident

    // ---- fp32 small-batch scratch, allocated on the first matvec_mat() ----
    bool   mc_alloc  = false;
    int    g1y_mc    = 0;   // ⌈W32 / GMC_CHUNK1⌉
    int    g2y_cap   = 0;   // Zmcp tile capacity
    float* Xmc       = nullptr;   // Npad × 8, column-major, tail rows always 0
    float* Wmc       = nullptr;   // M × 8, marker-major [j*NC+c]
    float* Ymcp      = nullptr;   // g1y_mc × 8 × M
    float* Zmcp      = nullptr;   // g2y_cap × 8 × Npad
    float* Zmc       = nullptr;   // N × 8, column-major
    float* Smc       = nullptr;   // 8
    float* Cmc       = nullptr;   // 8
};

std::size_t need_bytes(int N, int M)
{
    const Geom g = geom_of(N, M);
    const std::size_t f = sizeof(float);
    return (std::size_t)M * g.W32 * 4            // Ap
         + 2 * (std::size_t)M * f                // freq, invStd
         + 2 * (std::size_t)g.Npad * f           // xv, Zv
         + (std::size_t)g.g1y * M * f            // Ypart
         + (std::size_t)g.g2y * g.Npad * f       // Zpart
         + 2 * (std::size_t)M * f                // Yv, wv
         + 2 * f;                                // Sd, Cd
}

int         ctx_N(const Ctx* c)     { return c ? c->N : 0; }
int         ctx_M(const Ctx* c)     { return c ? c->M : 0; }
std::size_t ctx_bytes(const Ctx* c) { return c ? c->bytes : 0; }

void destroy(Ctx* c)
{
    if (!c) return;
    if (c->Ap)     cudaFree(c->Ap);
    if (c->freq)   cudaFree(c->freq);
    if (c->invStd) cudaFree(c->invStd);
    if (c->xv)     cudaFree(c->xv);
    if (c->Zv)     cudaFree(c->Zv);
    if (c->Ypart)  cudaFree(c->Ypart);
    if (c->Zpart)  cudaFree(c->Zpart);
    if (c->Yv)     cudaFree(c->Yv);
    if (c->wv)     cudaFree(c->wv);
    if (c->Sd)     cudaFree(c->Sd);
    if (c->Cd)     cudaFree(c->Cd);
    if (c->Xmc)    cudaFree(c->Xmc);
    if (c->Wmc)    cudaFree(c->Wmc);
    if (c->Ymcp)   cudaFree(c->Ymcp);
    if (c->Zmcp)   cudaFree(c->Zmcp);
    if (c->Zmc)    cudaFree(c->Zmc);
    if (c->Smc)    cudaFree(c->Smc);
    if (c->Cmc)    cudaFree(c->Cmc);
    delete c;
}

std::size_t mc_scratch_bytes(const Ctx* c)
{
    if (!c) return 0;
    const int g1y_mc = (c->W32 + GMC_CHUNK1 - 1)/GMC_CHUNK1;
    int g2y_cap = GMC_G2YCAP;
    const int maxy = (c->M + 255)/256;
    if (g2y_cap > maxy) g2y_cap = maxy;
    const std::size_t f = sizeof(float);
    return (std::size_t)c->Npad * GMC_NCMAX * f              // Xmc
         + (std::size_t)c->M    * GMC_NCMAX * f              // Wmc
         + (std::size_t)g1y_mc  * GMC_NCMAX * c->M * f       // Ymcp
         + (std::size_t)g2y_cap * GMC_NCMAX * c->Npad * f    // Zmcp
         + (std::size_t)c->N    * GMC_NCMAX * f              // Zmc
         + 2 * (std::size_t)GMC_NCMAX * f;                   // Smc, Cmc
}

// Shared by create() and create_rows(): everything except how the packed bytes
// are copied in. Exactly one of `flat` / `row_ptrs` is non-null.
static Ctx* create_impl(const unsigned char* flat,
                        const unsigned char* const* row_ptrs,
                        std::size_t stride_bytes,
                        int N, int M, const float* freq, const float* invstd)
{
    if ((!flat && !row_ptrs) || !freq || !invstd || N <= 0 || M <= 0) return nullptr;

    const Geom g = geom_of(N, M);
    // cudaMemcpy2D would silently truncate/overrun if a source row were wider
    // than a destination row. W32*4 = ⌈N/16⌉*4 ≥ ⌈N/4⌉ always, so this only
    // trips if the caller passes a padded PackedFlat stride.
    if (stride_bytes > (std::size_t)g.W32 * 4) {
        std::fprintf(stderr,
                     "[gemv2bit] stride_bytes=%zu exceeds row pitch %zu — refusing upload\n",
                     stride_bytes, (std::size_t)g.W32 * 4);
        return nullptr;
    }

    Ctx* c   = new Ctx();
    c->N     = N;
    c->M     = M;
    c->W32   = g.W32;
    c->Npad  = g.Npad;
    c->g1y   = g.g1y;
    c->g2y   = g.g2y;
    c->bytes = need_bytes(N, M);

    const std::size_t ap_bytes = (std::size_t)M * g.W32 * 4;

    struct { void** p; std::size_t bytes; const char* name; } bufs[] = {
        {(void**)&c->Ap,     ap_bytes,                                  "Ap"},
        {(void**)&c->freq,   (std::size_t)M * sizeof(float),            "freq"},
        {(void**)&c->invStd, (std::size_t)M * sizeof(float),            "invStd"},
        {(void**)&c->xv,     (std::size_t)g.Npad * sizeof(float),       "xv"},
        {(void**)&c->Zv,     (std::size_t)g.Npad * sizeof(float),       "Zv"},
        {(void**)&c->Ypart,  (std::size_t)g.g1y * M * sizeof(float),    "Ypart"},
        {(void**)&c->Zpart,  (std::size_t)g.g2y * g.Npad * sizeof(float),"Zpart"},
        {(void**)&c->Yv,     (std::size_t)M * sizeof(float),            "Yv"},
        {(void**)&c->wv,     (std::size_t)M * sizeof(float),            "wv"},
        {(void**)&c->Sd,     sizeof(float),                             "S"},
        {(void**)&c->Cd,     sizeof(float),                             "C"},
    };
    for (auto& b : bufs) {
        cudaError_t st = cudaMalloc(b.p, b.bytes);
        if (st != cudaSuccess) {
            std::fprintf(stderr, "[gemv2bit] cudaMalloc %s (%zu B) failed: %s\n",
                         b.name, b.bytes, cudaGetErrorString(st));
            destroy(c);
            return nullptr;
        }
    }

    // Fill with 0xFF first: a padding field of 0b11 decodes to g=0, so the
    // uint32s past the end of each row contribute nothing to either pass.
    // Whatever the host wrote into the samples past N inside the LAST real
    // byte is equally harmless — xv's tail is permanently 0 (pass 1 multiplies
    // it by zero) and pass 2's Zpart tail is dropped by gm_pass2_finish.
    cudaError_t st = cudaMemset(c->Ap, 0xFF, ap_bytes);
    if (st == cudaSuccess) {
        if (flat) {
            // One strided copy straight off the caller's buffer.
            st = cudaMemcpy2D(c->Ap, (std::size_t)g.W32 * 4,
                              flat, stride_bytes,
                              stride_bytes, (std::size_t)M,
                              cudaMemcpyHostToDevice);
        } else {
            // Scattered rows: M small copies. Still no host staging buffer.
            unsigned char* dst = reinterpret_cast<unsigned char*>(c->Ap);
            for (int m = 0; m < M && st == cudaSuccess; ++m)
                st = cudaMemcpy(dst + (std::size_t)m * g.W32 * 4, row_ptrs[m],
                                stride_bytes, cudaMemcpyHostToDevice);
        }
    }
    if (st == cudaSuccess)
        st = cudaMemset(c->xv, 0, (std::size_t)g.Npad * sizeof(float));
    if (st == cudaSuccess)
        st = cudaMemcpy(c->freq, freq, (std::size_t)M * sizeof(float),
                        cudaMemcpyHostToDevice);
    if (st == cudaSuccess)
        st = cudaMemcpy(c->invStd, invstd, (std::size_t)M * sizeof(float),
                        cudaMemcpyHostToDevice);
    if (st != cudaSuccess) {
        std::fprintf(stderr, "[gemv2bit] upload failed: %s\n", cudaGetErrorString(st));
        destroy(c);
        return nullptr;
    }
    return c;
}

Ctx* create(const unsigned char* packed, std::size_t stride_bytes,
            int N, int M, const float* freq, const float* invstd)
{
    return create_impl(packed, nullptr, stride_bytes, N, M, freq, invstd);
}

Ctx* create_rows(const unsigned char* const* row_ptrs, std::size_t stride_bytes,
                 int N, int M, const float* freq, const float* invstd)
{
    return create_impl(nullptr, row_ptrs, stride_bytes, N, M, freq, invstd);
}

bool matvec_range(Ctx* c, int j0, int jn, float inv_M,
                  const float* x, float* ret)
{
    if (!c || !x || !ret) return false;
    if (jn == 0) { std::memset(ret, 0, (std::size_t)c->N * sizeof(float)); return true; }
    if (j0 < 0 || jn < 0 || j0 + jn > c->M) {
        std::fprintf(stderr, "[gemv2bit] marker range [%d, %d) out of bounds M=%d\n",
                     j0, j0 + jn, c->M);
        return false;
    }

    G2B_CHECK_FALSE(cudaMemcpy(c->xv, x, (std::size_t)c->N * sizeof(float),
                               cudaMemcpyHostToDevice));

    const int M   = jn;
    const int g1x = (M + G_MTILE1 - 1)/G_MTILE1;
    const int g2x = (c->W32 + G_NTHREAD - 1)/G_NTHREAD;
    const int g2y = (M + G_MTILE2 - 1)/G_MTILE2;

    // All six launches stay on the default stream and S / C never come back to
    // the host: this path runs on every PCG iteration, and each extra
    // device→host sync is another few tens of microseconds of bubble.
    gm_sumx<<<1, G_NTHREAD>>>(c->xv, c->N, c->Sd);
    gm_pass1<<<dim3(g1x, c->g1y), G_NTHREAD>>>(c->Ap, c->W32, c->xv, c->Ypart,
                                               j0, M, c->M);
    gm_pass1_finish<<<(M + G_NTHREAD - 1)/G_NTHREAD, G_NTHREAD>>>(
        c->Ypart, c->g1y, M, c->M, c->freq + j0, c->invStd + j0, c->Sd, c->Yv, c->wv);
    gm_calcC<<<1, G_NTHREAD>>>(c->freq + j0, c->wv, M, c->Cd);
    gm_pass2<<<dim3(g2x, g2y), G_NTHREAD>>>(c->Ap, c->W32, c->wv, c->Zpart,
                                            j0, M, c->Npad);
    gm_pass2_finish<<<(c->N + G_NTHREAD - 1)/G_NTHREAD, G_NTHREAD>>>(
        c->Zpart, g2y, c->N, c->Npad, c->Cd, inv_M, c->Zv);

    G2B_CHECK_FALSE(cudaGetLastError());   // launch config errors surface here
    G2B_CHECK_FALSE(cudaMemcpy(ret, c->Zv, (std::size_t)c->N * sizeof(float),
                               cudaMemcpyDeviceToHost));
    return true;
}

// Z = inv_M · A_std (A_stdᵀ X), X/ret column-major N × ncol, fp32 throughout.
// ncol == 1 forwards to the single-column kernels; ncol > 8 is walked in
// chunks of 8; a chunk is rounded up to {2,4,8} with ZERO columns rather than
// split (5 → 4+1): every call re-reads the whole 2-bit matrix (≈7 ms/pass for
// big2's 5.7 GB), so padding only costs arithmetic — and arithmetic is not the
// wall below ncol=8. A zero column contributes 0 to S, to C and to both
// passes, so the real columns come out bit-identical to the unpadded case.
bool matvec_mat(Ctx* c, int ncol, float inv_M, const float* X, float* ret)
{
    if (!c || !X || !ret || ncol < 0) return false;
    if (ncol == 0) return true;
    if (ncol == 1) return matvec_range(c, 0, c->M, inv_M, X, ret);

    const int N = c->N, M = c->M;

    if (!c->mc_alloc) {
        c->g1y_mc  = (c->W32 + GMC_CHUNK1 - 1)/GMC_CHUNK1;
        c->g2y_cap = GMC_G2YCAP;
        const int maxy = (M + 255)/256;
        if (c->g2y_cap > maxy) c->g2y_cap = maxy;
        struct { float** p; std::size_t bytes; const char* name; } bufs[] = {
            {&c->Xmc,  (std::size_t)c->Npad*GMC_NCMAX*sizeof(float),             "Xmc"},
            {&c->Wmc,  (std::size_t)M*GMC_NCMAX*sizeof(float),                   "Wmc"},
            {&c->Ymcp, (std::size_t)c->g1y_mc*GMC_NCMAX*M*sizeof(float),         "Ymcp"},
            {&c->Zmcp, (std::size_t)c->g2y_cap*GMC_NCMAX*c->Npad*sizeof(float),  "Zmcp"},
            {&c->Zmc,  (std::size_t)N*GMC_NCMAX*sizeof(float),                   "Zmc"},
            {&c->Smc,  GMC_NCMAX*sizeof(float),                                  "Smc"},
            {&c->Cmc,  GMC_NCMAX*sizeof(float),                                  "Cmc"},
        };
        for (auto& b : bufs) {
            cudaError_t st = cudaMalloc((void**)b.p, b.bytes);
            if (st != cudaSuccess) {
                std::fprintf(stderr, "[gemv2bit] cudaMalloc %s (mc, %zu B) failed: %s\n",
                             b.name, b.bytes, cudaGetErrorString(st));
                return false;   // caller falls back; destroy() frees the rest
            }
        }
        // Zero X once: later chunks only overwrite the first N rows of the
        // first nb columns, so the padded tail rows stay 0 forever. Padding
        // COLUMNS still need an explicit clear below (they can hold the
        // previous chunk's data).
        G2B_CHECK_FALSE(cudaMemset(c->Xmc, 0,
                                   (std::size_t)c->Npad*GMC_NCMAX*sizeof(float)));
        c->mc_alloc = true;
        std::printf("[gemv2bit] fp32 batch path: %zu MB of scratch allocated\n",
                    mc_scratch_bytes(c) >> 20);
    }

    for (int c0 = 0; c0 < ncol; c0 += GMC_NCMAX) {
        const int nb = (c0 + GMC_NCMAX <= ncol) ? GMC_NCMAX : (ncol - c0);
        const int NC = (nb <= 2) ? 2 : (nb <= 4) ? 4 : 8;

        G2B_CHECK_FALSE(cudaMemcpy2D(c->Xmc, (std::size_t)c->Npad*sizeof(float),
                                     X + (std::size_t)c0*N, (std::size_t)N*sizeof(float),
                                     (std::size_t)N*sizeof(float), (std::size_t)nb,
                                     cudaMemcpyHostToDevice));
        if (nb < NC)
            G2B_CHECK_FALSE(cudaMemset(c->Xmc + (std::size_t)nb*c->Npad, 0,
                                       (std::size_t)(NC - nb)*c->Npad*sizeof(float)));

        const int T      = NC/2;
        const int g1x    = (M + G_NWARP*GMC_MPW - 1)/(G_NWARP*GMC_MPW);
        const int g2x    = (int)(((long long)c->W32*T + G_NTHREAD - 1)/G_NTHREAD);
        const int g2y    = gmc_pick_g2y(c->W32, T, M);
        const int mtile2 = (M + g2y - 1)/g2y;

        gm_colsum_mc<<<NC, G_NTHREAD>>>(c->Xmc, c->Npad, c->Smc);
        switch (NC) {
        case 2:
            gm_pass1_mc<2><<<dim3(g1x, c->g1y_mc), G_NTHREAD>>>(c->Ap, c->W32, c->Xmc, c->Npad, c->Ymcp, 0, M, M);
            gm_pass1_finish_mc<2><<<(M*2 + 255)/256, 256>>>(c->Ymcp, c->g1y_mc, M, M, c->freq, c->invStd, c->Smc, c->Wmc);
            gm_calcC_mc<2><<<2, G_NTHREAD>>>(c->freq, c->Wmc, M, c->Cmc);
            gm_pass2_mc<2><<<dim3(g2x, g2y), G_NTHREAD>>>(c->Ap, c->W32, c->Wmc, c->Zmcp, 0, M, c->Npad, mtile2);
            gm_pass2_finish_mc<2><<<(int)(((std::size_t)N*2 + 255)/256), 256>>>(c->Zmcp, g2y, N, c->Npad, c->Cmc, inv_M, c->Zmc);
            break;
        case 4:
            gm_pass1_mc<4><<<dim3(g1x, c->g1y_mc), G_NTHREAD>>>(c->Ap, c->W32, c->Xmc, c->Npad, c->Ymcp, 0, M, M);
            gm_pass1_finish_mc<4><<<(M*4 + 255)/256, 256>>>(c->Ymcp, c->g1y_mc, M, M, c->freq, c->invStd, c->Smc, c->Wmc);
            gm_calcC_mc<4><<<4, G_NTHREAD>>>(c->freq, c->Wmc, M, c->Cmc);
            gm_pass2_mc<4><<<dim3(g2x, g2y), G_NTHREAD>>>(c->Ap, c->W32, c->Wmc, c->Zmcp, 0, M, c->Npad, mtile2);
            gm_pass2_finish_mc<4><<<(int)(((std::size_t)N*4 + 255)/256), 256>>>(c->Zmcp, g2y, N, c->Npad, c->Cmc, inv_M, c->Zmc);
            break;
        default:
            gm_pass1_mc<8><<<dim3(g1x, c->g1y_mc), G_NTHREAD>>>(c->Ap, c->W32, c->Xmc, c->Npad, c->Ymcp, 0, M, M);
            gm_pass1_finish_mc<8><<<(M*8 + 255)/256, 256>>>(c->Ymcp, c->g1y_mc, M, M, c->freq, c->invStd, c->Smc, c->Wmc);
            gm_calcC_mc<8><<<8, G_NTHREAD>>>(c->freq, c->Wmc, M, c->Cmc);
            gm_pass2_mc<8><<<dim3(g2x, g2y), G_NTHREAD>>>(c->Ap, c->W32, c->Wmc, c->Zmcp, 0, M, c->Npad, mtile2);
            gm_pass2_finish_mc<8><<<(int)(((std::size_t)N*8 + 255)/256), 256>>>(c->Zmcp, g2y, N, c->Npad, c->Cmc, inv_M, c->Zmc);
            break;
        }
        G2B_CHECK_FALSE(cudaGetLastError());

        // The first nb columns are contiguous in Zmc (column-major, ld = N).
        G2B_CHECK_FALSE(cudaMemcpy(ret + (std::size_t)c0*N, c->Zmc,
                                   (std::size_t)nb*N*sizeof(float),
                                   cudaMemcpyDeviceToHost));
    }
    return true;
}

}  // namespace saige::gpu::g2b
