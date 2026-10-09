// gpu_scan_lowp.cuh — device kernels of the scan stage's fp32 and int8 modes
// (gpu_step2.cu, CreateArgs::precision). Included by gpu_step2.cu only; the
// fp64 path does not use anything here.
//
// FP32 ("chunked fp32, shifted and centred, fp64 accumulation")
// --------------------------------------------------------------
// The fp32 error of C = G^T B is mostly accumulation: one SGEMM sums a K = N
// long dot product in fp32, so the error scales with eps32 * sum_i |g_i b_ik|,
// far above |C| for a sum of nearly cancelling terms (a residual column). Three
// rewrites, each undone exactly in fp64, shrink it while the GEMMs stay fp32:
//   * K chunking: the samples are cut into kChunkF32 pieces (C2: kChunkF32Sq),
//     one SGEMM each into an fp32 scratch, the pieces summed into the fp64
//     result in fixed order. fp32 then only ever sums one chunk.
//   * integer shift: the decode writes g_i - r, r = rint(mean_i g_i) per
//     marker (an integer, so g - r is exact in fp32 for every hard call).
//   * column centring: a B column whose entries cluster around a nonzero mean
//     (|mean| > sd: mu2, W X's intercept) goes up as b - mean. Other columns
//     are not centred: that would cost their small entries -- all a rare
//     marker reads -- their relative precision.
//   sum_i g_i b_ik = sum_i (g_i - r)(b_ik - mu_k) + r colsum_k + mu_k sum_i (g_i - r),
// with colsum_k, mu_k from the fp64 B and sum_i g_i from the row's code
// counts, all in fp64. Measured on bingpu_test bt full against fp64: Tstat
// error 9e-6 -> 1.1e-6 of its sd, var 2.7e-6 -> 3e-8, SE 2.5e-4 -> 8e-7
// (relative); see PRECISION in gpu_step2.cu for what is left.
//
// INT8 (Ozaki-style split)
// ------------------------
// B's column k is written exactly-to-2^-(7k) as 2^(E_k-6) sum_s d_s 2^(-7s),
// |d_s| <= 64, d_s int8 (E_k from frexp of max |B_ik|, so |B 2^(6-E_k)| < 64;
// every step of the split is exact in fp64). The marker side is integer: each
// 2-bit code's table value v_c is written n_c + delta_c with n_c = rint(v_c);
// n_c goes into the int8 column, and every code with delta_c != 0 (in practice
// only the missing call's mean-imputed 2*altFreq) gets one extra 0/1 indicator
// column. C = sum_s 2^(E_k-6-7s) (G_int' d_s + sum_x delta_x I_x' d_s), the
// int32 products recombined in fp64. Overflow: |A| <= 4 (squared table), |d|
// <= 64, so |any partial sum| <= 256 N < 2^31 for N < 8,388,608; create()
// refuses larger N and reduce() refuses a table with |n_c| > 4.
#pragma once

#include <cstdint>

namespace saige {
namespace gpu2 {
namespace lowp {

constexpr int kChunkF32   = 4096;   // fp32 samples per partial sum of C1 (the last chunk takes < 2x)
constexpr int kChunkF32Sq = 16384;  // the same for C2 = (G % G)' B2: its terms are all >= 0
                                    // (mu2, diag Sigma^-1), so there is no cancellation to
                                    // protect, and narrow-K GEMMs at B2's few columns run
                                    // at poor efficiency
constexpr int kI8MaxA   = 4;      // largest |n_c| the int8 overflow bound allows

// ---- FP32 -----------------------------------------------------------------

// Code counts of one row (samples [0, N)), block-wide; result in sc[0..3].
__device__ __forceinline__ void rowCounts(const uint8_t* __restrict__ row, int N, unsigned* sc)
{
    if (threadIdx.x < 4) sc[threadIdx.x] = 0u;
    __syncthreads();
    unsigned c[4] = {0u, 0u, 0u, 0u};
    const int nq = N >> 2;
    for (int b = threadIdx.x; b < nq; b += blockDim.x) {
        const unsigned p = row[b];
        c[p & 3u]++; c[(p >> 2) & 3u]++; c[(p >> 4) & 3u]++; c[(p >> 6) & 3u]++;
    }
    for (int i = (nq << 2) + threadIdx.x; i < N; i += blockDim.x)
        c[(row[i >> 2] >> ((i & 3) * 2)) & 3u]++;
    for (int k = 0; k < 4; k++) if (c[k]) atomicAdd(&sc[k], c[k]);
    __syncthreads();
}

// One block per marker: shift r = rint(mean of the table over the row's
// codes), then g_i - r as float. shiftOut[2m] = r (an integer),
// shiftOut[2m+1] = sum_i (g_i - r), both in fp64 from the code counts.
__global__ void __launch_bounds__(256)
decode_f32_shift(const uint8_t* __restrict__ packed, std::size_t bpv, int N,
                 const double* __restrict__ lut, float* __restrict__ dG,
                 double* __restrict__ shiftOut)
{
    __shared__ unsigned sc[4];
    __shared__ double sh;
    const int m = blockIdx.x;
    const uint8_t* __restrict__ row = packed + (std::size_t)m * bpv;
    const double l0 = lut[4 * m], l1 = lut[4 * m + 1], l2 = lut[4 * m + 2], l3 = lut[4 * m + 3];
    rowCounts(row, N, sc);
    if (threadIdx.x == 0) {
        const double mean = ((double)sc[0] * l0 + (double)sc[1] * l1 +
                             (double)sc[2] * l2 + (double)sc[3] * l3) / (double)N;
        sh = rint(mean);
        const double gsum = (double)sc[0] * l0 + (double)sc[1] * l1 + (double)sc[2] * l2 + (double)sc[3] * l3;
        shiftOut[2 * m] = sh;
        shiftOut[2 * m + 1] = gsum - (double)N * sh;      // sum_i (g_i - r)
    }
    __syncthreads();
    const double r = sh;
    const float v0 = (float)(l0 - r), v1 = (float)(l1 - r), v2 = (float)(l2 - r), v3 = (float)(l3 - r);
    float* __restrict__ o = dG + (std::size_t)m * N;
    auto pk = [&](unsigned c) -> float {
        const float a = (c & 1u) ? v1 : v0;
        const float b = (c & 1u) ? v3 : v2;
        return (c & 2u) ? b : a;
    };
    if ((N & 3) == 0) {
        const int nq = N >> 2;
        for (int b = threadIdx.x; b < nq; b += blockDim.x) {
            const unsigned p = row[b];
            *reinterpret_cast<float4*>(o + 4 * b) =
                make_float4(pk(p & 3u), pk((p >> 2) & 3u), pk((p >> 4) & 3u), pk((p >> 6) & 3u));
        }
    } else {
        for (int i = threadIdx.x; i < N; i += blockDim.x)
            o[i] = pk((unsigned)((row[i >> 2] >> ((i & 3) * 2)) & 3u));
    }
}

// C[m + k*ldc] (=|+=) part[m + k*sc]; the first chunk also adds
// r_m colsum_k + (sum_i (g_i - r_m)) mean_k, shift = (r_m, sum_i (g_i - r_m))
// per slot and colsum = (colsum_k, mean_k): sum_i g_i b_ik =
// sum_i (g_i - r)(b_ik - mean_k) + r colsum_k + mean_k sum_i (g_i - r).
// Grid: (x over m, y over k); see accGrid in gpu_step2.cu.
__global__ void acc_f32(const float* __restrict__ part, int sc, int K,
                        double* __restrict__ C, std::size_t ldc, int first,
                        const double* __restrict__ shift, const double* __restrict__ colsum)
{
    for (int k = blockIdx.y; k < K; k += gridDim.y) {
        const float* __restrict__ pc = part + (std::size_t)k * sc;
        double* __restrict__ cc = C + (std::size_t)k * ldc;
        const double cs = colsum[k], mu = colsum[K + k];
        for (int m = blockIdx.x * blockDim.x + threadIdx.x; m < sc; m += gridDim.x * blockDim.x) {
            const double v = (double)pc[m];
            if (first) cc[m] = shift[2 * m] * cs + shift[2 * m + 1] * mu + v;
            else       cc[m] += v;
        }
    }
}

// ---- INT8 -----------------------------------------------------------------

// One block per A column (Np int8 values, samples >= N are 0). Columns
// [0, sc): slot j through rint(lut); [sc, sc+nx): indicator of code xc[x] on
// slot xs[x]; [sc+nx, sc+nx+sc2): slot j-sc-nx through rint(lut2) (sc2 = 0 or sc).
__global__ void __launch_bounds__(256)
decode_i8(const uint8_t* __restrict__ packed, std::size_t bpv, int N, int Np,
          const double* __restrict__ lut, const double* __restrict__ lut2,
          int sc, int nx, const int* __restrict__ xs, const int* __restrict__ xc,
          int8_t* __restrict__ A)
{
    const int j = blockIdx.x;
    int slot; int t[4];
    if (j < sc) {
        slot = j;
        for (int c = 0; c < 4; c++) t[c] = __double2int_rn(lut[4 * slot + c]);
    } else if (j < sc + nx) {
        const int x = j - sc;
        slot = xs[x];
        for (int c = 0; c < 4; c++) t[c] = (c == xc[x]) ? 1 : 0;
    } else {
        slot = j - sc - nx;
        for (int c = 0; c < 4; c++) t[c] = __double2int_rn(lut2[4 * slot + c]);
    }
    const uint8_t* __restrict__ row = packed + (std::size_t)slot * bpv;
    int8_t* __restrict__ o = A + (std::size_t)j * Np;
    const int nq = Np >> 2;
    for (int b = threadIdx.x; b < nq; b += blockDim.x) {
        int v[4];
        const unsigned p = (4 * b < N) ? row[b] : 0u;
        for (int q = 0; q < 4; q++) v[q] = (4 * b + q < N) ? t[(p >> (2 * q)) & 3u] : 0;
        *reinterpret_cast<char4*>(o + 4 * b) = make_char4((char)v[0], (char)v[1], (char)v[2], (char)v[3]);
    }
}

// One slice s of the recombination: rows [rowOff, rowOff+sc) of the int32
// product are the slots, rows [xOff, xOff+nx) the indicator columns.
//   C[m + k*ldc] (=|+=) 2^(E[k]-6-7s) * (P[slot m] + sum_x delta[x] P[x])
__global__ void acc_i8(const int* __restrict__ P, int rowsP, int rowOff, int xOff,
                       int sc, int K, const int* __restrict__ E, int s,
                       double* __restrict__ C, std::size_t ldc, int first,
                       const int* __restrict__ xFirst, const int* __restrict__ xN,
                       const double* __restrict__ delta)
{
    for (int k = blockIdx.y; k < K; k += gridDim.y) {
        const int* __restrict__ pc = P + (std::size_t)k * rowsP;
        double* __restrict__ cc = C + (std::size_t)k * ldc;
        const int ex = E[k] - 6 - 7 * s;
        for (int m = blockIdx.x * blockDim.x + threadIdx.x; m < sc; m += gridDim.x * blockDim.x) {
            double v = (double)pc[rowOff + m];
            const int x0 = xFirst[m], nxm = xN[m];
            for (int x = x0; x < x0 + nxm; x++) v += delta[x] * (double)pc[xOff + x];
            v = ldexp(v, ex);
            if (first) cc[m] = v;
            else       cc[m] += v;
        }
    }
}

}  // namespace lowp
}  // namespace gpu2
}  // namespace saige
