// gemv2bit.hpp — internal interface for the tier-4 packed-2-bit kernels.
//
// NOT the public facade. gpu_matvec.hpp stays the only thing SAIGE_step1_fast
// includes; this header exists so gpu_matvec.cu (tier dispatch, Handle
// lifetime) and gemv2bit.cu (the kernels) can be separate translation units.
//
// Provenance: ported from the SAIGE R package's gpu-opt branch,
// src/gpuSymMatMult.{cu,hpp} (class gpuSymMatMult, 2-bit single-column path).
// Only the fp32 single-column family is here — the fp16/wmma batch path and the
// dense-fp32 cuBLAS path were deliberately left behind (see §"NOT ported").
//
// What the kernels compute, given the packed genotype matrix A (N samples ×
// M markers, 2-bit, variant-major) and per-marker (freq, invStd):
//
//   A_std[i,j] = (g[i,j] − 2 freq[j]) · invStd[j]
//   ret        = inv_M · A_std (A_stdᵀ x)
//
// The standardization is lifted OUT of the inner product by a rank-one
// decomposition, so the GEMM body only touches the raw integers {0,1,2}:
//   pass 1   Y[j] = invStd[j] · ( Σ_i g[i,j] x[i] − 2 freq[j] · S ),  S = Σ_i x[i]
//   pass 2   Z[i] = Σ_j g[i,j] w[j] − C,   w[j] = invStd[j] Y[j],
//                                          C = Σ_j 2 freq[j] w[j]
// Both correction terms are O(M) / O(1), so compression costs zero extra error.
//
// Encoding contract (must match the host packer byte-for-byte):
//   0b00 HOM_REF → g=2,  0b10 HET → g=1,  0b11 HOM_ALT → g=0,  0b01 MISSING
//   never appears (marker_decoder.cpp fills missing with round(2·altFreq)
//   BEFORE packing). The kernel computes g = 2 − popc(field), which would map
//   0b01 to g=1 — silently disagreeing with the CPU path's round(2f). See the
//   assertion note in create().
//
// Reductions use partial buffers, never atomicAdd: float addition order is
// therefore fixed and two runs are bit-identical.
#pragma once

#include <cstddef>

namespace saige::gpu::g2b {

struct Ctx;   // opaque; defined in gemv2bit.cu

// Device bytes tier 4 needs for this problem size, for the tier decision in
// gpu_matvec.cu. Pure arithmetic, touches no CUDA state.
std::size_t need_bytes(int N, int M);

// Upload the packed bytes + per-marker vectors and allocate all device
// scratch. `packed` is M rows × stride_bytes bytes, variant-major — pass
// PackedFlat::raw() / PackedFlat::nbyte() directly, no gather buffer.
// Inputs are NOT retained: the caller may free them on return.
// Returns nullptr on any CUDA failure (caller falls back to a lower tier/CPU).
Ctx* create(const unsigned char* packed, std::size_t stride_bytes,
            int N, int M, const float* freq, const float* invstd);

// Same, but the marker rows are scattered: row_ptrs[m] points at marker m's
// `stride_bytes` packed bytes. Uploads row by row (M small H2D copies, once)
// rather than gathering into a staging buffer first — so the host-side peak
// stays exactly the size of the data that is already there.
// Needed by SAIGE's legacy genoVecofPointers storage, which is one heap
// allocation per marker and has no contiguous buffer to hand over.
Ctx* create_rows(const unsigned char* const* row_ptrs, std::size_t stride_bytes,
                 int N, int M, const float* freq, const float* invstd);

// ret = inv_M · B_std (B_stdᵀ x) over the marker range [j0, j0+jn).
// x and ret are host float32, length N. inv_M is applied by the last kernel,
// so no extra host-side pass over N.
// jn == M and j0 == 0 is the full-GRM case the PCG iteration walks.
bool matvec_range(Ctx* c, int j0, int jn, float inv_M,
                  const float* x, float* ret);

// Multi-RHS analogue: ret = inv_M · A_std (A_stdᵀ X) over ALL markers.
// X and ret are host float32, column-major N × ncol (an arma::fmat's memptr()
// drops straight in). ncol is unbounded — the implementation walks it in
// chunks of 8 and rounds each chunk up to {2,4,8} with zero columns. ncol == 1
// forwards to matvec_range.
// Device scratch for this path is allocated lazily on the first call, so a run
// that never batches pays nothing for it.
bool matvec_mat(Ctx* c, int ncol, float inv_M, const float* X, float* ret);

// Device bytes matvec_mat() will lazily allocate on first use, for logging.
std::size_t mc_scratch_bytes(const Ctx* c);

void destroy(Ctx* c);

// Reflection for logging / the memory report.
int  ctx_N(const Ctx* c);
int  ctx_M(const Ctx* c);
std::size_t ctx_bytes(const Ctx* c);

}  // namespace saige::gpu::g2b
