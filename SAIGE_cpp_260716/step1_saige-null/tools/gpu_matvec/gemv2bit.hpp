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
// therefore fixed and two runs are bit-identical. The scheme-C fill-correction
// scatter-adds below preserve that — see the TraitBind comment.
//
// ---------------------------------------------------------------------------
// Scheme C (optimization/missing_mt/SCHEME_C_DESIGN.md §3.3): the matrix and
// the standardization are separate objects. One Ctx holds the packed genotypes
// of the UNION of several phenotypes' sample sets; one TraitBind per phenotype
// holds that phenotype's freq / invStd / 1/M_t, the union-local rows it does
// NOT own, and the cells where its missing-value fill differs from the union's.
// Switching traits is a pointer swap — the matrix is never re-uploaded.
//
// With a bind whose mask and correction lists are both empty the kernels run
// exactly the arithmetic above; that degenerate case is bit-identical to the
// pre-scheme-C code (acceptance B1, baseline in b1_baseline_manifest.txt).
#pragma once

#include <cstddef>

namespace saige::gpu::g2b {

struct Ctx;         // opaque; defined in gemv2bit.cu
struct TraitBind;   // opaque; defined in gemv2bit.cu

// Device bytes tier 4 needs for this problem size, for the tier decision in
// gpu_matvec.cu. Pure arithmetic, touches no CUDA state.
// Includes 2·M floats that used to be Ctx's freq/invStd and now belong to the
// first TraitBind — the figure is unchanged on purpose so the tier threshold
// does not move, and it stays an honest estimate for the single-trait case.
// Each ADDITIONAL bind costs bind_bytes() on top.
std::size_t need_bytes(int N, int M);

// Device bytes one bind_trait() allocates. n_mask/n_corr are the list lengths.
std::size_t bind_bytes(int N, int M, int n_mask, int n_corr);

// Upload the packed bytes and allocate all device scratch. `packed` is M rows
// × stride_bytes bytes, variant-major — pass PackedFlat::raw() /
// PackedFlat::nbyte() directly, no gather buffer.
// Inputs are NOT retained: the caller may free them on return.
// Returns nullptr on any CUDA failure (caller falls back to a lower tier/CPU).
// No freq/invStd here any more: see bind_trait().
Ctx* create(const unsigned char* packed, std::size_t stride_bytes,
            int N, int M);

// Same, but the marker rows are scattered: row_ptrs[m] points at marker m's
// `stride_bytes` packed bytes. Uploads row by row (M small H2D copies, once)
// rather than gathering into a staging buffer first — so the host-side peak
// stays exactly the size of the data that is already there.
// Needed by SAIGE's legacy genoVecofPointers storage, which is one heap
// allocation per marker and has no contiguous buffer to hand over.
Ctx* create_rows(const unsigned char* const* row_ptrs, std::size_t stride_bytes,
                 int N, int M);

// Bind one phenotype's standardization + row mask + fill corrections to `c`.
// Several binds may coexist on one Ctx; switching between them costs nothing.
// All inputs are copied to the device and may be freed on return.
//
//   freq, invstd  length M (the Ctx's full marker count). invstd[j] == 0 marks
//                 a marker this phenotype's QC dropped: it contributes nothing
//                 to either pass, exactly as §1 requires.
//   mask_rows     n_mask union-local sample indices in [0,N) that this
//                 phenotype does NOT own. The kernels zero x on those rows ON
//                 THE DEVICE after upload (the host's PCG vectors are shared
//                 across the group, so the caller must not have to pre-mask),
//                 and zero `ret` on them before copying back — see below.
//   corr_*        n_corr triplets (row, col, delta): at union-local sample
//                 `row`, marker `col`, this phenotype's fill differs from the
//                 union's by `delta` ∈ {−2,−1,1,2}.
//                 ORDER AND UNIQUENESS ARE REQUIRED, not advisory: the list
//                 must be STRICTLY increasing in (col, row). Strictly, because
//                 Δ = fill_t − fill_U is one value per cell — a repeated
//                 (row,col) would be scatter-added twice and is a builder bug.
//                 bind_trait checks this on the host and refuses, rather than
//                 letting a bad segmentation produce quiet wrong numbers.
//
// NOTE what a bind deliberately does NOT hold: 1/M_t. The marker count that
// normalizes the result is per CALL, not per phenotype — a LOCO slice's M_t is
// the number of markers in [j0, j0+jn) passing THIS phenotype's QC, which the
// host already knows. Keeping it out means one bind per phenotype serves every
// chromosome instead of 23 binds each carrying a duplicate copy of the CSRs.
//
// Returns nullptr on bad arguments or any CUDA failure.
//
// Determinism: the corrections are turned into two CSR views at bind time —
// one indexed by marker (pass 1), one by sample (pass 2) — so each scatter-add
// is summed by a single thread walking a contiguous, sorted segment. No
// atomicAdd anywhere, and two runs stay bit-identical.
//
// ret ON MASKED ROWS IS ZERO — a REQUIREMENT of the contract, not an
// implementation detail, and callers may rely on it. §1 of the design defines
// Z[i] for every i of the union; on a row this phenotype does not own that is a
// finite but meaningless number (that row of A dotted with w). Handing it back
// would put junk into the caller's SHARED PCG vectors, where the damage is done
// by the host-side reductions — r·r, r·z — before the next matvec's re-mask of
// x can undo anything. Zeroing is what makes the union result substitutable for
// the single-phenotype one.
TraitBind* bind_trait(Ctx* c,
                      const float* freq, const float* invstd,
                      const int* mask_rows, int n_mask,
                      const int* corr_row, const int* corr_col,
                      const float* corr_delta, int n_corr);

void unbind_trait(TraitBind* t);

// ret = inv_M · B_std (B_stdᵀ x) over the marker range [j0, j0+jn), under the
// standardization/mask/corrections of `t`. x and ret are host float32, length
// N (the UNION's N — x is masked on the device, ret comes back zeroed on the
// masked rows). jn == M and j0 == 0 is the full-GRM case the PCG walks.
//
// inv_M is 1/M_t for THIS range: the number of markers in [j0, j0+jn) that pass
// this phenotype's QC. It is a call argument and not part of the bind exactly
// so a LOCO run can reuse one bind for all 23 slices. Applied by the last
// kernel, so no extra host pass over N.
//
// Corrections whose marker falls outside [j0, j0+jn) are SKIPPED — in both
// passes, and by construction: pass 1 gets the marker CSR offset to j0 so out
// of range segments are never visited, and pass 2 tests (col − j0) ∈ [0, jn)
// because its CSR stores a global marker index while w is range-local. A LOCO
// slice therefore agrees with the corresponding part of the full call.
bool matvec_range(Ctx* c, TraitBind* t, int j0, int jn, float inv_M,
                  const float* x, float* ret);

// Multi-RHS analogue: ret = inv_M · A_std (A_stdᵀ X) over ALL markers.
// X and ret are host float32, column-major N × ncol (an arma::fmat's memptr()
// drops straight in). ncol is unbounded — the implementation walks it in
// chunks of 8 and rounds each chunk up to {2,4,8} with zero columns. ncol == 1
// forwards to matvec_range.
// Every column is masked independently on the device; every output column
// comes back zeroed on the masked rows. inv_M is 1/M_t over ALL markers, same
// per-call rule as matvec_range.
// Device scratch for this path is allocated lazily on the first call, so a run
// that never batches pays nothing for it.
bool matvec_mat(Ctx* c, TraitBind* t, int ncol, float inv_M,
                const float* X, float* ret);

// Device bytes matvec_mat() will lazily allocate on first use, for logging.
std::size_t mc_scratch_bytes(const Ctx* c);

void destroy(Ctx* c);

// Reflection for logging / the memory report.
int  ctx_N(const Ctx* c);
int  ctx_M(const Ctx* c);
std::size_t ctx_bytes(const Ctx* c);

}  // namespace saige::gpu::g2b
