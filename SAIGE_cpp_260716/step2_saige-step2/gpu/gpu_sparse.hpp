// gpu_sparse.hpp — the within-block cross terms of g' Sigma^-1 g on the device,
// for the sparse-GRM variance of binary traits on the step-2 GPU path. Config
// key gpuSparse (needs gpuBinary); default off. No CUDA types cross this
// boundary; the no-CUDA build links gpu_sparse_stub.cpp and spqCreate()
// returns nullptr.
//
// ---------------------------------------------------------------------------
// What it computes
// ---------------------------------------------------------------------------
// With a sparse GRM, step 2's exact variance is var2 = g~' Sigma^-1 g~ with
// g~ = g - X z, z = XVX_inv_XV g (SAIGEClass::scoreTest). Sigma is block
// diagonal after permuting by the connected components of the GRM, so with
// B = Sigma^-1 (block diagonal, held explicitly) the variance expands into
//
//     var2 = g'Bg - 2 z'(BX)'g + z'(X'BX)z
//
// exactly as the dense kernel expands g~' W g~ (saige_mt.cpp). (BX)'g is one
// more block of columns in the reducer's first GEMM, z'(X'BX)z is a p x p
// contraction on the host, and g'Bg splits into
//
//     sum_i B_ii g_i^2                 one more column of the second GEMM
//                                      ((G % G)' Bdiag)
//   + sum_{i<j, same block} 2 B_ij g_i g_j      this module
//
// The second term is what is left of the per-marker Sigma solve: a quadratic
// form over the within-block pairs (i, j), one pass per marker column, all
// traits at once (the pairs are the partition's, shared by every trait; the
// weights 2 B_ij are per trait). One CUDA block per marker slot reads the
// slot's 2-bit column from the reducer's resident rows (gpu_step2.hpp
// devicePacked / deviceLut), so g is the very dosage the GEMMs decoded. A pair
// whose g_i g_j is zero -- nearly all of them on a rare marker -- costs two
// byte loads and nothing else. fp64; the sum order is fixed by (thread, pair)
// and a fixed tree, so a run repeats itself bit for bit.
#pragma once

#include <cstddef>
#include <cstdint>

namespace saige {
namespace gpu2 {

struct Reducer;
struct SpQuad;   // opaque; defined in gpu_sparse.cu

struct SpQuadCreateArgs {
    int device = 0;
    int N = 0;
    int nTraits = 0;                 // weight rows; output columns per slot
    long long nPairs = 0;            // within-block pairs i < j (may be 0)
    const int* pi = nullptr;         // nPairs sample indices
    const int* pj = nullptr;         // nPairs sample indices
    const double* w = nullptr;       // nTraits x nPairs, trait-major: w[t*nPairs + k]
    int maxSlots = 0;                // per spqRun() call
    // Pinned result sets (gpuOverlap): spqRun(..., t_set) writes set t_set and
    // reads the reducer's device set t_set. 1 = the original single buffer.
    int outSets = 1;
};

// nullptr on any failure; the caller then keeps the sparse variance on the CPU.
SpQuad* spqCreate(const SpQuadCreateArgs& t_args);
void    spqDestroy(SpQuad* t_q);

// out[slot * nTraits + t] = sum_k w[t][k] g_slot[pi[k]] g_slot[pj[k]] for slots
// [0, t_nSlots) of the reducer's last reduce(). Synchronous; the result is in
// pinned host memory (spqOut) until the next call. false on a CUDA error.
bool          spqRun(SpQuad* t_q, const Reducer* t_r, int t_nSlots, int t_set = -1);

// Different sample lists (gpuOwnSampleSets): the same sums on each trait's OWN
// genotype vector, g_t[i] = t_lut[(slot * nTraits + t) * 4 + code_i] -- the
// trait's code -> dosage table (its own imputation and flip). The pair list is
// over union indices and a trait's weights are zero on pairs outside its
// samples, so samples outside the trait contribute nothing. t_lut is host
// memory, nSlots x nTraits x 4, copied in by the call.
bool          spqRunOwn(SpQuad* t_q, const Reducer* t_r, int t_nSlots, const double* t_lut, int t_set = -1);
// t_set -1 = the set of the last call.
const double* spqOut(const SpQuad* t_q, int t_set = -1);

// Cumulative device seconds in the kernel, and slots processed.
void        spqTimings(const SpQuad* t_q, double* t_kernel, long long* t_slots);
std::size_t spqDeviceBytes(const SpQuad* t_q);
const char* spqLastError();

}  // namespace gpu2
}  // namespace saige
