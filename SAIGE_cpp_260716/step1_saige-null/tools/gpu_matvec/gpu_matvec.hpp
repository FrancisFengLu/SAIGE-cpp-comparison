// gpu_matvec.hpp — opaque facade for CUDA-accelerated K·u.
//
// Lifecycle:
//   Handle* h = saige::gpu::create(packed, freq, invstd, N);
//   if (!h)  // not built with CUDA, no device, or OOM → CPU fallback caller-side.
//   for (each matvec) saige::gpu::matvec(h, u, Au);
//   saige::gpu::destroy(h);
//
// The cu implementation lives in gpu_matvec.cu (compiled with nvcc).
// The cpu_stub.cpp version below returns nullptr/false for everything so
// the same translation units link cleanly when CUDA isn't available.
//
// Mathematical contract (same for both backends):
//   Given N samples × M_pass markers of standardized genotype data A,
//   matvec() computes out_Au = (A · Aᵀ / M_pass) · u = K · u.
//   out_Au and u are length N, float32, host memory.
//
// Numerics: single precision throughout, matches CPU parallelCrossProd's
// final output. All tiers reduce through partial buffers, never atomicAdd, so
// for a fixed tier two runs of the same binary are bit-identical. Switching
// tiers changes the summation order and therefore the last few ulps.
#pragma once

#include <cstddef>
#include <vector>

namespace saige { class PackedFlat; }

namespace saige::gpu {

struct Handle;   // opaque

// True if CUDA runtime is present AND at least one usable device exists.
bool available();

// Upload packed 2-bit genotypes + per-marker (freq, invstd) arrays to the
// device and preconfigure the matvec kernels.
//   packed        : PackedFlat from setGenoObj (n_stored × nbyte-per-marker)
//   freq          : length packed.n_stored() — alt allele frequency per marker
//   invstd        : length packed.n_stored() — 1/√(2 f (1−f)), 0 where f ∈ {0,1}
//   N             : number of phenotyped samples (packed.nbyte() ≥ ⌈N/4⌉)
//   tier_override : 0 (auto — prefers 4, then 3, then 1/2)
//                   1 = stream fp32 blocks from host per matvec (cuBLAS)
//                   3 = packed-resident, standardize per element in-kernel
//                   4 = packed-resident, rank-one standardization (gemv2bit.cu)
//                   (2 is not selectable: it is what tier 1 becomes when the
//                    whole fp32 A happens to fit in one resident block.)
// Returns nullptr on failure (no CUDA, OOM, device error). Inputs are NOT
// retained after return; caller can free them — EXCEPT for tiers 1/2, which
// keep reading `packed` on every matvec (tiers 3/4 do not).
Handle* create(const saige::PackedFlat& packed,
               const std::vector<float>& freq,
               const std::vector<float>& invstd,
               int                       N,
               int                       tier_override = 0);

// Compute out_Au = K · u. Returns false on error — caller falls back to CPU.
// Safe to call concurrently with other handles; NOT safe to call two matvecs
// on the same Handle concurrently from different threads.
bool matvec(Handle* h, const float* u, float* out_Au);

// Release all device buffers held by h.
void destroy(Handle* h);

// Reflection: which tier the handle actually runs (1/2/3/4, see create()).
// 0 means "not created".
int  tier(const Handle* h);

}  // namespace saige::gpu
