// gpu_step2.hpp — the ONLY header main.cpp includes for the step-2 GPU path.
//
// No CUDA types cross this boundary, so main.cpp compiles identically whether
// or not the build has CUDA. With USE_CUDA unset the implementation is
// gpu_step2_stub.cpp, whose available() always says "built without CUDA" and
// whose create() always returns nullptr — every caller then takes the CPU path
// it would have taken anyway.
//
// ---------------------------------------------------------------------------
// What this computes, and why exactly this
// ---------------------------------------------------------------------------
// SAIGE::scoreTestBatchMT (saige_mt.cpp) reduces a block of markers against the
// stacked per-trait constants. For a run in which EVERY trait is quantitative
// and every trait's sample list is the union's, the only things it touches in
// sample space are
//
//     Zall  = Astack^T G      (sumP x B)
//     GWqnt = Xstack^T G      (sumP x B)
//     GR    = G^T RES         (B x P)
//     Gsq   = colsum(G % G)   (B)
//
// and everything after that is O(p^2) and O(P) per marker. So one GEMM
//
//     C = G^T B,   B = [ Astack | Xstack | RES ]   (N x K, K = 2*sumP + P)
//
// is the entire sample-space cost, and C's columns split back into the three
// blocks above. Gsq is NOT computed here: the host already holds the marker's
// 2-bit code counts, so sum_i g_i^2 = sum_c count[c]*fd[c]^2 exactly, in
// double, for free (see mainMarkerMT). Keeping it off the GPU keeps it out of
// the fp32 error budget.
//
// ---------------------------------------------------------------------------
// The genotype matrix never becomes fp32 on the host
// ---------------------------------------------------------------------------
// The caller hands over the marker's PLINK 2-bit codes verbatim (packed(), 4
// samples per byte, sample i in bits 2*(i&3) of byte i>>2) plus a 4-entry
// code -> dosage table (lut(), from PLINK::finalizeFusedStats::fd). The decode
// kernel applies that table, so allele flip, missing-genotype imputation and
// the MAC-gated .clean() are all already baked in and the GPU reproduces the
// CPU's dosages exactly, cell for cell, with no separate missing-value path.
// That is the whole trick: 1/16 of the PCIe traffic of shipping fp32 dosages,
// and imputation is a per-marker constant rather than a per-cell special case.
//
// ---------------------------------------------------------------------------
// Streaming, not resident
// ---------------------------------------------------------------------------
// Step 2 reads every marker exactly once, so there is nothing for a resident
// matrix to amortise: uploading 12.5 GB once (3.7 s measured, V100 pinned-less)
// and then computing costs strictly more than streaming the same bytes chunk by
// chunk behind the compute. This reducer therefore always streams, in
// slots-per-pass batches the caller chooses, double-buffered across two streams
// so the H2D of batch k+1 overlaps the decode+GEMM of batch k. The
// "does the packed matrix fit in device memory" question does not arise.
#pragma once

#include <cstddef>
#include <string>

namespace saige {
namespace gpu2 {

struct Reducer;   // opaque; defined in gpu_step2.cu

// True if a usable device exists. `why` (optional) gets a one-line reason when
// the answer is false, for the startup line main() prints.
bool available(int t_device, std::string* t_why);

// One-line device description ("Tesla V100-SXM2-16GB, 16384 MiB, sm_70"), or
// "" when unavailable. For the startup banner only.
std::string describe(int t_device);

// Build the reducer.
//   t_N         samples
//   t_K         columns of B
//   t_B         N x K column-major DOUBLE (the caller's stacked constants as
//               they already are); narrowed here when t_fp64 is false. Copied
//               to the device, so it may be freed on return.
//   t_maxSlots  markers per reduce() call, at most
//   t_fp64      run the decode and the GEMM in double instead of float.
//               On a V100 this costs about 40% more GEMM time and about one
//               extra pass of device bandwidth in the decode -- a few percent
//               of a real step-2 run -- and it removes the precision question
//               rather than bounding it. See PRECISION in gpu_step2.cu.
// Returns nullptr on ANY failure (no device, allocation refused, ...). The
// caller must fall back to the CPU path; nothing is printed here.
Reducer* create(int t_device, int t_N, int t_K, const double* t_B, int t_maxSlots,
                bool t_fp64);
void     destroy(Reducer* t_r);

// Pinned staging buffers owned by the reducer. The caller writes marker data
// straight into them, so no copy stands between the decode and the H2D.
//   packed()  t_maxSlots * bpv bytes, bpv = (N+3)/4; slot j at offset j*bpv
//   lut()     t_maxSlots * 4 floats;  slot j's code->dosage table at 4*j
unsigned char* packed(Reducer* t_r);
float*         lut(Reducer* t_r);
std::size_t    bytesPerSlot(const Reducer* t_r);

// Reduce slots [0, t_nSlots). Returns false on any CUDA error, in which case
// the results are undefined and the caller must fall back for this batch.
bool reduce(Reducer* t_r, int t_nSlots);

// Results of the last reduce(): C is t_maxSlots x K column-major, so element
// (slot, k) sits at index (std::size_t)k * ldC() + slot. Pinned, owned here.
// Exactly one of the two accessors is non-null, per isFp64().
bool          isFp64(const Reducer* t_r);
const float*  outCf(const Reducer* t_r);
const double* outCd(const Reducer* t_r);
std::size_t   ldC(const Reducer* t_r);         // == maxSlots

// Cumulative device-side timings in seconds since create(), for the end-to-end
// breakdown. Measured with events on the compute stream, so they are the GPU's
// own numbers and not wall clock around the call.
void timings(const Reducer* t_r, double* t_h2d, double* t_decode, double* t_gemm,
             double* t_d2h);

// Peak device bytes allocated by this reducer.
std::size_t deviceBytes(const Reducer* t_r);

}  // namespace gpu2
}  // namespace saige
