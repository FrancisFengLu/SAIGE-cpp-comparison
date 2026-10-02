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
// stacked per-trait constants. For traits whose sample list is the union's,
// the only things it touches in sample space are
//
//     Zall  = Astack^T G      (sumP x B)       every trait
//     GWqnt = Xstack^T G      (sumPqnt x B)    quantitative traits
//     GWbin = WXstack^T G     (sumPbin x B)    binary traits  (WXstack = mu2 % X)
//     GR    = G^T RES         (B x P)          every trait
//     Gsq   = colsum(G % G)   (B)              quantitative traits
//     G2Mu2 = (G % G)^T MU2bin (B x nBin)      binary traits
//
// and everything after that is O(p^2) and O(P) per marker. So two GEMMs
//
//     C1 = G^T B1,         B1 = [ Astack | WXstack | Xstack_qnt | RES ]   (N x K1)
//     C2 = (G % G)^T B2,   B2 = MU2bin                                   (N x K2)
//
// are the entire sample-space cost; C1's columns split back into the first
// four blocks above and C2 is G2Mu2. Gsq is NOT computed here: the host already
// holds the marker's 2-bit code counts, so sum_i g_i^2 = sum_c count[c]*fd[c]^2
// exactly, in double, for free (see mainMarkerMTGpu). G2Mu2 cannot be had that
// way (it is weighted by mu2 per sample), so for binary traits the marker is
// decoded a second time through the SQUARED dosage table -- fd[c]*fd[c] formed
// on the host in double, which is cell for cell the number the CPU kernel's
// Gb % Gb holds -- and multiplied against MU2bin. A run with no binary trait
// passes K2 = 0 and this second GEMM does not exist.
//
// Binary traits also need, per (marker, trait), the number of cases carrying
// each 2-bit code (AF_case / AF_ctrl, popcount_af.hpp). Given one 2-bit mask
// per trait with 11 in every case sample's field, that is three popcounts per
// 64-bit word of (column AND mask); the reducer does them on the device for
// every slot and every mask and hands back the four counts, so the host never
// touches the column for it. Where the host decides the counts formula does
// not reproduce the sequential sum (a mean-imputed missing call), it replays
// the sum from the same packed column, which it still has in the pinned
// staging buffer.
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
// chunk behind the compute. This reducer therefore always streams, one
// superblock (maxSlots markers) at a time: the superblock's packed codes go up
// in one copy and stay resident until the next reduce() -- the SPA kernel
// (gpu_spa.hpp) reads flagged markers' columns straight from there -- while the
// decode+GEMM runs in slots-per-pass batches on two streams. The "does the
// packed matrix fit in device memory" question does not arise.
#pragma once

#include <cstddef>
#include <cstdint>
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

struct CreateArgs {
    int device = 0;
    int N = 0;                     // samples
    // Right operand of C1 = G^T B1: N x K1 column-major DOUBLE (the caller's
    // stacked constants as they already are); narrowed here when fp64 is
    // false. Copied to the device, so it may be freed on return.
    int K1 = 0; const double* B1 = nullptr;
    // Right operand of C2 = (G % G)^T B2, N x K2. K2 = 0: no second GEMM.
    int K2 = 0; const double* B2 = nullptr;
    int maxSlots = 0;              // markers per reduce() call, at most
    // Run the decode and the GEMMs in double instead of float. On a V100 this
    // costs about 40% more GEMM time and about one extra pass of device
    // bandwidth in the decode -- a few percent of a real step-2 run -- and it
    // removes the precision question rather than bounding it. See PRECISION in
    // gpu_step2.cu.
    bool fp64 = true;
    // Code-count masks: nMask rows of maskWords() 64-bit words each (the
    // layout popcount_af.hpp's pcBuildMask writes), 11 in every selected
    // sample's 2-bit field. nMask = 0: no counts.
    int nMask = 0; const uint64_t* masks = nullptr;
};

// Returns nullptr on ANY failure (no device, allocation refused, ...). The
// caller must fall back to the CPU path; nothing is printed here.
Reducer* create(const CreateArgs& t_args);
void     destroy(Reducer* t_r);

// Words per packed row, = popcount_af.hpp's pcWords(N) = (N + 31) / 32.
int maskWords(int t_N);

// Pinned staging buffers owned by the reducer. The caller writes marker data
// straight into them, so no copy stands between the decode and the H2D.
//   packed()  t_maxSlots * bytesPerSlot() bytes; slot j at offset
//             j*bytesPerSlot(). bytesPerSlot() is (N+3)/4 rounded UP to a
//             multiple of 8 so every row is a whole number of 64-bit words;
//             the caller writes (N+3)/4 bytes and leaves the rest, which the
//             reducer zeroed once and the kernels never see as a sample.
//   lut()     t_maxSlots * 4 doubles; slot j's code->dosage table at 4*j.
//             DOUBLE, not float, in both precision modes: three of its entries
//             are exact small integers but the fourth is the imputed mean
//             2*altFreq, and narrowing that on the host would put a 6e-8
//             relative error into every missing cell before the reduction even
//             starts -- which is a difference from the CPU's INPUT, not from
//             its arithmetic. In fp32 mode the kernel narrows it itself, so
//             the two modes still see the same table.
unsigned char* packed(Reducer* t_r);
double*        lut(Reducer* t_r);
std::size_t    bytesPerSlot(const Reducer* t_r);

// Reduce slots [0, t_nSlots). Returns false on any CUDA error, in which case
// the results are undefined and the caller must fall back for this batch.
bool reduce(Reducer* t_r, int t_nSlots);

// Results of the last reduce(): C1 is t_maxSlots x K1 column-major, so element
// (slot, k) sits at index (std::size_t)k * ldC() + slot; C2 likewise with K2
// columns. Pinned, owned here. Exactly one of each pair of accessors is
// non-null, per isFp64().
bool          isFp64(const Reducer* t_r);
const float*  outCf(const Reducer* t_r);
const double* outCd(const Reducer* t_r);
const float*  outC2f(const Reducer* t_r);
const double* outC2d(const Reducer* t_r);
std::size_t   ldC(const Reducer* t_r);         // == maxSlots
// Code counts of the last reduce(): for slot s and mask m, the four values at
// outCounts()[((std::size_t)s * nMask + m) * 4 + c] are the number of masked
// samples whose 2-bit code is c (c indexes PLINK codes 0..3, so c = 1 is the
// missing call). nullptr when nMask was 0.
const uint32_t* outCounts(const Reducer* t_r);

// Cumulative device-side timings in seconds since create(), for the end-to-end
// breakdown. Measured with events on the compute stream, so they are the GPU's
// own numbers and not wall clock around the call.
void timings(const Reducer* t_r, double* t_h2d, double* t_decode, double* t_gemm,
             double* t_d2h, double* t_popc = nullptr);

// Peak device bytes allocated by this reducer.
std::size_t deviceBytes(const Reducer* t_r);

// For gpu_spa.hpp: the device address of the resident packed rows and tables
// of the last reduce() (row stride bytesPerSlot()), and the stream they were
// written on. Opaque pointers; only the SPA module dereferences them.
const void* devicePacked(const Reducer* t_r);
const void* deviceLut(const Reducer* t_r);
void*       deviceStream(const Reducer* t_r);

}  // namespace gpu2
}  // namespace saige
