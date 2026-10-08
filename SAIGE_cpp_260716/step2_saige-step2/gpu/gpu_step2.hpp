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

#include "gpu_precision.hpp"

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
    // Arithmetic of the decode and the GEMMs (config key gpuPrecisionScan,
    // gpu_precision.hpp). FP64 (default): double throughout. On a V100 this
    // costs about 40% more GEMM time than FP32 and about one extra pass of
    // device bandwidth in the decode -- a few percent of a real step-2 run --
    // and it removes the precision question rather than bounding it. FP32:
    // float decode + SGEMM, float results (outCf / outC2f). INT8: the
    // Ozaki-style split (int8 genotypes, B1 / B2 split into int8Slices int8
    // slices with per-column scaling, int32 GEMMs, recombined in fp64); double
    // results (outCd / outC2d). Only the modes scanSupports() accepts may be
    // passed; create() returns nullptr for any other. See PRECISION in
    // gpu_step2.cu.
    Prec precision = Prec::FP64;
    // INT8 only: slices of the trait-side operand, kInt8SlicesMin..Max.
    int int8Slices = kInt8SlicesDefault;
    // Code-count masks: nMask rows of maskWords() 64-bit words each (the
    // layout popcount_af.hpp's pcBuildMask writes), 11 in every selected
    // sample's 2-bit field. nMask = 0: no counts.
    int nMask = 0; const uint64_t* masks = nullptr;
    // Pinned staging sets (packed + lut), each maxSlots wide. 1 is the
    // original single set. More than one lets the caller fill set s+1 on the
    // host while reduce(set s) and the host work after it run (the read /
    // compute overlap of config key gpuPrefetch, S2_PIPELINE.md). Nothing on
    // the device changes: every set uploads into the same device buffers.
    int stagingSets = 1;
    // fp64 decode with one 16-byte store per thread, lanes contiguous, instead
    // of two 16-byte stores per thread 32 bytes apart (default false; only
    // takes effect with fp64 and N % 4 == 0). Same output; on a V100 it
    // removes the read-modify-write that half-written sectors cost on ECC'd
    // HBM2 (S2_KERNEL_ROOFLINE.md).
    bool decodeX2 = false;
    // Device sets: resident packed rows + tables on the device, and pinned
    // host result buffers (C1, C2, counts), each maxSlots wide. 1 is the
    // original single set. More than one lets the caller keep superblock k's
    // device rows (for the device SPA / Firth) and host results (for the host
    // tail) alive while reduce() of a later superblock runs into another set
    // (config key gpuOverlap, S2_OVERLAP.md). The GEMM operands, the decode
    // buffers, the launch shapes and the cuBLAS calls are the same whichever
    // set is used, so the numbers are too.
    int deviceSets = 1;
};

// True when create() accepts CreateArgs::precision = t_p. main.cpp asks
// before create() and stops the run on false ("useGPU: scan precision <mode>
// is not implemented yet"); it never falls back to another precision.
bool scanSupports(Prec t_p);

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
//   t_set     which staging set (0 .. stagingSets()-1); the single set of a
//             reducer created with stagingSets = 1 is set 0.
unsigned char* packed(Reducer* t_r, int t_set = 0);
double*        lut(Reducer* t_r, int t_set = 0);
std::size_t    bytesPerSlot(const Reducer* t_r);
int            stagingSets(const Reducer* t_r);

// Reduce slots [0, t_nSlots) of staging set t_set. Returns false on any CUDA
// error, in which case the results are undefined and the caller must fall
// back for this batch. Synchronous: on return the host result buffers are
// filled and staging set t_set is no longer read by the device.
// t_devSet: the device set (0 .. deviceSets-1) the rows are uploaded into and
// the results come back into; -1 = 0. Only one thread may call reduce() at a
// time; the accessors below may be read from another thread for a set no
// reduce() in flight is writing.
bool reduce(Reducer* t_r, int t_nSlots, int t_set = 0, int t_devSet = -1);

// Results of the last reduce(): C1 is t_maxSlots x K1 column-major, so element
// (slot, k) sits at index (std::size_t)k * ldC() + slot; C2 likewise with K2
// columns. Pinned, owned here. Exactly one of each pair of accessors is
// non-null, per isFp64(), which says whether the RESULT buffers are double:
// true for precision FP64 and INT8 (recombined in fp64), false for FP32.
bool          isFp64(const Reducer* t_r);
const float*  outCf(const Reducer* t_r, int t_devSet = -1);
const double* outCd(const Reducer* t_r, int t_devSet = -1);
const float*  outC2f(const Reducer* t_r, int t_devSet = -1);
const double* outC2d(const Reducer* t_r, int t_devSet = -1);
std::size_t   ldC(const Reducer* t_r);         // == maxSlots
// Code counts of the last reduce(): for slot s and mask m, the four values at
// outCounts()[((std::size_t)s * nMask + m) * 4 + c] are the number of masked
// samples whose 2-bit code is c (c indexes PLINK codes 0..3, so c = 1 is the
// missing call). nullptr when nMask was 0.
const uint32_t* outCounts(const Reducer* t_r, int t_devSet = -1);

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
// t_devSet -1 = the set of the last reduce().
const void* devicePacked(const Reducer* t_r, int t_devSet = -1);
const void* deviceLut(const Reducer* t_r, int t_devSet = -1);
void*       deviceStream(const Reducer* t_r);
int         deviceSets(const Reducer* t_r);

// cudaSetDevice for a host thread other than the one that called create()
// (the gpuOverlap workers). false without CUDA or on error.
bool bindDevice(int t_device);
// Ask the driver to block, not spin, a host thread waiting on the device
// (cudaDeviceScheduleBlockingSync), so the gpuOverlap workers' waits leave
// their cores to the host tail. Returns the CUDA status as text ("" = ok).
std::string setBlockingSync(int t_device);

}  // namespace gpu2
}  // namespace saige
