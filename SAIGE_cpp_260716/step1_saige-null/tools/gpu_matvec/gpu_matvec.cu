// gpu_matvec.cu — tier-1 cuBLAS streamed sgemv backend.
// Implements the opaque facade declared in gpu_matvec.hpp.
//
// Design:
//   · At create() we standardize the packed 2-bit genotypes on the HOST into
//     a contiguous float32 A of shape (N × M_pass), column-major. We size
//     blocks to fit in the device and pin the host staging buffer so PCIe
//     transfers hit the fast path.
//   · matvec() does K·u in two cublasSgemv calls *per block*, accumulating
//     partials into a single length-N output on the device. Final result is
//     copied back to host and divided by M_pass.
//   · Everything is float32 end-to-end — matches CPU parallelCrossProd.
//
// Memory at 100 K × 250 K (N=100 k, M=250 k):
//     full A float32   = 100 GB  ← too big for any single GPU
//     one block        = MblkBytes = min(free_vram / 2 − handle_overhead, …)
//     e.g. 20 GB block = 50 k markers × 100 k samples × 4 bytes
//     → 5 blocks per matvec, ~0.8 s PCIe per block → ~4 s total
//
// NOTE: this tier streams blocks FROM HOST per matvec. That's fine for
// correctness and a ~4× gain, but it hits PCIe every matvec. G5 (tier-3)
// keeps the packed form resident on the device.
// Tier 4 (G6) lives in gemv2bit.cu: the packed bytes stay on the device AND
// the standardization is lifted out of the inner product by a rank-one
// decomposition, so the GEMM body only touches raw {0,1,2}. It supersedes
// tier 3 whenever it fits; tier 3 is kept as the per-element-standardizing
// cross-check baseline.
#include "gpu_matvec.hpp"
#include "gemv2bit.hpp"
#include "packed_store.hpp"

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace saige::gpu {

// ---------------------------------------------------------------- helpers
namespace {

#define CUDA_CHECK(expr)                                                       \
  do {                                                                         \
    cudaError_t _s = (expr);                                                   \
    if (_s != cudaSuccess) {                                                   \
      std::fprintf(stderr,                                                     \
                   "[gpu_matvec] CUDA error at %s:%d: %s\n",                   \
                   __FILE__, __LINE__, cudaGetErrorString(_s));                \
      return false;                                                            \
    }                                                                          \
  } while (0)

#define CUBLAS_CHECK(expr)                                                     \
  do {                                                                         \
    cublasStatus_t _s = (expr);                                                \
    if (_s != CUBLAS_STATUS_SUCCESS) {                                         \
      std::fprintf(stderr,                                                     \
                   "[gpu_matvec] cuBLAS error at %s:%d: %d\n",                 \
                   __FILE__, __LINE__, static_cast<int>(_s));                  \
      return false;                                                            \
    }                                                                          \
  } while (0)

// Standardize a range of markers into `A_blk` (column-major, N rows × n_col).
// A_blk[i + j*N] = (bufferGeno(i, m0+j) - 2*freq[m0+j]) * invstd[m0+j], where
// bufferGeno == 2 for hom A1, 1 for het, 0 for hom A2, and fill_missing for
// missing bytes. This matches Get_OneSNP_Geno_atBeginning's post-fill values.
//
// We re-derive the 2-bit BED codes from the packed layout: each byte encodes
// samples [4j, 4j+1, 4j+2, 4j+3] in lo-to-hi order, with code meanings
// identical to BED:
//     lo,hi = 00 → bufferGeno=2 (hom A1)
//           = 01 → missing
//           = 10 → bufferGeno=1 (het)
//           = 11 → bufferGeno=0 (hom A2)
void expand_block_cpu(const saige::PackedFlat& packed,
                      const std::vector<float>& freq,
                      const std::vector<float>& invstd,
                      int N,
                      int m0, int n_col,
                      float* A_blk) {
  const int nbyte = static_cast<int>(packed.nbyte());
  for (int j = 0; j < n_col; ++j) {
    const int m = m0 + j;
    const float f   = freq[m];
    const float s   = invstd[m];
    const int   fill = static_cast<int>(std::round(2.0f * f));  // post-fill value
    const float neg2f = 2.0f * f;
    const unsigned char* row = packed.raw() + static_cast<size_t>(m) * nbyte;
    float* col = A_blk + static_cast<size_t>(j) * N;
    int i = 0;
    for (int b = 0; b < nbyte && i < N; ++b) {
      unsigned char byte = row[b];
      for (int k = 0; k < 4 && i < N; ++k, byte >>= 2) {
        const int code = byte & 3;
        int g;
        switch (code) {
          case 0x0: g = 2; break;   // 0b00 hom A1 (stored HOM_REF in SAIGE's repack)
          case 0x2: g = 1; break;   // 0b10 het  (HET)
          case 0x3: g = 0; break;   // 0b11 hom A2 (HOM_ALT)
          default:  g = fill; break;// 0b01 never lands here after repack, kept for safety
        }
        col[i++] = (static_cast<float>(g) - neg2f) * s;
      }
    }
  }
}

} // namespace

// ---------------------------------------------------------------- Handle
struct Handle {
  // Config
  int    N          = 0;
  int    M          = 0;   // pass-QC markers
  int    nbyte      = 0;   // ⌈N/4⌉
  int    n_blocks   = 0;
  int    Mblk_max   = 0;   // max markers per block we sized for
  int    tier       = 1;   // 1=streamed cuBLAS, 2=resident float A,
                           // 3=packed custom kernel (per-element standardize),
                           // 4=packed + rank-one decomposition (gemv2bit.cu)

  // Device resources (tier 1/2: cuBLAS float A)
  cublasHandle_t cublas  = nullptr;
  float*         d_Ablk  = nullptr;   // Mblk_max × N  float32 (column-major)
  float*         d_u     = nullptr;   // N
  float*         d_y     = nullptr;   // Mblk_max (tier 1/2) or M (tier 3)
  float*         d_Au    = nullptr;   // N

  // Host staging (pinned)
  float*         h_Ablk_pinned = nullptr;

  // Device resources (tier 3: packed 2-bit bytes + per-marker freq/invstd)
  uint8_t*       d_packed  = nullptr;   // M × nbyte 2-bit bytes
  float*         d_freq    = nullptr;   // M
  float*         d_invstd  = nullptr;   // M

  // Source data (kept alive for streaming; we hold a const ref-like pointer)
  const saige::PackedFlat* packed  = nullptr;
  std::vector<float>       freq;     // owned copies — small (M floats)
  std::vector<float>       invstd;

  // Tier-2 shortcut: when the whole A fits in a single block, we expand +
  // upload it once at create() and set this flag so matvec() can skip the
  // per-call host-side expansion.
  bool           A_resident = false;

  // Tier-4 context (gemv2bit.cu). Owns its own device buffers.
  saige::gpu::g2b::Ctx* g2b = nullptr;
};

// -------- G5: custom kernels operating on packed 2-bit bytes --------------
//
// SAIGE's packed encoding (from setGenotype in SAIGE_step1_fast.cpp:196):
//   bit pair value → bufferGeno → standardized (g − 2f)/√(2f(1−f))
//     0b00 (0x0) = HOM_REF  → bufferGeno=2
//     0b10 (0x2) = HET      → bufferGeno=1
//     0b11 (0x3) = HOM_ALT  → bufferGeno=0
//     0b01 (0x1) = MISSING  → fill = round(2·freq); never appears after repack
//
// These kernels bypass the host-side expansion entirely.

__device__ __forceinline__
float decode_and_standardize(uint8_t byte, int bit_pos, float f, float invstd, float fill_g) {
  const int code = (byte >> bit_pos) & 0x3;
  float g;
  switch (code) {
    case 0x0: g = 2.0f; break;
    case 0x2: g = 1.0f; break;
    case 0x3: g = 0.0f; break;
    default:  g = fill_g; break;      // 0x1 = missing (post-fill value)
  }
  return (g - 2.0f * f) * invstd;
}

// Kernel T: y[m] = ⟨A[m,:], u⟩     grid: M blocks, block: 256 threads.
// Each block reduces across N samples for one marker.
__global__ void packed_sgemv_T(const uint8_t* __restrict__ packed,
                                const float*   __restrict__ freq,
                                const float*   __restrict__ invstd,
                                const float*   __restrict__ u,
                                float*                        y,
                                int N, int M, int nbyte) {
  const int m = blockIdx.x;
  if (m >= M) return;
  const uint8_t* row = packed + static_cast<size_t>(m) * nbyte;
  const float f      = freq[m];
  const float s      = invstd[m];
  const float fill_g = roundf(2.0f * f);   // match SAIGE host-side round(2·freq)

  float acc = 0.0f;
  for (int i = threadIdx.x; i < N; i += blockDim.x) {
    const int byte_idx = i >> 2;
    const int bit_pos  = (i & 3) << 1;
    const float x = decode_and_standardize(row[byte_idx], bit_pos, f, s, fill_g);
    acc += x * u[i];
  }

  // Warp-level then block-level reduce.
  for (int off = 16; off > 0; off >>= 1)
    acc += __shfl_down_sync(0xffffffff, acc, off);

  __shared__ float s_warps[32];
  const int lane = threadIdx.x & 31;
  const int warp = threadIdx.x >> 5;
  if (lane == 0) s_warps[warp] = acc;
  __syncthreads();

  if (warp == 0) {
    acc = (threadIdx.x < (blockDim.x + 31) / 32) ? s_warps[lane] : 0.0f;
    for (int off = 16; off > 0; off >>= 1)
      acc += __shfl_down_sync(0xffffffff, acc, off);
    if (threadIdx.x == 0) y[m] = acc;
  }
}

// Kernel N: Au[i] = Σ_m A[m,i] · y[m] / M     grid: ⌈N/256⌉, block: 256.
// Each thread handles one sample, loops over all markers.
// Memory access is strided by nbyte per marker — inefficient per-thread, but
// coalesced within a warp (32 consecutive samples hit 32 consecutive offsets
// within the same marker row, which lands 32 bytes = 2 128-byte cache lines).
__global__ void packed_sgemv_N(const uint8_t* __restrict__ packed,
                                const float*   __restrict__ freq,
                                const float*   __restrict__ invstd,
                                const float*   __restrict__ y,
                                float*                        Au,
                                int N, int M, int nbyte, float inv_M) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= N) return;
  const int byte_idx = i >> 2;
  const int bit_pos  = (i & 3) << 1;

  float acc = 0.0f;
  for (int m = 0; m < M; ++m) {
    const float f      = freq[m];
    const float s      = invstd[m];
    const float fill_g = roundf(2.0f * f);   // match SAIGE host-side round(2·freq)
    const uint8_t byte = packed[static_cast<size_t>(m) * nbyte + byte_idx];
    acc += decode_and_standardize(byte, bit_pos, f, s, fill_g) * y[m];
  }
  Au[i] = acc * inv_M;
}

static bool init_handle(Handle* h,
                        const saige::PackedFlat& packed,
                        const std::vector<float>& freq,
                        const std::vector<float>& invstd,
                        int N, int tier_override) {
  h->N       = N;
  h->M       = static_cast<int>(packed.n_stored());
  h->nbyte   = static_cast<int>(packed.nbyte());
  h->packed  = &packed;
  h->freq    = freq;
  h->invstd  = invstd;

  // Bind the device BEFORE the first cudaMalloc/cudaMemGetInfo, otherwise the
  // runtime picks device 0 implicitly and SAIGE_GPU_DEVICE has no effect.
  // (Ported from bind_device_for_rank() in the R package's gpuSymMatMult.cu,
  // minus the MPI rank: this binary is single-process.)
  {
    int dev = 0;
    if (const char* dv = std::getenv("SAIGE_GPU_DEVICE")) dev = std::atoi(dv);
    int count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&count));
    if (count <= 0) return false;
    if (dev < 0 || dev >= count) dev = 0;
    CUDA_CHECK(cudaSetDevice(dev));
  }

  size_t free_b = 0, total_b = 0;
  CUDA_CHECK(cudaMemGetInfo(&free_b, &total_b));

  // ---- decide tier ----
  //   tier 3 budget: packed (M*nbyte) + freq + invstd + u + Au + y + overhead
  const size_t per_col_float = static_cast<size_t>(h->N) * sizeof(float);
  const size_t packed_bytes  = static_cast<size_t>(h->M) * h->nbyte;
  const size_t aux_bytes     = 2 * static_cast<size_t>(h->M) * sizeof(float)  // freq,invstd
                              + 2 * per_col_float                              // d_u, d_Au
                              + static_cast<size_t>(h->M) * sizeof(float);     // d_y
  const size_t tier3_need = packed_bytes + aux_bytes + 256 * 1024 * 1024;
  // Tier 4 carries the same packed bytes plus the two partial-reduction
  // buffers (Ypart ≈ ⌈W32/128⌉·M·4, Zpart ≈ ⌈M/256⌉·Npad·4); at UKB scale
  // Zpart is ~690 MB, so it is NOT a rounding term — ask gemv2bit for the
  // exact figure rather than re-deriving it here.
  const size_t tier4_need =
      saige::gpu::g2b::need_bytes(h->N, h->M) + 256 * 1024 * 1024;

  int chosen_tier = 1;
  if (tier_override == 4 || (tier_override == 0 && tier4_need <= free_b * 90 / 100)) {
    chosen_tier = 4;
  } else if (tier_override == 3 || (tier_override == 0 && tier3_need <= free_b * 90 / 100)) {
    chosen_tier = 3;
  }

  if (chosen_tier == 4) {
    // PackedFlat is one contiguous buffer with row stride == nbyte(), so the
    // upload is a single cudaMemcpy2D straight off raw() — no gather into a
    // staging buffer, no host-side fp32 expansion. (The R package had to
    // gather because its store was a vector of per-block pointers; that copy
    // is what made host RSS 23.2/29 GiB at UKB scale.)
    h->g2b = saige::gpu::g2b::create(packed.raw(), packed.nbyte(),
                                     h->N, h->M, freq.data(), invstd.data());
    if (!h->g2b) {
      std::fprintf(stderr, "[gpu_matvec] tier-4 create failed"
                           " (need %zu MB, free %zu MB)\n",
                   tier4_need >> 20, free_b >> 20);
      return false;
    }
    h->tier       = 4;
    h->n_blocks   = 1;
    h->Mblk_max   = h->M;
    h->A_resident = false;

    std::printf("[gpu_matvec] tier=4 N=%d  M=%d  nbyte=%d  device=%zu MB "
                "(packed %zu MB + scratch)  VRAM free=%zu MB  "
                "(2-BIT RESIDENT, rank-one standardization)\n",
                h->N, h->M, h->nbyte,
                saige::gpu::g2b::ctx_bytes(h->g2b) >> 20,
                packed_bytes >> 20, free_b >> 20);
    return true;
  }

  CUBLAS_CHECK(cublasCreate(&h->cublas));

  if (chosen_tier == 3) {
    CUDA_CHECK(cudaMalloc(&h->d_packed, packed_bytes));
    CUDA_CHECK(cudaMalloc(&h->d_freq,   h->M * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&h->d_invstd, h->M * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&h->d_u,      per_col_float));
    CUDA_CHECK(cudaMalloc(&h->d_y,      h->M * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&h->d_Au,     per_col_float));

    CUDA_CHECK(cudaMemcpy(h->d_packed, packed.raw(), packed_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(h->d_freq,   freq.data(),  h->M * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(h->d_invstd, invstd.data(),h->M * sizeof(float), cudaMemcpyHostToDevice));

    h->tier       = 3;
    h->n_blocks   = 1;
    h->Mblk_max   = h->M;
    h->A_resident = false;

    std::printf("[gpu_matvec] tier=3 N=%d  M=%d  nbyte=%d  packed=%zu MB "
                "VRAM free=%zu MB  (PACKED-RESIDENT, custom kernel)\n",
                h->N, h->M, h->nbyte, packed_bytes >> 20, free_b >> 20);
    return true;
  }

  // -------- tier 1/2 fallback --------
  const size_t overhead = size_t(4) * h->N * sizeof(float) + size_t(256) * 1024 * 1024;
  const size_t budget = (free_b > overhead) ? (free_b - overhead) * 80 / 100 : size_t(0);
  const size_t per_col = per_col_float;
  size_t Mblk = (per_col > 0) ? (budget / per_col) : 0;
  if (Mblk == 0) {
    std::fprintf(stderr,
                 "[gpu_matvec] VRAM too small: free=%zu MB, need ≥ %zu MB\n",
                 free_b >> 20, (overhead + per_col) >> 20);
    return false;
  }
  Mblk = std::min<size_t>(Mblk, h->M);
  h->Mblk_max = static_cast<int>(Mblk);
  h->n_blocks = (h->M + h->Mblk_max - 1) / h->Mblk_max;

  CUDA_CHECK(cudaMalloc(&h->d_Ablk, per_col * h->Mblk_max));
  CUDA_CHECK(cudaMalloc(&h->d_u,    per_col));
  CUDA_CHECK(cudaMalloc(&h->d_y,    h->Mblk_max * sizeof(float)));
  CUDA_CHECK(cudaMalloc(&h->d_Au,   per_col));
  CUDA_CHECK(cudaMallocHost(&h->h_Ablk_pinned, per_col * h->Mblk_max));

  if (h->n_blocks == 1) {
    expand_block_cpu(*h->packed, h->freq, h->invstd, h->N, 0, h->M, h->h_Ablk_pinned);
    CUDA_CHECK(cudaMemcpy(h->d_Ablk, h->h_Ablk_pinned, per_col * h->M, cudaMemcpyHostToDevice));
    h->A_resident = true;
  }
  h->tier = (h->A_resident ? 2 : 1);

  std::printf("[gpu_matvec] tier=%d N=%d  M=%d  Mblk_max=%d  n_blocks=%d  "
              "VRAM free=%zu MB  A_resident=%s\n",
              h->tier, h->N, h->M, h->Mblk_max, h->n_blocks, free_b >> 20,
              h->A_resident ? "yes" : "no");
  return true;
}

static void free_handle(Handle* h) {
  if (h->h_Ablk_pinned) cudaFreeHost(h->h_Ablk_pinned);
  if (h->d_Ablk)        cudaFree(h->d_Ablk);
  if (h->d_u)           cudaFree(h->d_u);
  if (h->d_y)           cudaFree(h->d_y);
  if (h->d_Au)          cudaFree(h->d_Au);
  if (h->d_packed)      cudaFree(h->d_packed);
  if (h->d_freq)        cudaFree(h->d_freq);
  if (h->d_invstd)      cudaFree(h->d_invstd);
  if (h->cublas)        cublasDestroy(h->cublas);
  if (h->g2b)           saige::gpu::g2b::destroy(h->g2b);
  // No memset(h, 0, sizeof(*h)) here: Handle holds two std::vector<float>
  // members, and zeroing them out from under the destructor leaked their
  // buffers (2 × M × 4 bytes per handle). Null the raw pointers explicitly so
  // a double destroy() is still safe.
  h->h_Ablk_pinned = nullptr;
  h->d_Ablk = h->d_u = h->d_y = h->d_Au = nullptr;
  h->d_freq = h->d_invstd = nullptr;
  h->d_packed = nullptr;
  h->cublas   = nullptr;
  h->g2b      = nullptr;
  h->packed   = nullptr;
  h->tier     = 0;
}

// ---------------------------------------------------------------- API
bool available() {
  int count = 0;
  if (cudaGetDeviceCount(&count) != cudaSuccess) return false;
  return count > 0;
}

Handle* create(const saige::PackedFlat& packed,
               const std::vector<float>& freq,
               const std::vector<float>& invstd,
               int N,
               int tier_override) {
  if (!available()) return nullptr;
  if (packed.n_stored() == 0 || N <= 0) return nullptr;
  if (static_cast<int>(freq.size())   < static_cast<int>(packed.n_stored())) return nullptr;
  if (static_cast<int>(invstd.size()) < static_cast<int>(packed.n_stored())) return nullptr;

  auto* h = new Handle();
  if (!init_handle(h, packed, freq, invstd, N, tier_override)) {
    free_handle(h);
    delete h;
    return nullptr;
  }
  return h;
}

void destroy(Handle* h) {
  if (!h) return;
  free_handle(h);
  delete h;
}

int tier(const Handle* h) { return h ? h->tier : 0; }

bool matvec(Handle* h, const float* u, float* out_Au) {
  if (!h || !u || !out_Au) return false;

  // ---- tier 4: 2-bit resident + rank-one standardization ----
  // Whole-GRM range [0, M); gemv2bit uploads u and applies 1/M itself.
  if (h->tier == 4) {
    return saige::gpu::g2b::matvec_range(h->g2b, 0, h->M,
                                         1.0f / static_cast<float>(h->M),
                                         u, out_Au);
  }

  // Push u once to the device.
  const size_t per_col = static_cast<size_t>(h->N) * sizeof(float);
  CUDA_CHECK(cudaMemcpy(h->d_u, u, per_col, cudaMemcpyHostToDevice));

  // ---- tier 3: custom packed kernels ----
  if (h->tier == 3) {
    const int tpb = 256;
    packed_sgemv_T<<<h->M, tpb>>>(h->d_packed, h->d_freq, h->d_invstd,
                                   h->d_u, h->d_y,
                                   h->N, h->M, h->nbyte);
    CUDA_CHECK(cudaGetLastError());
    const int blocks = (h->N + tpb - 1) / tpb;
    const float inv_M = 1.0f / static_cast<float>(h->M);
    packed_sgemv_N<<<blocks, tpb>>>(h->d_packed, h->d_freq, h->d_invstd,
                                     h->d_y, h->d_Au,
                                     h->N, h->M, h->nbyte, inv_M);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpy(out_Au, h->d_Au, per_col, cudaMemcpyDeviceToHost));
    return true;
  }

  // ---- tier 1/2: cuBLAS sgemv (streamed or A-resident) ----
  CUDA_CHECK(cudaMemsetAsync(h->d_Au, 0, per_col));

  const float alpha = 1.0f;
  const float beta_y = 0.0f;  // y ← αAᵀu + βy
  const float beta_z = 1.0f;  // Au ← αAy + βAu   (accumulate)

  for (int b = 0; b < h->n_blocks; ++b) {
    const int m0  = b * h->Mblk_max;
    const int n_c = std::min(h->Mblk_max, h->M - m0);

    // Tier-2 fast path: A is already on the device from create(); skip the
    // re-expand + re-upload. Tier-1 streaming path (A doesn't fit) still does
    // both per matvec.
    if (!h->A_resident) {
      expand_block_cpu(*h->packed, h->freq, h->invstd, h->N, m0, n_c,
                       h->h_Ablk_pinned);
      CUDA_CHECK(cudaMemcpy(h->d_Ablk, h->h_Ablk_pinned,
                            static_cast<size_t>(h->N) * n_c * sizeof(float),
                            cudaMemcpyHostToDevice));
    }

    // y ← Aᵀ · u                      (N × n_c)' · (N) → (n_c)
    CUBLAS_CHECK(cublasSgemv(
        h->cublas, CUBLAS_OP_T,
        /*m=*/h->N, /*n=*/n_c,
        &alpha,
        h->d_Ablk, /*lda=*/h->N,
        h->d_u, /*incx=*/1,
        &beta_y,
        h->d_y, /*incy=*/1));
    // Au ← A · y + Au                  (N × n_c) · (n_c) → (N)
    CUBLAS_CHECK(cublasSgemv(
        h->cublas, CUBLAS_OP_N,
        /*m=*/h->N, /*n=*/n_c,
        &alpha,
        h->d_Ablk, /*lda=*/h->N,
        h->d_y, /*incx=*/1,
        &beta_z,
        h->d_Au, /*incy=*/1));
  }

  // Divide by M — do it on the device, then copy back.
  const float invM = 1.0f / static_cast<float>(h->M);
  CUBLAS_CHECK(cublasSscal(h->cublas, h->N, &invM, h->d_Au, 1));
  CUDA_CHECK(cudaMemcpy(out_Au, h->d_Au, per_col, cudaMemcpyDeviceToHost));
  return true;
}

} // namespace saige::gpu
