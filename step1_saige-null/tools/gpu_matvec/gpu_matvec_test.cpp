// gpu_matvec_test.cpp — synthetic correctness harness for the GPU facade.
//
// G1 scope: just verify the link works. With HAS_NVCC=no the stub returns
// false/nullptr for every call; we report that and exit 0. G2 swaps in the
// real tier-1 implementation and this test grows an Au_cpu ≡ Au_gpu diff.
// gpu_matvec_test.cpp — G2 correctness harness.
// Builds a synthetic packed BED of N samples × M markers, computes the same
// K·u both ways (CPU rank-1 accumulation + GPU cuBLAS streamed sgemv), and
// requires ‖Au_cpu − Au_gpu‖_∞ < 1e-3 on float32 (rel ≤ 1e-4 vs max|Au|).
#include "gpu_matvec.hpp"
#include "packed_store.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

namespace {

// Expand one marker's packed bytes into a length-N float32 standardized
// column (same math as gpu_matvec.cu::expand_block_cpu but here as a
// standalone helper so the test doesn't depend on the cu file's internal
// symbols).
void expand_col(const saige::PackedFlat& packed,
                int m, int N, float f, float s,
                std::vector<float>& col) {
  const int nbyte = static_cast<int>(packed.nbyte());
  const int fill  = static_cast<int>(std::round(2.0f * f));
  const float neg2f = 2.0f * f;
  col.assign(N, 0.0f);
  const unsigned char* row = packed.raw() + static_cast<std::size_t>(m) * nbyte;
  int i = 0;
  for (int b = 0; b < nbyte && i < N; ++b) {
    unsigned char byte = row[b];
    for (int k = 0; k < 4 && i < N; ++k, byte >>= 2) {
      const int code = byte & 3;
      int g;
      switch (code) {
        case 0x0: g = 2; break;
        case 0x2: g = 1; break;
        case 0x3: g = 0; break;
        default:  g = fill; break;
      }
      col[i++] = (static_cast<float>(g) - neg2f) * s;
    }
  }
}

// Reference K·u on the CPU: Au = (sum over m of (A_m · u) * A_m) / M.
void matvec_cpu(const saige::PackedFlat& packed,
                const std::vector<float>& freq,
                const std::vector<float>& invstd,
                int N, int M,
                const std::vector<float>& u,
                std::vector<float>& Au) {
  Au.assign(N, 0.0f);
  std::vector<float> col;
  for (int m = 0; m < M; ++m) {
    expand_col(packed, m, N, freq[m], invstd[m], col);
    float val = 0.0f;
    for (int i = 0; i < N; ++i) val += col[i] * u[i];
    for (int i = 0; i < N; ++i) Au[i] += val * col[i];
  }
  const float invM = 1.0f / static_cast<float>(M);
  for (float& v : Au) v *= invM;
}

} // namespace

int main() {
  const bool have_gpu = saige::gpu::available();
  std::printf("saige::gpu::available() = %s\n",
              have_gpu ? "true" : "false (CPU fallback)");

  // Bigger than G1 to actually exercise the streaming logic but small
  // enough to finish in seconds.
  const int N = 4096;
  const int M = 2048;
  const int nbyte = (N + 3) / 4;

  saige::PackedFlat packed;
  packed.init(M, nbyte);
  std::vector<float> freq(M), invstd(M);

  std::mt19937 rng(202604);
  std::vector<unsigned char> row(nbyte);
  for (int m = 0; m < M; ++m) {
    for (int b = 0; b < nbyte; ++b) row[b] = static_cast<unsigned char>(rng());
    packed.write(m, row.data());
    // Alt allele freq between 0.03 and 0.4.
    freq[m]   = 0.03f + (m % 100) * 0.0037f;
    invstd[m] = 1.0f / std::sqrt(2.0f * freq[m] * (1.0f - freq[m]));
  }
  packed.set_n_stored(M);

  std::vector<float> u(N);
  std::uniform_real_distribution<float> d(-1, 1);
  for (int i = 0; i < N; ++i) u[i] = d(rng);

  // CPU reference
  std::vector<float> Au_cpu;
  matvec_cpu(packed, freq, invstd, N, M, u, Au_cpu);
  float max_abs_cpu = 0.0f;
  for (float v : Au_cpu) max_abs_cpu = std::max(max_abs_cpu, std::fabs(v));
  std::printf("CPU reference: max|Au|=%.6g  Au_cpu[0:5]= %.6g %.6g %.6g %.6g %.6g\n",
              max_abs_cpu, Au_cpu[0], Au_cpu[1], Au_cpu[2], Au_cpu[3], Au_cpu[4]);

  if (!have_gpu) {
    std::printf("No GPU — scaffold test exits OK\n");
    return 0;
  }

  auto* h = saige::gpu::create(packed, freq, invstd, N);
  if (!h) { std::fprintf(stderr, "gpu::create failed\n"); return 1; }
  std::printf("GPU handle created (tier=%d)\n", saige::gpu::tier(h));

  std::vector<float> Au_gpu(N, 0.0f);
  if (!saige::gpu::matvec(h, u.data(), Au_gpu.data())) {
    std::fprintf(stderr, "gpu::matvec failed\n");
    saige::gpu::destroy(h);
    return 1;
  }
  std::printf("GPU result    : Au_gpu[0:5]= %.6g %.6g %.6g %.6g %.6g\n",
              Au_gpu[0], Au_gpu[1], Au_gpu[2], Au_gpu[3], Au_gpu[4]);

  // Diff
  float max_abs_diff = 0.0f;
  for (int i = 0; i < N; ++i)
    max_abs_diff = std::max(max_abs_diff, std::fabs(Au_cpu[i] - Au_gpu[i]));
  const float rel = max_abs_diff / std::max(1e-12f, max_abs_cpu);
  std::printf("Δ (CPU vs GPU): max|d|=%.4g  max|d|/max|Au_cpu|=%.4g\n",
              max_abs_diff, rel);

  // Inter-run determinism check (GPU vs GPU on two matvec calls with same u)
  std::vector<float> Au_gpu2(N, 0.0f);
  saige::gpu::matvec(h, u.data(), Au_gpu2.data());
  float inter = 0.0f;
  for (int i = 0; i < N; ++i) inter = std::max(inter, std::fabs(Au_gpu[i] - Au_gpu2[i]));
  std::printf("GPU inter-run: max|Δ|=%.4g  %s\n",
              inter, (inter < 1e-6f) ? "deterministic" : "nondet");

  saige::gpu::destroy(h);

  const float tol_rel = 1e-4f;
  if (rel > tol_rel) {
    std::fprintf(stderr, "FAIL: rel diff %.4g > %.4g\n", rel, tol_rel);
    return 1;
  }
  std::printf("OK (rel diff %.4g < %.4g)\n", rel, tol_rel);
  return 0;
}
