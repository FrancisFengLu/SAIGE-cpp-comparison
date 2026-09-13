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
  // Only emit the three codes the real packer can produce: 0b00 HOM_REF,
  // 0b10 HET, 0b11 HOM_ALT. 0b01 (MISSING) is filled in with round(2·altFreq)
  // by marker_decoder BEFORE packing, so it never reaches the device — and the
  // two backends disagree on it by construction (tier 1/3 substitute
  // round(2f), tier 4's g = 2 − popc(field) would call it a het). Feeding
  // random bytes here would therefore test a state that cannot occur.
  const unsigned char kCodes[3] = {0x0, 0x2, 0x3};
  for (int m = 0; m < M; ++m) {
    for (int b = 0; b < nbyte; ++b) {
      unsigned char byte = 0;
      for (int k = 0; k < 4; ++k) byte |= (unsigned char)(kCodes[rng() % 3] << (2*k));
      row[b] = byte;
    }
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

  // Exercise every tier the facade can be forced into. 2 is not selectable
  // (it is what 1 becomes when A happens to fit resident), so 1 covers both.
  const int   tiers[]  = {1, 3, 4};
  const float tol_rel  = 1e-4f;
  int         failures = 0;

  for (int t : tiers) {
    auto* h = saige::gpu::create(packed, freq, invstd, N, t);
    if (!h) {
      std::fprintf(stderr, "gpu::create(tier_override=%d) failed\n", t);
      ++failures;
      continue;
    }
    const int got = saige::gpu::tier(h);
    std::printf("\n-- tier_override=%d -> tier=%d --\n", t, got);

    std::vector<float> Au_gpu(N, 0.0f);
    if (!saige::gpu::matvec(h, u.data(), Au_gpu.data())) {
      std::fprintf(stderr, "gpu::matvec failed (tier %d)\n", got);
      saige::gpu::destroy(h);
      ++failures;
      continue;
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
    int bitdiff = 0;
    for (int i = 0; i < N; ++i)
      if (std::memcmp(&Au_gpu[i], &Au_gpu2[i], sizeof(float)) != 0) ++bitdiff;
    std::printf("GPU inter-run: %d/%d words differ  %s\n",
                bitdiff, N, bitdiff == 0 ? "bit-identical" : "NONDETERMINISTIC");
    if (bitdiff != 0) ++failures;

    saige::gpu::destroy(h);

    if (rel > tol_rel) {
      std::fprintf(stderr, "FAIL tier %d: rel diff %.4g > %.4g\n", got, rel, tol_rel);
      ++failures;
    } else {
      std::printf("OK tier %d (rel diff %.4g < %.4g)\n", got, rel, tol_rel);
    }
  }

  // ---- create_rows(): scattered per-marker pointers, no host gather -------
  // PackedFlat's rows happen to be contiguous, so we can point at them and
  // still exercise the row-by-row upload used for genoVecofPointers.
  {
    std::vector<const unsigned char*> rows(M);
    for (int m = 0; m < M; ++m)
      rows[m] = packed.raw() + static_cast<std::size_t>(m) * nbyte;
    auto* h = saige::gpu::create_rows(rows.data(), nbyte, M, freq, invstd, N, 4);
    if (!h) {
      std::fprintf(stderr, "gpu::create_rows failed\n");
      ++failures;
    } else {
      std::printf("\n-- create_rows -> tier=%d --\n", saige::gpu::tier(h));
      std::vector<float> Au(N, 0.0f);
      if (!saige::gpu::matvec(h, u.data(), Au.data())) {
        std::fprintf(stderr, "matvec after create_rows failed\n");
        ++failures;
      } else {
        float d = 0.0f;
        for (int i = 0; i < N; ++i) d = std::max(d, std::fabs(Au_cpu[i] - Au[i]));
        const float rel = d / std::max(1e-12f, max_abs_cpu);
        std::printf("Δ (CPU vs GPU): max|d|/max|Au_cpu|=%.4g  %s\n",
                    rel, rel <= tol_rel ? "OK" : "FAIL");
        if (rel > tol_rel) ++failures;
      }
      saige::gpu::destroy(h);
    }
  }

  // ---- matvec_mat(): batch K·U vs k separate matvecs ----------------------
  // k=5  → one chunk padded from 5 to NC=8 (zero columns must not leak)
  // k=12 → two chunks, 8 then 4
  for (int k : {5, 12}) {
    auto* h = saige::gpu::create(packed, freq, invstd, N, 4);
    if (!h) { std::fprintf(stderr, "create(tier 4) failed\n"); ++failures; continue; }
    if (!saige::gpu::matvec_mat_available(h)) {
      std::fprintf(stderr, "tier 4 reports no batch kernel\n");
      ++failures; saige::gpu::destroy(h); continue;
    }
    std::printf("\n-- matvec_mat k=%d --\n", k);

    std::vector<float> U(static_cast<std::size_t>(N) * k);
    for (auto& v : U) v = d(rng);

    std::vector<float> KU(static_cast<std::size_t>(N) * k, 0.0f);
    if (!saige::gpu::matvec_mat(h, U.data(), k, KU.data())) {
      std::fprintf(stderr, "matvec_mat failed\n"); ++failures;
      saige::gpu::destroy(h); continue;
    }
    // Reference: the same handle, one column at a time.
    float max_rel = 0.0f;
    std::vector<float> col(N, 0.0f);
    for (int c = 0; c < k; ++c) {
      saige::gpu::matvec(h, U.data() + static_cast<std::size_t>(c) * N, col.data());
      float num = 0.0f, den = 0.0f;
      for (int i = 0; i < N; ++i) {
        const float dv = KU[static_cast<std::size_t>(c) * N + i] - col[i];
        num += dv * dv;
        den += col[i] * col[i];
      }
      max_rel = std::max(max_rel, std::sqrt(num) / std::max(1e-20f, std::sqrt(den)));
    }
    std::printf("max per-column rel-L2 vs single-column kernel = %.4g  %s\n",
                max_rel, max_rel <= 1e-5f ? "OK" : "FAIL");
    if (max_rel > 1e-5f) ++failures;

    // Determinism
    std::vector<float> KU2(static_cast<std::size_t>(N) * k, 0.0f);
    saige::gpu::matvec_mat(h, U.data(), k, KU2.data());
    std::printf("batch inter-run: %s\n",
                std::memcmp(KU.data(), KU2.data(), KU.size() * sizeof(float)) == 0
                    ? "bit-identical" : "NONDETERMINISTIC");
    if (std::memcmp(KU.data(), KU2.data(), KU.size() * sizeof(float)) != 0) ++failures;

    saige::gpu::destroy(h);
  }

  if (failures) {
    std::fprintf(stderr, "\n%d check(s) FAILED\n", failures);
    return 1;
  }
  std::printf("\nAll checks OK\n");
  return 0;
}
