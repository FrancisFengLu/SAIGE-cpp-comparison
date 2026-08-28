// Step-1 microbenchmark for the AVX2 fused 2-bit decode kernels.
// Standalone: g++ -O3 -march=native -std=c++17 avx2_bench.cpp -o /tmp/avx2_bench
//
// Validates against the exact scalar decode used by genoClass
// (Get_OneSNP_StdGeno: g = 2-(a+b), std value = (g-2f)*invStd) and against
// an fp64 reference for the dot/axpy numerics. Times scalar vs AVX2 for
// pass1 (sum g*x) and pass2 (b += s*g), reporting cycles/genotype.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>
#include <chrono>
#include <x86intrin.h>

#include "avx2_kernel.hpp"

using std::size_t;
using std::uint8_t;

static const int SCALAR_LUT_G[4] = {2, 1, 1, 0};  // code -> g, same as 2-(a+b)

// --- scalar reference paths (replicating SAIGE_step1_fast.cpp decode) -------

// Get_OneSNP_StdGeno-equivalent: byte-serial decode into std values, then
// fp32 dot + axpy exactly like CorssProd::operator() does.
static void scalar_marker(const uint8_t* row, size_t N, float freq,
                          float invStd, const float* x, float* bout,
                          float* val1_out) {
  const float lut[3] = {(0 - 2 * freq) * invStd, (1 - 2 * freq) * invStd,
                        (2 - 2 * freq) * invStd};
  std::vector<float> vec(N);
  size_t ind = 0;
  const size_t nbyte = (N + 3) / 4;
  for (size_t i = 0; i < nbyte && ind < N; ++i) {
    unsigned geno1 = row[i];
    for (int j = 0; j < 4 && ind < N; ++j) {
      int b = geno1 & 1; geno1 >>= 1;
      int a = geno1 & 1; geno1 >>= 1;
      vec[ind++] = lut[2 - (a + b)];
    }
  }
  float val1 = 0.f;
  for (size_t i = 0; i < N; ++i) val1 += vec[i] * x[i];
  for (size_t i = 0; i < N; ++i) bout[i] += val1 * vec[i];
  *val1_out = val1;
}

// fp64 reference for numerics.
static void ref64_marker(const uint8_t* row, size_t N, double freq,
                         double invStd, const float* x, double* bout,
                         double* val1_out) {
  double val1 = 0.0;
  for (size_t i = 0; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    const double g = (SCALAR_LUT_G[code] - 2 * freq) * invStd;
    val1 += g * (double)x[i];
  }
  for (size_t i = 0; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    const double g = (SCALAR_LUT_G[code] - 2 * freq) * invStd;
    bout[i] += val1 * g;
  }
  *val1_out = val1;
}

int main(int argc, char** argv) {
  size_t N = 50000, M = 2000;
  if (argc > 1) N = strtoull(argv[1], nullptr, 10);
  if (argc > 2) M = strtoull(argv[2], nullptr, 10);
  const size_t nbyte = (N + 3) / 4;

  printf("N=%zu  M=%zu  nbyte=%zu  (N mod 128 = %zu)\n", N, M, nbyte, N % 128);
  if (!saige_avx2::available()) { printf("AVX2 not available!\n"); return 1; }

  // --- build data ----------------------------------------------------------
  std::mt19937 gen(20260828);
  std::uniform_real_distribution<float> uf(0.02f, 0.5f);
  std::normal_distribution<float> nd(0.f, 1.f);

  std::vector<uint8_t> rows(M * nbyte);
  std::vector<float> freq(M), invstd(M);
  // per-marker genotype distribution driven by its freq; codes only from
  // {0b00(g2), 0b10(g1), 0b11(g0)} like real packed storage, but also inject
  // some 0b01 codes in a few markers to prove LUT covers all 4 codes.
  for (size_t m = 0; m < M; ++m) {
    float f = uf(gen);
    freq[m] = f;
    invstd[m] = 1.0f / std::sqrt(2.f * f * (1.f - f));
    std::discrete_distribution<int> gd{(1 - f) * (1 - f), 2 * f * (1 - f),
                                       f * f};
    static const unsigned code_of_g[3] = {0b11, 0b10, 0b00};  // g=0,1,2
    for (size_t i = 0; i < N; ++i) {
      unsigned g = gd(gen);
      unsigned code = code_of_g[g];
      if (m % 97 == 0 && i % 37 == 0) code = 0b01;  // exercise LUT[1]
      rows[m * nbyte + (i >> 2)] |= code << (2 * (i & 3));
    }
    // dirty the padding bits of the last byte (must be ignored by kernels)
    if (N % 4) rows[m * nbyte + nbyte - 1] |= 0xC0u;
  }
  std::vector<float> x(N);
  for (auto& v : x) v = nd(gen);

  // --- 1) decode bit-exactness --------------------------------------------
  {
    std::vector<uint8_t> g_scalar(N), g_perm(N), g_perm_ref(N);
    size_t bad = 0;
    for (size_t m = 0; m < M; ++m) {
      const uint8_t* row = rows.data() + m * nbyte;
      for (size_t i = 0; i < N; ++i)
        g_scalar[i] =
            (uint8_t)SCALAR_LUT_G[(row[i >> 2] >> (2 * (i & 3))) & 3u];
      saige_avx2::decode_row_perm(row, g_perm.data(), N);
      // build expected perm order from scalar decode
      const size_t S = N / 128;
      for (size_t s = 0; s < S; ++s)
        for (size_t k = 0; k < 4; ++k)
          for (size_t i = 0; i < 32; ++i)
            g_perm_ref[128 * s + 32 * k + i] = g_scalar[128 * s + 4 * i + k];
      for (size_t i = 128 * S; i < N; ++i) g_perm_ref[i] = g_scalar[i];
      bad += (size_t)(memcmp(g_perm.data(), g_perm_ref.data(), N) != 0);
    }
    printf("[decode] bit-exact vs scalar over %zu markers: %s (%zu bad)\n", M,
           bad ? "FAIL" : "PASS", bad);
    if (bad) return 1;
  }

  // --- 2) numerics: val1 + bout vs fp64 reference --------------------------
  // Per-marker relative error is meaningless when val1 ~ 0 (the scalar fp32
  // path itself shows rel errors up to ~6e-3 there). Metric: absolute error
  // normalized by RMS(val1) over markers — for both the existing scalar fp32
  // path and the AVX2 rank-one path. PASS line: AVX2 error <= max(2e-5,
  // 4x scalar), i.e. no worse than the fp32 arithmetic SAIGE already does.
  {
    std::vector<float> xp(N);
    saige_avx2::permute_fwd(x.data(), xp.data(), N);
    double Sx64 = 0;
    for (size_t i = 0; i < N; ++i) Sx64 += x[i];
    const float Sx = (float)Sx64;  // integration computes Sx in fp64 once/psi-v

    std::vector<float>  bout_sc(N, 0.f), bp(N, 0.f), bout_avx(N);
    std::vector<double> bout64(N, 0.0);
    float coffset = 0.f;

    double max_abs_avx = 0, max_abs_sc = 0, rms = 0;
    for (size_t m = 0; m < M; ++m) {
      const uint8_t* row = rows.data() + m * nbyte;
      float v_sc; double v64;
      scalar_marker(row, N, freq[m], invstd[m], x.data(), bout_sc.data(), &v_sc);
      ref64_marker(row, N, freq[m], invstd[m], x.data(), bout64.data(), &v64);
      const float raw = saige_avx2::pass1_sum_gx(row, xp.data(), N);
      const float v_avx = invstd[m] * (raw - 2.f * freq[m] * Sx);
      const float sval = v_avx * invstd[m];
      saige_avx2::pass2_axpy(row, sval, bp.data(), N);
      coffset += 2.f * freq[m] * sval;

      rms += v64 * v64;
      max_abs_avx = std::max(max_abs_avx, std::fabs((double)v_avx - v64));
      max_abs_sc  = std::max(max_abs_sc, std::fabs((double)v_sc - v64));
    }
    rms = std::sqrt(rms / M);
    saige_avx2::unpermute_sub(bp.data(), bout_avx.data(), N, coffset);

    double num = 0, den = 0, num_sc = 0;
    for (size_t i = 0; i < N; ++i) {
      num += (bout_avx[i] - bout64[i]) * (bout_avx[i] - bout64[i]);
      num_sc += (bout_sc[i] - bout64[i]) * (bout_sc[i] - bout64[i]);
      den += bout64[i] * bout64[i];
    }
    const double e_avx = max_abs_avx / rms, e_sc = max_abs_sc / rms;
    const bool val1_pass = e_avx <= std::max(2e-5, 4.0 * e_sc);
    printf("[val1 ] max|err|/RMS(val1) vs fp64:  avx2 = %.3e   scalar fp32 = "
           "%.3e -> %s\n", e_avx, e_sc, val1_pass ? "PASS" : "FAIL");
    printf("[bout ] rel-L2 vs fp64:  avx2 = %.3e   scalar fp32 = %.3e -> %s\n",
           std::sqrt(num / den), std::sqrt(num_sc / den),
           std::sqrt(num / den) < 1e-5 ? "PASS" : "FAIL");
    if (!val1_pass) return 1;
  }

  // --- 3) timing -----------------------------------------------------------
  {
    std::vector<float> xp(N), bp(N, 0.f), bout(N, 0.f);
    saige_avx2::permute_fwd(x.data(), xp.data(), N);
    float Sx = 0; for (size_t i = 0; i < N; ++i) Sx += x[i];

    auto time_it = [&](auto&& fn, int reps) {
      double best = 1e300;
      for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        const uint64_t c0 = __rdtsc();
        fn();
        const uint64_t c1 = __rdtsc();
        const auto t1 = std::chrono::steady_clock::now();
        (void)c0; (void)c1;
        best = std::min(best,
                        std::chrono::duration<double>(t1 - t0).count());
      }
      return best;
    };
    // TSC frequency for cycle conversion
    const uint64_t ct0 = __rdtsc();
    const auto tt0 = std::chrono::steady_clock::now();
    while (std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         tt0).count() < 0.1) {}
    const double tsc_hz =
        (__rdtsc() - ct0) /
        std::chrono::duration<double>(std::chrono::steady_clock::now() - tt0)
            .count();

    volatile float sink = 0;
    const double t_sc = time_it([&] {
      float acc = 0;
      for (size_t m = 0; m < M; ++m) {
        float v;
        scalar_marker(rows.data() + m * nbyte, N, freq[m], invstd[m], x.data(),
                      bout.data(), &v);
        acc += v;
      }
      sink = acc;
    }, 3);
    const double t_p1 = time_it([&] {
      float acc = 0;
      for (size_t m = 0; m < M; ++m)
        acc += saige_avx2::pass1_sum_gx(rows.data() + m * nbyte, xp.data(), N);
      sink = acc;
    }, 5);
    const double t_p2 = time_it([&] {
      for (size_t m = 0; m < M; ++m)
        saige_avx2::pass2_axpy(rows.data() + m * nbyte, 0.001f, bp.data(), N);
    }, 5);
    const double t_fused = time_it([&] {
      float coff = 0;
      for (size_t m = 0; m < M; ++m) {
        const uint8_t* row = rows.data() + m * nbyte;
        const float raw = saige_avx2::pass1_sum_gx(row, xp.data(), N);
        const float v = invstd[m] * (raw - 2.f * freq[m] * Sx);
        const float sval = v * invstd[m];
        saige_avx2::pass2_axpy(row, sval, bp.data(), N);
        coff += 2.f * freq[m] * sval;
      }
      sink = coff;
    }, 5);
    (void)sink;

    const double G = (double)N * M;
    printf("\nTSC ~= %.2f GHz (nominal; cycles below use this)\n",
           tsc_hz / 1e9);
    printf("%-34s %10s %14s %12s %10s\n", "kernel", "ms", "cyc/genotype",
           "genos/cyc", "GB/s(row)");
    auto rep = [&](const char* name, double t, double passes) {
      printf("%-34s %10.2f %14.4f %12.2f %10.2f\n", name, t * 1e3,
             t * tsc_hz / (G * passes), (G * passes) / (t * tsc_hz),
             passes * M * nbyte / t / 1e9);
    };
    rep("scalar decode+dot+axpy (2 passes)", t_sc, 2.0);
    rep("avx2 pass1 (sum g*x)", t_p1, 1.0);
    rep("avx2 pass2 (b += s*g)", t_p2, 1.0);
    rep("avx2 fused pass1+pass2", t_fused, 2.0);
    printf("\nspeedup scalar/(avx2 fused): %.1fx  (single core)\n",
           t_sc / t_fused);
  }

  // --- 4) Phase-2 multi-RHS kernels ---------------------------------------
  // Correctness: one multi-column call vs k single-column calls (same xp
  // data). Different accumulator trees => not bit-exact; acceptance is
  // rel-L2 < 1e-6 on the val1 vector across markers and on each bout column.
  // Timing: full fused sweep over M markers, multi(k) vs k x single.
  {
    printf("\n=== multi-RHS (Phase 2) ===\n");
    auto time_it = [&](auto&& fn, int reps) {
      double best = 1e300;
      for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        fn();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best,
                        std::chrono::duration<double>(t1 - t0).count());
      }
      return best;
    };
    volatile float sink = 0;
    (void)sink;
    const int ks[] = {2, 4, 8, 16, 30};
    const size_t kmax = 30;
    std::mt19937 gen2(4711);
    std::normal_distribution<float> nd2(0.f, 1.f);
    std::vector<float> X(N * kmax), Xp(N * kmax);
    for (auto& v : X) v = nd2(gen2);
    for (size_t j = 0; j < kmax; ++j)
      saige_avx2::permute_fwd(X.data() + j * N, Xp.data() + j * N, N);
    std::vector<float> sval(kmax);
    for (size_t j = 0; j < kmax; ++j) sval[j] = 0.001f * (float)(j + 1);

    // correctness at kmax (covers the <8 tail path via 30 = 8*3+6)
    {
      std::vector<float> raw_m(kmax), raw_s(kmax);
      double num_v = 0, den_v = 0;
      std::vector<float> bp_m(N * kmax, 0.f), bp_s(N * kmax, 0.f);
      for (size_t m = 0; m < M; ++m) {
        const uint8_t* row = rows.data() + m * nbyte;
        saige_avx2::pass1_sum_gx_multi(row, Xp.data(), N, kmax, N,
                                       raw_m.data());
        saige_avx2::pass2_axpy_multi(row, sval.data(), bp_m.data(), N, kmax,
                                     N);
        for (size_t j = 0; j < kmax; ++j) {
          raw_s[j] = saige_avx2::pass1_sum_gx(row, Xp.data() + j * N, N);
          saige_avx2::pass2_axpy(row, sval[j], bp_s.data() + j * N, N);
          num_v += (double)(raw_m[j] - raw_s[j]) * (raw_m[j] - raw_s[j]);
          den_v += (double)raw_s[j] * raw_s[j];
        }
      }
      double num_b = 0, den_b = 0;
      for (size_t i = 0; i < N * kmax; ++i) {
        num_b += (double)(bp_m[i] - bp_s[i]) * (bp_m[i] - bp_s[i]);
        den_b += (double)bp_s[i] * bp_s[i];
      }
      const double e_v = std::sqrt(num_v / den_v);
      const double e_b = std::sqrt(num_b / den_b);
      printf("[multi] k=%zu vs %zu single calls: rel-L2 pass1 = %.3e, "
             "pass2 = %.3e -> %s\n", kmax, kmax, e_v, e_b,
             (e_v < 1e-6 && e_b < 1e-6) ? "PASS" : "FAIL");
    }

    // timing curve. The multi path replicates the app's CorssProdMat worker:
    // sweep 1 (dots) and sweep 2 (axpys) are sample-blocked so the k-column
    // RHS/bout slice stays L2-resident across the marker range.
    printf("%-6s %14s %16s %12s %14s\n", "k", "multi(ms)", "k x single(ms)",
           "speedup", "ms/column");
    std::vector<float> raw_t(kmax);
    std::vector<float> bp_t(N * kmax, 0.f);
    std::vector<float> raws_t(M * kmax);
    for (int k : ks) {
      size_t sbs = ((size_t)(24576 / k)) & ~(size_t)127;
      if (sbs < 512) sbs = 512;
      const double t_multi = time_it([&] {
        for (size_t i = 0; i < M * (size_t)k; ++i) raws_t[i] = 0.f;
        for (size_t sb = 0; sb < N; sb += sbs) {
          const size_t len = std::min(sbs, N - sb);
          for (size_t m = 0; m < M; ++m) {
            saige_avx2::pass1_sum_gx_multi(rows.data() + m * nbyte + sb / 4,
                                           Xp.data() + sb, N, (size_t)k, len,
                                           raw_t.data());
            for (int j = 0; j < k; ++j) raws_t[m * k + j] += raw_t[j];
          }
        }
        for (size_t sb = 0; sb < N; sb += sbs) {
          const size_t len = std::min(sbs, N - sb);
          for (size_t m = 0; m < M; ++m)
            saige_avx2::pass2_axpy_multi(rows.data() + m * nbyte + sb / 4,
                                         raws_t.data() + m * k,
                                         bp_t.data() + sb, N, (size_t)k, len);
        }
        sink = raws_t[0];
      }, 5);
      const double t_single = time_it([&] {
        float acc = 0;
        for (int j = 0; j < k; ++j)
          for (size_t m = 0; m < M; ++m) {
            const uint8_t* row = rows.data() + m * nbyte;
            acc += saige_avx2::pass1_sum_gx(row, Xp.data() + j * N, N);
            saige_avx2::pass2_axpy(row, sval[j], bp_t.data() + j * N, N);
          }
        sink = acc;
      }, 3);
      printf("%-6d %14.2f %16.2f %11.2fx %14.3f\n", k, t_multi * 1e3,
             t_single * 1e3, t_single / t_multi, t_multi * 1e3 / k);
    }
  }
  return 0;
}
