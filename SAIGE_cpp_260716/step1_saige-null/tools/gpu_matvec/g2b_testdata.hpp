// g2b_testdata.hpp — deterministic synthetic inputs shared by the tier-4
// kernel tests and by the B1 reference-capture tool.
//
// Why a header instead of ad-hoc generation inside each test: B1 compares the
// raw float bits of today's kernel output (captured from the code at
// ce9bacd6, BEFORE scheme C's bind_trait split) against tomorrow's. That
// comparison is only meaningful if both binaries see byte-identical inputs, so
// the generator has to live in one place and must not depend on anything whose
// implementation could drift (hence a hand-rolled splitmix64 rather than
// <random>'s distributions, whose output is not standardized).
//
// Layout produced here matches what tools/bed_reader hands the GPU:
//   packed  : M rows × nbyte(N) bytes, variant-major, 4 samples per byte,
//             sample i in bits [2*(i%4), 2*(i%4)+1] of byte i/4.
//   codes   : 0b00 HOM_REF → g=2, 0b10 HET → g=1, 0b11 HOM_ALT → g=0.
//             0b01 (MISSING) is never emitted — marker_decoder fills missing
//             cells with round(2·altFreq) before packing, and the tier-4
//             kernel's g = 2 − popc(field) has no encoding for it.
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace g2btest {

// splitmix64. Self-contained and bit-stable across compilers/libstdc++.
struct Rng {
  std::uint64_t s;
  explicit Rng(std::uint64_t seed) : s(seed) {}
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  std::uint32_t u32() { return static_cast<std::uint32_t>(next() >> 32); }
  int below(int n) { return static_cast<int>(u32() % static_cast<std::uint32_t>(n)); }
  // Uniform in [-1, 1) on a 2^-24 grid: exactly representable in fp32, so the
  // host never rounds when it stores these.
  float pm1() {
    return static_cast<float>(static_cast<std::int32_t>(u32() >> 8) - (1 << 15)) /
           static_cast<float>(1 << 15);
  }
};

inline std::size_t nbyte_of(int N) { return static_cast<std::size_t>((N + 3) / 4); }

// Genotype of sample i, marker m, straight off the packed bytes.
inline int geno_at(const unsigned char* packed, std::size_t stride, int m, int i) {
  const unsigned char byte = packed[static_cast<std::size_t>(m) * stride + (i >> 2)];
  switch ((byte >> (2 * (i & 3))) & 3) {
    case 0x0: return 2;
    case 0x2: return 1;
    case 0x3: return 0;
    default:  return -1;   // 0b01: cannot occur, see the header comment
  }
}

inline void set_geno(unsigned char* packed, std::size_t stride, int m, int i, int g) {
  const unsigned char code = (g == 2) ? 0x0 : (g == 1) ? 0x2 : 0x3;
  unsigned char& byte = packed[static_cast<std::size_t>(m) * stride + (i >> 2)];
  const int sh = 2 * (i & 3);
  byte = static_cast<unsigned char>((byte & ~(3u << sh)) | (code << sh));
}

// Packed matrix + the per-marker freq/invstd a single-trait run would compute.
// Genotypes are drawn uniformly from {0,1,2} rather than from HWE so that no
// marker is degenerate; freq is a fixed sweep over [0.03, 0.4) so the rank-one
// correction terms are never negligible.
inline void gen_matrix(int N, int M, std::uint64_t seed,
                       std::vector<unsigned char>& packed, std::size_t& stride,
                       std::vector<float>& freq, std::vector<float>& invstd) {
  stride = nbyte_of(N);
  packed.assign(static_cast<std::size_t>(M) * stride, 0xFF);   // pad fields → g=0
  freq.resize(M);
  invstd.resize(M);
  Rng rng(seed);
  static const int kG[3] = {2, 1, 0};
  for (int m = 0; m < M; ++m) {
    for (int i = 0; i < N; ++i) set_geno(packed.data(), stride, m, i, kG[rng.below(3)]);
    freq[m]   = 0.03f + static_cast<float>(m % 100) * 0.0037f;
    invstd[m] = 1.0f / std::sqrt(2.0f * freq[m] * (1.0f - freq[m]));
  }
}

inline void gen_vec(std::size_t n, std::uint64_t seed, std::vector<float>& v) {
  v.resize(n);
  Rng rng(seed);
  for (std::size_t i = 0; i < n; ++i) v[i] = rng.pm1();
}

// FNV-1a over raw bytes — used to pin inputs and outputs in the B1 manifest.
inline std::uint64_t fnv1a(const void* p, std::size_t n) {
  const unsigned char* b = static_cast<const unsigned char*>(p);
  std::uint64_t h = 1469598103934665603ull;
  for (std::size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
  return h;
}

// The shapes B1 sweeps. Every N is ≫ 2048 = one pass-1 x tile (G_TILE_W·16),
// so the multi-tile partial-buffer reduction is actually exercised:
// g1y = ⌈⌈N/16⌉/128⌉ is 3 / 10 / 17 for the three cases.
struct Case { int N; int M; std::uint64_t seed; };
inline const Case* cases(int* n) {
  static const Case kC[] = {
      { 5000,  700, 0x51A6E1ull},
      {20000, 3000, 0x51A6E2ull},
      {33000, 1500, 0x51A6E3ull},
  };
  *n = 3;
  return kC;
}

}  // namespace g2btest
