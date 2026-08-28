// AVX2 fused 2-bit decode kernels for the SAIGE step-1 ψv hot path.
// Design: AVX2_KERNEL_PLAN.md (vpshufb 4-phase decode, amortized interleave
// permutation, two passes per marker with re-decode instead of materialize).
//
// Packed layout (per marker "row" of nbyte = (N+3)/4 bytes):
//   byte i holds samples 4i..4i+3, 2 bits each (bit 2j..2j+1 for sample 4i+j).
//   code c = (byte >> 2j) & 3;  genotype g = LUT4[c] with LUT4 = {2,1,1,0}
//   (identical to the scalar decode g = 2 - (bit1 + bit0)).
//
// Permuted vector layout ("perm order"): for each full superblock s of 128
// samples, phase-major:  xp[128s + 32k + i] = x[128s + 4i + k]  (k=0..3,
// i=0..31). The tail (N mod 128 samples) stays in natural order — the scalar
// tail loop handles it and naturally skips the padding bits of the last byte.
//
// All loads are unaligned (rows are nbyte-strided, nbyte is not a multiple
// of 32). No AVX-512/VL instructions — AVX2 + FMA only.

#pragma once

#include <cstddef>
#include <cstdint>

#if defined(__AVX2__) && defined(__FMA__)
#define SAIGE_AVX2_KERNEL_AVAILABLE 1
#include <immintrin.h>
#else
#define SAIGE_AVX2_KERNEL_AVAILABLE 0
#endif

namespace saige_avx2 {

constexpr bool available() {
#if SAIGE_AVX2_KERNEL_AVAILABLE
  return true;
#else
  return false;
#endif
}

// Scalar 2-bit code -> genotype value {2,1,1,0}. code 0b01 (PLINK MISSING)
// never occurs in QC'd storage (missing is mean-imputed at pack time), but
// LUT4[1] = 1 matches the scalar decode 2-(a+b) anyway.
inline int decode_code(unsigned code) {
  static const int LUT4[4] = {2, 1, 1, 0};
  return LUT4[code & 3u];
}

// ---------------------------------------------------------------------------
// Permutation helpers (O(N), done once per ψv call — amortized over M rows)
// ---------------------------------------------------------------------------

// x (natural order, length N) -> xp (perm order, length N)
inline void permute_fwd(const float* x, float* xp, std::size_t N) {
  const std::size_t S = N / 128;
  for (std::size_t s = 0; s < S; ++s) {
    const float* xs = x + 128 * s;
    float*       ps = xp + 128 * s;
    for (std::size_t k = 0; k < 4; ++k)
      for (std::size_t i = 0; i < 32; ++i)
        ps[32 * k + i] = xs[4 * i + k];
  }
  for (std::size_t i = 128 * S; i < N; ++i) xp[i] = x[i];
}

// bp (perm order) -> b (natural order), applying b = bp_unperm - coffset
inline void unpermute_sub(const float* bp, float* b, std::size_t N,
                          float coffset) {
  const std::size_t S = N / 128;
  for (std::size_t s = 0; s < S; ++s) {
    const float* ps = bp + 128 * s;
    float*       bs = b + 128 * s;
    for (std::size_t k = 0; k < 4; ++k)
      for (std::size_t i = 0; i < 32; ++i)
        bs[4 * i + k] = ps[32 * k + i] - coffset;
  }
  for (std::size_t i = 128 * S; i < N; ++i) b[i] = bp[i] - coffset;
}

#if SAIGE_AVX2_KERNEL_AVAILABLE

namespace detail {

// 16-entry in-lane LUT for vpshufb: indices 0..3 -> {2,1,1,0}, rest unused.
inline __m256i geno_lut() {
  return _mm256_setr_epi8(2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

// Extract phase k (k = 0..3) of 32 packed bytes: result byte i = genotype of
// sample 4i+k within this 128-sample superblock. Cross-byte bits that vpsrlw
// smears into bits 2..7 are masked off before the table lookup.
inline __m256i decode_phase(__m256i B, int k, __m256i lut, __m256i m3) {
  __m256i sh;
  switch (k) {
    case 0: sh = B; break;
    case 1: sh = _mm256_srli_epi16(B, 2); break;
    case 2: sh = _mm256_srli_epi16(B, 4); break;
    default: sh = _mm256_srli_epi16(B, 6); break;
  }
  return _mm256_shuffle_epi8(lut, _mm256_and_si256(sh, m3));
}

// Widen 8 genotype bytes (starting at byte j*8 of gk) to 8 fp32 lanes.
inline __m256 widen8(__m128i g8) {
  return _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(g8));
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Pass 1: raw = sum_i g[i] * xp[i]  (xp in perm order).
// Caller applies the rank-one identity: val1 = invStd * (raw - 2*freq*Sx).
// ---------------------------------------------------------------------------
inline float pass1_sum_gx(const std::uint8_t* row, const float* xp,
                          std::size_t N) {
  const __m256i lut = detail::geno_lut();
  const __m256i m3  = _mm256_set1_epi8(0x03);
  const std::size_t S = N / 128;

  __m256 acc0 = _mm256_setzero_ps(), acc1 = _mm256_setzero_ps();
  __m256 acc2 = _mm256_setzero_ps(), acc3 = _mm256_setzero_ps();

  for (std::size_t s = 0; s < S; ++s) {
    const __m256i B =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + 32 * s));
    const float* xs = xp + 128 * s;
#define SAIGE_P1_PHASE(K, ACC)                                              \
    {                                                                       \
      const __m256i gk = detail::decode_phase(B, K, lut, m3);               \
      const __m128i lo = _mm256_castsi256_si128(gk);                        \
      const __m128i hi = _mm256_extracti128_si256(gk, 1);                   \
      const float* xk = xs + 32 * (K);                                      \
      ACC = _mm256_fmadd_ps(detail::widen8(lo), _mm256_loadu_ps(xk), ACC);  \
      ACC = _mm256_fmadd_ps(detail::widen8(_mm_srli_si128(lo, 8)),          \
                            _mm256_loadu_ps(xk + 8), ACC);                  \
      ACC = _mm256_fmadd_ps(detail::widen8(hi),                             \
                            _mm256_loadu_ps(xk + 16), ACC);                 \
      ACC = _mm256_fmadd_ps(detail::widen8(_mm_srli_si128(hi, 8)),          \
                            _mm256_loadu_ps(xk + 24), ACC);                 \
    }
    SAIGE_P1_PHASE(0, acc0)
    SAIGE_P1_PHASE(1, acc1)
    SAIGE_P1_PHASE(2, acc2)
    SAIGE_P1_PHASE(3, acc3)
#undef SAIGE_P1_PHASE
  }

  // Horizontal sum of the 4 accumulators.
  __m256 acc = _mm256_add_ps(_mm256_add_ps(acc0, acc1),
                             _mm256_add_ps(acc2, acc3));
  __m128 v4  = _mm_add_ps(_mm256_castps256_ps128(acc),
                          _mm256_extractf128_ps(acc, 1));
  v4 = _mm_add_ps(v4, _mm_movehl_ps(v4, v4));
  v4 = _mm_add_ss(v4, _mm_shuffle_ps(v4, v4, 1));
  float raw = _mm_cvtss_f32(v4);

  // Scalar tail (natural order in xp; skips last-byte padding bits).
  for (std::size_t i = 128 * S; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    raw += static_cast<float>(decode_code(code)) * xp[i];
  }
  return raw;
}

// ---------------------------------------------------------------------------
// Pass 2: bp[i] += sval * g[i]  (bp in perm order; caller subtracts the
// uniform offset 2*freq*sval once per marker via the accumulated Coffset).
// ---------------------------------------------------------------------------
inline void pass2_axpy(const std::uint8_t* row, float sval, float* bp,
                       std::size_t N) {
  const __m256i lut = detail::geno_lut();
  const __m256i m3  = _mm256_set1_epi8(0x03);
  const __m256  sv  = _mm256_set1_ps(sval);
  const std::size_t S = N / 128;

  for (std::size_t s = 0; s < S; ++s) {
    const __m256i B =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + 32 * s));
    float* bs = bp + 128 * s;
#define SAIGE_P2_PHASE(K)                                                   \
    {                                                                       \
      const __m256i gk = detail::decode_phase(B, K, lut, m3);               \
      const __m128i lo = _mm256_castsi256_si128(gk);                        \
      const __m128i hi = _mm256_extracti128_si256(gk, 1);                   \
      float* bk = bs + 32 * (K);                                            \
      _mm256_storeu_ps(bk, _mm256_fmadd_ps(detail::widen8(lo), sv,          \
                                           _mm256_loadu_ps(bk)));           \
      _mm256_storeu_ps(bk + 8,                                              \
          _mm256_fmadd_ps(detail::widen8(_mm_srli_si128(lo, 8)), sv,        \
                          _mm256_loadu_ps(bk + 8)));                        \
      _mm256_storeu_ps(bk + 16, _mm256_fmadd_ps(detail::widen8(hi), sv,     \
                                                _mm256_loadu_ps(bk + 16))); \
      _mm256_storeu_ps(bk + 24,                                             \
          _mm256_fmadd_ps(detail::widen8(_mm_srli_si128(hi, 8)), sv,        \
                          _mm256_loadu_ps(bk + 24)));                       \
    }
    SAIGE_P2_PHASE(0)
    SAIGE_P2_PHASE(1)
    SAIGE_P2_PHASE(2)
    SAIGE_P2_PHASE(3)
#undef SAIGE_P2_PHASE
  }

  for (std::size_t i = 128 * S; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    bp[i] += sval * static_cast<float>(decode_code(code));
  }
}

// ---------------------------------------------------------------------------
// Phase-2 multi-RHS kernels (block-PCG). k columns share ONE decode of the
// packed row. Per 128-sample superblock the row is decoded once into an
// L1-resident fp32 buffer gbuf[128] (in perm order — decode_phase's byte
// order IS the perm order), then a light per-column loop runs 16 FMA groups
// against gbuf. This keeps register pressure flat in k (an early variant
// with k accumulators live in the inner loop spilled and ran SLOWER than k
// single-column calls).
//
// Numerics: pass1 uses 4 accumulator chains keyed by phase (groups 0-3 →
// chain0, 4-7 → chain1, ...) and the same ((a0+a1)+(a2+a3)) horizontal
// reduction as pass1_sum_gx — the per-column result is BIT-EXACT equal to a
// single-column pass1_sum_gx call. pass2 is elementwise, also bit-exact.
//
// Columns live in a column-major matrix (arma::fmat layout): column j starts
// at xp + j*ldx, each column independently in perm order (permute_fwd/col).
// ---------------------------------------------------------------------------

namespace detail {

// Decode one 32-byte superblock (128 samples) of the packed row into fp32
// genotypes in perm order. gbuf must be 32-byte aligned.
inline void decode_block_ps(__m256i B, __m256i lut, __m256i m3, float* gbuf) {
#define SAIGE_MC_DECODE_PHASE(KPH)                                          \
  {                                                                         \
    const __m256i gk = decode_phase(B, KPH, lut, m3);                       \
    const __m128i lo = _mm256_castsi256_si128(gk);                          \
    const __m128i hi = _mm256_extracti128_si256(gk, 1);                     \
    float* gb = gbuf + 32 * (KPH);                                          \
    _mm256_store_ps(gb, widen8(lo));                                        \
    _mm256_store_ps(gb + 8, widen8(_mm_srli_si128(lo, 8)));                 \
    _mm256_store_ps(gb + 16, widen8(hi));                                   \
    _mm256_store_ps(gb + 24, widen8(_mm_srli_si128(hi, 8)));                \
  }
  SAIGE_MC_DECODE_PHASE(0)
  SAIGE_MC_DECODE_PHASE(1)
  SAIGE_MC_DECODE_PHASE(2)
  SAIGE_MC_DECODE_PHASE(3)
#undef SAIGE_MC_DECODE_PHASE
}

}  // namespace detail

// Max columns per call (bounds the on-stack accumulator arrays: 64 cols x 4
// chains x 32 B = 8 KB). Callers with more RHS chunk at this granularity.
constexpr std::size_t kMaxMultiCols = 64;

// Multi-column pass 1: raw[j] = sum_i g[i] * xp[i + j*ldx], j = 0..k-1.
// Per column bit-exact equal to pass1_sum_gx(row, xp + j*ldx, N).
inline void pass1_sum_gx_multi(const std::uint8_t* row, const float* xp,
                               std::size_t ldx, std::size_t k, std::size_t N,
                               float* raw) {
  while (k > kMaxMultiCols) {  // chunk very wide batches (row stays cached)
    pass1_sum_gx_multi(row, xp, ldx, kMaxMultiCols, N, raw);
    xp += kMaxMultiCols * ldx;
    raw += kMaxMultiCols;
    k -= kMaxMultiCols;
  }
  const __m256i lut = detail::geno_lut();
  const __m256i m3  = _mm256_set1_epi8(0x03);
  const std::size_t S = N / 128;
  alignas(32) float gbuf[128];

  __m256 acc0[kMaxMultiCols], acc1[kMaxMultiCols];
  __m256 acc2[kMaxMultiCols], acc3[kMaxMultiCols];
  const __m256 z = _mm256_setzero_ps();
  for (std::size_t j = 0; j < k; ++j) {
    acc0[j] = z; acc1[j] = z; acc2[j] = z; acc3[j] = z;
  }

  for (std::size_t s = 0; s < S; ++s) {
    const __m256i B =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + 32 * s));
    detail::decode_block_ps(B, lut, m3, gbuf);
    for (std::size_t j = 0; j < k; ++j) {
      const float* xj = xp + j * ldx + 128 * s;
      __m256 a0 = acc0[j], a1 = acc1[j], a2 = acc2[j], a3 = acc3[j];
#define SAIGE_MC_P1_GROUP(T, ACC)                                           \
      ACC = _mm256_fmadd_ps(_mm256_load_ps(gbuf + 8 * (T)),                 \
                            _mm256_loadu_ps(xj + 8 * (T)), ACC);
      SAIGE_MC_P1_GROUP(0, a0)  SAIGE_MC_P1_GROUP(1, a0)
      SAIGE_MC_P1_GROUP(2, a0)  SAIGE_MC_P1_GROUP(3, a0)
      SAIGE_MC_P1_GROUP(4, a1)  SAIGE_MC_P1_GROUP(5, a1)
      SAIGE_MC_P1_GROUP(6, a1)  SAIGE_MC_P1_GROUP(7, a1)
      SAIGE_MC_P1_GROUP(8, a2)  SAIGE_MC_P1_GROUP(9, a2)
      SAIGE_MC_P1_GROUP(10, a2) SAIGE_MC_P1_GROUP(11, a2)
      SAIGE_MC_P1_GROUP(12, a3) SAIGE_MC_P1_GROUP(13, a3)
      SAIGE_MC_P1_GROUP(14, a3) SAIGE_MC_P1_GROUP(15, a3)
#undef SAIGE_MC_P1_GROUP
      acc0[j] = a0; acc1[j] = a1; acc2[j] = a2; acc3[j] = a3;
    }
  }

  for (std::size_t j = 0; j < k; ++j) {
    // Same reduction tree as pass1_sum_gx: (a0+a1)+(a2+a3), then hsum.
    __m256 acc = _mm256_add_ps(_mm256_add_ps(acc0[j], acc1[j]),
                               _mm256_add_ps(acc2[j], acc3[j]));
    __m128 v4  = _mm_add_ps(_mm256_castps256_ps128(acc),
                            _mm256_extractf128_ps(acc, 1));
    v4 = _mm_add_ps(v4, _mm_movehl_ps(v4, v4));
    v4 = _mm_add_ss(v4, _mm_shuffle_ps(v4, v4, 1));
    raw[j] = _mm_cvtss_f32(v4);
  }

  // Scalar tail (natural order in every column; skips padding bits).
  for (std::size_t i = 128 * S; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    const float g = static_cast<float>(decode_code(code));
    for (std::size_t j = 0; j < k; ++j) raw[j] += g * xp[j * ldx + i];
  }
}

// Multi-column pass 2: bp[i + j*ldb] += sval[j] * g[i] (columns perm order).
// Per column bit-exact equal to pass2_axpy(row, sval[j], bp + j*ldb, N).
inline void pass2_axpy_multi(const std::uint8_t* row, const float* sval,
                             float* bp, std::size_t ldb, std::size_t k,
                             std::size_t N) {
  while (k > kMaxMultiCols) {
    pass2_axpy_multi(row, sval, bp, ldb, kMaxMultiCols, N);
    sval += kMaxMultiCols;
    bp += kMaxMultiCols * ldb;
    k -= kMaxMultiCols;
  }
  const __m256i lut = detail::geno_lut();
  const __m256i m3  = _mm256_set1_epi8(0x03);
  const std::size_t S = N / 128;
  alignas(32) float gbuf[128];

  for (std::size_t s = 0; s < S; ++s) {
    const __m256i B =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + 32 * s));
    detail::decode_block_ps(B, lut, m3, gbuf);
    for (std::size_t j = 0; j < k; ++j) {
      float* bj = bp + j * ldb + 128 * s;
      const __m256 sv = _mm256_set1_ps(sval[j]);
#define SAIGE_MC_P2_GROUP(T)                                                \
      _mm256_storeu_ps(bj + 8 * (T),                                        \
          _mm256_fmadd_ps(_mm256_load_ps(gbuf + 8 * (T)), sv,               \
                          _mm256_loadu_ps(bj + 8 * (T))));
      SAIGE_MC_P2_GROUP(0)  SAIGE_MC_P2_GROUP(1)  SAIGE_MC_P2_GROUP(2)
      SAIGE_MC_P2_GROUP(3)  SAIGE_MC_P2_GROUP(4)  SAIGE_MC_P2_GROUP(5)
      SAIGE_MC_P2_GROUP(6)  SAIGE_MC_P2_GROUP(7)  SAIGE_MC_P2_GROUP(8)
      SAIGE_MC_P2_GROUP(9)  SAIGE_MC_P2_GROUP(10) SAIGE_MC_P2_GROUP(11)
      SAIGE_MC_P2_GROUP(12) SAIGE_MC_P2_GROUP(13) SAIGE_MC_P2_GROUP(14)
      SAIGE_MC_P2_GROUP(15)
#undef SAIGE_MC_P2_GROUP
    }
  }

  for (std::size_t i = 128 * S; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    const float g = static_cast<float>(decode_code(code));
    for (std::size_t j = 0; j < k; ++j) bp[j * ldb + i] += sval[j] * g;
  }
}

// Bench/verify helper: decode a full row into genotype bytes in PERM order
// (g_perm[128s+32k+i] = g[sample 128s+4i+k]; tail natural). Bit-exactness of
// the vpshufb decode is validated against the scalar decode through this.
inline void decode_row_perm(const std::uint8_t* row, std::uint8_t* g_perm,
                            std::size_t N) {
  const __m256i lut = detail::geno_lut();
  const __m256i m3  = _mm256_set1_epi8(0x03);
  const std::size_t S = N / 128;
  for (std::size_t s = 0; s < S; ++s) {
    const __m256i B =
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(row + 32 * s));
    for (int k = 0; k < 4; ++k) {
      const __m256i gk = detail::decode_phase(B, k, lut, m3);
      _mm256_storeu_si256(
          reinterpret_cast<__m256i*>(g_perm + 128 * s + 32 * k), gk);
    }
  }
  for (std::size_t i = 128 * S; i < N; ++i) {
    const unsigned code = (row[i >> 2] >> (2 * (i & 3))) & 3u;
    g_perm[i] = static_cast<std::uint8_t>(decode_code(code));
  }
}

#endif  // SAIGE_AVX2_KERNEL_AVAILABLE

}  // namespace saige_avx2
