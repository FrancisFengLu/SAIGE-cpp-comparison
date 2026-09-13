#include "marker_decoder.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace saige {

namespace {
// Packed-output 2-bit codes — mirror SAIGE_step1_fast.cpp:38-41.
constexpr unsigned char HOM_REF = 0x0;  // 0b00 — bufferGeno == 2
constexpr unsigned char HET     = 0x2;  // 0b10 — bufferGeno == 1
constexpr unsigned char HOM_ALT = 0x3;  // 0b11 — bufferGeno == 0
// MISSING is never written; we fill-in before packing.
} // namespace

void build_bed_lookup(BedLut& out) {
  // Mirror of SAIGE_step1_fast.cpp lines 836-848:
  //   lo,hi from bits [2j, 2j+1]:
  //     01 → 3 (missing), 00 → 2 (hom A1), 10 → 1 (het), 11 → 0 (hom A2)
  for (int byte = 0; byte < 256; ++byte) {
    int b = byte;
    for (int j = 0; j < 4; ++j) {
      const int lo = b & 1; b >>= 1;
      const int hi = b & 1; b >>= 1;
      if      (lo == 1 && hi == 0) out[byte][j] = 3;  // missing
      else if (lo == 0 && hi == 0) out[byte][j] = 2;  // hom A1
      else if (lo == 0 && hi == 1) out[byte][j] = 1;  // het
      else                         out[byte][j] = 0;  // hom A2
    }
  }
}

void decode_marker(const unsigned char* raw,
                   std::size_t N,
                   const int* ptrsub, std::size_t Nnomissing,
                   float min_maf, float max_miss,
                   const BedLut& lut,
                   MarkerStats& stats,
                   unsigned char* packed_out) {
  const VarRatioRule no_vr;
  bool passVR = false;
  decode_marker(raw, N, ptrsub, Nnomissing, min_maf, max_miss, lut,
                no_vr, /*vr_drawn=*/false, stats, passVR, packed_out);
}

void decode_marker(const unsigned char* raw,
                   std::size_t N,
                   const int* ptrsub, std::size_t Nnomissing,
                   float min_maf, float max_miss,
                   const BedLut& lut,
                   const VarRatioRule& vr, bool vr_drawn,
                   MarkerStats& stats,
                   bool& passVR,
                   unsigned char* packed_out) {
  const std::size_t nbyte_in  = (N + 3) / 4;
  const std::size_t nbyte_out = (Nnomissing + 3) / 4;

  // Pass 1: decode raw BED byte stream into a length-N per-FAM-sample
  // bufferGeno vector on the stack/heap. We only SUM + COUNT over the
  // phenotyped samples, using ptrsub (which is 1-based FAM index per GRM slot).
  //
  // Using a fixed-size 16KB stack buffer when possible; otherwise heap.
  // On UKB N=408970 → 409 KB, need heap.
  // (SAIGE's current code uses a thread_local N-length vector.)
  static thread_local std::vector<int> genoAll;
  if (genoAll.size() < N) genoAll.assign(N, 0);

  for (std::size_t i = 0; i < nbyte_in; ++i) {
    const unsigned char byte = raw[i];
    const int* l = lut[byte];
    const std::size_t base = i * 4;
    for (int j = 0; j < 4; ++j) {
      const std::size_t fam_idx = base + j;
      if (fam_idx >= N) break;
      genoAll[fam_idx] = l[j];
    }
  }

  // Tally over phenotyped samples (GRM slots 0..Nnomissing-1).
  int alleleCount = 0;
  int numMissing  = 0;
  for (std::size_t k = 0; k < Nnomissing; ++k) {
    const int g = genoAll[ptrsub[k] - 1];
    if (g == 3) {
      ++numMissing;
    } else {
      alleleCount += g;
    }
  }

  // First-pass altFreq over *non-missing* phenotyped samples.
  const int    n_nonmiss = static_cast<int>(Nnomissing) - numMissing;
  const float  altFreq_pre = (n_nonmiss > 0)
      ? alleleCount / static_cast<float>(n_nonmiss * 2)
      : 0.0f;
  const float  missingRate = (Nnomissing > 0)
      ? numMissing / static_cast<float>(Nnomissing)
      : 0.0f;

  // Fill-in: missings become round(2 * altFreq_pre). Update alleleCount to
  // reflect the post-fill tally, then recompute altFreq over all Nnomissing.
  const int fillin = static_cast<int>(std::round(2.0f * altFreq_pre));
  if (numMissing > 0) alleleCount += fillin * numMissing;
  const float altFreq = (Nnomissing > 0)
      ? alleleCount / static_cast<float>(Nnomissing * 2)
      : 0.0f;
  const float maf = std::min(altFreq, 1.0f - altFreq);
  const int   mac = std::min(alleleCount,
                             static_cast<int>(Nnomissing) * 2 - alleleCount);
  bool passQC = (maf >= min_maf) && (missingRate <= max_miss);

  // Variance-ratio claim — byte-for-byte the rule in
  // SAIGE_step1_fast.cpp:518-571. Note `mac` here is the same post-fill MAC
  // the serial path tests, so the two paths select the same markers.
  passVR = false;
  if (vr.enabled) {
    if (vr.max_mac != -1.0f) {                       // categorical VR bins
      if (mac >= vr.min_mac && mac < vr.max_mac) {
        passVR = true;
      } else if (mac >= vr.max_mac) {
        passVR = vr_drawn;
      }
    } else {                                         // single common-MAC bin
      if (mac >= vr.min_mac) passVR = vr_drawn;
    }
    // A VR marker never contributes to the GRM.
    if (passVR) passQC = false;
  }

  stats.altFreq     = altFreq;
  stats.missingRate = missingRate;
  stats.alleleCount = alleleCount;
  stats.mac         = mac;
  stats.numMissing  = numMissing;
  stats.passQC      = passQC;

  if (!passQC && !passVR) {
    // Don't waste cycles packing; caller won't use packed_out.
    std::memset(packed_out, 0, nbyte_out);
    return;
  }

  // Pass 2: re-pack Nnomissing samples (with missings filled in) into ⌈N'/4⌉
  // bytes, using SAIGE's HOM_ALT/HET/HOM_REF 2-bit convention.
  std::memset(packed_out, 0, nbyte_out);
  unsigned char geno2 = 0;
  for (std::size_t idx = 0; idx < Nnomissing; ++idx) {
    const int u = static_cast<int>(idx & 3);
    int buf = genoAll[ptrsub[idx] - 1];
    if (buf == 3) buf = fillin;
    unsigned char code = 0;
    switch (buf) {
      case 0: code = HOM_ALT; break;
      case 1: code = HET;     break;
      case 2: code = HOM_REF; break;
      default: /* shouldn't happen */ break;
    }
    geno2 |= static_cast<unsigned char>(code << (u << 1));

    if (u == 3 || idx == Nnomissing - 1) {
      packed_out[idx >> 2] = geno2;
      geno2 = 0;
    }
  }
}

} // namespace saige
