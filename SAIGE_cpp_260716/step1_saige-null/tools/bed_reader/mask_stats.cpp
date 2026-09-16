#include "mask_stats.hpp"

#include <cmath>
#include <stdexcept>

namespace saige {

namespace {

// SAIGE_step1_fast.cpp:1123 — the only place invstd is born. Reproduced here
// rather than referenced, because the parallel loader computes it inline.
inline float inv_std_from_freq(float altFreq) {
  const float Std = std::sqrt(2.0f * altFreq * (1.0f - altFreq));
  return (Std == 0.0f) ? 0.0f : 1.0f / Std;
}

} // namespace

void compute_trait_stats(const MarkerStats* union_stats, std::size_t M,
                         const TraitTally* tally, int P, int t,
                         int n_t,
                         float min_maf, float max_miss,
                         const VarRatioRule& vr, const unsigned char* vr_drawn,
                         TraitMarkerStats& out) {
  if (union_stats == nullptr) throw std::runtime_error("compute_trait_stats: null union_stats");
  if (tally != nullptr && (t < 0 || t >= P))
    throw std::runtime_error("compute_trait_stats: trait index out of range");
  if (n_t <= 0) throw std::runtime_error("compute_trait_stats: n_t must be > 0");

  out.freq.resize(M);
  out.invstd.resize(M);
  out.missingRate.resize(M);
  out.passQC.resize(M);
  out.passVR.resize(M);
  out.fill.resize(M);
  out.mac.resize(M);
  out.numMissing.resize(M);
  out.alleleCount.resize(M);
  out.M_t = 0;

  const std::size_t stride = static_cast<std::size_t>(P);
  for (std::size_t j = 0; j < M; ++j) {
    const MarkerStats& u = union_stats[j];
    int a_excl = 0, m_excl = 0;
    if (tally != nullptr) {
      const TraitTally& d = tally[j * stride + static_cast<std::size_t>(t)];
      a_excl = d.alleleRawExcl;
      m_excl = d.numMissingExcl;
    }

    MarkerStats s;
    bool        pvr = false;
    const bool  drawn = (vr.enabled && vr_drawn != nullptr) ? (vr_drawn[j] != 0) : false;
    marker_stats_from_counts(u.alleleRaw - a_excl, u.numMissing - m_excl,
                             static_cast<std::size_t>(n_t),
                             min_maf, max_miss, vr, drawn, s, pvr);

    out.freq[j]        = s.altFreq;
    out.missingRate[j] = s.missingRate;
    out.passQC[j]      = s.passQC ? 1 : 0;
    out.passVR[j]      = pvr ? 1 : 0;
    out.fill[j]        = s.fillin;
    out.mac[j]         = s.mac;
    out.numMissing[j]  = s.numMissing;
    out.alleleCount[j] = s.alleleCount;
    // §1 step 3: a marker this trait fails contributes nothing to its GRM.
    out.invstd[j]      = s.passQC ? inv_std_from_freq(s.altFreq) : 0.0f;
    if (s.passQC) ++out.M_t;
  }
}

void build_fill_corrections(const int* union_fill,
                            const int* trait_fill,
                            std::size_t M,
                            const std::vector<std::vector<int>>& missing_cells,
                            const unsigned char* in_trait,
                            FillCorrection& out) {
  out.row.clear();
  out.col.clear();
  out.delta.clear();
  if (missing_cells.size() < M)
    throw std::runtime_error("build_fill_corrections: missing_cells shorter than M");

  for (std::size_t j = 0; j < M; ++j) {
    const int d = trait_fill[j] - union_fill[j];
    if (d == 0) continue;                       // whole column needs nothing
    const std::vector<int>& rows = missing_cells[j];
    for (int i : rows) {
      if (in_trait[i] == 0) continue;           // masked out anyway, x_i == 0
      out.row.push_back(i);
      out.col.push_back(static_cast<int>(j));
      out.delta.push_back(static_cast<float>(d));
    }
  }
}

void union_keep_flags(const MarkerStats* union_stats, std::size_t M,
                      const TraitTally* tally, int P, const int* n_t,
                      float min_maf, float max_miss,
                      const VarRatioRule& vr, const unsigned char* vr_drawn,
                      std::vector<char>& keep_out) {
  keep_out.assign(M, 0);
  const std::size_t stride = static_cast<std::size_t>(P);
  for (std::size_t j = 0; j < M; ++j) {
    const MarkerStats& u = union_stats[j];
    const bool drawn = (vr.enabled && vr_drawn != nullptr) ? (vr_drawn[j] != 0) : false;
    for (int t = 0; t < P; ++t) {
      const TraitTally& d = tally[j * stride + static_cast<std::size_t>(t)];
      MarkerStats s;
      bool        pvr = false;
      marker_stats_from_counts(u.alleleRaw - d.alleleRawExcl,
                               u.numMissing - d.numMissingExcl,
                               static_cast<std::size_t>(n_t[t]),
                               min_maf, max_miss, vr, drawn, s, pvr);
      if (s.passQC) { keep_out[j] = 1; break; }
    }
  }
}

} // namespace saige
