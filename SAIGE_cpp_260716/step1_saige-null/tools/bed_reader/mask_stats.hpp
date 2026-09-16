// mask_stats.hpp — §3.2 of optimization/missing_mt/SCHEME_C_DESIGN.md.
//
// Scheme C loads the genotype matrix once over the UNION U of the traits'
// sample sets and masks rows per trait, instead of re-loading once per distinct
// sample set. Everything a trait needs that is *not* the matrix itself is
// rebuilt here from the union's per-marker counters plus the per-trait
// deductions the same decode pass produced (marker_decoder.hpp: TraitTally):
//
//   alleleRaw_t[j]  = alleleRaw_U[j]  - tally[j][t].alleleRawExcl
//   numMissing_t[j] = numMissing_U[j] - tally[j][t].numMissingExcl
//
// Both are integer subtractions, so they are exact, and the stats that follow
// are computed by marker_stats_from_counts() — literally the decoder's own
// block. That is what makes gate A1 (bit-identical freq / invstd / fill /
// passQC / mac / M_t vs decoding S_t on its own) achievable at all: the design
// doc §4 demands bit-identity here, and only shared code can promise it.
//
// The one quantity that is NOT an exact subtraction is the *fill* value: a
// trait's round(2*f_pre) can differ from the union's. The union matrix stores
// fill_U in every missing cell, so each trait needs a correction list — see
// build_fill_corrections().
#pragma once

#include "marker_decoder.hpp"

#include <cstddef>
#include <vector>

namespace saige {

// Per-trait, per-marker stats over the union's marker order (length M).
struct TraitMarkerStats {
  std::vector<float>         freq;         // altFreq_t (post-fill)
  std::vector<float>         invstd;       // 1/sqrt(2 f (1-f)); 0 when !passQC
  std::vector<float>         missingRate;
  std::vector<unsigned char> passQC;
  std::vector<unsigned char> passVR;
  std::vector<int>           fill;         // 0/1/2 — round(2*f_pre_t)
  std::vector<int>           mac;
  std::vector<int>           numMissing;
  std::vector<int>           alleleCount;  // post-fill
  int                        M_t = 0;      // #{j : passQC_t[j]}
};

// Rebuild trait `t`'s per-marker stats from the union's stats + the tallies.
//
//   union_stats — length M, as parallel_decode_bed produced over U
//   tally       — length M*P, marker-major (tally[j*P + t]); may be null when
//                 P == 1 and the trait IS the union (no deduction)
//   vr_drawn    — length M 0/1, or null; only read when vr.enabled
//
// invstd mirrors SAIGE_step1_fast.cpp:1123 exactly:
//   Std = sqrt(2.0f * freq * (1.0f - freq)); invStd = Std == 0 ? 0 : 1/Std
void compute_trait_stats(const MarkerStats* union_stats, std::size_t M,
                         const TraitTally* tally, int P, int t,
                         int n_t,
                         float min_maf, float max_miss,
                         const VarRatioRule& vr, const unsigned char* vr_drawn,
                         TraitMarkerStats& out);

// §1 step 4: the union matrix holds fill_U[j] in every missing cell; trait t
// needs fill_t[j] there. One entry per (union row i, marker j) with
//   i is missing at j,  i ∈ S_t,  fill_t[j] != fill_U[j]
// and delta = fill_t[j] - fill_U[j] ∈ {-2,-1,1,2}.
// Sorted by marker ascending, and within a marker by row ascending.
struct FillCorrection {
  std::vector<int>   row;    // union-local sample index
  std::vector<int>   col;    // marker index (union marker order)
  std::vector<float> delta;
  std::size_t size() const { return col.size(); }
};

// `missing_cells` is parallel_decode_bed's per-marker list of union-local rows
// that are missing (length M; empty entries for markers with no missing).
// `in_trait` is length n_union: non-zero iff that union row belongs to S_t.
void build_fill_corrections(const int* union_fill,
                            const int* trait_fill,
                            std::size_t M,
                            const std::vector<std::vector<int>>& missing_cells,
                            const unsigned char* in_trait,
                            FillCorrection& out);

// Convenience: the same rule decode_marker's DecodeAux::keep_any_trait applies
// (§2), recomputed from the tallies. Provided so a caller that already has the
// tallies can re-derive the flag without a second decode; parallel_decode_bed
// already returns it in `keep_union`.
void union_keep_flags(const MarkerStats* union_stats, std::size_t M,
                      const TraitTally* tally, int P, const int* n_t,
                      float min_maf, float max_miss,
                      const VarRatioRule& vr, const unsigned char* vr_drawn,
                      std::vector<char>& keep_out);

} // namespace saige
