// marker_decoder.hpp — PR-2 of PARALLEL_BED_PLAN_2026-04-22.md.
//
// Purpose: take the ⌈N/4⌉ raw BED bytes produced by BedReaderPool and produce
//   (a) MAF / MAC / missing-rate / passQC for one marker,
//   (b) the Nnomissing-sample re-packed bytes suitable for downstream GRM math.
//
// This is a faithful port of SAIGE_step1_fast.cpp::Get_OneSNP_Geno_atBeginning,
// kept in the tools/bed_reader/ sandbox so the production build is untouched.
// Once PR-2 passes its byte-level regression it becomes the core of the new
// parallel `setGenoObj`.
//
// BED 2-bit decoding conventions (identical to SAIGE's):
//   BED bit pair (lo,hi)   → bufferGeno   meaning
//       00                → 2            hom A1     (2 alt-ref alleles)
//       01                → 3            missing
//       10                → 1            het
//       11                → 0            hom A2
// The `altFreq` computed here is what SAIGE's current code calls `altFreq`,
// i.e., mean(bufferGeno) / 2 after filling missings with round(2*altFreq_hat).
//
// Re-packing convention (matches SAIGE's setGenotype() + HOM_ALT/HET/HOM_REF
// constants in SAIGE_step1_fast.cpp:38-41):
//   bufferGeno = 0  → HOM_ALT = 0b11 = 0x3
//   bufferGeno = 1  → HET     = 0b10 = 0x2
//   bufferGeno = 2  → HOM_REF = 0b00 = 0x0
//   (missing gets filled first, then mapped like above)
// Each packed byte holds 4 samples, sample j in bits [2j, 2j+1].
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace saige {

struct MarkerStats {
  float altFreq     = 0.0f;   // post-fill allele frequency over Nnomissing samples
  float missingRate = 0.0f;   // numMissing / Nnomissing (pre-fill)
  int   alleleCount = 0;      // post-fill sum of bufferGeno (0/1/2) over Nnomissing
  int   mac         = 0;      // min(alleleCount, 2*Nnomissing - alleleCount)
  int   numMissing  = 0;      // raw missings before fill-in
  bool  passQC      = false;  // altFreq-based MAF ≥ min_maf  AND  missingRate ≤ max_miss
  // Scheme C (SCHEME_C_DESIGN.md §1) needs the two pre-fill quantities that the
  // original code kept only as locals: the raw (pre-fill) allele sum, so a
  // subset's raw sum is an integer subtraction away, and the fill value, so the
  // per-trait fill correction table can be built.
  int   alleleRaw   = 0;      // pre-fill sum of bufferGeno over NON-MISSING samples
  int   fillin      = 0;      // round(2*altFreq_pre), the value missings took
};

// Variance-ratio marker rule — mirrors SAIGE_step1_fast.cpp:518-571
// (Get_OneSNP_Geno_atBeginning, the isVarRatio block).
//
// A marker is claimed for the variance-ratio pool when
//   max_mac != -1 :  min_mac <= mac <  max_mac                       (categorical bin)
//               or:  mac >= max_mac  AND  the marker was drawn into the random pool
//   max_mac == -1 :  mac >= min_mac  AND  the marker was drawn into the random pool
// A claimed marker is REMOVED from the GRM (passQC forced false) so the GRM
// and the VR pool never share a marker.
struct VarRatioRule {
  bool  enabled = false;
  float min_mac = 0.0f;    // genoClass::g_minMACVarRatio
  float max_mac = -1.0f;   // genoClass::g_maxMACVarRatio (-1 == non-categorical)
};

// ---------------------------------------------------------------------------
// Scheme C (optimization/missing_mt/SCHEME_C_DESIGN.md) — optional extra
// outputs of the same decode pass. All of this is inert unless the caller
// passes a DecodeAux; with aux == nullptr decode_marker does exactly what it
// did before, instruction for instruction.
// ---------------------------------------------------------------------------

// Per-trait rows of the union that the trait does NOT own (§3.1).
// `excl[t]` holds UNION-LOCAL indices, i.e. slots into `ptrsub` (0-based,
// 0..Nnomissing-1), ascending. `n_t[t]` is |S_t| = Nnomissing - excl[t].size().
struct ExclusionSets {
  int P = 0;
  const std::vector<int>* excl = nullptr;   // P lists
  const int*              n_t  = nullptr;   // P counts
};

// Per (marker, trait) deduction — what has to come off the union's counters to
// get the trait's own counters. Integers, so the subtraction is exact.
struct TraitTally {
  int alleleRawExcl  = 0;   // Σ g over excluded, non-missing rows
  int numMissingExcl = 0;   // # excluded rows that are missing
};

struct DecodeAux {
  // in — when non-null, per-trait tallies (and the union keep flag) are built.
  const ExclusionSets* excl = nullptr;
  // out — length excl->P, marker-major slice for this marker. Optional.
  TraitTally*          tally = nullptr;
  // out — UNION-LOCAL row indices where this marker is missing, ascending.
  //       Cleared then filled on every call; empty when numMissing == 0.
  std::vector<int>*    missing_rows = nullptr;
  // out — §2's union-pack rule: true iff at least one trait's passQC_t is true.
  //       Requires `excl`. Optional.
  bool*                keep_any_trait = nullptr;
  // in — when true, `packed_out` is written iff the §2 keep rule fires, instead
  //      of the default (union's own passQC || passVR). Requires `excl`.
  //      Note this decouples "was packed" from stats.passQC: with pack_if_keep
  //      the caller must key storage off keep_any_trait, not passQC.
  bool                 pack_if_keep = false;
};

// The per-marker arithmetic that turns raw counters into stats + QC + VR.
// This is decode_marker's own block (lines 93-130 of the pre-scheme-C file),
// lifted verbatim so that the scheme-C per-trait path and the decoder share one
// copy of it. SCHEME_C_DESIGN.md §4 requires the two to agree BIT FOR BIT, and
// sharing the code is the only way to keep that true under later edits.
//
//   alleleRaw   — Σ bufferGeno over NON-MISSING samples of the set
//   numMissing  — # missing samples of the set
//   n_samples   — |set|
// Writes every field of `stats` (including alleleRaw/fillin) and `passVR`.
void marker_stats_from_counts(int alleleRaw, int numMissing,
                              std::size_t n_samples,
                              float min_maf, float max_miss,
                              const VarRatioRule& vr, bool vr_drawn,
                              MarkerStats& stats, bool& passVR);

// Lookup table built once per process: one BED byte → 4 bufferGeno values.
// Exposed so the pipeline test can verify its contents byte-wise against
// SAIGE_step1_fast.cpp's version. Values are {0,1,2,3} as documented above.
using BedLut = int[256][4];
void build_bed_lookup(BedLut& out);

// Decode one marker's raw BED bytes, compute stats, re-pack into `packed_out`.
//
//   raw         — ⌈N/4⌉ bytes, exactly as BedReaderPool::read_marker() returns
//   N           — number of FAM samples (== raw covers ceil(N/4) bytes)
//   ptrsub      — length Nnomissing, 1-based FAM index per GRM slot
//                 (i.e., GRM sample k takes its geno from FAM row ptrsub[k]-1)
//   Nnomissing  — number of phenotyped samples == length of ptrsub == nbyte_new*4 (padded)
//   min_maf     — QC threshold on min(altFreq, 1-altFreq)
//   max_miss    — QC threshold on missingRate
//   lut         — prebuilt BED lookup table
//
// Outputs (always populated):
//   stats.altFreq, alleleCount, missingRate, numMissing, mac, passQC
//   packed_out  — ⌈Nnomissing/4⌉ bytes, valid only when stats.passQC is true.
//                 Bytes beyond the last-written position are zero-initialized.
void decode_marker(const unsigned char* raw,
                   std::size_t N,
                   const int* ptrsub, std::size_t Nnomissing,
                   float min_maf, float max_miss,
                   const BedLut& lut,
                   MarkerStats& stats,
                   unsigned char* packed_out);

// Variance-ratio aware overload. Identical to the above except that it also
// applies `vr` (see VarRatioRule):
//   in      vr_drawn — true when this marker index is in the random VR draw
//                      (genoClass::g_randMarkerIndforVR). Ignored when
//                      vr.enabled is false.
//   out     passVR   — marker claimed for the VR pool
//   stats.passQC     — forced false when passVR is true
//   packed_out       — written whenever (stats.passQC || passVR), i.e. a VR
//                      marker that fails the GRM MAF filter still gets packed
void decode_marker(const unsigned char* raw,
                   std::size_t N,
                   const int* ptrsub, std::size_t Nnomissing,
                   float min_maf, float max_miss,
                   const BedLut& lut,
                   const VarRatioRule& vr, bool vr_drawn,
                   MarkerStats& stats,
                   bool& passVR,
                   unsigned char* packed_out);

// Scheme-C overload. `aux == nullptr` is byte-for-byte the overload above.
void decode_marker(const unsigned char* raw,
                   std::size_t N,
                   const int* ptrsub, std::size_t Nnomissing,
                   float min_maf, float max_miss,
                   const BedLut& lut,
                   const VarRatioRule& vr, bool vr_drawn,
                   MarkerStats& stats,
                   bool& passVR,
                   unsigned char* packed_out,
                   const DecodeAux* aux);

} // namespace saige
