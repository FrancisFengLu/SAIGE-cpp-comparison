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

namespace saige {

struct MarkerStats {
  float altFreq     = 0.0f;   // post-fill allele frequency over Nnomissing samples
  float missingRate = 0.0f;   // numMissing / Nnomissing (pre-fill)
  int   alleleCount = 0;      // post-fill sum of bufferGeno (0/1/2) over Nnomissing
  int   mac         = 0;      // min(alleleCount, 2*Nnomissing - alleleCount)
  int   numMissing  = 0;      // raw missings before fill-in
  bool  passQC      = false;  // altFreq-based MAF ≥ min_maf  AND  missingRate ≤ max_miss
};

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

} // namespace saige
