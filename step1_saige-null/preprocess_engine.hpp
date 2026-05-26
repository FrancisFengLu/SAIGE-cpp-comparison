#pragma once
#include "saige_null.hpp"
#include <string>

namespace saige {

struct PreOut {
  Design      design;     // numeric, IID-aligned to genotype order
  LocoRanges  chr;        // LOCO start/end per chr (0-based, inclusive)
  FitNullConfig cfg;      // possibly sanitized
};

void apply_row_subset(Design& d, const std::vector<size_t>& keep);

class PreprocessEngine {
public:
  PreprocessEngine(const Paths& paths, const FitNullConfig& cfg);

  // Run preprocessing. Assumes `design_in` is already numeric (pre-encoded).
  // The engine will:
  //  - intersect/reorder rows by genotype (FAM or sparse GRM IDs)
  //  - apply inverse normalization for quantitative if requested
  //  - bin survival event time if requested
  //  - compute LOCO chr ranges from BIM (if LOCO enabled and PLINK available)
  PreOut run(const Design& design_in);
  static void apply_sex_filter_if_requested(Design& design, const FitNullConfig& cfg);

private:
  const Paths paths_;
  FitNullConfig cfg_;

  // helpers
  std::vector<std::string> read_fam_iids_() const;
  std::vector<std::string> read_sparse_iids_() const;
  void intersect_and_reorder_(Design& d, const std::vector<std::string>& geno_iids) const;
  void inv_normalize_quant_(Design& d) const;
  void bin_event_time_if_needed_(Design& d) const;
  LocoRanges compute_chr_ranges_from_bim_() const;

  // utils
  static std::vector<std::string> split_ws_(const std::string& s);
  static bool is_number_(const std::string& s);
};

} // namespace saige
