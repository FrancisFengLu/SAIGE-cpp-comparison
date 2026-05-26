#pragma once
#include "saige_null.hpp"
#include <string>

namespace saige {

// Real variance ratio computation matching R's extractVarianceRatio()
void compute_variance_ratio(const Paths& paths,
                            const FitNullConfig& cfg,
                            const LocoRanges& chr,
                            const FitNullResult& fit,
                            const Design& design,
                            std::string& out_vr_path,
                            std::string& out_marker_results_path);

} // namespace saige
