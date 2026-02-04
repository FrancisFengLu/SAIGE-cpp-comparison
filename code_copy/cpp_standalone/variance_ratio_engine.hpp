#pragma once
#include "saige_null.hpp"
#include <string>

namespace saige {

class VarianceRatioEngine {
public:
  VarianceRatioEngine(const Paths& paths, const FitNullConfig& cfg, const LocoRanges& chr);

  // Run VR extraction; returns an updated FitNullResult with vr_path and marker-out path filled.
  // Requires that a VR runner hook be registered (see register_vr_runner below).
  FitNullResult run(const FitNullResult& in, const Design& design);

  // Optional: validate the written variance-ratio file
  static void validate_vr_file(const std::string& vr_path, bool expect_categorical);

private:
  const Paths paths_;
  const FitNullConfig cfg_;
  const LocoRanges chr_;
  bool should_skip_writing_vr_(const std::string& outfile,
                                const FitNullConfig& cfg) const;
};

// ---- Hook registration for your native VR implementation ----
using VRRunnerFn = void (*)(const Paths&,
                            const FitNullConfig&,
                            const LocoRanges&,
                            const FitNullResult&,   // fitted theta/alpha/offset
                            const Design&,          // design (n,p,iid)
                            std::string& /*out_vr_path*/,
                            std::string& /*out_marker_results_path*/);

void register_vr_runner(VRRunnerFn fn);

} // namespace saige
