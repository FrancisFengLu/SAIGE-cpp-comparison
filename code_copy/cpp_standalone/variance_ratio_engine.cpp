#include "variance_ratio_engine.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cctype>
#include <filesystem>
#include "SAIGE_step1_fast.hpp"
namespace saige {

static inline bool file_exists_(const std::string& p) {
   return std::filesystem::exists(std::filesystem::path(p));
}

bool VarianceRatioEngine::should_skip_writing_vr_(const std::string& outfile,
                                                  const FitNullConfig& cfg) const
{
    // Expect a boolean like cfg.overwrite_vr (R: IsOverwriteVarianceRatioFile)
    if (!file_exists_(outfile)) return false;
    if (cfg.overwrite_vr)       return false;

    std::cerr << "[VR] File exists and overwrite disabled: " << outfile
              << " — skipping write.\n";
    return true;
}

static VRRunnerFn g_vr_runner = nullptr;
void register_vr_runner(VRRunnerFn fn) { g_vr_runner = fn; }

VarianceRatioEngine::VarianceRatioEngine(const Paths& paths,
                                         const FitNullConfig& cfg,
                                         const LocoRanges& chr)
  : paths_(paths), cfg_(cfg), chr_(chr) {}

FitNullResult VarianceRatioEngine::run(const FitNullResult& in, const Design& design) {
  if (!g_vr_runner)
    throw std::runtime_error("Variance ratio runner not registered. Call register_vr_runner().");

  // ---- Overwrite guard (R: IsOverwriteVarianceRatioFile) ----
  // Convention: if the runner doesn't override, it will write to:
  //   <paths_.out_prefix_vr>.varianceRatio.txt   (optionally with per-chr suffix inside runner)
  auto default_vr_path = paths_.out_prefix_vr + ".varianceRatio.txt";
  auto file_exists = [](const std::string& p) {
    return std::filesystem::exists(std::filesystem::path(p));
  };

  // If your runner writes per-chromosome files (e.g., add ".chr" + chr_),
  // you can derive that here instead of default_vr_path — or keep the check
  // inside the runner. This pre-check covers the common single-file case.
  if (!cfg_.overwrite_vr && file_exists(default_vr_path)) {
    std::cerr << "[VR] File exists and overwrite_vr=false: " << default_vr_path
              << " — skipping variance ratio computation.\n";
    FitNullResult out = in;
    out.vr_path = default_vr_path;
    // markers_out_path left as-is; if you also persist SPA markers, consider
    // guarding them similarly (e.g., <out_prefix_vr>.#markers.SPAOut.txt).
    return out;
  }

  // ---- Run the registered VR implementation ----
  FitNullResult out = in; // copy, we’ll just add paths
  std::string vr_path, marker_out_path;
  g_vr_runner(paths_, cfg_, chr_, in, design, vr_path, marker_out_path);
  out.vr_path = vr_path.empty() ? default_vr_path : vr_path;
  out.markers_out_path = marker_out_path;
  return out;
}

// ---- simple VR file validator (optional) ----
static std::vector<std::string> split_csv_(const std::string& s) {
  // naive split on commas; your writer likely uses tabs—adapt if needed
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == ',') { out.push_back(cur); cur.clear(); }
    else          { cur.push_back(c); }
  }
  out.push_back(cur);
  return out;
}

void VarianceRatioEngine::validate_vr_file(const std::string& vr_path, bool expect_categorical) {
  std::ifstream in(vr_path);
  if (!in) throw std::runtime_error("Cannot open VR file: " + vr_path);

  std::string header;
  if (!std::getline(in, header)) throw std::runtime_error("Empty VR file: " + vr_path);
  auto cols = split_csv_(header);
  // Accept common variants: Category,MACBin,VarianceRatio,NMarkers,...
  auto has = [&](const std::string& name){
    for (auto& c : cols) if (c == name) return true;
    return false;
  };
  if (!has("VarianceRatio"))
    throw std::runtime_error("VR file missing VarianceRatio column");

  // Optional: ensure positivity and basic schema
  std::string line;
  int n_rows = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto toks = split_csv_(line);
    // find VarianceRatio column index
    int vr_idx = -1;
    for (int i = 0; i < (int)cols.size(); ++i) if (cols[i] == "VarianceRatio") { vr_idx = i; break; }
    if (vr_idx < 0 || vr_idx >= (int)toks.size())
      throw std::runtime_error("Malformed VR row (VarianceRatio missing): " + line);
    double vr = std::stod(toks[vr_idx]);
    if (!(vr > 0.0)) throw std::runtime_error("Non-positive VarianceRatio encountered");
    ++n_rows;
  }
  if (n_rows == 0) throw std::runtime_error("VR file has no data rows");

  // If categorical expected, ensure at least 2 rows
  if (expect_categorical && n_rows < 2)
    throw std::runtime_error("Expected multiple MAC categories in VR file");
}

} // namespace saige
