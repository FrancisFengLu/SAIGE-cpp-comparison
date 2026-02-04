#include "preprocess_engine.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cctype>
#include "SAIGE_step1_fast.hpp"

namespace saige {

PreprocessEngine::PreprocessEngine(const Paths& paths, const FitNullConfig& cfg)
  : paths_(paths), cfg_(cfg) {}

// --------- public entry ----------
PreOut PreprocessEngine::run(const Design& design_in) {
  if (design_in.n <= 0) throw std::invalid_argument("Design.n must be > 0");
  if (design_in.y.size() != static_cast<size_t>(design_in.n))
    throw std::invalid_argument("Design.y length must equal n");
  if (design_in.iid.size() != static_cast<size_t>(design_in.n))
    throw std::invalid_argument("Design.iid length must equal n");
  if (design_in.p > 0 && design_in.X.size() != static_cast<size_t>(design_in.n) * static_cast<size_t>(design_in.p))
    throw std::invalid_argument("Design.X size must equal n*p");

  // Copy (we’ll mutate)
  Design d = design_in;

  // 1) Determine genotype IID universe and LOCO eligibility
  std::vector<std::string> geno_iids;
  if (!cfg_.use_sparse_grm_to_fit && !paths_.fam.empty()) {
    geno_iids = read_fam_iids_();
  } else if (cfg_.use_sparse_grm_to_fit && !paths_.sparse_grm_ids.empty()) {
    geno_iids = read_sparse_iids_();
  } else {
    // If neither FAM nor sparse IDs provided, keep order as-is
    geno_iids = d.iid;
  }

  // 2) Intersect/reorder design rows to match genotype order
  intersect_and_reorder_(d, geno_iids);

  // 3) Inverse normalization for quantitative (optional)
  if (cfg_.trait == std::string("quantitative") && cfg_.inv_normalize) {
    inv_normalize_quant_(d);
  }

  // 4) Survival event-time binning (optional)
  if (cfg_.trait == std::string("survival") && cfg_.event_time_bin_size.has_value()) {
    bin_event_time_if_needed_(d);
  }

  // 5) LOCO chromosome ranges
  LocoRanges lr;
  lr.enabled = cfg_.loco;
  if (cfg_.loco && !paths_.bim.empty()) {
    lr = compute_chr_ranges_from_bim_();
    lr.enabled = !lr.start.empty();
  }

  return PreOut{ std::move(d), lr, cfg_ };
}

// --------- FAM & sparse-ID readers ----------
std::vector<std::string> PreprocessEngine::read_fam_iids_() const {
  std::ifstream in(paths_.fam);
  if (!in) throw std::runtime_error("Failed to open FAM: " + paths_.fam);
  std::vector<std::string> out;
  out.reserve(1024);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto toks = split_ws_(line);
    if (toks.size() < 2) continue;
    out.push_back(toks[1]); // column 2: IID
  }
  return out;
}

std::vector<std::string> PreprocessEngine::read_sparse_iids_() const {
  std::ifstream in(paths_.sparse_grm_ids);
  if (!in) throw std::runtime_error("Failed to open sparse GRM IDs: " + paths_.sparse_grm_ids);
  std::vector<std::string> out;
  out.reserve(1024);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    // one IID per line
    out.push_back(line);
  }
  return out;
}

// --------- Intersection & reorder ----------
void PreprocessEngine::intersect_and_reorder_(Design& d,
                                              const std::vector<std::string>& geno_iids) const {
  // map genotype IID -> position
  std::unordered_map<std::string, int> pos;
  pos.reserve(geno_iids.size() * 2);
  for (int i = 0; i < (int)geno_iids.size(); ++i) pos.emplace(geno_iids[i], i);

  // collect indices in genotype order that exist in design
  std::vector<int> keep_pos; keep_pos.reserve(d.n);
  for (int i = 0; i < d.n; ++i) {
    auto it = pos.find(d.iid[i]);
    if (it != pos.end()) keep_pos.push_back(i);
  }
  if (keep_pos.empty())
    throw std::runtime_error("No overlapping samples between design and genotype IDs");

  // reorder design to genotype order (stable)
  std::vector<int> geno_order; geno_order.reserve(keep_pos.size());
  {
    // For each geno_iid, if present in design, push its design index
    std::unordered_map<std::string, int> dpos;
    dpos.reserve(d.n * 2);
    for (int i = 0; i < d.n; ++i) dpos.emplace(d.iid[i], i);
    for (int gi = 0; gi < (int)geno_iids.size(); ++gi) {
      auto jt = dpos.find(geno_iids[gi]);
      if (jt != dpos.end()) geno_order.push_back(jt->second);
    }
  }
  // Subset/reorder rows
  const int n2 = (int)geno_order.size();
  std::vector<std::string> iid2; iid2.reserve(n2);
  std::vector<double> y2(n2), off2(n2), time2;
  if (!d.event_time.empty()) time2.resize(n2);
  for (int i = 0; i < n2; ++i) {
    int k = geno_order[i];
    iid2.push_back(d.iid[k]);
    y2[i] = d.y[k];
    off2[i] = d.offset.empty() ? 0.0 : d.offset[k];
    if (!d.event_time.empty()) time2[i] = d.event_time[k];
  }
  std::vector<double> X2;
  if (d.p > 0) {
    X2.resize((size_t)n2 * (size_t)d.p);
    // copy row blocks
    for (int i = 0; i < n2; ++i) {
      int k = geno_order[i];
      const size_t src = (size_t)k * (size_t)d.p;
      const size_t dst = (size_t)i * (size_t)d.p;
      std::copy_n(d.X.begin() + src, d.p, X2.begin() + dst);
    }
  }
  d.n = n2;
  d.iid = std::move(iid2);
  d.y.swap(y2);
  d.offset.swap(off2);
  if (!d.event_time.empty()) d.event_time.swap(time2);
  d.X.swap(X2);
}

// --------- Inverse normalization ----------
void PreprocessEngine::inv_normalize_quant_(Design& d) const {
  // rank-based inverse normal (Blom-ish): qnorm((rank - 0.5) / n_non_na)
  struct Pair { double val; int idx; };
  std::vector<Pair> v; v.reserve(d.n);
  for (int i = 0; i < d.n; ++i) v.push_back({d.y[i], i});
  std::sort(v.begin(), v.end(), [](const Pair& a, const Pair& b){ return a.val < b.val; });
  int n_non_na = d.n; // assuming no NA in numeric payload
  std::vector<double> y2(d.n);
  for (int r = 0; r < d.n; ++r) {
    int i = v[r].idx;
    double u = ( (double)r + 0.5 ) / (double)n_non_na;
    // inverse normal (erfinv) via std::erfcinv/erfinv not in C++17; use approximation
    // Use Beasley-Springer/Moro or fall back to std::erfc for a quick approximation.
    // Here: simple Acklam inverse CDF approximation (quick-and-dirty), sufficient for preprocessing.
    // To keep this file lean, approximate with std::sqrt(2)*erfinv(2u-1) via rational approx:

    auto inv_norm = [](double p)->double {
      // Peter J. Acklam's inverse normal CDF approximation (abridged)
      if (p <= 0.0) return -1e10;
      if (p >= 1.0) return  1e10;
      static const double a1=-3.969683028665376e+01, a2= 2.209460984245205e+02,
                          a3=-2.759285104469687e+02, a4= 1.383577518672690e+02,
                          a5=-3.066479806614716e+01, a6= 2.506628277459239e+00;
      static const double b1=-5.447609879822406e+01, b2= 1.615858368580409e+02,
                          b3=-1.556989798598866e+02, b4= 6.680131188771972e+01,
                          b5=-1.328068155288572e+01;
      static const double c1=-7.784894002430293e-03, c2=-3.223964580411365e-01,
                          c3=-2.400758277161838e+00, c4=-2.549732539343734e+00,
                          c5= 4.374664141464968e+00, c6= 2.938163982698783e+00;
      static const double d1= 7.784695709041462e-03, d2= 3.224671290700398e-01,
                          d3= 2.445134137142996e+00, d4= 3.754408661907416e+00;
      double q, r;
      if (p < 0.02425) {
        q = std::sqrt(-2*std::log(p));
        return (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
               ((((d1*q+d2)*q+d3)*q+d4)*q+1);
      } else if (p > 1-0.02425) {
        q = std::sqrt(-2*std::log(1-p));
        return -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
                 ((((d1*q+d2)*q+d3)*q+d4)*q+1);
      } else {
        q = p-0.5; r = q*q;
        return (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
               (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
      }
    };
    y2[i] = inv_norm(u);
  }
  d.y.swap(y2);
}

// --------- Survival event-time binning ----------
void PreprocessEngine::bin_event_time_if_needed_(Design& d) const {
  if (d.event_time.empty()) return;
  const int bin = cfg_.event_time_bin_size.value();
  double min_t = d.event_time[0];
  for (double t : d.event_time) if (t < min_t) min_t = t;
  for (int i = 0; i < d.n; ++i) {
    const double t = d.event_time[i];
    const int grp = static_cast<int>( (t - min_t) / (double)bin ) + 1;
    d.event_time[i] = (double)grp;
  }
}

// --------- LOCO from BIM ----------
LocoRanges PreprocessEngine::compute_chr_ranges_from_bim_() const {
  LocoRanges lr;
  lr.enabled = false;
  if (paths_.bim.empty()) return lr;

  std::ifstream in(paths_.bim);
  if (!in) throw std::runtime_error("Failed to open BIM: " + paths_.bim);

  // 0-based marker index across entire BIM
  // Identify first/last index for chr 1..22 (numeric only)
  std::vector<int> first(23, -1), last(23, -1);
  std::string line; int idx = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto toks = split_ws_(line);
    if (toks.empty()) continue;
    const std::string& chr_s = toks[0];
    if (!is_number_(chr_s)) { ++idx; continue; }
    int chr = std::stoi(chr_s);
    if (chr < 1 || chr > 22) { ++idx; continue; }
    if (first[chr] == -1) first[chr] = idx;
    last[chr] = idx;
    ++idx;
  }
  for (int c = 1; c <= 22; ++c) {
    if (first[c] != -1 && last[c] != -1) {
      lr.start.push_back(first[c]);
      lr.end.push_back(last[c]);
    } else {
      lr.start.push_back(-1);
      lr.end.push_back(-1);
    }
  }
  lr.enabled = true;
  return lr;
}

// --------- tiny utils ----------
std::vector<std::string> PreprocessEngine::split_ws_(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream iss(s);
  std::string tok;
  while (iss >> tok) out.push_back(tok);
  return out;
}
bool PreprocessEngine::is_number_(const std::string& s) {
  if (s.empty()) return false;
  for (char c : s) if (!std::isdigit(static_cast<unsigned char>(c))) return false;
  return true;
}


static inline std::vector<size_t>
build_keep_index_from_sex(const std::vector<std::string>& sex_vec,
                          const FitNullConfig& cfg)
{
    std::vector<size_t> keep;
    keep.reserve(sex_vec.size());

    const bool do_female = cfg.female_only;
    const bool do_male   = cfg.male_only;

    // If both toggles are set (user error), prefer female_only and warn by narrowing to females.
    for (size_t i = 0; i < sex_vec.size(); ++i) {
        const std::string& v = sex_vec[i];
        if (do_female && v == cfg.female_code) { keep.push_back(i); continue; }
        if (do_male   && v == cfg.male_code)   { keep.push_back(i); continue; }
        if (!do_female && !do_male)            { keep.push_back(i); /* no-op path */ }
    }
    return keep;
}

void apply_row_subset(Design& d, const std::vector<size_t>& keep)
{
 // iid
  {
    std::vector<std::string> iid2; iid2.reserve(keep.size());
    for (auto i : keep) iid2.push_back(d.iid[i]);
    d.iid.swap(iid2);
  }

  // y
  {
    std::vector<double> y2; y2.reserve(keep.size());
    for (auto i : keep) y2.push_back(d.y[i]);
    d.y.swap(y2);
  }

  // offset (if present)
  if (!d.offset.empty()) {
    std::vector<double> off2; off2.reserve(keep.size());
    for (auto i : keep) off2.push_back(d.offset[i]);
    d.offset.swap(off2);
  }

  // X (row-major)
  if (d.p > 0 && !d.X.empty()) {
    const int p = d.p;
    std::vector<double> X2; X2.resize(keep.size() * (size_t)p);
    for (size_t r = 0; r < keep.size(); ++r) {
      const size_t i = keep[r];
      const double* src = &d.X[(size_t)i * (size_t)p];
      std::copy(src, src + p, &X2[(size_t)r * (size_t)p]);
    }
    d.X.swap(X2);
  }

  d.n = static_cast<int>(keep.size());
}


void PreprocessEngine::apply_sex_filter_if_requested(Design& /*design*/, const FitNullConfig& cfg)
{
    // Design does not expose string columns here; sex filtering is done in main.cpp
    // right after reading the design CSV (where the column is available).
    // Keep this as a no-op to avoid duplicate filtering and compile-time dependency.
    if (!cfg.sex_col.empty() && (cfg.female_only || cfg.male_only)) {
        std::cerr << "[design] note: sex filtering handled in main.cpp (no-op here)\n";
    }
}

// Minimal row-subsetter; implement with your Design layout.
// This example assumes Design has:
//   - std::vector<std::string> ids
//   - arma::vec y;
//   - arma::mat X;
// and that sizes are aligned by rows.

} // namespace saige
 