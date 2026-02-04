// main.cpp
// ------------------------------------------------------------------
// CLI orchestration for SAIGE null fitting + optional Variance Ratio,
// with categorical covariates, covariate-offset, inverse-normalization,
// LOCO ranges from BIM, sparse-GRM reuse/build, IID whitelist.
//
// Build:
//   - yaml-cpp, cxxopts, Armadillo
//   - Assumes SAIGE kernels expose sparse-GRM hooks (see calls below)
//   - Optionally SAIGE_step1_fast.hpp for genoClass (kept optional)
//
// Usage:
//   saige-null -c config.yaml -d design.tsv
//   saige-null -c config.yaml -o fit.nthreads=32 -o paths.out_prefix=out/run2
// ------------------------------------------------------------------

#include "saige_null.hpp"     // Paths, FitNullConfig, Design, FitNullResult, register_default_solvers()
#include "glmm.hpp"
#include "SAIGE_step1_fast.hpp"   // (optional) genoClass decl — comment out if not available
#include "preprocess_engine.hpp"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <armadillo>
#include <iomanip>
#include <cmath>

namespace fs = std::filesystem;
using saige::FitNullConfig;
using saige::Paths;
using saige::Design;
using saige::FitNullResult;

// ------------------ small helpers ------------------
static inline void ensure_parent_dir(const std::string& path) {
  fs::path p(path);
  auto dir = p.parent_path();
  if (!dir.empty()) fs::create_directories(dir);
}
static bool ieq(const std::string& a, const std::string& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (std::tolower(static_cast<unsigned char>(a[i])) != std::tolower(static_cast<unsigned char>(b[i])))
      return false;
  return true;
}
static std::string trim(const std::string& s) {
  size_t i = 0, j = s.size();
  while (i < j && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
  return s.substr(i, j - i);
}
static std::vector<std::string> split_simple(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == delim) { out.push_back(cur); cur.clear(); }
    else            { cur.push_back(c); }
  }
  out.push_back(cur);
  return out;
}
static YAML::Node parse_scalar_to_yaml(const std::string& v) {
  YAML::Node out;
  if (ieq(v, "true"))  { out = true;  return out; }
  if (ieq(v, "false")) { out = false; return out; }
  if (ieq(v, "null"))  { out = YAML::Node(); return out; }
  char* end = nullptr;
  long val_i = std::strtol(v.c_str(), &end, 10);
  if (end && *end == '\0') { out = static_cast<int>(val_i); return out; }
  end = nullptr;
  double val_d = std::strtod(v.c_str(), &end);
  if (end && *end == '\0') { out = val_d; return out; }
  out = v;
  return out;
}
static void yaml_set_dotted(YAML::Node& root,
                            const std::string& dotted,
                            const YAML::Node& value) {
  // Split on '.' but allow '\.' to mean a literal dot in a key
  auto split_dotted = [](const std::string& s) {
    std::vector<std::string> parts;
    std::string cur; cur.reserve(s.size());
    bool esc = false;
    for (char c : s) {
      if (esc) { cur.push_back(c); esc = false; continue; }
      if (c == '\\') { esc = true; continue; }
      if (c == '.') { parts.push_back(cur); cur.clear(); }
      else          { cur.push_back(c); }
    }
    parts.push_back(cur);
    return parts;
  };

  const auto parts = split_dotted(dotted);
  if (parts.empty()) return;

  // Disallow clobbering an entire map with a non-map value (e.g., -o paths=foo)
  if (parts.size() == 1 && !value.IsMap()) {
    // If the target currently is a map (or expected to be a map), refuse
    YAML::Node existing = root[parts[0]];
    if (existing && existing.IsMap()) {
      throw std::runtime_error(
        "Refusing to replace map '" + parts[0] +
        "' with a scalar via override '" + dotted +
        "'. Use '-o " + parts[0] + ".someKey=...'");
    }
  }

  // Walk/create intermediate maps
  YAML::Node node = root;
  for (size_t i = 0; i + 1 < parts.size(); ++i) {
    const std::string& k = parts[i];
    YAML::Node next = node[k];
    if (!next || next.IsNull()) {
      node[k] = YAML::Node(YAML::NodeType::Map);
      next = node[k];
    } else if (!next.IsMap()) {
      // Don't smash existing non-map nodes when asked to go deeper
      throw std::runtime_error(
        "Override '" + dotted + "' expects '" + k +
        "' to be a map, but it is " +
        (next.IsScalar() ? "a scalar" :
         next.IsSequence()? "a sequence" : "unknown") + ".");
    }
    node = next;
  }

  // Set the leaf (this updates only the designated item)
  node[parts.back()] = value;
}

static std::vector<std::string>
read_column_from_csv(const std::string& csv_path, const std::string& col_name) {
  std::ifstream f(csv_path);
  if (!f) throw std::runtime_error("Failed to open design CSV: " + csv_path);

  std::string line;
  if (!std::getline(f, line)) throw std::runtime_error("Empty design CSV");

  // header -> find column index
  std::vector<std::string> header;
  { std::istringstream iss(line); std::string tok;
    while (std::getline(iss, tok, (line.find('\t') != std::string::npos ? '\t' : ','))) header.push_back(tok);
  }
  int idx = -1;
  for (int i=0;i<(int)header.size();++i) if (header[i] == col_name) { idx = i; break; }
  if (idx < 0) throw std::runtime_error("Column not found in design CSV: " + col_name);

  // read the values
  std::vector<std::string> out; out.reserve(1024);
  char sep = (line.find('\t') != std::string::npos ? '\t' : ',');
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::vector<std::string> cells; cells.reserve(header.size());
    std::string cell; std::istringstream iss2(line);
    while (std::getline(iss2, cell, sep)) cells.push_back(cell);
    if ((int)cells.size() <= idx) continue;
    out.push_back(cells[idx]);
  }
  return out;
}

static std::vector<size_t>
build_keep_index_from_sex(const std::vector<std::string>& sex_vec, const saige::FitNullConfig& cfg) {
  std::vector<size_t> keep; keep.reserve(sex_vec.size());
  const bool do_female = cfg.female_only;
  const bool do_male   = cfg.male_only;

  for (size_t i=0; i<sex_vec.size(); ++i) {
    const auto& v = sex_vec[i];
    if (do_female && v == cfg.female_code) { keep.push_back(i); continue; }
    if (do_male   && v == cfg.male_code)   { keep.push_back(i); continue; }
    if (!do_female && !do_male) keep.push_back(i); // nothing requested -> keep all
  }
  return keep;
}

// ------------------ YAML loaders ------------------
static FitNullConfig load_cfg(const YAML::Node& y) {
  FitNullConfig c;
  const auto f = y["fit"];
  if (!f) return c;

  auto get = [&](const char* k){ return f[k]; };
  if (get("trait")) c.trait = get("trait").as<std::string>();
  if (get("loco")) c.loco = get("loco").as<bool>();
  if (get("lowmem_loco")) c.lowmem_loco = get("lowmem_loco").as<bool>();
  if (get("use_sparse_grm_to_fit")) c.use_sparse_grm_to_fit = get("use_sparse_grm_to_fit").as<bool>();
  if (get("use_sparse_grm_for_vr")) c.use_sparse_grm_for_vr = get("use_sparse_grm_for_vr").as<bool>();
  if (get("covariate_qr")) c.covariate_qr = get("covariate_qr").as<bool>();
  if (get("covariate_offset")) c.covariate_offset = get("covariate_offset").as<bool>();
  if (get("inv_normalize")) c.inv_normalize = get("inv_normalize").as<bool>();
  if (get("include_nonauto_for_vr")) c.include_nonauto_for_vr = get("include_nonauto_for_vr").as<bool>();

  if (get("tol")) c.tol = get("tol").as<double>();
  if (get("maxiter")) c.maxiter = get("maxiter").as<int>();
  if (get("tolPCG")) c.tolPCG = get("tolPCG").as<double>();
  if (get("maxiterPCG")) c.maxiterPCG = get("maxiterPCG").as<int>();
  if (get("nrun")) c.nrun = get("nrun").as<int>();
  if (get("nthreads")) c.nthreads = get("nthreads").as<int>();
  if (get("traceCVcutoff")) c.traceCVcutoff = get("traceCVcutoff").as<double>();
  if (get("ratio_cv_cutoff")) c.ratio_cv_cutoff = get("ratio_cv_cutoff").as<double>();
  if (get("min_maf_grm")) c.min_maf_grm = get("min_maf_grm").as<double>();
  if (get("max_miss_grm")) c.max_miss_grm = get("max_miss_grm").as<double>();
  if (get("num_markers_for_vr")) c.num_markers_for_vr = get("num_markers_for_vr").as<int>();
  if (get("event_time_bin_size") && !get("event_time_bin_size").IsNull())
    c.event_time_bin_size = get("event_time_bin_size").as<int>();
  if (get("relatedness_cutoff")) c.relatedness_cutoff = get("relatedness_cutoff").as<double>();
  if (get("make_sparse_grm_only")) c.make_sparse_grm_only = get("make_sparse_grm_only").as<bool>();
  if (get("memory_chunk_gb")) c.memory_chunk_gb = get("memory_chunk_gb").as<double>();
  if (get("vr_min_mac")) c.vr_min_mac = get("vr_min_mac").as<int>();
  if (get("vr_max_mac")) c.vr_max_mac = get("vr_max_mac").as<int>();
  if (get("diag_one")) c.isDiagofKinSetAsOne = get("diag_one").as<bool>();
  return c;
}

// --- helpers (put near your other helpers) ---
static inline std::string trim_copy(std::string s) {
  auto issp = [](unsigned char c){ return std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](char c){ return !issp((unsigned char)c); }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [&](char c){ return !issp((unsigned char)c); }).base(), s.end());
  // strip optional surrounding quotes
  if (s.size() >= 2 && ((s.front()=='"' && s.back()=='"') || (s.front()=='\'' && s.back()=='\'')))
    s = s.substr(1, s.size()-2);
  return s;
}

// more robust ext-stripper: handles .bed / .bim / .fam and .*.gz (case-insensitive)
static inline std::string strip_plink_ext_if_any(std::string s) {
  auto lower = s;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  auto try_drop = [&](const std::string& ext){
    if (lower.size() >= ext.size() && lower.compare(lower.size()-ext.size(), ext.size(), ext) == 0) {
      s.erase(s.size()-ext.size()); lower.erase(lower.size()-ext.size()); return true;
    }
    return false;
  };
  // drop .gz first if present
  if (try_drop(".gz")) {
    // after dropping .gz, fall through to try base ext
  }
  (void)(try_drop(".bed") | try_drop(".bim") | try_drop(".fam")); // bitwise | to evaluate all
  return s;
}


// If you already have a "rebase" helper, use that instead.
static inline std::string rebase_to_yaml_dir(const std::string& p,
                                             const std::string& yaml_dir) {
  namespace fs = std::filesystem;
  if (p.empty() || yaml_dir.empty()) return p;
  fs::path P(p);
  if (P.is_absolute()) return p;
  return (fs::path(yaml_dir) / P).string();
}

static inline const char* node_type(const YAML::Node& n) {
  if (!n) return "Undefined";
  if (n.IsNull()) return "Null";
  if (n.IsScalar()) return "Scalar";
  if (n.IsSequence()) return "Sequence";
  if (n.IsMap()) return "Map";
  return "Unknown";
}

static inline std::string norm_key(std::string s) {
  auto issp=[](unsigned char c){ return std::isspace(c); };
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](char c){ return !issp((unsigned char)c); }));
  s.erase(std::find_if(s.rbegin(), s.rend(), [&](char c){ return !issp((unsigned char)c); }).base(), s.end());
  std::string out; out.reserve(s.size());
  for (unsigned char c: s) if (c!='_') out.push_back(std::tolower(c));
  return out;
}

// Find a 'paths' map, even if root isn't a map or the key has weird casing/spaces
static YAML::Node find_paths_node(const YAML::Node& root) {
  // Case A: root is a map
  if (root && root.IsMap()) {
    // exact first
    YAML::Node r = root["paths"];
    if (r && r.IsMap()) return r;
    // tolerant scan
    for (auto it : root) {
      std::string k; try { k = it.first.as<std::string>(); } catch (...) { continue; }
      if (norm_key(k) == "paths" && it.second.IsMap()) return it.second;
    }
  }
  // Case B: root is a sequence (multi-doc or list root)
  if (root && root.IsSequence()) {
    for (std::size_t i=0;i<root.size();++i) {
      YAML::Node elem = root[i];
      if (elem && elem.IsMap()) {
        // exact, then tolerant
        YAML::Node r = elem["paths"];
        if (r && r.IsMap()) return r;
        for (auto it : elem) {
          std::string k; try { k = it.first.as<std::string>(); } catch (...) { continue; }
          if (norm_key(k) == "paths" && it.second.IsMap()) return it.second;
        }
      }
    }
  }
  return YAML::Node(); // Undefined
}


// Pass yaml_dir = directory of the loaded YAML file ("" if unknown)
static Paths load_paths_v2(const YAML::Node& y, const std::string& yaml_dir = "") {
  namespace fs = std::filesystem;

  Paths p;

  YAML::Node r = find_paths_node(y);

  if (!r) {
    std::ostringstream oss;
    oss << "Could not find a 'paths' map in the YAML root.\n"
        << "Root type: " << node_type(y) << "\n"
        << "Hint: ensure your config has a top-level 'paths:' map, "
          "and avoid '-o paths=...'; use '-o paths.plinkFile=...'\n";
    throw std::runtime_error(oss.str());
  }  

  auto get = [&](const char* k) -> YAML::Node { return r[k]; };
  auto as_str = [&](const char* k) -> std::string {
    auto n = get(k);
    return n ? n.as<std::string>() : std::string();
  };

  // 1) Read explicit files if present
  p.bed = as_str("bed");
  p.bim = as_str("bim");
  p.fam = as_str("fam");


  // 2) Read plink prefix (support both styles)
  std::string plink_prefix = as_str("plinkFile");
  // debug

  std::cout << "plink_prefix: " << plink_prefix << std::endl; 
  //
  if (plink_prefix.empty()) plink_prefix = as_str("plinkfile");

  // 3) Other paths
  p.sparse_grm     = as_str("sparse_grm");
  p.sparse_grm_ids = as_str("sparse_grm_ids");
  p.out_prefix     = as_str("out_prefix");
  p.out_prefix_vr  = as_str("out_prefix_vr");
  // optional extras in your YAML:
  // p.pheno          = as_str("pheno");          // if Paths has it
  // p.include_sample = as_str("include_sample"); // if Paths has it

  plink_prefix = trim_copy(plink_prefix);

  // If a prefix exists, ALWAYS synthesize the trio (fill empties; warn on overwrites).
  // This avoids any weirdness with empty-string values in YAML.
  if (!plink_prefix.empty()) {
    static bool once=false; if (!once) { std::cerr << "[load_paths] synthesizing from plinkFile/plinkfile\n"; once=true; }

    std::string prefix = strip_plink_ext_if_any(plink_prefix);

    // If caller explicitly set any of the trio non-empty, keep it; otherwise synthesize.
    if (p.bed.empty()) p.bed = prefix + ".bed";
    if (p.bim.empty()) p.bim = prefix + ".bim";
    if (p.fam.empty()) p.fam = prefix + ".fam";
  }

  // 5) Rebase relative paths to YAML directory (so config is portable)
  auto rebase = [&](std::string& s) {
    s = rebase_to_yaml_dir(s, yaml_dir);
  };
  rebase(p.bed);
  rebase(p.bim);
  rebase(p.fam);
  rebase(p.sparse_grm);
  rebase(p.sparse_grm_ids);
  rebase(p.out_prefix);
  rebase(p.out_prefix_vr);
  // rebase(p.pheno);          // if you have it
  // rebase(p.include_sample); // if you have it

  // 6) Default out_prefix_vr to out_prefix if empty
  if (p.out_prefix_vr.empty()) p.out_prefix_vr = p.out_prefix;

  return p;
}

// -------- MatrixMarket (COO) helpers for sparse GRM --------
static void write_matrix_market_coo(const arma::umat& loc,
                                    const arma::vec&  val,
                                    int n,
                                    const std::string& path)
{
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to write " + path);
  out.setf(std::ios::fixed); out << std::setprecision(10);
  out << "%%MatrixMarket matrix coordinate real general\n%\n";
  out << n << " " << n << " " << val.n_elem << "\n";
  for (arma::uword k = 0; k < val.n_elem; ++k) {
    out << (loc(0,k) + 1) << " " << (loc(1,k) + 1) << " " << val(k) << "\n";
  }
}
static void write_id_list(const std::vector<std::string>& ids,
                          const std::string& path)
{
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Failed to write " + path);
  for (const auto& s : ids) out << s << "\n";
}
static void load_matrix_market_coo(const std::string& path,
                                   arma::umat& loc,
                                   arma::vec&  val,
                                   int& n_out)
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Failed to open " + path);
  std::string header_line;
  if (!std::getline(in, header_line)) throw std::runtime_error("Empty MM file: " + path);

  // Check if symmetric (header contains "symmetric")
  bool is_symmetric = (header_line.find("symmetric") != std::string::npos);
  if (is_symmetric) {
    std::cout << "[load_matrix_market] Detected SYMMETRIC MatrixMarket format" << std::endl;
  }

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '%') continue;
    std::istringstream iss(line);
    int nr, nc, nnz;
    if (!(iss >> nr >> nc >> nnz)) throw std::runtime_error("Bad size line in " + path);
    if (nr != nc) throw std::runtime_error("Non-square matrix in " + path);
    n_out = nr;

    // For symmetric matrices, we need to expand off-diagonal entries
    std::vector<std::tuple<arma::uword, arma::uword, double>> entries;
    entries.reserve(is_symmetric ? 2*nnz : nnz);

    int i, j; double v; int k = 0;
    while (in >> i >> j >> v) {
      if (k >= nnz) break;
      arma::uword row = static_cast<arma::uword>(i - 1);
      arma::uword col = static_cast<arma::uword>(j - 1);
      entries.push_back({row, col, v});
      // For symmetric matrices, add transpose entry for off-diagonal
      if (is_symmetric && row != col) {
        entries.push_back({col, row, v});
      }
      ++k;
    }
    if (k != nnz) throw std::runtime_error("Unexpected EOF while reading entries from " + path);

    int actual_nnz = static_cast<int>(entries.size());
    loc.set_size(2, actual_nnz);
    val.set_size(actual_nnz);
    for (int idx = 0; idx < actual_nnz; ++idx) {
      loc(0, idx) = std::get<0>(entries[idx]);
      loc(1, idx) = std::get<1>(entries[idx]);
      val(idx) = std::get<2>(entries[idx]);
    }

    if (is_symmetric) {
      std::cout << "[load_matrix_market] Expanded " << nnz << " entries to " << actual_nnz << " for symmetric matrix" << std::endl;
    }
    break;
  }
}

// ------------------ Probit approx (Acklam) ------------------
static double probit(double p) {
  // clamp
  if (p <= 0.0) return -std::numeric_limits<double>::infinity();
  if (p >= 1.0) return  std::numeric_limits<double>::infinity();
  // coefficients
  const double a1=-3.969683028665376e+01, a2= 2.209460984245205e+02,
               a3=-2.759285104469687e+02, a4= 1.383577518672690e+02,
               a5=-3.066479806614716e+01, a6= 2.506628277459239e+00;
  const double b1=-5.447609879822406e+01, b2= 1.615858368580409e+02,
               b3=-1.556989798598866e+02, b4= 6.680131188771972e+01,
               b5=-1.328068155288572e+01;
  const double c1=-7.784894002430293e-03, c2=-3.223964580411365e-01,
               c3=-2.400758277161838e+00, c4=-2.549732539343734e+00,
               c5= 4.374664141464968e+00, c6= 2.938163982698783e+00;
  const double d1= 7.784695709041462e-03, d2= 3.224671290700398e-01,
               d3= 2.445134137142996e+00, d4= 3.754408661907416e+00;
  const double pl=0.02425, pu=1.0-pl;
  double q, r;
  if (p < pl) {
    q = std::sqrt(-2*std::log(p));
    return (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
           ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  } else if (p > pu) {
    q = std::sqrt(-2*std::log(1-p));
    return -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) /
             ((((d1*q+d2)*q+d3)*q+d4)*q+1);
  } else {
    q = p-0.5; r = q*q;
    return (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q /
           (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1);
  }
}

// ------------------ LOCO: scan BIM for per-chr ranges ------------------
static std::vector<std::pair<size_t,size_t>> scan_bim_chr_ranges(const std::string& bim) {
  std::ifstream in(bim);
  if (!in) throw std::runtime_error("Failed to open BIM: " + bim);
  std::string chr,id,c3,pos,a1,a2;
  struct R { size_t lo=SIZE_MAX, hi=0; bool any=false; };
  std::unordered_map<std::string,R> map;
  size_t idx=0;
  while (in >> chr >> id >> c3 >> pos >> a1 >> a2) {
    auto &r = map[chr];
    r.any=true; r.lo = std::min(r.lo, idx); r.hi = std::max(r.hi, idx);
    ++idx;
  }
  std::vector<std::pair<size_t,size_t>> out; out.reserve(map.size());
  auto push = [&](const std::string& c){
    auto it=map.find(c); if (it!=map.end() && it->second.any) out.emplace_back(it->second.lo, it->second.hi);
  };
  for (int c=1;c<=22;++c) push(std::to_string(c));
  push("X"); push("Y"); push("MT");
  return out;
}

// ------------------ Design helpers ------------------
static bool is_missing(const std::string& s){
  return s.empty() || s=="NA" || s=="NaN" || s=="nan" || s=="NULL";
}
static bool looks_numeric(const std::string& s){
  if (is_missing(s)) return true;
  char* e=nullptr; std::strtod(s.c_str(), &e);
  return e && *e=='\0';
}
struct CategoricalPlan {
  // for each X col: either numeric (levels empty), or list of kept levels (reference dropped)
  std::vector<bool> is_num;
  std::vector<std::vector<std::string>> kept_levels;
  std::vector<std::string> ref_level;
  int out_p=0;
};

// two-pass plan: detect & choose reference (most frequent; tie → lexicographically smallest)
static CategoricalPlan plan_categoricals(const std::vector<std::vector<std::string>>& rows,
                                         const std::vector<int>& x_idx,
                                         const std::vector<std::string>& header,
                                         bool drop_reference=true)
{
  const size_t n = rows.size();
  const size_t p = x_idx.size();
  CategoricalPlan plan;
  plan.is_num.assign(p,true);
  plan.ref_level.assign(p,"");
  plan.kept_levels.resize(p);

  std::vector<std::unordered_map<std::string,size_t>> counts(p);

  for (size_t i=0;i<n;++i){
    const auto& r = rows[i];
    for (size_t j=0;j<p;++j){
      const std::string &s = r[x_idx[j]];
      if (plan.is_num[j] && !looks_numeric(s)) plan.is_num[j]=false;
      if (!plan.is_num[j] && !is_missing(s)) ++counts[j][s];
    }
  }

  plan.out_p = 0;
  for (size_t j=0;j<p;++j){
    if (plan.is_num[j]) { ++plan.out_p; continue; }
    // choose reference
    size_t best_ct=0; std::string best;
    for (auto& kv: counts[j]){
      if (kv.second>best_ct || (kv.second==best_ct && (best.empty() || kv.first<best)))
        { best_ct=kv.second; best=kv.first; }
    }
    plan.ref_level[j]=best;
    // build kept levels, sorted
    std::vector<std::string> lv; lv.reserve(counts[j].size());
    for (auto& kv: counts[j]) lv.push_back(kv.first);
    std::sort(lv.begin(), lv.end());
    for (auto& v: lv){
      if (drop_reference && v==plan.ref_level[j]) continue;
      plan.kept_levels[j].push_back(v);
    }
    plan.out_p += (int)plan.kept_levels[j].size();
  }
  return plan;
}

static void drop_low_count_binaries_in_place(Design& d,
                                             const std::vector<std::string>& xnames,
                                             int minc)
{
  if (minc<=0 || d.p<=0) return;
  arma::mat X(d.n, d.p);
  for (int i=0;i<d.n;++i)
    for (int j=0;j<d.p;++j)
      X(i,j) = d.X[(size_t)i*(size_t)d.p + (size_t)j];

  std::vector<size_t> keep;
  std::vector<std::string> newn;
  for (int j=0;j<d.p;++j){
    arma::vec col = X.col(j);
    arma::uvec fin = arma::find_finite(col);
    double ones = arma::sum(col.elem(fin));
    double zeros = fin.n_elem - ones;
    if (std::min(ones,zeros) >= (double)minc) { keep.push_back(j); newn.push_back(xnames[j]); }
  }
  if ((int)keep.size()==d.p) return;
  arma::mat Xk(d.n, keep.size());
  for (int i=0;i<d.n;++i)
    for (size_t k=0;k<keep.size();++k)
      Xk(i,k) = X(i, keep[k]);

  d.X.assign(Xk.begin(), Xk.end());
  d.p = (int)keep.size();
  // (optional) you can store xnames in Design if you have a slot
}

// ------------------ Design CSV/TSV/space parser + categorical encoding ------------------
// Expected header columns (case-insensitive): IID, y, [offset], [time|event_time|eventTime], X...
static Design load_design_csv(const std::string& path,
                              int min_covariate_count,
                              bool categorical_drop_reference,
                              const std::vector<std::string>& covar_col_names = {})
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Failed to open design file: " + path);

  std::string header;
  if (!std::getline(in, header)) throw std::runtime_error("Empty design file: " + path);

  char delim;
  if (header.find('\t')!=std::string::npos)      delim = '\t';
  else if (header.find(' ')!=std::string::npos)  delim = ' ';
  else                                           delim = ',';

  auto cols = split_simple(header, delim);
  for (auto& c: cols) c = trim(c);

  auto find_col = [&](std::initializer_list<const char*> names) -> int {
    for (int i = 0; i < (int)cols.size(); ++i)
      for (auto n : names) if (ieq(cols[i], n)) return i;
    return -1;
  };

  int idx_iid    = find_col({"IID"});
  int idx_y      = find_col({"y"});
  int idx_offset = find_col({"offset", "covoffset"});
  int idx_time   = find_col({"time","event_time","eventTime"});

  if (idx_iid < 0 || idx_y < 0)
    throw std::runtime_error("Design file must contain IID and y columns");

  // FIX: Only use columns specified in covar_col_names (not all numeric columns!)
  std::vector<int> x_idx;
  std::vector<std::string> x_names;
  if (!covar_col_names.empty()) {
    // Use only the specified covariate columns
    for (const auto& cov_name : covar_col_names) {
      int idx = find_col({cov_name.c_str()});
      if (idx < 0) {
        throw std::runtime_error("Covariate column not found: " + cov_name);
      }
      x_idx.push_back(idx);
      x_names.push_back(cov_name);
    }
    std::cout << "[design] Using specified covariates: ";
    for (const auto& n : x_names) std::cout << n << " ";
    std::cout << std::endl;
  } else {
    // covar_col_names is empty -> NO covariates (this is the key fix!)
    std::cout << "[design] covar_cols is empty -> using NO covariates\n";
  }

  // Read all rows as strings
  std::vector<std::vector<std::string>> rows;
  rows.reserve(1024);
  std::string line;
  while (std::getline(in, line)){
    if (line.empty()) continue;
    auto toks = split_simple(line, delim);
    // pad short rows
    if ((int)toks.size() < (int)cols.size()) toks.resize(cols.size(), "");
    for (auto& t: toks) t = trim(t);
    rows.push_back(std::move(toks));
  }
  const int n = (int)rows.size();

  // Build Design core vectors
  Design d;
  d.n = n;
  d.iid.resize(n);
  d.y.resize(n);
  if (idx_offset>=0) d.offset.assign(n, 0.0);
  if (idx_time>=0)   d.event_time.assign(n, 0.0);

  for (int i=0;i<n;++i){
    d.iid[i] = rows[i][idx_iid];
    d.y[i]   = std::stod(rows[i][idx_y]);
    if (idx_offset>=0 && !rows[i][idx_offset].empty())
      d.offset[i] = std::stod(rows[i][idx_offset]);
    if (idx_time>=0 && !rows[i][idx_time].empty())
      d.event_time[i] = std::stod(rows[i][idx_time]);
  }

  // If no covariates:
  if (x_idx.empty()) { d.p=0; d.X.clear(); return d; }

  // Plan categorical encoding for X columns
  auto plan = plan_categoricals(rows, x_idx, cols, /*drop_reference=*/categorical_drop_reference);

  // Allocate numeric X and fill
  d.p = plan.out_p;
  d.X.assign((size_t)n*(size_t)d.p, 0.0);

  size_t col_out = 0;
  for (size_t j=0;j<x_idx.size();++j){
    if (plan.is_num[j]){
      for (int i=0;i<n;++i){
        const std::string& s = rows[i][x_idx[j]];
        double v = s.empty() ? std::numeric_limits<double>::quiet_NaN() : std::strtod(s.c_str(), nullptr);
        d.X[(size_t)i*(size_t)d.p + col_out] = v;
      }
      ++col_out;
    } else {
      // one-hot for kept levels (reference dropped)
      const auto& kept = plan.kept_levels[j];
      for (const auto& lvl : kept){
        for (int i=0;i<n;++i){
          const std::string& s = rows[i][x_idx[j]];
          double v = (!s.empty() && s==lvl) ? 1.0 : 0.0; // missing -> 0 (acts like reference)
          d.X[(size_t)i*(size_t)d.p + col_out] = v;
        }
        ++col_out;
      }
    }
  }

  // Optional: drop low-count dummies
  if (min_covariate_count > 0) {
    drop_low_count_binaries_in_place(d, x_names, min_covariate_count);
  }
  return d;
}

static bool design_has_intercept(const saige::Design& d) {
  if (d.p <= 0) return false;
  for (int j = 0; j < d.p; ++j) {
    bool all_one = true;
    for (int i = 0; i < d.n; ++i) {
      double v = d.X[(size_t)i*(size_t)d.p + (size_t)j];
      if (!std::isfinite(v) || std::fabs(v - 1.0) > 1e-12) { all_one = false; break; }
    }
    if (all_one) return true;  // found an all-ones column => intercept already present
  }
  return false;
}

static void add_intercept_if_missing(saige::Design& d) {
  if (d.n <= 0) return;
  if (design_has_intercept(d)) return;

  std::vector<double> X2;
  X2.resize((size_t)d.n * (size_t)(d.p + 1));

  for (int i = 0; i < d.n; ++i) {
    // new col 0 is intercept
    X2[(size_t)i*(size_t)(d.p + 1) + 0] = 1.0;
    // shift existing X to the right by 1 column
    for (int j = 0; j < d.p; ++j) {
      X2[(size_t)i*(size_t)(d.p + 1) + (size_t)(j + 1)] =
          d.X[(size_t)i*(size_t)d.p + (size_t)j];
    }
  }
  d.X.swap(X2);
  d.p += 1;
  std::cout << "[design] added intercept column, new p=" << d.p << "\n";
}

// ------------------ FAM IID reader ------------------
static std::vector<std::string> read_fam_iids(const std::string& fam_path) {
  std::ifstream in(fam_path);
  if (!in) throw std::runtime_error("Failed to open FAM: " + fam_path);
  std::vector<std::string> ids;
  std::string fid, iid, p1, p2, sex, pheno;
  ids.reserve(1024);
  while (in >> fid >> iid >> p1 >> p2 >> sex >> pheno) {
    ids.push_back(iid);
  }
  return ids;
}

// Simple slicer if you need to subset Design rows
static void design_take_rows(Design& d, const std::vector<size_t>& keep) {
  const int n2 = (int)keep.size();
  auto take_vec = [&](std::vector<double>& v){
    if (v.empty()) return;
    std::vector<double> out; out.reserve(n2);
    for (auto i: keep) out.push_back(v[i]);
    v.swap(out);
  };
  auto take_str = [&](std::vector<std::string>& v){
    std::vector<std::string> out; out.reserve(n2);
    for (auto i: keep) out.push_back(v[i]);
    v.swap(out);
  };
  // y, offset, time, iid
  take_vec(d.y);
  take_vec(d.offset);
  take_vec(d.event_time);
  take_str(d.iid);
  // X
  if (d.p>0){
    arma::mat X(d.n, d.p);
    for (int i=0;i<d.n;++i)
      for (int j=0;j<d.p;++j)
        X(i,j) = d.X[(size_t)i*(size_t)d.p + (size_t)j];
    arma::mat Xk(n2, d.p);
    for (int r=0;r<n2;++r) Xk.row(r) = X.row(keep[r]);
    d.X.assign(Xk.begin(), Xk.end());
  }
  d.n = n2;
}
 
// ------------------ main ------------------
int main(int argc, char** argv) {
  cxxopts::Options opts("saige-null", "Null GLMM fitting with LOCO/VR (genoClass-integrated)");
  opts.add_options()
    ("c,config",   "YAML config path", cxxopts::value<std::string>())
    ("d,design",   "Override design CSV/TSV path", cxxopts::value<std::string>()->default_value(""))
    ("o,override", "YAML dot-override, e.g., fit.nthreads=32", cxxopts::value<std::vector<std::string>>()->default_value({}))
    ("v,verbose",  "Verbose", cxxopts::value<bool>()->default_value("false"))
    ("h,help",     "Show help");

  auto res = opts.parse(argc, argv);
  if (res.count("help") || !res.count("config")) {
    std::cout << opts.help() << "\n";
    return 0;
  }

  const std::string cfg_path = res["config"].as<std::string>();
  YAML::Node y = YAML::LoadFile(cfg_path);

// debug
  // std::cerr << "[yaml] top keys: ";
  // for (auto it : y) std::cerr << it.first.as<std::string>() << " ";
  // std::cerr << "\n";
  // if (!y["paths"]) {
  //   std::cerr << "[yaml] 'paths' is missing or undefined\n";
  // } else {
  //   std::cerr << "[yaml] paths node type: " 
  //             << (y["paths"].IsMap() ? "Map" : y["paths"].IsScalar() ? "Scalar" : "Other")
  //             << "\n";
  //   std::cerr << "[yaml] paths content:\n" << YAML::Dump(y["paths"]) << "\n";
  // }

  // Apply dot overrides
  if (res.count("override")) {
    for (const auto& kv : res["override"].as<std::vector<std::string>>()) {
      auto pos = kv.find('=');
      if (pos == std::string::npos) {
        std::cerr << "Ignoring override without '=': " << kv << "\n";
        continue;
      }
      auto key = kv.substr(0, pos);
      auto val = kv.substr(pos + 1);
      yaml_set_dotted(y, key, parse_scalar_to_yaml(val));
    }
  }

  // Load config and paths
  FitNullConfig cfg = load_cfg(y);

  std::string yaml_dir;
  try {
    yaml_dir = std::filesystem::path(cfg_path).parent_path().string();
  } catch (...) {
    yaml_dir.clear();
  }

  Paths paths = load_paths_v2(y, yaml_dir);

  if (paths.out_prefix_vr.empty()) paths.out_prefix_vr = paths.out_prefix;

  // Sparse-GRM → force nthreads=1, disable LOCO for null fit (R behavior)
  if (cfg.use_sparse_grm_to_fit) {
    if (cfg.nthreads != 1) std::cerr << "[note] use_sparse_grm_to_fit=true → forcing nthreads=1\n";
    cfg.nthreads = 1;
    cfg.loco = false;
  }

  saige::register_default_solvers();

  // Load design path (CLI > YAML)
  std::string design_csv;
  if (res.count("design") && !res["design"].as<std::string>().empty()) {
    design_csv = res["design"].as<std::string>();
  } else if (y["design"] && y["design"]["csv"]) {
    design_csv = y["design"]["csv"].as<std::string>();
  } else {
    std::cerr << "No design CSV provided. Use -d or set design.csv in YAML.\n";
    return 1;
  }

  // Design knobs
  const int  min_cov_ct = (y["design"] && y["design"]["min_covariate_count"]) ? y["design"]["min_covariate_count"].as<int>() : -1;
  const bool drop_ref   = true; // can expose via YAML if desired

  // FIX: Read covar_cols from config (previously ignored!)
  std::vector<std::string> covar_col_names;
  if (y["design"] && y["design"]["covar_cols"]) {
    const auto& covar_node = y["design"]["covar_cols"];
    if (covar_node.IsSequence()) {
      for (const auto& item : covar_node) {
        covar_col_names.push_back(item.as<std::string>());
      }
    }
  }
  std::cout << "[config] covar_cols=[";
  for (size_t i = 0; i < covar_col_names.size(); ++i) {
    std::cout << covar_col_names[i];
    if (i < covar_col_names.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << (covar_col_names.empty() ? " (NO COVARIATES)" : "") << std::endl;

  // Parse design (with categoricals) - NOW using covar_col_names from config!
  Design design = load_design_csv(design_csv, min_cov_ct, drop_ref, covar_col_names);
  add_intercept_if_missing(design);

  // DEBUG: Verify covariate count
  std::cout << "============================================" << std::endl;
  std::cout << "[DEBUG COVARIATE CHECK]" << std::endl;
  std::cout << "  design.p (number of columns in X) = " << design.p << std::endl;
  std::cout << "  design.n (number of samples) = " << design.n << std::endl;
  std::cout << "  design.X.size() = " << design.X.size() << std::endl;
  std::cout << "  Expected X size (n*p) = " << (design.n * design.p) << std::endl;
  if (design.p == 1) {
    std::cout << "  => INTERCEPT ONLY (no x1, x2 covariates)" << std::endl;
  } else {
    std::cout << "  => WARNING: p > 1, covariates ARE present!" << std::endl;
  }
  std::cout << "============================================" << std::endl;

  if (!cfg.sex_col.empty() && (cfg.female_only || cfg.male_only)) {
    auto sex_vec = read_column_from_csv(design_csv, cfg.sex_col);
    if ((int)sex_vec.size() != design.n) {
      throw std::runtime_error("sex column length != design.n ("
                              + std::to_string(sex_vec.size()) + " vs "
                              + std::to_string(design.n) + ")");
    }
    auto keep = build_keep_index_from_sex(sex_vec, cfg);
    apply_row_subset(design, keep);  // <<— the helper from 2a, make it visible or move it here
    std::cout << "[design] after sex filter: n=" << design.n << "\n";
  }

  // IID whitelist (optional)
  if (y["design"] && y["design"]["whitelist_ids"]) {
    std::ifstream w(y["design"]["whitelist_ids"].as<std::string>());
    if (!w) throw std::runtime_error("Failed to open whitelist: " + y["design"]["whitelist_ids"].as<std::string>());
    std::unordered_set<std::string> ids;
    std::string s; while (std::getline(w, s)) if (!s.empty()) ids.insert(s);
    if (!ids.empty()) {
      std::vector<size_t> keep; keep.reserve(design.n);
      for (size_t i=0;i<design.iid.size();++i) if (ids.count(design.iid[i])) keep.push_back(i);
      design_take_rows(design, keep);
      std::cout << "[design] after whitelist: n=" << design.n << "\n";
    }
  }
  // Sex-stratified filtering (if requested via YAML: design.sex_col + female_only/male_only)
  try {
    saige::PreprocessEngine::apply_sex_filter_if_requested(design, cfg);
    std::cout << "[design] after sex filter: n=" << design.n << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[design] sex filter skipped: " << e.what() << "\n";
  }

  std::cout << "PATH" << paths.bed << "\n";
  // Ensure paths exist (unless sparse-only make)
  auto must_exist = [&](const std::string& p, const char* what){
    if (p.empty() || !fs::exists(p)) {
      std::ostringstream oss; oss << "ERROR: " << what << " not found: " << p;
      throw std::runtime_error(oss.str());
    }
  };
  if (!cfg.use_sparse_grm_to_fit || !cfg.make_sparse_grm_only) {
    must_exist(paths.bed, "BED");
    must_exist(paths.bim, "BIM");
    must_exist(paths.fam, "FAM");
  }

  // Inverse-normalize for quantitative (optional)
  if (cfg.inv_normalize && ieq(cfg.trait,"quantitative")) {
    arma::vec yv(design.n);
    for (int i=0;i<design.n;++i) yv(i)=design.y[i];
    arma::uvec fin = arma::find_finite(yv);
    arma::vec sub = yv.elem(fin);
    arma::uvec ord = arma::sort_index(sub);
    arma::uvec ranks(ord.n_elem);
    for (size_t k=0;k<ord.n_elem;++k) ranks(ord(k)) = k+1;
    for (size_t t=0;t<fin.n_elem;++t) {
      double p = (ranks(t)-0.5) / double(fin.n_elem);
      yv(fin(t)) = probit(p);
    }
    for (int i=0;i<design.n;++i) design.y[i]=yv(i);
    std::cout << "[design] inverse-normalized phenotype\n";
  }

  // Covariate offset path: fit β once, offset=Xβ, drop X
  if (cfg.covariate_offset && design.p>0) {
    arma::mat X(design.n, design.p);
    arma::vec yv(design.n);
    for (int i=0;i<design.n;++i) yv(i)=design.y[i];
    for (int i=0;i<design.n;++i)
      for (int j=0;j<design.p;++j)
        X(i,j)=design.X[(size_t)i*(size_t)design.p + (size_t)j];

    // Very small ridge to avoid pathologies; use Gaussian closed form if trait is quant,
    // otherwise 5-6 IRLS steps for logistic.
    arma::vec beta(design.p, arma::fill::zeros);
    if (ieq(cfg.trait,"quantitative")) {
      double lam = 1e-8;
      beta = arma::solve(X.t()*X + lam*arma::eye(design.p,design.p), X.t()*yv);
    } else {
      arma::vec b(design.p, arma::fill::zeros);
      arma::vec mu, w, z;
      for (int it=0; it<8; ++it) {
        mu = 1.0 / (1.0 + arma::exp(-X*b));
        w  = mu % (1.0 - mu) + 1e-8;
        z  = X*b + (yv - mu) / w;
        arma::mat XtW = X.t() * arma::diagmat(w);
        b = arma::solve(XtW * X, XtW * z);
      }
      beta = b;
    }
    arma::vec off = X * beta;
    if (design.offset.size() != (size_t)design.n) design.offset.assign(design.n, 0.0);
    for (int i=0;i<design.n;++i) design.offset[i] += off(i);
    design.X.clear(); design.p=0;
    std::cout << "[design] covariate_offset=true → added Xβ to offset and dropped X\n";
  }

  // Ensure output dirs exist
  ensure_parent_dir(paths.out_prefix + ".touch");
  ensure_parent_dir(paths.out_prefix_vr + ".touch");

  // LOCO: precompute ranges (pass into solver/VR if your engines accept it)
  std::vector<std::pair<size_t,size_t>> loco_ranges;
  if (cfg.loco) {
    loco_ranges = scan_bim_chr_ranges(paths.bim);
    std::cout << "[loco] ranges=" << loco_ranges.size() << "\n";
    // TODO: plumb loco_ranges to your engine via cfg or a setter.
  }

  // FAM alignment (IID->1-based index)
  // FIXED: indicatorWithPheno should have N elements (FAM size), not design.n elements
  // This matches R's: indicatorGenoSamplesWithPheno = (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)
  std::vector<int> subSampleInGeno;
  std::vector<bool> indicatorWithPheno;
  {
    auto fam_iids = read_fam_iids(paths.fam);
    int N_fam = static_cast<int>(fam_iids.size());

    // Create mapping from IID to FAM index (1-based)
    std::unordered_map<std::string,int> fam_pos; fam_pos.reserve(N_fam*2);
    for (int i=0;i<N_fam;++i) fam_pos.emplace(fam_iids[i], i+1); // 1-based

    // Initialize indicator with N elements, all false
    indicatorWithPheno.resize(N_fam, false);

    // For each phenotype sample, mark its FAM position in the indicator
    subSampleInGeno.reserve(design.n);
    for (const auto& id : design.iid) {
      auto it = fam_pos.find(id);
      if (it == fam_pos.end()) throw std::runtime_error("IID in design not found in FAM: " + id);
      int fam_idx_1based = it->second;
      subSampleInGeno.push_back(fam_idx_1based);
      indicatorWithPheno[fam_idx_1based - 1] = true;  // Convert to 0-based index
    }

    std::cout << "[FAM] N_fam=" << N_fam << ", design.n=" << design.n
              << ", indicatorWithPheno.size()=" << indicatorWithPheno.size() << std::endl;
  }

  // Initialize genotype data BEFORE sparse GRM section (needed for build_sparse_grm_in_place)
  init_global_geno(paths.bed, paths.bim, paths.fam, subSampleInGeno, indicatorWithPheno, cfg.isDiagofKinSetAsOne, cfg.min_maf_grm, cfg.max_miss_grm);

  // ============ CHECKPOINT 1: After genotype loading ============
  {
    static const std::string CP_DIR = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/checkpoints";
    int CP1_M = gettotalMarker();
    int CP1_N = (int)subSampleInGeno.size();  // N samples with phenotype
    int CP1_Nnomissing = getNnomissingOut();
    arma::fvec CP1_alleleFreq = getAlleleFreqVec();

    std::cout << "\n=== C++ CHECKPOINT 1: After genotype loading ===" << std::endl;
    std::cout << "CP1: M (markers) = " << CP1_M << std::endl;
    std::cout << "CP1: N (samples) = " << CP1_N << std::endl;
    std::cout << "CP1: Nnomissing = " << CP1_Nnomissing << std::endl;
    std::cout << "CP1: alleleFreq[0:10] = ";
    for (int i = 0; i < std::min(10, (int)CP1_alleleFreq.n_elem); ++i) std::cout << CP1_alleleFreq[i] << " ";
    std::cout << std::endl;

    // Save to files
    std::ofstream f_cp1(CP_DIR + "/CPP_CP1_genotype.csv");
    f_cp1 << "M," << CP1_M << "\n";
    f_cp1 << "N," << CP1_N << "\n";
    f_cp1 << "Nnomissing," << CP1_Nnomissing << "\n";
    f_cp1.close();

    CP1_alleleFreq.save(CP_DIR + "/CPP_CP1_alleleFreq.csv", arma::csv_ascii);
    std::cout << "Saved checkpoint files to " << CP_DIR << std::endl;
  }

  // -------- Sparse GRM build/reuse (+enforcement) --------
  if (cfg.use_sparse_grm_to_fit) {
    const bool have_files =
      !paths.sparse_grm.empty()     && fs::exists(paths.sparse_grm) &&
      !paths.sparse_grm_ids.empty() && fs::exists(paths.sparse_grm_ids);

    if (have_files) {
      arma::umat loc; arma::vec val; int n_mtx=0;
      load_matrix_market_coo(paths.sparse_grm, loc, val, n_mtx);
      setupSparseGRM(n_mtx, loc, val);
      setisUseSparseSigmaforInitTau(true);
      setisUseSparseSigmaforNullModelFitting(true);
      std::cout << "[sparse] Reusing GRM: " << paths.sparse_grm
                << "  n=" << n_mtx << "  nnz=" << val.n_elem << "\n";
    } else {
      double rc = (cfg.relatedness_cutoff > 0.0 ? cfg.relatedness_cutoff : 0.05);
      build_sparse_grm_in_place(rc, cfg.min_maf_grm, cfg.max_miss_grm);
      auto loc = export_sparse_grm_locations();
      auto val = export_sparse_grm_values();
      int  n   = export_sparse_grm_dim();
      if (!paths.sparse_grm.empty()) {
        write_matrix_market_coo(loc, val, n, paths.sparse_grm);
        std::cout << "[sparse] Saved GRM to " << paths.sparse_grm
                  << "  n=" << n << "  nnz=" << val.n_elem << "\n";
      }
      if (!paths.sparse_grm_ids.empty()) {
        auto fam_iids = read_fam_iids(paths.fam);
        std::vector<std::string> id_out; id_out.reserve(subSampleInGeno.size());
        for (int pos1b : subSampleInGeno) id_out.push_back(fam_iids[pos1b-1]);
        if ((int)id_out.size() != n) {
          std::cerr << "[warn] ID list size (" << id_out.size()
                    << ") differs from GRM n (" << n << ").\n";
        }
        write_id_list(id_out, paths.sparse_grm_ids);
      }
      setisUseSparseSigmaforInitTau(true);
      setisUseSparseSigmaforNullModelFitting(true);
    }
    if (!get_isUseSparseSigmaforModelFitting()) {
      throw std::runtime_error("Sparse GRM enabled, but sparse-Sigma flag is off; aborting.");
    }
  }

  // Early exit: construct-only
  if (cfg.make_sparse_grm_only) {
    std::cout << "[ok] make_sparse_grm_only=true: exiting before null model fit.\n";
    return 0;
  }

  // Variance-ratio overwrite guard
  if (cfg.num_markers_for_vr > 0) {
    std::string vr_txt = paths.out_prefix_vr + ".varianceRatio.txt";
    bool allow_overwrite = (y["paths"] && y["paths"]["overwrite_varratio"]) ? y["paths"]["overwrite_varratio"].as<bool>() : false;
    if (!allow_overwrite && fs::exists(vr_txt)) {
      std::ostringstream oss;
      oss << "Refusing to overwrite existing variance-ratio file: " << vr_txt
          << " (set paths.overwrite_varratio=true to allow).";
      throw std::runtime_error(oss.str());
    }
  }

  // ------------------ genoClass integration (optional) ------------------
  // genoClass geno;
  // geno.isVarRatio = (cfg.num_markers_for_vr > 0 || cfg.use_sparse_grm_for_vr);
  // geno.g_minMACVarRatio = static_cast<float>(cfg.vr_min_mac > 0 ? cfg.vr_min_mac : 1);
  // geno.g_maxMACVarRatio = static_cast<float>(cfg.vr_max_mac != 0 ? cfg.vr_max_mac : -1);
  // geno.setGenoObj(paths.bed, paths.bim, paths.fam,
  //                 subSampleInGeno, indicatorWithPheno,
  //                 static_cast<float>(cfg.memory_chunk_gb > 0 ? cfg.memory_chunk_gb : 1.0f),
  //                 /* isDiagofKinSetAsOne */ false);

  // NOTE: init_global_geno moved earlier (before sparse GRM section) for make_sparse_grm_only support

  // ------------------ Run null fit ------------------
  // If you have the overload with genoClass&, call:
  // FitNullResult out = saige::fit_null(cfg, paths, design, geno);

  FitNullResult out = saige::fit_null(cfg, paths, design);

  // ------------------ Output GRM diagonal (after fit_null, same as R version) ------------------
  output_grm_diagonal(paths.out_prefix + ".grm_diag.txt");

  // ------------------ Report artifacts ------------------
  std::cout << "== SAIGE Null Fit Completed ==\n";
  std::cout << "Model artifact: " << out.model_rda_path << "\n";
  if (!out.vr_path.empty())           std::cout << "Variance ratio: " << out.vr_path << "\n";
  if (!out.markers_out_path.empty())  std::cout << "Marker results: " << out.markers_out_path << "\n";
  std::cout << "LOCO: " << (out.loco ? "on" : "off")
            << "  LowMem: " << (out.lowmem_loco ? "yes" : "no") << "\n";
  return 0;
}
