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
#include "covariate_offset.hpp"
#include "glmm.hpp"
#include "loco_engine.hpp"
#include <dlfcn.h>
#include <cstdlib>
#include <thread>
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "SAIGE_step1_fast.hpp"   // (optional) genoClass decl — comment out if not available
#include "preprocess_engine.hpp"

#include <yaml-cpp/yaml.h>
#include <cxxopts.hpp>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <RcppArmadillo.h>
#include <Rembedded.h>
#include <iomanip>
#include <chrono>
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
      next.reset(node[k]);
    } else if (!next.IsMap()) {
      // Don't smash existing non-map nodes when asked to go deeper
      throw std::runtime_error(
        "Override '" + dotted + "' expects '" + k +
        "' to be a map, but it is " +
        (next.IsScalar() ? "a scalar" :
         next.IsSequence()? "a sequence" : "unknown") + ".");
    }
    // yaml-cpp: Node::operator= copies the value INTO the node this handle
    // refers to, while reset() rebinds the handle. With '=' the first step
    // wrote the child map over root itself, so any nested override such as
    // '-o fit.nthreads=1' replaced the whole config with the 'fit' map.
    node.reset(next);
  }

  // Set the leaf (this updates only the designated item)
  node[parts.back()] = value;
}

// ------------------ YAML loaders ------------------
static FitNullConfig load_cfg(const YAML::Node& y) {
  FitNullConfig c;

  // Sex-specific fit (R: sexCol / FemaleOnly / MaleOnly / FemaleCode /
  // MaleCode). These were declared in FitNullConfig but never parsed, so a
  // config that set them silently fitted both sexes. Accepted under design:
  // (next to the other column names and row filters) or under fit:; setting
  // the same key in both places is an error rather than a silent pick.
  {
    const auto sd = y["design"];
    const auto sf = y["fit"];
    auto where = [&](const char* k) -> int {
      const bool in_d = sd && sd.IsMap() && sd[k];
      const bool in_f = sf && sf.IsMap() && sf[k];
      if (in_d && in_f)
        throw std::runtime_error(std::string("config: '") + k +
                                 "' is set under both design: and fit:; set it in one place");
      return in_d ? 1 : (in_f ? 2 : 0);
    };
    auto str_of = [&](const char* k, std::string& dst) {
      if (int w = where(k)) dst = trim((w == 1 ? sd[k] : sf[k]).as<std::string>());
    };
    auto bool_of = [&](const char* k, bool& dst) {
      if (int w = where(k)) dst = (w == 1 ? sd[k] : sf[k]).as<bool>();
    };
    str_of("sex_col", c.sex_col);
    bool_of("female_only", c.female_only);
    bool_of("male_only", c.male_only);
    str_of("female_code", c.female_code);
    str_of("male_code", c.male_code);
  }

  const auto f = y["fit"];
  if (!f) return c;

  auto get = [&](const char* k){ return f[k]; };
  if (get("trait")) c.trait = get("trait").as<std::string>();
  if (get("loco")) c.loco = get("loco").as<bool>();
  if (get("lowmem_loco")) c.lowmem_loco = get("lowmem_loco").as<bool>();
  if (get("use_sparse_grm_to_fit")) c.use_sparse_grm_to_fit = get("use_sparse_grm_to_fit").as<bool>();
  if (get("use_gpu")) c.use_gpu = get("use_gpu").as<bool>();
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
  if (get("trace_seed")) c.trace_seed = get("trace_seed").as<int>();
  if (get("nthreads")) c.nthreads = get("nthreads").as<int>();
  if (get("traceCVcutoff")) c.traceCVcutoff = get("traceCVcutoff").as<double>();
  if (get("ratio_cv_cutoff")) c.ratio_cv_cutoff = get("ratio_cv_cutoff").as<double>();
  if (get("min_maf_grm")) c.min_maf_grm = get("min_maf_grm").as<double>();
  if (get("max_miss_grm")) c.max_miss_grm = get("max_miss_grm").as<double>();
  if (get("num_markers_for_vr")) c.num_markers_for_vr = get("num_markers_for_vr").as<int>();

  // step-2 knobs carried through nullmodel.json. Defaults keep the previously
  // hardcoded values; note fast_test and impute_method do NOT match the R CLI
  // defaults (R: is_fastTest=FALSE, impute_method=best_guess), so set them
  // explicitly when comparing against R.
  if (get("fast_test"))          c.fast_test          = get("fast_test").as<bool>();
  if (get("impute_method"))      c.impute_method      = get("impute_method").as<std::string>();
  if (get("spa_cutoff"))         c.spa_cutoff         = get("spa_cutoff").as<double>();
  if (get("p_cutoff_for_firth")) c.p_cutoff_for_firth = get("p_cutoff_for_firth").as<double>();
  if (get("firth_beta"))         c.firth_beta         = get("firth_beta").as<bool>() ? 1 : 0;
  // Categorical variance ratio (R: --isCateVarianceRatio + MAC bin vectors)
  if (get("isCateVarianceRatio")) c.isCateVarianceRatio = get("isCateVarianceRatio").as<bool>();
  if (get("cateVarRatioMinMACVecExclude") && get("cateVarRatioMinMACVecExclude").IsSequence()) {
    c.cateVarRatioMinMACVecExclude.clear();
    for (const auto& v : get("cateVarRatioMinMACVecExclude")) c.cateVarRatioMinMACVecExclude.push_back(v.as<double>());
  }
  if (get("cateVarRatioMaxMACVecInclude") && get("cateVarRatioMaxMACVecInclude").IsSequence()) {
    c.cateVarRatioMaxMACVecInclude.clear();
    for (const auto& v : get("cateVarRatioMaxMACVecInclude")) c.cateVarRatioMaxMACVecInclude.push_back(v.as<double>());
  }
  if (get("cateVarRatioIndexVec") && get("cateVarRatioIndexVec").IsSequence()) {
    c.cateVarRatioIndexVec.clear();
    for (const auto& v : get("cateVarRatioIndexVec")) c.cateVarRatioIndexVec.push_back(v.as<int>());
  }
  if (get("event_time_bin_size") && !get("event_time_bin_size").IsNull())
    c.event_time_bin_size = get("event_time_bin_size").as<int>();
  if (get("relatedness_cutoff")) c.relatedness_cutoff = get("relatedness_cutoff").as<double>();
  if (get("make_sparse_grm_only")) c.make_sparse_grm_only = get("make_sparse_grm_only").as<bool>();
  if (get("memory_chunk_gb")) c.memory_chunk_gb = get("memory_chunk_gb").as<double>();
  if (get("vr_min_mac")) c.vr_min_mac = get("vr_min_mac").as<int>();
  if (get("vr_max_mac")) c.vr_max_mac = get("vr_max_mac").as<int>();
  if (get("diag_one")) c.isDiagofKinSetAsOne = get("diag_one").as<bool>();
  if (get("use_pcg_with_sparse_grm")) c.use_pcg_with_sparse_grm = get("use_pcg_with_sparse_grm").as<bool>();
  if (get("multi_lockstep")) c.multi_lockstep = get("multi_lockstep").as<bool>();
  if (get("mask_missing")) c.mask_missing = get("mask_missing").as<bool>();
  if (get("mask_min_coverage")) c.mask_min_coverage = get("mask_min_coverage").as<double>();
  if (get("scheme_c_break")) c.scheme_c_break = get("scheme_c_break").as<std::string>();
  if (get("use_blocked_gemv")) c.use_blocked_gemv = get("use_blocked_gemv").as<bool>();
  if (get("gemv_block_size")) c.gemv_block_size = get("gemv_block_size").as<int>();
  if (get("gemv_verify")) c.gemv_verify = get("gemv_verify").as<bool>();
  if (get("overwrite_vr")) c.overwrite_vr = get("overwrite_vr").as<bool>();
  if (get("skip_model_fitting")) c.skip_model_fitting = get("skip_model_fitting").as<bool>();
  if (get("model_file")) c.model_file = get("model_file").as<std::string>();
  if (get("dry_run")) c.dry_run = get("dry_run").as<bool>();

  // Parse q_covar_cols from YAML (categorical covariate names)
  const auto d = y["design"];
  if (d && d["q_covar_cols"]) {
    const auto& qnode = d["q_covar_cols"];
    if (qnode.IsSequence()) {
      for (const auto& item : qnode)
        c.q_covar_cols.push_back(item.as<std::string>());
    }
  }
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

// NOTE: scan_bim_chr_ranges() used to live here. It scanned the BIM for raw
// per-chromosome marker ranges, printed a count and threw the result away, and
// its indices were in raw-BIM space rather than the compacted post-QC space the
// genotype object uses. LOCO ranges are now computed once in
// PreprocessEngine::compute_chr_ranges_from_bim_() and flow through PreOut::chr.

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

// ------------------ Sex-specific row filter ------------------
// R (SAIGE_fitGLMM_fast.R:1189-1212): after complete.cases() on the phenotype /
// covariate / ID columns, FemaleOnly keeps data[which(data[, sexCol] == FemaleCode), ]
// (MaleOnly likewise), and stops if nothing is left. A missing sex value
// compares NA and is dropped by which().
struct SexRowFilter {
  std::string col;    // design.sex_col
  std::string code;   // female_code or male_code
  std::string label;  // "female_only" / "male_only", for messages
};

// R compares with `==`: fread reads an all-numeric sex column as numbers, so a
// file cell "1.0" matches the default code "1". Match on the exact string, or
// numerically when both the cell and the code are plain numbers.
static bool sex_code_matches(const std::string& v, const std::string& code) {
  if (is_missing(v)) return false;
  if (v == code) return true;
  char* ev = nullptr; char* ec = nullptr;
  const double a = std::strtod(v.c_str(), &ev);
  const double b = std::strtod(code.c_str(), &ec);
  return ev != v.c_str() && *ev == '\0' && ec != code.c_str() && *ec == '\0'
         && std::isfinite(a) && a == b;
}

// ------------------ Design CSV/TSV/space parser + categorical encoding ------------------
// Expected header columns (case-insensitive): <iid_col>, <y_col>, [offset], [time|event_time|eventTime], X...
// iid_col and y_col default to "IID" and "y" for backward compatibility.
// Rows where the phenotype cell is empty, "NA", or "NaN" are silently dropped.
// With `sex` set, rows whose own sex cell does not match are dropped at the
// same point, before categorical levels, min_covariate_count and the intercept
// check are worked out, so a sex-specific run sees exactly what a design file
// pre-filtered to that sex would give.
static Design load_design_csv(const std::string& path,
                              int min_covariate_count,
                              bool categorical_drop_reference,
                              const std::vector<std::string>& covar_col_names = {},
                              const std::string& iid_col = "IID",
                              const std::string& y_col   = "y",
                              const SexRowFilter* sex = nullptr)
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

  int idx_iid    = find_col({iid_col.c_str()});
  int idx_y      = find_col({y_col.c_str()});
  int idx_offset = find_col({"offset", "covoffset"});
  int idx_time   = find_col({"time","event_time","eventTime"});

  if (idx_iid < 0)
    throw std::runtime_error("Design file: sample ID column '" + iid_col + "' not found. "
                             "Set design.iid_col in YAML if your file uses a different name.");
  if (idx_y < 0)
    throw std::runtime_error("Design file: phenotype column '" + y_col + "' not found. "
                             "Set design.y_col in YAML if your file uses a different name.");
  int idx_sex = -1;
  if (sex) {
    idx_sex = find_col({sex->col.c_str()});
    if (idx_sex < 0)
      throw std::runtime_error("ERROR: column for sex '" + sex->col +
                               "' (design.sex_col) does not exist in the design file " + path);
  }

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

  // ===== Step 12: Drop rows with any missing value (R line 1430: complete.cases) =====
  // R: data = data[complete.cases(data),,drop=F]
  // Checks phenotype AND all covariate columns for empty/NA/NaN.
  {
    int n_na = 0;
    std::vector<std::vector<std::string>> valid_rows;
    valid_rows.reserve(rows.size());
    for (auto& row : rows) {
      bool any_missing = false;
      // Check phenotype
      const std::string& yval = row[idx_y];
      if (yval.empty() || ieq(yval, "NA") || ieq(yval, "NaN")) {
        any_missing = true;
      }
      // Check all covariate columns
      if (!any_missing) {
        for (int j : x_idx) {
          const std::string& cv = row[j];
          if (cv.empty() || ieq(cv, "NA") || ieq(cv, "NaN")) {
            any_missing = true;
            break;
          }
        }
      }
      if (any_missing) {
        ++n_na;
      } else {
        valid_rows.push_back(std::move(row));
      }
    }
    if (n_na > 0)
      std::cout << "[design] dropped " << n_na
                << " row(s) with missing phenotype or covariates (complete.cases)\n";
    rows = std::move(valid_rows);
  }

  // Sex-specific fit: keep rows whose sex cell matches the requested code.
  if (sex) {
    std::vector<std::vector<std::string>> kept;
    kept.reserve(rows.size());
    int n_other = 0, n_miss = 0;
    const int before = (int)rows.size();
    for (auto& row : rows) {
      const std::string& v = row[idx_sex];
      if (sex_code_matches(v, sex->code)) kept.push_back(std::move(row));
      else if (is_missing(v))             ++n_miss;
      else                                ++n_other;
    }
    rows = std::move(kept);
    std::cout << "[design] " << sex->label << ": " << sex->col << " == " << sex->code
              << " kept " << rows.size() << " of " << before << " row(s) ("
              << n_other << " other value, " << n_miss << " missing sex)\n";
    if (rows.empty())
      throw std::runtime_error("ERROR: no samples in the phenotype are coded as " +
                               sex->code + " in the column " + sex->col);
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
  // X (row-major: d.X[i*p + j]).
  // BUGFIX: this used to round-trip through an arma::mat and write the result
  // back with d.X.assign(Xk.begin(), Xk.end()). Armadillo iterators walk
  // COLUMN-major, so the row-major buffer every other reader assumes was left
  // transposed whenever p > 1 — silently corrupting the covariates on the two
  // paths that call this (duplicate-IID removal and the design.whitelist_ids
  // filter), and on the multi-phenotype intersection added later. Copy rows
  // directly instead; apply_row_subset (preprocess_engine.cpp:379) already did.
  if (d.p>0 && !d.X.empty()){
    std::vector<double> X2((size_t)n2 * (size_t)d.p);
    for (int r=0;r<n2;++r) {
      const double* src = &d.X[(size_t)keep[r]*(size_t)d.p];
      std::copy(src, src + d.p, &X2[(size_t)r*(size_t)d.p]);
    }
    d.X.swap(X2);
  }
  d.n = n2;
}
 
// ------------------ main ------------------
int main(int argc, char** argv) {
  // Initialize R's embedded runtime so we can use R's RNG (Mersenne Twister)
  // This makes set_seed() and Rf_rbinom() work in standalone mode
  {
    // R_HOME for the embedded R runtime. Honor an existing R_HOME; otherwise use
    // SAIGE_R_HOME if provided (set it to your R install, e.g. $CONDA_PREFIX/lib/R).
    if (!std::getenv("R_HOME")) { if (const char* rh = std::getenv("SAIGE_R_HOME")) setenv("R_HOME", rh, 0); }
    char* r_argv[] = { (char*)"saige-null", (char*)"--vanilla", (char*)"--no-readline", (char*)"--silent", nullptr };
    Rf_initialize_R(4, r_argv);
    extern uintptr_t R_CStackLimit;
    R_CStackLimit = (uintptr_t)-1;  // disable stack checking
    setup_Rmainloop();
    std::cout << "[R] Embedded R runtime initialized" << std::endl;
  }

  cxxopts::Options opts("saige-null", "Null GLMM fitting with LOCO/VR (genoClass-integrated)");
  opts.add_options()
    ("c,config",   "YAML config path", cxxopts::value<std::string>())
    ("d,design",   "Override design CSV/TSV path", cxxopts::value<std::string>()->default_value(""))
    ("o,override", "YAML dot-override, e.g., fit.nthreads=32", cxxopts::value<std::vector<std::string>>()->default_value({}))
    ("t,threads",  "Override fit.nthreads (shortcut for -o fit.nthreads=N; 0 = hw concurrency)", cxxopts::value<int>()->default_value("-1"))
    ("gpu",        "Enable GPU (cuBLAS) K·u acceleration; falls back to CPU if unavailable", cxxopts::value<bool>()->default_value("false"))
    ("blocked-gemv","Enable CPU blocked-GEMV K·u (parallelCrossProd_blocked); off by default", cxxopts::value<bool>()->default_value("false"))
    ("gemv-block",  "Block size for blocked-GEMV K·u (e.g. 64, 128, 256)", cxxopts::value<int>()->default_value("-1"))
    ("gemv-verify", "Run both K·u paths and print rel/abs error (debug)", cxxopts::value<bool>()->default_value("false"))
    ("v,verbose",  "Verbose", cxxopts::value<bool>()->default_value("false"))
    ("dry-run",    "Validate inputs only (no genotype loading or solver)", cxxopts::value<bool>()->default_value("false"))
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

  // CLI --dry-run overrides YAML
  if (res["dry-run"].as<bool>()) cfg.dry_run = true;

  // CLI --threads N overrides cfg.nthreads
  {
    int t_cli = res["threads"].as<int>();
    if (t_cli >= 0) cfg.nthreads = t_cli;  // 0 = auto (hardware_concurrency resolved in configure_threads)
  }

  // CLI --gpu overrides cfg.use_gpu
  if (res["gpu"].as<bool>()) cfg.use_gpu = true;
  // Propagate to setGenoObj via env var (same pattern as RCPP_PARALLEL_NUM_THREADS)
  setenv("SAIGE_USE_GPU", cfg.use_gpu ? "1" : "0", 1);

  // CLI overrides for blocked-GEMV K·u path
  if (res["blocked-gemv"].as<bool>()) cfg.use_blocked_gemv = true;
  if (res["gemv-verify"].as<bool>()) cfg.gemv_verify = true;
  {
    int bs_cli = res["gemv-block"].as<int>();
    if (bs_cli > 0) cfg.gemv_block_size = bs_cli;
  }
  setenv("SAIGE_USE_BLOCKED_GEMV", cfg.use_blocked_gemv ? "1" : "0", 1);
  setenv("SAIGE_GEMV_BLOCK_SIZE",  std::to_string(cfg.gemv_block_size).c_str(), 1);
  setenv("SAIGE_GEMV_VERIFY",      cfg.gemv_verify ? "1" : "0", 1);

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

  // Configure thread pools (RcppParallel / OpenMP / OpenBLAS).
  // Design:
  //   - cfg.nthreads drives RcppParallel (TBB) for marker-level parallelism
  //     (the parallelCrossProd worker) via the RCPP_PARALLEL_NUM_THREADS env var.
  //   - OpenMP for any #pragma omp regions.
  //   - OpenBLAS is pinned to 1 thread to avoid over-subscription when
  //     marker-level TBB workers call into BLAS.
  //   - If the user pre-set OMP_NUM_THREADS / OPENBLAS_NUM_THREADS / RCPP_PARALLEL_NUM_THREADS
  //     in the environment, respect it (don't overwrite).
  {
    int n = cfg.nthreads;
    if (n <= 0) n = static_cast<int>(std::thread::hardware_concurrency());
    if (n <= 0) n = 1;
    auto set_if_unset = [](const char* k, const std::string& v) {
      if (std::getenv(k) == nullptr) setenv(k, v.c_str(), 1);
    };
    set_if_unset("RCPP_PARALLEL_NUM_THREADS", std::to_string(n));
    set_if_unset("OMP_NUM_THREADS", std::to_string(n));
    // Pin BLAS to 1 when we're doing outer parallelism to avoid over-subscription
    set_if_unset("OPENBLAS_NUM_THREADS", std::string(n > 1 ? "1" : "1"));
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
    // Best-effort runtime OpenBLAS thread control (symbol may be absent)
    if (auto sym = dlsym(RTLD_DEFAULT, "openblas_set_num_threads")) {
      reinterpret_cast<void(*)(int)>(sym)(n > 1 ? 1 : 1);
    }
    std::cout << "[threads] cfg.nthreads=" << cfg.nthreads
              << " effective=" << n
              << " RCPP_PARALLEL_NUM_THREADS=" << (std::getenv("RCPP_PARALLEL_NUM_THREADS") ?: "(unset)")
              << " OMP_NUM_THREADS=" << (std::getenv("OMP_NUM_THREADS") ?: "(unset)")
              << " OPENBLAS_NUM_THREADS=" << (std::getenv("OPENBLAS_NUM_THREADS") ?: "(unset)")
              << "\n";
  }

  saige::register_default_solvers();
  saige::register_default_loco_batch();

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
  std::string iid_col_name = "IID";
  std::string y_col_name   = "y";
  if (y["design"] && y["design"]["iid_col"])
    iid_col_name = y["design"]["iid_col"].as<std::string>();
  if (y["design"] && y["design"]["y_col"])
    y_col_name = y["design"]["y_col"].as<std::string>();

  // ===== Tier-1 multi-phenotype: model specs =====
  // Three mutually exclusive ways to name the phenotype(s):
  //   design.y_col: y                     (single, legacy — unchanged behaviour)
  //   design.y_cols: [q1, q2, ...]        (shorthand; per-trait paths derived)
  //   models: [{y_col:, out_prefix:, out_prefix_vr:}, ...]   (explicit paths)
  // The derived-path rule for y_cols mirrors how the two prefixes are consumed
  // downstream: out_prefix is used as a DIRECTORY (nullmodel.json + *.arma live
  // inside it) so the trait becomes a subdirectory, while out_prefix_vr is a
  // FILE prefix (+ ".varianceRatio.txt") so the trait is suffixed.
  struct ModelSpec { std::string y_col, out_prefix, out_prefix_vr; };
  std::vector<ModelSpec> models;
  {
    const bool has_models = (bool)y["models"];
    const bool has_ycols  = (bool)(y["design"] && y["design"]["y_cols"]);
    const bool has_ycol   = (bool)(y["design"] && y["design"]["y_col"]);
    if (has_models && (has_ycols || has_ycol))
      throw std::runtime_error("config: `models:` is mutually exclusive with design.y_col / design.y_cols");
    if (has_models && has_ycols)
      throw std::runtime_error("config: `models:` is mutually exclusive with design.y_cols");
    if (has_ycols && has_ycol)
      throw std::runtime_error("config: design.y_cols is mutually exclusive with design.y_col");

    if (has_models) {
      const auto& mn = y["models"];
      if (!mn.IsSequence() || mn.size() == 0)
        throw std::runtime_error("config: `models:` must be a non-empty sequence");
      for (const auto& m : mn) {
        ModelSpec ms;
        if (!m["y_col"]) throw std::runtime_error("config: every models[] entry needs y_col");
        ms.y_col         = m["y_col"].as<std::string>();
        ms.out_prefix    = m["out_prefix"]    ? m["out_prefix"].as<std::string>()    : (paths.out_prefix + "/" + ms.y_col);
        ms.out_prefix_vr = m["out_prefix_vr"] ? m["out_prefix_vr"].as<std::string>() : (paths.out_prefix_vr + "_" + ms.y_col);
        models.push_back(std::move(ms));
      }
    } else if (has_ycols) {
      const auto& cn = y["design"]["y_cols"];
      if (!cn.IsSequence() || cn.size() == 0)
        throw std::runtime_error("config: design.y_cols must be a non-empty sequence");
      for (const auto& c : cn) {
        ModelSpec ms;
        ms.y_col         = c.as<std::string>();
        ms.out_prefix    = paths.out_prefix + "/" + ms.y_col;
        ms.out_prefix_vr = paths.out_prefix_vr + "_" + ms.y_col;
        models.push_back(std::move(ms));
      }
    } else {
      // Single-phenotype: byte-for-byte the legacy path (prefixes untouched).
      models.push_back(ModelSpec{y_col_name, paths.out_prefix, paths.out_prefix_vr});
    }
    {
      std::set<std::string> seen_y, seen_p;
      for (const auto& m : models) {
        if (!seen_y.insert(m.y_col).second)
          throw std::runtime_error("config: duplicate phenotype column '" + m.y_col + "'");
        if (!seen_p.insert(m.out_prefix).second)
          throw std::runtime_error("config: duplicate out_prefix '" + m.out_prefix + "' across models");
      }
    }
    if (models.size() > 1) {
      std::cout << "[multi-pheno] P=" << models.size() << " phenotypes (genotype load shared per sample-set group):\n";
      for (const auto& m : models)
        std::cout << "    " << m.y_col << "  ->  " << m.out_prefix
                  << "   vr: " << m.out_prefix_vr << "\n";
    }
  }

  // ===== Step 14: Validate q_covar_cols subset of covar_cols (R lines 1446-1454) =====
  // R: if(!all(qCovarCol %in% covarColList)) stop("ERROR! all covariates in qCovarCol must be in covarColList")
  if (!cfg.q_covar_cols.empty()) {
    std::unordered_set<std::string> covar_set(covar_col_names.begin(), covar_col_names.end());
    for (const auto& qc : cfg.q_covar_cols) {
      if (covar_set.find(qc) == covar_set.end()) {
        throw std::runtime_error(
            "ERROR: categorical covariate '" + qc
            + "' in q_covar_cols is not in covar_cols. "
            "All q_covar_cols must be a subset of covar_cols.");
      }
    }
    std::cout << "[config] q_covar_cols (categorical): ";
    for (const auto& qc : cfg.q_covar_cols) std::cout << qc << " ";
    std::cout << "\n";
  }

  std::cout << "[config] covar_cols=[";
  for (size_t i = 0; i < covar_col_names.size(); ++i) {
    std::cout << covar_col_names[i];
    if (i < covar_col_names.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << (covar_col_names.empty() ? " (NO COVARIATES)" : "") << std::endl;
  std::cout << "[config] iid_col=" << iid_col_name << "  y_col=" << y_col_name << "\n";

  // ===== Step 15: sex-specific fit (R SAIGE_fitGLMM_fast.R:1060-1072, 1189-1212) =====
  // The filter itself runs inside load_design_csv on each row's own sex cell
  // (stage A), so it drops the same samples for every phenotype and a trait's
  // row set is (non-missing for that trait) AND (sex matches).
  std::optional<SexRowFilter> sex_filter;
  {
    // R: if (FemaleOnly & MaleOnly) stop("Both FemaleOnly and MaleOnly are TRUE...")
    if (cfg.female_only && cfg.male_only)
      throw std::runtime_error(
          "ERROR: Both female_only and male_only are true. "
          "Please specify only one to run a sex-specific job.");
    if ((cfg.female_only || cfg.male_only) && cfg.sex_col.empty())
      throw std::runtime_error(
          "ERROR: female_only or male_only is true but sex_col is not specified in config.");
    if (cfg.female_only || cfg.male_only) {
      sex_filter = SexRowFilter{cfg.sex_col,
                                cfg.female_only ? cfg.female_code : cfg.male_code,
                                cfg.female_only ? "female_only" : "male_only"};
      std::cout << "[config] " << (cfg.female_only ? "Female" : "Male")
                << "-specific model will be fitted: only samples coded as "
                << sex_filter->code << " in column " << cfg.sex_col << " are included\n";
    } else if (!cfg.sex_col.empty()) {
      std::cout << "[config] sex_col=" << cfg.sex_col
                << " is set but neither female_only nor male_only is true -> no sex filter\n";
    }
  }

  // Parse design (with categoricals) - using configurable column names
  auto T0 = std::chrono::steady_clock::now();
  // Everything from the CSV parse through the covariate-offset GLM is
  // per-phenotype and is run in the SAME order as the single-trait path, on the
  // CSV row order (before the FAM reorder below). Keeping that order is what
  // makes each trait of a P>1 run byte-identical to running it alone.
  // Stage A: parse + de-duplicate. This is what FIXES the row set (rows with a
  // missing value in THIS phenotype or in a covariate, and with sex_filter set
  // rows of the other / missing sex, are dropped by load_design_csv), so it has
  // to run for every trait before the sample sets can be reconciled.
  auto load_design_stageA = [&](const std::string& y_col_name) -> Design {
  Design design = load_design_csv(design_csv, min_cov_ct, drop_ref, covar_col_names,
                                  iid_col_name, y_col_name,
                                  sex_filter ? &*sex_filter : nullptr);
  add_intercept_if_missing(design);

  // ===== Step 7: Duplicate sample ID removal (R line 1437) =====
  // R: sampleIDInclude[!duplicated(sampleIDInclude)]
  {
    std::unordered_set<std::string> seen;
    std::vector<size_t> keep;
    keep.reserve(design.n);
    int n_dup = 0;
    for (size_t i = 0; i < (size_t)design.n; ++i) {
      if (seen.insert(design.iid[i]).second) {
        keep.push_back(i);
      } else {
        ++n_dup;
      }
    }
    if (n_dup > 0) {
      std::cerr << "[warning] removed " << n_dup
                << " duplicate sample ID(s), keeping first occurrence\n";
      design_take_rows(design, keep);
    }
  }
  return design;
  };  // end load_design_stageA

  // Stage B: everything downstream that reads y (validity checks, optional
  // inverse-normalisation, the covariate-offset GLM) plus the y-independent
  // whitelist filter, which drops the same rows for every trait. Runs
  // AFTER the sample sets are reconciled, in the original order, on the CSV row
  // order — the covariate GLM accumulates over rows, so reordering first would
  // perturb beta in the last fp digits.
  auto finish_design_stageB = [&](Design& design, const Paths& paths) {

  // ===== Step 1: Binary phenotype must be 0 or 1 (R lines 1754-1757) =====
  // R: uniqPheno = sort(unique(y)); if (uniqPheno[1] != 0 | uniqPheno[2] != 1) stop(...)
  if (ieq(cfg.trait, "binary")) {
    std::set<double> unique_y;
    for (int i = 0; i < design.n; ++i) {
      double yval = design.y[i];
      if (yval != 0.0 && yval != 1.0) {
        throw std::runtime_error(
            "ERROR: binary phenotype value must be 0 or 1, found: "
            + std::to_string(yval)
            + " at sample " + design.iid[i]);
      }
      unique_y.insert(yval);
    }
    if (unique_y.size() < 2) {
      std::cerr << "[warning] binary phenotype has only one unique level ("
                << *unique_y.begin() << "), model fitting may be degenerate\n";
    }
  }

  // ===== Step 2: Phenotype variance check for quantitative (R lines 879-880) =====
  // R: if (abs(var(Y)) < 0.1) stop("WARNING: variance of the phenotype is much smaller than 1...")
  if (ieq(cfg.trait, "quantitative")) {
    double sum_y = 0.0, sum_y2 = 0.0;
    for (int i = 0; i < design.n; ++i) {
      sum_y  += design.y[i];
      sum_y2 += design.y[i] * design.y[i];
    }
    double mean_y = sum_y / design.n;
    double var_y  = sum_y2 / design.n - mean_y * mean_y;
    if (std::fabs(var_y) < 0.1) {
      throw std::runtime_error(
          "ERROR: variance of the phenotype (" + std::to_string(var_y)
          + ") is much smaller than 1. Please consider setting inv_normalize: true in config.");
    }
  }

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
  // Delegates to saige::fit_covariate_offset (covariate_offset.cpp) so the initial
  // GLM can be unit-tested against R's glm.fit across configs without rebuilding the
  // full benchmark pipeline. The function dumps <paths.out_prefix>_cov_{input,beta,fit,trace}.csv
  // which a companion R harness ingests to diff against glm.fit on identical (X, y).
  if (cfg.covariate_offset && design.p>0) {
    auto cor = saige::fit_covariate_offset(design, cfg.trait,
                                            /*max_iter=*/25,
                                            /*tol=*/1e-8,
                                            paths.out_prefix);
    if (design.offset.size() != (size_t)design.n) design.offset.assign(design.n, 0.0);
    for (int i = 0; i < design.n; ++i) design.offset[i] += cor.offset[i];
    // Snapshot the full-covariate X before collapsing — Step-2 needs it for the
    // marker-genotype projection g_tilde = g − X(X'X)⁻¹X'g.  The GLMM solver
    // itself runs on the collapsed (intercept-only + offset) design; we only
    // restore the full X at the .arma-save block in null_model_engine.cpp.
    design.X_full = design.X;
    design.p_full = design.p;
    // Keep the full-length initial-GLM beta so nullmodel.json can report a
    // p-length alpha (the GLMM below only re-fits the intercept).
    design.beta_full = cor.beta;
    design.X.assign(design.n, 1.0);
    design.p = 1;
    std::cout << "[design] covariate_offset=true → added Xβ(covariates) to offset, "
              << "X -> intercept only for fit (p=1), snapshot p_full=" << design.p_full
              << " preserved for Step-2 export\n";
  }

  };  // end finish_design_stageB

  // Build one Design per phenotype.
  std::vector<Design> designs;
  designs.reserve(models.size());
  for (const auto& m : models) {
    if (models.size() > 1)
      std::cout << "\n[multi-pheno] ===== design for phenotype '" << m.y_col << "' =====\n";
    designs.push_back(load_design_stageA(m.y_col));
  }

  // Traits whose sample sets differ (different missingness per phenotype).
  // Rows are still in CSV order here, so "same sample set" is just
  // iid-sequence equality.
  //
  //   - default: GROUP the traits by sample set (below, after stage B). Each
  //     group gets its own genotype load, GRM and GPU upload, so every trait
  //     keeps the property that a P>1 run reproduces its solo run exactly.
  //   - design.intersect_samples: true: keep the common samples only, one
  //     group. This is what the R tier-1 prototype does. It makes every
  //     trait's fit differ from its solo run, because the solo run would have
  //     used more samples — so it is opt-in and says so loudly.
  if (designs.size() > 1) {
    bool same = true;
    for (size_t k = 1; k < designs.size() && same; ++k)
      same = (designs[k].n == designs[0].n && designs[k].iid == designs[0].iid);

    const bool do_intersect =
        (y["design"] && y["design"]["intersect_samples"])
          ? y["design"]["intersect_samples"].as<bool>() : false;
    if (!same && do_intersect) {
      // Intersection: an IID kept by every trait. Each design's rows are a
      // subsequence of the CSV, so filtering each to the intersection leaves
      // all of them in the same order — no sorting needed.
      std::unordered_map<std::string,int> cnt;
      for (const auto& d : designs)
        for (const auto& id : d.iid) ++cnt[id];
      const int P_ = (int)designs.size();
      for (size_t k = 0; k < designs.size(); ++k) {
        std::vector<size_t> keep;
        keep.reserve(designs[k].n);
        for (size_t i = 0; i < (size_t)designs[k].n; ++i)
          if (cnt[designs[k].iid[i]] == P_) keep.push_back(i);
        const int before = designs[k].n;
        if ((int)keep.size() != before) design_take_rows(designs[k], keep);
        std::cout << "[multi-pheno] " << models[k].y_col << ": " << before
                  << " -> " << designs[k].n << " samples after intersection\n";
      }
      for (size_t k = 1; k < designs.size(); ++k)
        if (designs[k].iid != designs[0].iid)
          throw std::runtime_error("multi-phenotype: intersection did not align "
                                   "sample lists (internal error)");
      std::cerr << "[warning] design.intersect_samples=true: all "
                << designs.size() << " traits were fitted on the "
                << designs[0].n << " samples common to every phenotype. These "
                   "results are NOT comparable to single-trait runs, which "
                   "would each use more samples.\n";
    }
  }

  // Now the y-dependent work, per trait, in the original order.
  for (size_t mi = 0; mi < models.size(); ++mi) {
    if (models.size() > 1)
      std::cout << "\n[multi-pheno] ===== covariate fit for '" << models[mi].y_col << "' =====\n";
    Paths mp = paths;
    mp.out_prefix    = models[mi].out_prefix;
    mp.out_prefix_vr = models[mi].out_prefix_vr;
    ensure_parent_dir(mp.out_prefix + ".touch");
    ensure_parent_dir(mp.out_prefix_vr + ".touch");
    finish_design_stageB(designs[mi], mp);
  }

  // ===== Sample-set groups =====
  // Traits are grouped by their FINAL row set (after stage B, whose whitelist
  // filter is the last thing that can drop rows), still in CSV
  // order, so equal iid sequences <=> equal sample sets. Groups are ordered by
  // first appearance and list their traits in config order. With P=1, with
  // identical sample sets, or after intersect_samples this is one group and
  // everything below is the single-load path it always was.
  std::vector<std::vector<size_t>> groups;
  for (size_t k = 0; k < designs.size(); ++k) {
    size_t g = 0;
    for (; g < groups.size(); ++g) {
      const Design& d0 = designs[groups[g][0]];
      if (designs[k].n == d0.n && designs[k].iid == d0.iid) break;
    }
    if (g == groups.size()) groups.emplace_back();
    groups[g].push_back(k);
  }
  // Not const: fit.mask_missing may merge these base groups into fewer
  // "mask groups", each of which is one genotype load (see below).
  size_t G = groups.size();
  auto group_traits = [&](size_t g) {
    std::string s;
    for (size_t k : groups[g]) { if (!s.empty()) s += " "; s += models[k].y_col; }
    return s;
  };
  if (designs.size() > 1) {
    if (G == 1) {
      std::cout << "[multi-pheno] all " << designs.size()
                << " phenotypes share the same " << designs[0].n << " samples\n";
    } else {
      std::cout << "\n[multi-pheno] " << designs.size() << " phenotypes fall into "
                << G << " sample-set groups; each group gets its own genotype "
                   "load, GRM and GPU upload:\n";
      for (size_t g = 0; g < G; ++g)
        std::cout << "    group " << (g + 1) << "/" << G << ": n="
                  << designs[groups[g][0]].n << "  traits(" << groups[g].size()
                  << "): " << group_traits(g) << "\n";
    }
  }

  // LOCO ranges are computed inside PreprocessEngine::compute_chr_ranges_from_bim_()
  // (post-QC/compacted marker index space, which is what the genotype object
  // indexes in) and flow to NullModelEngine via PreOut::chr. There used to be a
  // duplicate raw-BIM scan here whose result was discarded; it has been removed.

  // FAM alignment (IID->1-based index), once per group.
  // FIXED: indicatorWithPheno should have N elements (FAM size), not design.n elements
  // This matches R's: indicatorGenoSamplesWithPheno = (sampleListwithGeno$IndexGeno %in% dataMerge_sort$IndexGeno)
  std::vector<std::vector<int>>  group_subSampleInGeno(G);
  std::vector<std::vector<bool>> group_indicatorWithPheno(G);
  int N_fam_total = 0;
  {
    auto fam_iids = read_fam_iids(paths.fam);
    int N_fam = static_cast<int>(fam_iids.size());
    N_fam_total = N_fam;

    // Create mapping from IID to FAM index (1-based)
    std::unordered_map<std::string,int> fam_pos; fam_pos.reserve(N_fam*2);
    for (int i=0;i<N_fam;++i) fam_pos.emplace(fam_iids[i], i+1); // 1-based

    for (size_t g = 0; g < G; ++g) {
    Design& design = designs[groups[g][0]];
    std::vector<int>&  subSampleInGeno    = group_subSampleInGeno[g];
    std::vector<bool>& indicatorWithPheno = group_indicatorWithPheno[g];

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

    // Reorder design + subSampleInGeno to match R's convention (ascending FAM row).
    // R's dataMerge_sort sorts the merged pheno/FAM table by IndexGeno. Without this,
    // design rows and genotype rows use the same indices but different physical samples
    // than R, which breaks any R-written bypass (random vectors, QR) when read by C++.
    {
      std::vector<int> perm(design.n);
      std::iota(perm.begin(), perm.end(), 0);
      std::sort(perm.begin(), perm.end(),
                [&](int a, int b) { return subSampleInGeno[a] < subSampleInGeno[b]; });

      // The permutation is derived from subSampleInGeno, which is built from
      // design.iid — identical across the phenotypes of one group — so one
      // perm serves every trait in the group. Field-for-field identical to the
      // single-trait code.
      for (size_t k : groups[g]) {
        Design& d = designs[k];
        std::vector<std::string> new_iid(d.n);
        std::vector<double> new_y(d.n);
        std::vector<double> new_X(static_cast<size_t>(d.n) * d.p);
        std::vector<double> new_offset;
        std::vector<double> new_event;
        if (!d.offset.empty())     new_offset.resize(d.n);
        if (!d.event_time.empty()) new_event.resize(d.n);

        for (int i = 0; i < d.n; ++i) {
          int s = perm[i];
          new_iid[i] = d.iid[s];
          new_y[i]   = d.y[s];
          for (int j = 0; j < d.p; ++j)
            new_X[static_cast<size_t>(i) * d.p + j] = d.X[static_cast<size_t>(s) * d.p + j];
          if (!d.offset.empty())     new_offset[i] = d.offset[s];
          if (!d.event_time.empty()) new_event[i]  = d.event_time[s];
        }
        d.iid = std::move(new_iid);
        d.y   = std::move(new_y);
        d.X   = std::move(new_X);
        if (!new_offset.empty()) d.offset     = std::move(new_offset);
        if (!new_event.empty())  d.event_time = std::move(new_event);
      }
      {
        std::vector<int> new_sub(design.n);
        for (int i = 0; i < design.n; ++i) new_sub[i] = subSampleInGeno[perm[i]];
        subSampleInGeno = std::move(new_sub);
      }

      if (G > 1) std::cout << "[FAM] group " << (g + 1) << "/" << G << ":\n";
      std::cout << "[FAM] reordered design/subSampleInGeno to ascending FAM row (matches R)\n";
      std::cout << "[FAM] subSampleInGeno[0:5] after reorder: ";
      for (int i = 0; i < std::min(5, design.n); ++i) std::cout << subSampleInGeno[i] << " ";
      std::cout << std::endl;
    }

    std::cout << "[FAM] N_fam=" << N_fam << ", design.n=" << design.n
              << ", indicatorWithPheno.size()=" << indicatorWithPheno.size() << std::endl;
    }
  }

  // ===== Scheme C: mask groups (optimization/missing_mt/SCHEME_C_DESIGN.md §6) =====
  // The sample-set groups above are the unit of "identical rows". fit.mask_missing
  // relaxes that: several of them share ONE genotype load over their union, and
  // each phenotype zeroes the union rows it does not own inside the
  // multiplication. Worth it only when the sets overlap heavily, because every
  // multiplication then reads |U| rows instead of |S_t| — mask_min_coverage is
  // the whole cost model.
  //
  // A mask group holding a single sample set is kept, not unwound: it is the
  // degenerate case (empty exclusion list, empty correction list, per-trait
  // stats equal to the union's) and it must reproduce the non-masked run.
  struct MaskInfo {
    bool                          on = false;
    std::vector<std::vector<int>> excl;      // per sample set: union rows NOT owned
    std::vector<std::vector<int>> scatter;   // per sample set: union rows owned
    std::vector<std::string>      names;     // per sample set
    std::vector<int>              member_bind;  // per trait of the group
  };
  std::vector<MaskInfo> unit_mask(G);
  if (cfg.mask_missing) {
    std::string why;
    if (!cfg.use_gpu)
      why = "fit.use_gpu is off — SCHEME_C_DESIGN.md §7 keeps the CPU path on grouping";
    else if (cfg.loco)
      why = "fit.loco is on — per-chromosome M_t and diagonals are the second cut (§7)";
    else if (cfg.use_sparse_grm_to_fit)
      why = "fit.use_sparse_grm_to_fit is on — the sparse fit never goes through the psi kernel (§7)";
    else if (cfg.use_sparse_grm_for_vr || cfg.make_sparse_grm_only)
      // Not in §7's list, but it has to be: the subset sparse GRM below is
      // built ONCE per group from designs[members[0]], because until now every
      // trait of a group had the same rows. A mask group holds several sample
      // sets, so that subset would be the wrong matrix (and the wrong
      // dimension) for every trait but the first.
      why = "a sparse GRM is loaded (use_sparse_grm_for_vr / make_sparse_grm_only) — "
            "it is subset once per group from one trait's sample set, which a mask "
            "group does not have";
    else if (cfg.use_blocked_gemv || cfg.gemv_verify)
      // parallelCrossProd_blocked is a CPU path over the packed store; under
      // masking it would read the union's rows with the trait's marker count.
      why = "fit.use_blocked_gemv / fit.gemv_verify select a CPU K.u path, and §7 keeps "
            "the CPU on grouping";
    if (!why.empty()) {
      std::cout << "[mask] fit.mask_missing requested but NOT used: " << why
                << ". Falling back to sample-set grouping.\n";
    } else {
      // Greedy (§6): largest sample set first; a set joins an existing mask
      // group only if every member still covers >= mask_min_coverage of the
      // resulting union.
      std::vector<size_t> order(G);
      std::iota(order.begin(), order.end(), size_t{0});
      std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return group_subSampleInGeno[a].size() > group_subSampleInGeno[b].size();
      });
      std::vector<std::vector<size_t>> units;    // base-group indices
      std::vector<std::vector<int>>    uni;      // FAM 1-based union, ascending
      for (size_t gi : order) {
        const std::vector<int>& sg = group_subSampleInGeno[gi];
        bool placed = false;
        for (size_t u = 0; u < units.size() && !placed; ++u) {
          std::vector<int> cand;
          cand.reserve(uni[u].size() + sg.size());
          std::set_union(uni[u].begin(), uni[u].end(), sg.begin(), sg.end(),
                         std::back_inserter(cand));
          double mincov = (double)sg.size() / (double)cand.size();
          for (size_t bg : units[u])
            mincov = std::min(mincov,
                              (double)group_subSampleInGeno[bg].size() / (double)cand.size());
          if (mincov >= cfg.mask_min_coverage) {
            units[u].push_back(gi);
            uni[u] = std::move(cand);
            placed = true;
          }
        }
        if (!placed) { units.push_back({gi}); uni.push_back(sg); }
      }
      // Order the units by their earliest phenotype so the log reads in config
      // order; the fit itself does not depend on it.
      std::vector<size_t> uorder(units.size());
      std::iota(uorder.begin(), uorder.end(), size_t{0});
      std::stable_sort(uorder.begin(), uorder.end(), [&](size_t a, size_t b) {
        auto first_of = [&](size_t u) {
          size_t f = designs.size();
          for (size_t bg : units[u]) f = std::min(f, groups[bg].front());
          return f;
        };
        return first_of(a) < first_of(b);
      });

      std::vector<std::vector<size_t>> new_groups;
      std::vector<std::vector<int>>    new_sub;
      std::vector<std::vector<bool>>   new_ind;
      std::vector<MaskInfo>            new_mask;
      for (size_t uu : uorder) {
        // Base groups of this unit, in config order.
        std::vector<size_t> bgs = units[uu];
        std::stable_sort(bgs.begin(), bgs.end(), [&](size_t a, size_t b) {
          return groups[a].front() < groups[b].front();
        });
        const std::vector<int>& U = uni[uu];
        std::unordered_map<int,int> pos;               // FAM 1-based -> union row
        pos.reserve(U.size() * 2);
        for (size_t i = 0; i < U.size(); ++i) pos.emplace(U[i], (int)i);

        MaskInfo mi;
        mi.on = true;
        std::vector<size_t> members;
        for (size_t bi = 0; bi < bgs.size(); ++bi) {
          const size_t bg = bgs[bi];
          const std::vector<int>& sg = group_subSampleInGeno[bg];
          std::vector<int> sc; sc.reserve(sg.size());
          for (int fam1 : sg) sc.push_back(pos.at(fam1));
          std::vector<char> own(U.size(), 0);
          for (int r : sc) own[(size_t)r] = 1;
          std::vector<int> ex;
          ex.reserve(U.size() - sc.size());
          for (size_t i = 0; i < U.size(); ++i) if (!own[i]) ex.push_back((int)i);
          mi.scatter.push_back(std::move(sc));
          mi.excl.push_back(std::move(ex));
          mi.names.push_back(group_traits(bg));
          for (size_t k : groups[bg]) { members.push_back(k); mi.member_bind.push_back((int)bi); }
        }
        std::vector<bool> ind(N_fam_total, false);
        for (int fam1 : U) ind[(size_t)fam1 - 1] = true;
        new_groups.push_back(std::move(members));
        new_sub.push_back(U);
        new_ind.push_back(std::move(ind));
        new_mask.push_back(std::move(mi));
      }
      groups                 = std::move(new_groups);
      group_subSampleInGeno  = std::move(new_sub);
      group_indicatorWithPheno = std::move(new_ind);
      unit_mask              = std::move(new_mask);
      G                      = groups.size();
      std::cout << "\n[mask] fit.mask_missing: " << designs.size()
                << " phenotypes over " << G << " mask group(s), min coverage "
                << cfg.mask_min_coverage << ":\n";
      for (size_t g = 0; g < G; ++g) {
        std::cout << "    mask group " << (g + 1) << "/" << G << ": union n="
                  << group_subSampleInGeno[g].size() << ", "
                  << unit_mask[g].scatter.size() << " sample set(s):";
        for (size_t b = 0; b < unit_mask[g].scatter.size(); ++b)
          std::cout << "  [" << unit_mask[g].names[b] << " n="
                    << unit_mask[g].scatter[b].size() << " cov="
                    << (double)unit_mask[g].scatter[b].size() /
                       (double)group_subSampleInGeno[g].size() << "]";
        std::cout << "\n";
      }
    }
  }

  Design& design = designs[0];
  const std::vector<int>& subSampleInGeno = group_subSampleInGeno[0];

  // ===== Step 18: Dry-run exit (validate inputs only, no genotype loading) =====
  if (cfg.dry_run) {
    std::cout << "\n============================================\n";
    std::cout << "=== DRY RUN: Input Validation Summary ===\n";
    std::cout << "============================================\n";
    std::cout << "Config:\n";
    std::cout << "  trait:          " << cfg.trait << "\n";
    std::cout << "  tol:            " << cfg.tol << "\n";
    std::cout << "  maxiter:        " << cfg.maxiter << "\n";
    std::cout << "  nthreads:       " << cfg.nthreads << "\n";
    std::cout << "  loco:           " << (cfg.loco ? "true" : "false") << "\n";
    std::cout << "  covariate_qr:   " << (cfg.covariate_qr ? "true" : "false") << "\n";
    std::cout << "  covariate_offset: " << (cfg.covariate_offset ? "true" : "false") << "\n";
    std::cout << "  inv_normalize:  " << (cfg.inv_normalize ? "true" : "false") << "\n";
    std::cout << "  use_sparse_grm: " << (cfg.use_sparse_grm_to_fit ? "true" : "false") << "\n";
    std::cout << "  num_markers_vr: " << cfg.num_markers_for_vr << "\n";
    std::cout << "Paths:\n";
    std::cout << "  bed:            " << paths.bed << (fs::exists(paths.bed) ? " [OK]" : " [MISSING]") << "\n";
    std::cout << "  bim:            " << paths.bim << (fs::exists(paths.bim) ? " [OK]" : " [MISSING]") << "\n";
    std::cout << "  fam:            " << paths.fam << (fs::exists(paths.fam) ? " [OK]" : " [MISSING]") << "\n";
    if (!paths.sparse_grm.empty())
      std::cout << "  sparse_grm:     " << paths.sparse_grm << (fs::exists(paths.sparse_grm) ? " [OK]" : " [MISSING]") << "\n";
    if (!paths.sparse_grm_ids.empty())
      std::cout << "  sparse_grm_ids: " << paths.sparse_grm_ids << (fs::exists(paths.sparse_grm_ids) ? " [OK]" : " [MISSING]") << "\n";
    std::cout << "  out_prefix:     " << paths.out_prefix << "\n";
    std::cout << "Design:\n";
    std::cout << "  n (samples):    " << design.n << "\n";
    std::cout << "  p (covariates): " << design.p << "\n";
    if (design.n > 0) {
      // Phenotype summary
      double y_min = *std::min_element(design.y.begin(), design.y.end());
      double y_max = *std::max_element(design.y.begin(), design.y.end());
      double y_sum = 0.0;
      for (auto v : design.y) y_sum += v;
      std::cout << "  y range:        [" << y_min << ", " << y_max << "]  mean=" << (y_sum / design.n) << "\n";
      if (ieq(cfg.trait, "binary")) {
        int n0 = 0, n1 = 0;
        for (auto v : design.y) { if (v == 0.0) ++n0; else if (v == 1.0) ++n1; }
        std::cout << "  binary counts:  0=" << n0 << "  1=" << n1 << "\n";
      }
      // First few IIDs
      std::cout << "  first IIDs:     ";
      for (int i = 0; i < std::min(5, design.n); ++i) std::cout << design.iid[i] << " ";
      std::cout << "\n";
    }
    std::cout << "FAM alignment:\n";
    std::cout << "  samples matched: " << subSampleInGeno.size() << " / " << design.n << "\n";
    if (G > 1) {
      std::cout << "Sample-set groups: " << G << "\n";
      for (size_t g = 0; g < G; ++g)
        std::cout << "  group " << (g + 1) << ": n=" << designs[groups[g][0]].n
                  << "  traits: " << group_traits(g) << "\n";
    }
    std::cout << "============================================\n";
    std::cout << "DRY RUN PASSED: all input validations succeeded.\n";
    std::cout << "============================================\n";
    return 0;
  }

  // ===== Step 16: Skip model fitting (R lines 1252-1254) =====
  if (cfg.skip_model_fitting) {
    if (cfg.model_file.empty()) {
      throw std::runtime_error(
          "skip_model_fitting=true but no model_file specified in config. "
          "Set fit.model_file to the path of an existing model output.");
    }
    if (!fs::exists(cfg.model_file)) {
      throw std::runtime_error(
          "skip_model_fitting=true but model_file does not exist: " + cfg.model_file);
    }
    std::cout << "[skip_model_fitting] Using existing model: " << cfg.model_file << "\n";
    std::cout << "[skip_model_fitting] NOTE: Loading pre-fitted models is not yet fully implemented.\n";
    std::cout << "  The model file exists and is accessible. To proceed with VR estimation\n";
    std::cout << "  from a pre-fitted model, full deserialization support is needed.\n";
    return 0;
  }

  // Sparse GRM: whether it is read from files or built from the genotypes is
  // decided ONCE, before any group runs — a solo run decides it before it
  // writes anything, so an earlier group must not flip a later group's choice
  // by writing the files.
  const bool need_sparse = (cfg.use_sparse_grm_to_fit || cfg.use_sparse_grm_for_vr);
  const bool sparse_have_files =
    !paths.sparse_grm.empty()     && fs::exists(paths.sparse_grm) &&
    !paths.sparse_grm_ids.empty() && fs::exists(paths.sparse_grm_ids);
  if (G > 1 && cfg.make_sparse_grm_only)
    throw std::runtime_error(
        "make_sparse_grm_only=true with " + std::to_string(G) + " sample-set groups: "
        "the sparse GRM is a property of ONE sample set. Build it once on the full "
        "cohort (single phenotype, or no missingness) and pass it via "
        "paths.sparse_grm / paths.sparse_grm_ids.");
  if (G > 1 && need_sparse && !sparse_have_files)
    throw std::runtime_error(
        "sparse GRM requested with " + std::to_string(G) + " sample-set groups but "
        "paths.sparse_grm / paths.sparse_grm_ids do not name existing files, so each "
        "group would build (and write) its own GRM on its own samples. Build the "
        "sparse GRM once on the full cohort (fit.make_sparse_grm_only) and pass the "
        "files; every group is then subset from them.");

  // Variance-ratio overwrite guard (checked for EVERY trait before any fitting,
  // so a P-trait run cannot die half way through after writing some outputs —
  // with several sample-set groups that means before the first group loads).
  // A make_sparse_grm_only run writes no variance ratio, and used to exit before
  // this guard, so it stays exempt.
  if (cfg.num_markers_for_vr > 0 && !cfg.make_sparse_grm_only) {
    bool allow_overwrite = (y["paths"] && y["paths"]["overwrite_varratio"]) ? y["paths"]["overwrite_varratio"].as<bool>() : false;
    for (const auto& m : models) {
      std::string vr_txt = m.out_prefix_vr + ".varianceRatio.txt";
      if (!allow_overwrite && fs::exists(vr_txt)) {
        std::ostringstream oss;
        oss << "Refusing to overwrite existing variance-ratio file: " << vr_txt
            << " (set paths.overwrite_varratio=true to allow).";
        throw std::runtime_error(oss.str());
      }
    }
  }

  // Propagate the AI-REML trace-estimator seed override (fit.trace_seed) to the
  // GetTrace / GetTrace_q RNG. -1 keeps the builtin per-trait defaults (10/200).
  setTraceSeed(cfg.trace_seed);
  if (cfg.trace_seed >= 0)
    std::cout << "[config] trace_seed override = " << cfg.trace_seed << "\n";

  {
    auto t = std::chrono::steady_clock::now();
    printf("[TIMER-MAIN] %-40s %8.2fs\n", "Design + preprocessing",
           std::chrono::duration<double>(t - T0).count());
  }

  // The sparse GRM file is the same for every group; parse it once (lazily, on
  // the first group) and subset it per group below.
  bool sparse_parsed = false;
  arma::umat sparse_loc; arma::vec sparse_val; int sparse_n_mtx = 0;
  std::vector<std::string> sparse_grm_ids;

  // ================= one pass per sample-set group =================
  // Everything from here to the end of the loop body depends on the group's
  // sample set: the genotype object (subSampleInGeno, allele frequencies,
  // invstd, QC'd marker list, VR marker pool, GRM diagonal, LOCO diagonals),
  // the subset sparse GRM, the GPU-resident matrix and the psi*U trace cache.
  // reset_step1_state_for_new_sample_set() tears all of that down between
  // groups (see its comment in SAIGE_step1_fast.cpp for the list).
  for (size_t gi = 0; gi < G; ++gi) {
  const std::vector<size_t>& members = groups[gi];
  Design& design = designs[members[0]];
  const std::vector<int>& subSampleInGeno = group_subSampleInGeno[gi];

  if (gi > 0) reset_step1_state_for_new_sample_set();
  if (G > 1) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "=== sample-set group " << (gi + 1) << "/" << G << ": n="
              << (unit_mask[gi].on ? (int)subSampleInGeno.size() : design.n)
              << (unit_mask[gi].on ? "  (union)" : "")
              << "  traits(" << members.size() << "): " << group_traits(gi) << "\n";
    std::cout << std::string(70, '=') << "\n";
  }

  // Match R: set isVarRatio=true so genotype loading excludes VR markers from GRM
  if (cfg.num_markers_for_vr > 0) {
    setminMAC_VarianceRatio(20.0f, -1.0f, true);
  }

  // Initialize genotype data BEFORE sparse GRM section (needed for build_sparse_grm_in_place)
  auto T1 = std::chrono::steady_clock::now();
  {
    // init_global_geno takes non-const refs; it copies them into the genotype object.
    std::vector<int>  sub_copy = subSampleInGeno;
    std::vector<bool> ind_copy = group_indicatorWithPheno[gi];
    if (unit_mask[gi].on) {
      set_scheme_c_break(cfg.scheme_c_break);
      // Scheme C: one decode over the union, per-trait stats / corrections /
      // VR pools rebuilt from it (SCHEME_C_DESIGN.md §1-§3, §7).
      init_global_geno_masked(paths.bed, paths.bim, paths.fam, sub_copy, ind_copy,
                              cfg.isDiagofKinSetAsOne, cfg.min_maf_grm, cfg.max_miss_grm,
                              unit_mask[gi].excl, unit_mask[gi].scatter,
                              unit_mask[gi].names);
    } else {
      init_global_geno(paths.bed, paths.bim, paths.fam, sub_copy, ind_copy, cfg.isDiagofKinSetAsOne, cfg.min_maf_grm, cfg.max_miss_grm);
    }
  }
  {
    auto t = std::chrono::steady_clock::now();
    printf("[TIMER-MAIN] %-40s %8.2fs\n", "Genotype loading (init_global_geno)",
           std::chrono::duration<double>(t - T1).count());
  }

  // -------- Sparse GRM build/reuse (+enforcement) --------
  // Need sparse GRM in memory when EITHER:
  //   - use_sparse_grm_to_fit: GLMM null model uses sparse Σ for the random effect
  //   - use_sparse_grm_for_vr: variance-ratio computation needs Σ_sparse⁻¹·g
  //                            (SAIGE-GENE+ cate-VR "sparse" row).
  // The setisUseSparseSigma* flags downstream are still only flipped when
  // use_sparse_grm_to_fit=TRUE; otherwise the dense path runs the GLMM fit and
  // only VR computation uses the loaded sparse Σ.
  if (need_sparse) {
    if (sparse_have_files) {
      if (!sparse_parsed) {
      auto T2 = std::chrono::steady_clock::now();
      load_matrix_market_coo(paths.sparse_grm, sparse_loc, sparse_val, sparse_n_mtx);
      {
        auto t = std::chrono::steady_clock::now();
        printf("[TIMER-MAIN] %-40s %8.2fs\n", "Sparse GRM file parse (MTX)",
               std::chrono::duration<double>(t - T2).count());
      }
        // Read one ID per line from the sparse GRM sample IDs file
        std::ifstream id_in(paths.sparse_grm_ids);
        if (!id_in) throw std::runtime_error("Failed to open sparse GRM IDs: " + paths.sparse_grm_ids);
        sparse_grm_ids.reserve(sparse_n_mtx);
        std::string line;
        while (std::getline(id_in, line)) {
          // trim whitespace
          size_t s = line.find_first_not_of(" \t\r\n");
          size_t e = line.find_last_not_of(" \t\r\n");
          if (s != std::string::npos) sparse_grm_ids.push_back(line.substr(s, e - s + 1));
        }
        if ((int)sparse_grm_ids.size() != sparse_n_mtx) {
          std::cerr << "[warning] Sparse GRM IDs count (" << sparse_grm_ids.size()
                    << ") differs from GRM dimension (" << sparse_n_mtx << ").\n";
        }
        sparse_parsed = true;
      }
      const arma::umat& loc = sparse_loc;
      const arma::vec&  val = sparse_val;
      const int n_mtx = sparse_n_mtx;
      const std::vector<std::string>& grm_ids = sparse_grm_ids;

      // ===== Subset sparse GRM to phenotyped samples =====
      // The loaded GRM may cover all samples in the cohort (e.g. 488K),
      // but we only need the subset that overlaps with our phenotyped samples
      // (this group's samples, in this group's row order).
      // Build GRM-index -> design-index map.
      {
        // Build map: design IID -> design index (0-based)
        std::unordered_map<std::string, int> design_pos;
        design_pos.reserve(design.n * 2);
        for (int i = 0; i < design.n; ++i) {
          design_pos[design.iid[i]] = i;
        }

        // Build map: GRM 0-based index -> design 0-based index (-1 if not phenotyped)
        std::vector<int> grm_to_design(n_mtx, -1);
        int n_overlap = 0;
        for (int i = 0; i < (int)grm_ids.size(); ++i) {
          auto it = design_pos.find(grm_ids[i]);
          if (it != design_pos.end()) {
            grm_to_design[i] = it->second;
            ++n_overlap;
          }
        }
        std::cout << "[sparse] GRM samples: " << n_mtx
                  << ", phenotyped: " << design.n
                  << ", overlap: " << n_overlap << "\n";
        if (n_overlap == 0) {
          throw std::runtime_error("No overlap between sparse GRM sample IDs and phenotyped samples. "
                                   "Check that sparse_grm_ids file contains matching sample IDs.");
        }

        // Filter COO entries: keep only entries where both row and col are phenotyped
        int total_nnz = (int)val.n_elem;
        std::vector<arma::uword> new_rows, new_cols;
        std::vector<double> new_vals;
        new_rows.reserve(total_nnz);
        new_cols.reserve(total_nnz);
        new_vals.reserve(total_nnz);

        for (int k = 0; k < total_nnz; ++k) {
          int r_old = (int)loc(0, k);
          int c_old = (int)loc(1, k);
          if (r_old < n_mtx && c_old < n_mtx &&
              grm_to_design[r_old] >= 0 && grm_to_design[c_old] >= 0) {
            new_rows.push_back((arma::uword)grm_to_design[r_old]);
            new_cols.push_back((arma::uword)grm_to_design[c_old]);
            new_vals.push_back(val(k));
          }
        }

        // Check that all phenotyped samples have a diagonal entry
        std::vector<bool> has_diag(design.n, false);
        for (size_t k = 0; k < new_rows.size(); ++k) {
          if (new_rows[k] == new_cols[k]) has_diag[new_rows[k]] = true;
        }
        for (int i = 0; i < design.n; ++i) {
          if (!has_diag[i]) {
            // Add identity diagonal for samples with no GRM entry
            new_rows.push_back((arma::uword)i);
            new_cols.push_back((arma::uword)i);
            new_vals.push_back(1.0);
          }
        }

        int sub_nnz = (int)new_rows.size();
        arma::umat sub_loc(2, sub_nnz);
        arma::vec sub_val(sub_nnz);
        for (int k = 0; k < sub_nnz; ++k) {
          sub_loc(0, k) = new_rows[k];
          sub_loc(1, k) = new_cols[k];
          sub_val(k) = new_vals[k];
        }

        std::cout << "[sparse] Subsetted GRM: " << design.n << "x" << design.n
                  << "  nnz=" << sub_nnz << " (from " << total_nnz << ")\n";

        setupSparseGRM(design.n, sub_loc, sub_val);
      }
      // Only flip the GLMM-fit sparse path when explicitly requested. When the
      // sparse GRM was loaded only for VR ("use_sparse_grm_for_vr"), we keep
      // the dense GLMM fit; VR computation reads locationMat / valueVec directly.
      if (cfg.use_sparse_grm_to_fit) {
        setisUseSparseSigmaforInitTau(true);
        setisUseSparseSigmaforNullModelFitting(true);
        std::cout << "[sparse] GRM reused for fit + VR: " << paths.sparse_grm
                  << "  n=" << design.n << "\n";
      } else {
        std::cout << "[sparse] GRM reused for VR only: " << paths.sparse_grm
                  << "  n=" << design.n << "\n";
      }
    } else {
      // Only reachable with G == 1 (refused above otherwise).
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
      if (cfg.use_sparse_grm_to_fit) {
        setisUseSparseSigmaforInitTau(true);
        setisUseSparseSigmaforNullModelFitting(true);
      }
    }
    // Only an error when we INTENDED to fit with sparse Sigma. When the sparse
    // GRM was loaded for VR only (use_sparse_grm_for_vr, dense GLMM fit), the
    // fitting flag is deliberately left off — that is not a failure.
    if (cfg.use_sparse_grm_to_fit && !get_isUseSparseSigmaforModelFitting()) {
      throw std::runtime_error("Sparse GRM-to-fit requested, but sparse-Sigma flag is off; aborting.");
    }
    setisUsePCGwithSparseSigma(cfg.use_pcg_with_sparse_grm);
    std::cout << "[sparse] use_pcg_with_sparse_grm=" << (cfg.use_pcg_with_sparse_grm ? "true" : "false")
              << (cfg.use_pcg_with_sparse_grm ? " (PCG solver)" : " (direct sparse solve, R default)") << "\n";
  }

  // Early exit: construct-only (G == 1 here; refused above otherwise)
  if (cfg.make_sparse_grm_only) {
    std::cout << "[ok] make_sparse_grm_only=true: exiting before null model fit.\n";
    return 0;
  }

  // ------------------ Tier-1: fit each phenotype on the shared genotype ------
  // The genotype object, the 2-bit GRM (and its GPU-resident copy) were built
  // once above for this group; every trait of the group reuses them. Nothing
  // inside fit_null carries state across traits: the VR marker order is a
  // fresh std::mt19937(200) and the Hutchinson probe stream is re-seeded on
  // every GetTrace/GetTrace_q entry, so trait k's numbers do not depend on
  // traits 0..k-1 having run first.
  // ------------------ Tier-2: lockstep the group's AI-REML loops ------------
  // fit.multi_lockstep advances all the group's AI-REML iterations together so
  // the fixed-effect PCG solves of every still-active trait are issued as ONE
  // batched multi-Sigma solve. Off by default: it is not bit-identical to the
  // per-trait path (psi*B reduces in a different order), and tier-1's
  // "a P>1 run reproduces each solo run exactly" property is worth keeping as
  // the default. A group with one trait always takes the per-trait path.
  // §3.4: masking and lockstep are mutually exclusive in the first cut —
  // lockstep fits several phenotypes at once and one global activate_trait()
  // cannot serve them. Say so; never turn it off silently.
  if (cfg.multi_lockstep && unit_mask[gi].on && members.size() > 1)
    std::cout << "[mask] fit.multi_lockstep is OFF for this group: masking fits one "
                 "phenotype at a time (SCHEME_C_DESIGN.md §3.4).\n";
  const bool use_lockstep =
      (cfg.multi_lockstep && members.size() > 1 && !unit_mask[gi].on);
  std::vector<FitNullResult> lockstep_fits;
  if (use_lockstep) {
    std::cout << "\n" << std::string(70, '#') << "\n";
    std::cout << "### lockstep multi-phenotype fit: P=" << members.size();
    if (G > 1) std::cout << " (group " << (gi + 1) << "/" << G << ")";
    std::cout << "\n" << std::string(70, '#') << "\n";
    std::vector<Paths> mpaths_all(members.size(), paths);
    // The lockstep driver takes the group's designs as one vector. Moved, not
    // copied: under lockstep designs[mi] is not read again after the fit.
    std::vector<Design> group_designs;
    group_designs.reserve(members.size());
    for (size_t k = 0; k < members.size(); ++k) {
      mpaths_all[k].out_prefix    = models[members[k]].out_prefix;
      mpaths_all[k].out_prefix_vr = models[members[k]].out_prefix_vr;
      group_designs.push_back(std::move(designs[members[k]]));
    }
    auto T_all = std::chrono::steady_clock::now();
    lockstep_fits = saige::fit_null_multi(cfg, mpaths_all, group_designs);
    auto t = std::chrono::steady_clock::now();
    const double s = std::chrono::duration<double>(t - T_all).count();
    printf("[TIMER-MAIN] lockstep fit_null_multi P=%zu %8.2fs (%.2fs/trait)\n",
           members.size(), s, s / (double)members.size());
  } else if (cfg.multi_lockstep && models.size() == 1) {
    std::cout << "[multi-pheno] fit.multi_lockstep requested with P=1 — "
                 "nothing to lockstep, using the per-trait path.\n";
  } else if (cfg.multi_lockstep) {
    std::cout << "[multi-pheno] fit.multi_lockstep: group " << (gi + 1) << "/" << G
              << " has one trait — nothing to lockstep, using the per-trait path.\n";
  }

  for (size_t k = 0; k < members.size(); ++k) {
    const size_t mi = members[k];
    const auto& m = models[mi];
    Paths mpaths = paths;
    mpaths.out_prefix    = m.out_prefix;
    mpaths.out_prefix_vr = m.out_prefix_vr;

    if (models.size() > 1) {
      std::cout << "\n" << std::string(70, '#') << "\n";
      std::cout << "### phenotype " << (mi + 1) << "/" << models.size()
                << ": " << m.y_col;
      if (G > 1) std::cout << "  (group " << (gi + 1) << "/" << G << ")";
      std::cout << "\n" << std::string(70, '#') << "\n";
    }
    auto T_ph = std::chrono::steady_clock::now();

    // Scheme C §3.4: point the genotype object (and the GPU bind, the psi*U
    // cache, the GRM-diagonal caches) at this phenotype. No-op when masking is
    // off, which is what keeps the non-masked path untouched.
    if (unit_mask[gi].on) activate_trait_for_fit(unit_mask[gi].member_bind[k]);

    FitNullResult out = use_lockstep ? std::move(lockstep_fits[k])
                                     : saige::fit_null(cfg, mpaths, designs[mi]);

    // ------------------ Output GRM diagonal (after fit_null, same as R version) ------------------
    output_grm_diagonal(mpaths.out_prefix + ".grm_diag.txt");

    // ------------------ Report artifacts ------------------
    std::cout << "== SAIGE Null Fit Completed ==\n";
    if (models.size() > 1) std::cout << "Phenotype: " << m.y_col << "\n";
    std::cout << "Converged: " << (out.converged ? "yes" : "NO") << "\n";
    std::cout << "Iterations: " << out.iterations << "\n";
    std::cout << "Model artifact: " << out.model_rda_path << "\n";
    if (!out.vr_path.empty())           std::cout << "Variance ratio: " << out.vr_path << "\n";
    if (!out.markers_out_path.empty())  std::cout << "Marker results: " << out.markers_out_path << "\n";
    // out.loco is set by the engine only when the LOCO batch actually ran, so this
    // no longer claims "on" for a run that quietly skipped LOCO.
    std::cout << "LOCO: " << (out.loco ? "on" : "off")
              << "  LowMem: " << (out.lowmem_loco ? "yes" : "no");
    if (out.loco) {
      std::cout << "  chroms:";
      for (int c : out.loco_chroms) std::cout << " " << c;
    }
    std::cout << "\n";
    if (models.size() > 1) {
      auto t = std::chrono::steady_clock::now();
      // Under lockstep the fit already happened in the one fit_null_multi call
      // above, so this window covers only the post-fit work; say so rather than
      // letting it read as a per-trait fit time.
      printf("[TIMER-MAIN] phenotype %-28s %8.2fs%s\n", m.y_col.c_str(),
             std::chrono::duration<double>(t - T_ph).count(),
             use_lockstep ? "  (post-fit only; fit is in fit_null_multi)" : "");
    }

    // Free this trait's design as soon as it is fitted: at P=32 on a big cohort
    // the N*(p_full+3) doubles per trait are not negligible next to the GRM.
    designs[mi] = Design{};
  }
  }  // end of the per-group loop

  // Diagnostic only: exercise the tier-2 lockstep multi-Sigma primitives against
  // the real psi now that the genotype object (and the GPU handle) are live
  // (the last group's, when there are several).
  // SAIGE_MULTISIGMA_SELFTEST=<P>, optionally SAIGE_MULTISIGMA_SELFTEST_NRHS=<k>.
  if (const char* e = std::getenv("SAIGE_MULTISIGMA_SELFTEST")) {
    const int P_test = std::max(1, std::atoi(e));
    const char* e2 = std::getenv("SAIGE_MULTISIGMA_SELFTEST_NRHS");
    const int nrhs  = e2 ? std::max(1, std::atoi(e2)) : 2;
    runMultiSigmaSelfTest(P_test, nrhs, cfg.maxiterPCG,
                          static_cast<float>(cfg.tolPCG));
  }
  return 0;
}
