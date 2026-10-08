// saige_cli.cpp -- the command-line front end of the C++/GPU SAIGE.
//
//   <prog> step1   [flags]     null model fitting   (runs saige-null)
//   <prog> step2   [flags]     association tests    (runs saige-step2)
//   <prog> sgs2txt [args]      .sgs -> text         (runs sgs2txt)
//
// Flags follow R SAIGE 1.5.2 (extdata/step1_fitNULLGLMM.R, step2_SPAtests.R).
// The front end computes nothing. It (1) loads --config if given, (2) applies
// the flags on top of it, (3) writes the resolved config to
// <outDir>/step1.yaml or <outDir>/step2.yaml, and (4) execs the engine on that
// file, so a CLI run and a hand-run of the written config are the same run.
//
// The program name comes from the build (Makefile variable PROGRAM).

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

#ifndef SAIGE_CLI_NAME
#define SAIGE_CLI_NAME "saige-gpu-cpp"
#endif
#ifndef SAIGE_CLI_VERSION
#define SAIGE_CLI_VERSION "unknown"
#endif
#ifndef SAIGE_CLI_BUILD
#define SAIGE_CLI_BUILD ""
#endif

namespace fs = std::filesystem;

namespace {

const char* PROG = SAIGE_CLI_NAME;

// Engine executables, installed in the same directory as this program.
const char* ENGINE_STEP1 = "saige-null";
const char* ENGINE_STEP2 = "saige-step2";
const char* ENGINE_SGS2TXT = "sgs2txt";

struct UsageError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// ------------------------------------------------------------------ strings
std::string trim(const std::string& s) {
  size_t i = 0, j = s.size();
  while (i < j && std::isspace((unsigned char)s[i])) ++i;
  while (j > i && std::isspace((unsigned char)s[j - 1])) --j;
  return s.substr(i, j - i);
}
std::string lower(std::string s) {
  for (auto& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}
bool starts_with(const std::string& s, const std::string& p) {
  return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}
bool ends_with(const std::string& s, const std::string& p) {
  return s.size() >= p.size() && s.compare(s.size() - p.size(), p.size(), p) == 0;
}
std::vector<std::string> split(const std::string& s, char d) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == d) { out.push_back(cur); cur.clear(); }
    else cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}
std::string shell_quote(const std::string& s) {
  if (!s.empty() && s.find_first_not_of(
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-.,/:=+@%") ==
        std::string::npos)
    return s;
  std::string o = "'";
  for (char c : s) { if (c == '\'') o += "'\\''"; else o += c; }
  return o + "'";
}
size_t edit_distance(const std::string& a, const std::string& b) {
  std::vector<size_t> prev(b.size() + 1), cur(b.size() + 1);
  for (size_t j = 0; j <= b.size(); ++j) prev[j] = j;
  for (size_t i = 1; i <= a.size(); ++i) {
    cur[0] = i;
    for (size_t j = 1; j <= b.size(); ++j)
      cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1,
                         prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1)});
    std::swap(prev, cur);
  }
  return prev[b.size()];
}

// ------------------------------------------------------------------ values
// R's optparse logical: TRUE/FALSE/T/F (any case). Also true/false, 1/0, yes/no.
bool parse_bool_lit(const std::string& v, bool& out) {
  const std::string l = lower(trim(v));
  if (l == "true" || l == "t" || l == "1" || l == "yes") { out = true; return true; }
  if (l == "false" || l == "f" || l == "0" || l == "no") { out = false; return true; }
  return false;
}
bool is_number(const std::string& v) {
  const std::string t = trim(v);
  if (t.empty()) return false;
  char* end = nullptr;
  errno = 0;
  std::strtod(t.c_str(), &end);
  return end && *end == '\0' && errno != ERANGE;
}
bool is_integer(const std::string& v) {
  const std::string t = trim(v);
  if (t.empty()) return false;
  char* end = nullptr;
  errno = 0;
  long long x = std::strtoll(t.c_str(), &end, 10);
  (void)x;
  return end && *end == '\0' && errno != ERANGE;
}

std::string abs_path(const std::string& p, const fs::path& base = fs::current_path()) {
  if (p.empty()) return p;
  std::string q = p;
  if (q == "~" || starts_with(q, "~/")) {
    const char* h = std::getenv("HOME");
    if (h) q = std::string(h) + q.substr(1);
  }
  fs::path P(q);
  if (!P.is_absolute()) P = base / P;
  return P.lexically_normal().string();
}

// ------------------------------------------------------------------ YAML paths
// Dotted access ("fit.nthreads") into nested maps; step 2 keys have no dots.
YAML::Node get_path(const YAML::Node& root, const std::string& dotted) {
  const auto parts = split(dotted, '.');
  YAML::Node cur = YAML::Clone(root);  // read-only walk on a copy: no zombie inserts
  for (const auto& k : parts) {
    if (!cur || !cur.IsMap()) return YAML::Node(YAML::NodeType::Undefined);
    YAML::Node next = cur[k];
    if (!next) return YAML::Node(YAML::NodeType::Undefined);
    cur.reset(next);
  }
  return cur;
}
bool has_path(const YAML::Node& root, const std::string& dotted) {
  YAML::Node n = get_path(root, dotted);
  return n.IsDefined() && !n.IsNull();
}
void set_path(YAML::Node& root, const std::string& dotted, const YAML::Node& value) {
  const auto parts = split(dotted, '.');
  YAML::Node node = root;  // handle to the same tree
  for (size_t i = 0; i + 1 < parts.size(); ++i) {
    YAML::Node next = node[parts[i]];
    if (!next || !next.IsMap()) {
      node[parts[i]] = YAML::Node(YAML::NodeType::Map);
      next.reset(node[parts[i]]);
    }
    node.reset(next);
  }
  node[parts.back()] = value;
}
void erase_path(YAML::Node& root, const std::string& dotted) {
  const auto parts = split(dotted, '.');
  YAML::Node node = root;
  for (size_t i = 0; i + 1 < parts.size(); ++i) {
    if (!node.IsMap() || !node[parts[i]]) return;
    YAML::Node next = node[parts[i]];
    node.reset(next);
  }
  if (node.IsMap() && node[parts.back()]) node.remove(parts.back());
}
std::string get_str(const YAML::Node& root, const std::string& dotted) {
  YAML::Node n = get_path(root, dotted);
  if (!n.IsDefined() || n.IsNull() || !n.IsScalar()) return "";
  return n.as<std::string>();
}

// ------------------------------------------------------------------ flags
enum class Kind { Bool, Str, Num, Int, StrList, NumList, Path };

struct Ctx;
using Apply = std::function<void(Ctx&, const std::string&)>;

struct Flag {
  std::string name;   // without the leading --
  Kind kind;
  std::string meta;   // value placeholder for --help
  std::string def;    // default shown in --help
  std::string help;
  bool ext;           // not an R SAIGE flag
  Apply apply;
  std::string key;    // what it sets in the config (for --help-markdown)
  std::string dkey;   // engine config key whose default --help shows ("" = def is the text)
};
struct Section {
  std::string title;
  std::vector<Flag> flags;
};
struct Rejected {
  std::string name;
  std::string why;
  std::string accept_if;  // value for which the flag is a harmless no-op ("" = never)
};

struct Ctx {
  int step = 0;
  YAML::Node cfg{YAML::NodeType::Map};
  std::string configPath;
  bool dryRun = false;
  std::vector<std::string> argvAll;
  std::set<std::string> given;  // flag names seen on the command line

  // shared
  std::string outDir;
  std::vector<std::string> phenoCols;
  int gpuDevice = -1;

  // step 1
  std::string outputPrefix, outputPrefixVR;
  std::string plinkFile;
  bool plinkTrioFromFlags = false;

  // step 2 genotype
  std::string s2plink, s2bed, s2bim, s2fam;
  std::string pgenPrefix, pgenFile, pvarFile, psamFile;
  std::string bgenFile, bgenIndex, sampleFile;
  std::string vcfFile, vcfIndex;
  // step 2 models
  std::string step1Dir, gmmat, vrFile, outFile;
  bool overwriteOutput = true;

  std::vector<std::string> notes;
};

// Typed setter: validates the value and writes it at `key`.
YAML::Node typed_node(Kind k, const std::string& flag, const std::string& v) {
  switch (k) {
    case Kind::Bool: {
      bool b;
      if (!parse_bool_lit(v, b))
        throw UsageError("--" + flag + " expects TRUE or FALSE, got '" + v + "'");
      return YAML::Node(b);
    }
    case Kind::Num:
      if (!is_number(v)) throw UsageError("--" + flag + " expects a number, got '" + v + "'");
      return YAML::Node(trim(v));
    case Kind::Int:
      if (!is_integer(v)) throw UsageError("--" + flag + " expects an integer, got '" + v + "'");
      return YAML::Node(trim(v));
    case Kind::Str:
      return YAML::Node(v);
    case Kind::Path:
      if (trim(v).empty()) throw UsageError("--" + flag + " expects a path");
      return YAML::Node(abs_path(trim(v)));
    case Kind::StrList: {
      YAML::Node seq(YAML::NodeType::Sequence);
      if (trim(v).empty()) return seq;
      for (auto& e : split(v, ',')) {
        const std::string t = trim(e);
        if (t.empty()) throw UsageError("--" + flag + ": empty element in '" + v + "'");
        seq.push_back(t);
      }
      return seq;
    }
    case Kind::NumList: {
      YAML::Node seq(YAML::NodeType::Sequence);
      for (auto& e : split(v, ',')) {
        const std::string t = trim(e);
        if (!is_number(t))
          throw UsageError("--" + flag + " expects comma-separated numbers, got '" + v + "'");
        seq.push_back(YAML::Node(t));
      }
      return seq;
    }
  }
  return YAML::Node();
}

Flag F(const std::string& name, Kind kind, const std::string& meta, const std::string& def,
       const std::string& help, const std::string& key, bool ext = false) {
  Flag f{name, kind, meta, def, help, ext, nullptr, key, key};
  f.apply = [key, kind, name](Ctx& c, const std::string& v) {
    set_path(c.cfg, key, typed_node(kind, name, v));
  };
  return f;
}
Flag FX(const std::string& name, Kind kind, const std::string& meta, const std::string& def,
        const std::string& help, Apply a, bool ext = false) {
  return Flag{name, kind, meta, def, help, ext, std::move(a), "", ""};
}

// Config keys of the flags that are not a plain "set this key".
void annotate(std::vector<Section>& S, const std::map<std::string, std::string>& keys,
              const std::map<std::string, std::string>& dkeys = {}) {
  for (auto& sec : S)
    for (auto& f : sec.flags) {
      auto it = keys.find(f.name);
      if (it != keys.end()) f.key = it->second;
      auto jt = dkeys.find(f.name);
      if (jt != dkeys.end()) f.dkey = jt->second;
    }
}

std::string need_value_str(const std::string& flag, const std::string& v) {
  if (trim(v).empty()) throw UsageError("--" + flag + " expects a value");
  return trim(v);
}
bool bool_of(const std::string& flag, const std::string& v) {
  bool b;
  if (!parse_bool_lit(v, b))
    throw UsageError("--" + flag + " expects TRUE or FALSE, got '" + v + "'");
  return b;
}

// Flags common to both steps.
void common_flags(std::vector<Section>& S, int step) {
  Section run{"Run", {}};
  run.flags.push_back(FX("config", Kind::Path, "FILE", "none",
      "YAML config to start from; flags given here override its keys",
      [](Ctx&, const std::string&) {}, true));
  run.flags.push_back(FX("outDir", Kind::Path, "DIR", "none",
      step == 1 ? "output directory (created): models/<trait>/, vr_<trait>.varianceRatio.txt, step1.yaml"
                : "output directory (created): <trait>.txt (or <trait>.txt.sgs), step2.yaml",
      [](Ctx& c, const std::string& v) { c.outDir = abs_path(need_value_str("outDir", v)); }, true));
  run.flags.push_back(FX("useGPU", Kind::Bool, "", "",
      step == 1 ? "GPU for the GRM-vector products (falls back to the CPU when no device)"
                : "GPU for the marker scan, SPA, Firth, ER (falls back to the CPU when not possible)",
      [step](Ctx& c, const std::string& v) {
        set_path(c.cfg, step == 1 ? "fit.use_gpu" : "useGPU", YAML::Node(bool_of("useGPU", v)));
      }, true));
  run.flags.push_back(FX("gpuDevice", Kind::Int, "N", step == 1 ? "not set" : "",
      step == 1 ? "CUDA device index (passed as SAIGE_GPU_DEVICE)" : "CUDA device index",
      [step](Ctx& c, const std::string& v) {
        if (!is_integer(v) || std::stoi(v) < 0)
          throw UsageError("--gpuDevice expects a device index >= 0, got '" + v + "'");
        c.gpuDevice = std::stoi(v);
        if (step == 2) set_path(c.cfg, "gpuDevice", YAML::Node(trim(v)));
      }, true));
  run.flags.push_back(F("nThreads", Kind::Int, "N", "", "CPU threads",
      step == 1 ? "fit.nthreads" : "nThreads"));
  run.flags.push_back(FX("set", Kind::Str, "KEY=VALUE", "none",
      step == 1 ? "set any config key directly, e.g. fit.selective_geno_load=false (repeatable)"
                : "set any config key directly, e.g. gpuSpaOrder=marker (repeatable)",
      [](Ctx& c, const std::string& v) {
        auto eq = v.find('=');
        if (eq == std::string::npos || eq == 0)
          throw UsageError("--set expects KEY=VALUE, got '" + v + "'");
        const std::string key = trim(v.substr(0, eq));
        YAML::Node val;
        try { val = YAML::Load(v.substr(eq + 1)); }
        catch (const std::exception& e) {
          throw UsageError("--set " + key + ": cannot parse value: " + e.what());
        }
        set_path(c.cfg, key, val);
      }, true));
  run.flags.push_back(FX("dryRun", Kind::Bool, "", "FALSE",
      "write the config and print the engine command, do not run it",
      [](Ctx& c, const std::string& v) { c.dryRun = bool_of("dryRun", v); }, true));
  S.push_back(run);
}

// ------------------------------------------------------------------ step 1
std::vector<Section> step1_sections() {
  std::vector<Section> S;
  Section in{"Input", {}};
  in.flags.push_back(FX("plinkFile", Kind::Path, "PREFIX", "none",
      "PLINK .bed/.bim/.fam prefix (full GRM and variance-ratio markers)",
      [](Ctx& c, const std::string& v) { c.plinkFile = abs_path(need_value_str("plinkFile", v)); }));
  for (const char* x : {"bed", "bim", "fam"}) {
    const std::string nm = std::string(x) + "File";
    const std::string key = std::string("paths.") + x;
    in.flags.push_back(FX(nm, Kind::Path, "FILE", "none",
        std::string("PLINK .") + x + " file (instead of --plinkFile)",
        [key](Ctx& c, const std::string& v) {
          set_path(c.cfg, key, YAML::Node(abs_path(trim(v))));
          c.plinkTrioFromFlags = true;
        }));
  }
  in.flags.push_back(F("phenoFile", Kind::Path, "FILE", "",
      "phenotype/covariate table (tab, space or comma separated, header line)", "design.csv"));
  in.flags.push_back(FX("phenoCol", Kind::StrList, "COL[,COL...]", "none",
      "phenotype column(s); a comma list fits several traits in one run (R: one)",
      [](Ctx& c, const std::string& v) {
        c.phenoCols.clear();
        for (auto& e : split(v, ',')) {
          const std::string t = trim(e);
          if (t.empty()) throw UsageError("--phenoCol: empty column name in '" + v + "'");
          c.phenoCols.push_back(t);
        }
      }));
  in.flags.push_back(F("sampleIDColinphenoFile", Kind::Str, "COL", "",
      "sample ID column of --phenoFile", "design.iid_col"));
  in.flags.push_back(F("covarColList", Kind::StrList, "COL[,COL...]", "",
      "covariate columns", "design.covar_cols"));
  in.flags.push_back(F("qCovarColList", Kind::StrList, "COL[,COL...]", "",
      "categorical covariates (must also be in --covarColList)", "design.q_covar_cols"));
  in.flags.push_back(F("SampleIDIncludeFile", Kind::Path, "FILE", "",
      "only these sample IDs (one per line) are used", "design.whitelist_ids"));
  in.flags.push_back(F("minCovariateCount", Kind::Int, "N", "",
      "binary covariates with fewer than N samples in a level are dropped (-1 = off)",
      "design.min_covariate_count"));
  in.flags.push_back(F("sexCol", Kind::Str, "COL", "", "sex column for --FemaleOnly/--MaleOnly",
      "design.sex_col"));
  in.flags.push_back(F("FemaleOnly", Kind::Bool, "", "", "fit female samples only",
      "design.female_only"));
  in.flags.push_back(F("MaleOnly", Kind::Bool, "", "", "fit male samples only",
      "design.male_only"));
  in.flags.push_back(F("FemaleCode", Kind::Str, "CODE", "", "female value in --sexCol",
      "design.female_code"));
  in.flags.push_back(F("MaleCode", Kind::Str, "CODE", "", "male value in --sexCol",
      "design.male_code"));
  S.push_back(in);

  Section md{"Model", {}};
  md.flags.push_back(FX("traitType", Kind::Str, "TYPE", "",
      "binary, quantitative or survival",
      [](Ctx& c, const std::string& v) {
        const std::string t = lower(trim(v));
        if (t != "binary" && t != "quantitative" && t != "survival")
          throw UsageError("--traitType must be binary, quantitative or survival, got '" + v + "'");
        set_path(c.cfg, "fit.trait", YAML::Node(t));
      }));
  md.flags.push_back(F("invNormalize", Kind::Bool, "", "",
      "inverse-normal transform a quantitative phenotype first", "fit.inv_normalize"));
  md.flags.push_back(FX("eventTimeCol", Kind::Str, "COL", "time",
      "survival event time column; must be named time, event_time or eventTime",
      [](Ctx&, const std::string& v) {
        const std::string t = lower(trim(v));
        if (t != "time" && t != "event_time" && t != "eventtime")
          throw UsageError("--eventTimeCol: this step 1 reads the event time from a column named "
                           "time, event_time or eventTime; rename column '" + v + "' to one of them");
      }));
  md.flags.push_back(F("eventTimeBinSize", Kind::Int, "N", "",
      "bin survival event times by this size", "fit.event_time_bin_size"));
  md.flags.push_back(F("LOCO", Kind::Bool, "", "",
      "leave-one-chromosome-out models (needs >= 2 autosomes)", "fit.loco"));
  md.flags.push_back(F("isLowMemLOCO", Kind::Bool, "", "",
      "low-memory LOCO", "fit.lowmem_loco"));
  md.flags.push_back(F("isCovariateTransform", Kind::Bool, "", "",
      "QR-transform the covariates", "fit.covariate_qr"));
  md.flags.push_back(F("isCovariateOffset", Kind::Bool, "", "",
      "fit covariate effects once and use them as an offset", "fit.covariate_offset"));
  md.flags.push_back(F("tol", Kind::Num, "X", "", "convergence tolerance of tau", "fit.tol"));
  md.flags.push_back(F("maxiter", Kind::Int, "N", "", "max AI-REML iterations", "fit.maxiter"));
  md.flags.push_back(F("tolPCG", Kind::Num, "X", "", "PCG tolerance", "fit.tolPCG"));
  md.flags.push_back(F("maxiterPCG", Kind::Int, "N", "", "max PCG iterations", "fit.maxiterPCG"));
  md.flags.push_back(F("nrun", Kind::Int, "N", "", "random vectors for the trace estimate",
      "fit.nrun"));
  md.flags.push_back(F("traceCVcutoff", Kind::Num, "X", "",
      "trace-estimate CV threshold", "fit.traceCVcutoff"));
  S.push_back(md);

  Section grm{"GRM", {}};
  grm.flags.push_back(F("minMAFforGRM", Kind::Num, "X", "", "min MAF of GRM markers",
      "fit.min_maf_grm"));
  grm.flags.push_back(F("maxMissingRateforGRM", Kind::Num, "X", "",
      "max missing rate of GRM markers", "fit.max_miss_grm"));
  grm.flags.push_back(F("isDiagofKinSetAsOne", Kind::Bool, "", "",
      "set the GRM diagonal to 1", "fit.diag_one"));
  grm.flags.push_back(F("sparseGRMFile", Kind::Path, "FILE", "",
      "sparse GRM (MatrixMarket); written with --makeSparseGRMOnly", "paths.sparse_grm"));
  grm.flags.push_back(F("sparseGRMSampleIDFile", Kind::Path, "FILE", "",
      "sample IDs of the sparse GRM, one per line", "paths.sparse_grm_ids"));
  grm.flags.push_back(F("useSparseGRMtoFitNULL", Kind::Bool, "", "",
      "fit the null model on the sparse GRM", "fit.use_sparse_grm_to_fit"));
  grm.flags.push_back(F("usePCGwithSparseGRM", Kind::Bool, "", "",
      "PCG instead of the direct solve on the sparse GRM", "fit.use_pcg_with_sparse_grm"));
  grm.flags.push_back(F("useSparseGRMforVarRatio", Kind::Bool, "", "",
      "sparse GRM for the variance ratio", "fit.use_sparse_grm_for_vr"));
  grm.flags.push_back(F("relatednessCutoff", Kind::Num, "X", "(R: 0)",
      "building a sparse GRM: entries below this are dropped (0 is read as 0.05)",
      "fit.relatedness_cutoff"));
  grm.flags.push_back(FX("makeSparseGRMOnly", Kind::Bool, "", "",
      "only build the sparse GRM into --sparseGRMFile/--sparseGRMSampleIDFile, then stop",
      [](Ctx& c, const std::string& v) {
        const bool b = bool_of("makeSparseGRMOnly", v);
        set_path(c.cfg, "fit.make_sparse_grm_only", YAML::Node(b));
        if (b) set_path(c.cfg, "fit.use_sparse_grm_to_fit", YAML::Node(true));
      }, true));
  S.push_back(grm);

  Section vr{"Variance ratio", {}};
  vr.flags.push_back(F("numRandomMarkerforVarianceRatio", Kind::Int, "N", "",
      "markers for the variance ratio", "fit.num_markers_for_vr"));
  vr.flags.push_back(FX("skipVarianceRatioEstimation", Kind::Bool, "", "FALSE",
      "no variance ratio (step 2 needs one)",
      [](Ctx& c, const std::string& v) {
        if (bool_of("skipVarianceRatioEstimation", v))
          set_path(c.cfg, "fit.num_markers_for_vr", YAML::Node(0));
      }));
  vr.flags.push_back(F("ratioCVcutoff", Kind::Num, "X", "",
      "variance-ratio CV threshold", "fit.ratio_cv_cutoff"));
  vr.flags.push_back(F("isCateVarianceRatio", Kind::Bool, "", "",
      "variance ratios per MAC category", "fit.isCateVarianceRatio"));
  vr.flags.push_back(F("cateVarRatioMinMACVecExclude", Kind::NumList, "X,X,...", "",
      "MAC category lower bounds (exclusive)", "fit.cateVarRatioMinMACVecExclude"));
  vr.flags.push_back(F("cateVarRatioMaxMACVecInclude", Kind::NumList, "X,...", "",
      "MAC category upper bounds (inclusive)", "fit.cateVarRatioMaxMACVecInclude"));
  vr.flags.push_back(F("includeNonautoMarkersforVarRatio", Kind::Bool, "", "",
      "allow non-autosomal variance-ratio markers", "fit.include_nonauto_for_vr"));
  vr.flags.push_back(F("IsOverwriteVarianceRatioFile", Kind::Bool, "", "",
      "overwrite an existing variance-ratio file", "paths.overwrite_varratio"));
  S.push_back(vr);

  Section st2{"Step-2 settings stored in the model", {}};
  st2.flags.push_back(F("SPAcutoff", Kind::Num, "X", "",
      "SPA is applied when |z| > X", "fit.spa_cutoff"));
  st2.flags.push_back(F("is_Firth_beta", Kind::Bool, "", "",
      "Firth beta for binary traits with p < --pCutoffforFirth (R: a step-2 flag)",
      "fit.firth_beta", true));
  st2.flags.push_back(F("pCutoffforFirth", Kind::Num, "X", "",
      "p-value cutoff for Firth (R: a step-2 flag)", "fit.p_cutoff_for_firth", true));
  st2.flags.push_back(F("is_fastTest", Kind::Bool, "", "(R step 2: FALSE)",
      "fast test for sparse-GRM models (R: a step-2 flag)", "fit.fast_test", true));
  st2.flags.push_back(FX("impute_method", Kind::Str, "M", "(R step 2: best_guess)",
      "missing genotypes in step 2: best_guess, mean or minor (R: a step-2 flag)",
      [](Ctx& c, const std::string& v) {
        const std::string t = trim(v);
        if (t != "best_guess" && t != "mean" && t != "minor")
          throw UsageError("--impute_method must be best_guess, mean or minor, got '" + v + "'");
        set_path(c.cfg, "fit.impute_method", YAML::Node(t));
      }, true));
  S.push_back(st2);

  Section out{"Output", {}};
  out.flags.push_back(FX("outputPrefix", Kind::Path, "PREFIX", "none",
      "R-style output prefix (instead of --outDir): model in PREFIX/, PREFIX.varianceRatio.txt",
      [](Ctx& c, const std::string& v) { c.outputPrefix = abs_path(need_value_str("outputPrefix", v)); }));
  out.flags.push_back(FX("outputPrefix_varRatio", Kind::Path, "PREFIX", "--outputPrefix",
      "variance-ratio prefix (with --outputPrefix)",
      [](Ctx& c, const std::string& v) {
        c.outputPrefixVR = abs_path(need_value_str("outputPrefix_varRatio", v));
      }));
  S.push_back(out);

  common_flags(S, 1);
  annotate(S, {
    {"config", "(the file is loaded first)"},
    {"outDir", "paths.out_prefix = DIR/models, paths.out_prefix_vr = DIR/vr"},
    {"useGPU", "fit.use_gpu"},
    {"gpuDevice", "environment SAIGE_GPU_DEVICE"},
    {"set", "KEY (dotted)"},
    {"dryRun", "-"},
    {"plinkFile", "paths.plinkFile"},
    {"bedFile", "paths.bed"}, {"bimFile", "paths.bim"}, {"famFile", "paths.fam"},
    {"phenoCol", "design.y_cols (design.y_col with --outputPrefix and one trait)"},
    {"traitType", "fit.trait"},
    {"eventTimeCol", "- (the column name is fixed)"},
    {"makeSparseGRMOnly", "fit.make_sparse_grm_only, fit.use_sparse_grm_to_fit"},
    {"skipVarianceRatioEstimation", "fit.num_markers_for_vr = 0"},
    {"outputPrefix", "paths.out_prefix"},
    {"outputPrefix_varRatio", "paths.out_prefix_vr"},
    {"impute_method", "fit.impute_method"},
  }, {
    {"useGPU", "fit.use_gpu"}, {"traitType", "fit.trait"},
    {"makeSparseGRMOnly", "fit.make_sparse_grm_only"}, {"impute_method", "fit.impute_method"},
  });
  return S;
}

std::vector<Rejected> step1_rejected() {
  return {
    {"skipModelFitting", "loading an existing model is not implemented in this step 1", "FALSE"},
    {"tauInit", "initial tau values cannot be set; the fit starts from R's default 0,0", "0,0"},
    {"memoryChunk", "the genotype store is not chunked; remove the flag", ""},
    {"pcgforUhatforSurvAnalysis", "not implemented in this step 1", "FALSE"},
  };
}

// ------------------------------------------------------------------ step 2
std::vector<Section> step2_sections() {
  std::vector<Section> S;
  Section in{"Genotypes (one format)", {}};
  in.flags.push_back(FX("plinkFile", Kind::Path, "PREFIX", "none", "PLINK .bed/.bim/.fam prefix",
      [](Ctx& c, const std::string& v) { c.s2plink = abs_path(need_value_str("plinkFile", v)); }, true));
  in.flags.push_back(FX("bedFile", Kind::Path, "FILE", "none", "PLINK .bed (with --bimFile, --famFile; same prefix)",
      [](Ctx& c, const std::string& v) { c.s2bed = abs_path(need_value_str("bedFile", v)); }));
  in.flags.push_back(FX("bimFile", Kind::Path, "FILE", "none", "PLINK .bim",
      [](Ctx& c, const std::string& v) { c.s2bim = abs_path(need_value_str("bimFile", v)); }));
  in.flags.push_back(FX("famFile", Kind::Path, "FILE", "none", "PLINK .fam",
      [](Ctx& c, const std::string& v) { c.s2fam = abs_path(need_value_str("famFile", v)); }));
  in.flags.push_back(FX("pgenPrefix", Kind::Path, "PREFIX", "none", "PGEN .pgen/.pvar/.psam prefix",
      [](Ctx& c, const std::string& v) { c.pgenPrefix = abs_path(need_value_str("pgenPrefix", v)); }));
  in.flags.push_back(FX("pgenFile", Kind::Path, "FILE", "none",
      "PGEN .pgen (.pvar/.psam default to the same prefix)",
      [](Ctx& c, const std::string& v) { c.pgenFile = abs_path(need_value_str("pgenFile", v)); }, true));
  in.flags.push_back(FX("pvarFile", Kind::Path, "FILE", "<pgen prefix>.pvar", "PGEN .pvar",
      [](Ctx& c, const std::string& v) { c.pvarFile = abs_path(need_value_str("pvarFile", v)); }, true));
  in.flags.push_back(FX("psamFile", Kind::Path, "FILE", "<pgen prefix>.psam", "PGEN .psam",
      [](Ctx& c, const std::string& v) { c.psamFile = abs_path(need_value_str("psamFile", v)); }, true));
  in.flags.push_back(FX("bgenFile", Kind::Path, "FILE", "none", "BGEN 1.2 file (needs --sampleFile)",
      [](Ctx& c, const std::string& v) { c.bgenFile = abs_path(need_value_str("bgenFile", v)); }));
  in.flags.push_back(FX("bgenFileIndex", Kind::Path, "FILE", "<bgenFile>.bgi",
      "BGEN index; read from <bgenFile>.bgi only, so it must be that file",
      [](Ctx& c, const std::string& v) { c.bgenIndex = abs_path(need_value_str("bgenFileIndex", v)); }));
  in.flags.push_back(FX("sampleFile", Kind::Path, "FILE", "none", "BGEN .sample file",
      [](Ctx& c, const std::string& v) { c.sampleFile = abs_path(need_value_str("sampleFile", v)); }));
  in.flags.push_back(FX("vcfFile", Kind::Path, "FILE", "none", "VCF / VCF.GZ / BCF",
      [](Ctx& c, const std::string& v) { c.vcfFile = abs_path(need_value_str("vcfFile", v)); }));
  in.flags.push_back(FX("vcfFileIndex", Kind::Path, "FILE", "none",
      "accepted; not needed (the VCF is read from start to end)",
      [](Ctx& c, const std::string& v) { c.vcfIndex = abs_path(need_value_str("vcfFileIndex", v)); }));
  in.flags.push_back(FX("vcfField", Kind::Str, "DS|GT", "", "VCF field to read",
      [](Ctx& c, const std::string& v) {
        const std::string t = trim(v);
        if (t != "DS" && t != "GT") throw UsageError("--vcfField must be DS or GT, got '" + v + "'");
        set_path(c.cfg, "vcfField", YAML::Node(t));
      }));
  in.flags.push_back(FX("AlleleOrder", Kind::Str, "ORDER", "(PLINK; PGEN and BGEN: ref-first)",
      "alt-first or ref-first",
      [](Ctx& c, const std::string& v) {
        const std::string t = trim(v);
        if (t != "alt-first" && t != "ref-first")
          throw UsageError("--AlleleOrder must be alt-first or ref-first, got '" + v + "'");
        set_path(c.cfg, "AlleleOrder", YAML::Node(t));
      }));
  in.flags.push_back(F("chrom", Kind::Str, "CHR", "",
      "chromosome being tested (required with --LOCO=TRUE)", "chrom"));
  in.flags.push_back(F("is_imputed_data", Kind::Bool, "", "",
      "dosages are imputed (adds imputation columns)", "isImputation"));
  S.push_back(in);

  Section md{"Models and output", {}};
  md.flags.push_back(FX("step1Dir", Kind::Path, "DIR", "none",
      "a step-1 --outDir: tests every trait fitted there (needs --outDir)",
      [](Ctx& c, const std::string& v) { c.step1Dir = abs_path(need_value_str("step1Dir", v)); }, true));
  md.flags.push_back(FX("phenoCol", Kind::StrList, "T[,T...]", "all",
      "with --step1Dir: test only these traits",
      [](Ctx& c, const std::string& v) {
        c.phenoCols.clear();
        for (auto& e : split(v, ',')) {
          const std::string t = trim(e);
          if (t.empty()) throw UsageError("--phenoCol: empty trait name in '" + v + "'");
          c.phenoCols.push_back(t);
        }
      }, true));
  md.flags.push_back(FX("GMMATmodelFile", Kind::Path, "DIR", "none",
      "one step-1 model directory (instead of --step1Dir)",
      [](Ctx& c, const std::string& v) { c.gmmat = abs_path(need_value_str("GMMATmodelFile", v)); }));
  md.flags.push_back(FX("varianceRatioFile", Kind::Path, "FILE", "none",
      "its variance-ratio file (with --GMMATmodelFile)",
      [](Ctx& c, const std::string& v) { c.vrFile = abs_path(need_value_str("varianceRatioFile", v)); }));
  md.flags.push_back(FX("SAIGEOutputFile", Kind::Path, "FILE", "<outDir>/<model name>.txt",
      "result file (with --GMMATmodelFile)",
      [](Ctx& c, const std::string& v) { c.outFile = abs_path(need_value_str("SAIGEOutputFile", v)); }));
  md.flags.push_back(FX("outputFormat", Kind::Str, "text|sgs", "",
      "sgs: binary columnar output, convert with `sgs2txt` (multi-trait or --useGPU runs)",
      [](Ctx& c, const std::string& v) {
        const std::string t = trim(v);
        if (t != "text" && t != "sgs") throw UsageError("--outputFormat must be text or sgs, got '" + v + "'");
        set_path(c.cfg, "outputFormat", YAML::Node(t));
      }, true));
  md.flags.push_back(FX("sgsPrecision", Kind::Str, "fp64|fp32", "",
      "width of the .sgs floating-point columns (fp32 is not exact)",
      [](Ctx& c, const std::string& v) {
        const std::string t = trim(v);
        if (t != "fp64" && t != "fp32") throw UsageError("--sgsPrecision must be fp64 or fp32, got '" + v + "'");
        set_path(c.cfg, "sgsPrecision", YAML::Node(t));
      }, true));
  md.flags.push_back(F("is_output_moreDetails", Kind::Bool, "", "",
      "extra output columns", "isMoreOutput"));
  md.flags.push_back(FX("is_overwrite_output", Kind::Bool, "", "TRUE",
      "FALSE: refuse to run when a result file exists",
      [](Ctx& c, const std::string& v) { c.overwriteOutput = bool_of("is_overwrite_output", v); }));
  md.flags.push_back(F("markers_per_chunk", Kind::Int, "N", "",
      "markers per progress report", "marker_chunksize"));
  S.push_back(md);

  Section qc{"Marker QC", {}};
  qc.flags.push_back(F("minMAF", Kind::Num, "X", "", "min minor allele frequency", "minMAF"));
  qc.flags.push_back(F("minMAC", Kind::Num, "X", "", "min minor allele count", "minMAC"));
  qc.flags.push_back(F("maxMissing", Kind::Num, "X", "", "max missing rate", "maxMissRate"));
  qc.flags.push_back(F("minInfo", Kind::Num, "X", "", "min imputation info", "minINFO"));
  qc.flags.push_back(F("dosage_zerod_cutoff", Kind::Num, "X", "",
      "dosages <= X are set to 0 for markers with MAC <= --dosage_zerod_MAC_cutoff",
      "dosage_zerod_cutoff"));
  qc.flags.push_back(F("dosage_zerod_MAC_cutoff", Kind::Num, "X", "",
      "see --dosage_zerod_cutoff", "dosage_zerod_MAC_cutoff"));
  S.push_back(qc);

  Section ts{"Tests", {}};
  ts.flags.push_back(F("LOCO", Kind::Bool, "", "(R: TRUE)",
      "use the LOCO models of --chrom", "LOCO"));
  ts.flags.push_back(FX("is_Firth_beta", Kind::Bool, "", "",
      "Firth beta for binary traits with p < --pCutoffforFirth",
      [](Ctx& c, const std::string& v) {
        const bool b = bool_of("is_Firth_beta", v);
        set_path(c.cfg, "is_Firth_beta", YAML::Node(b));
        set_path(c.cfg, "isFirth", YAML::Node(b));
      }));
  ts.flags.push_back(F("pCutoffforFirth", Kind::Num, "X", "",
      "p-value cutoff for Firth", "pCutoffforFirth"));
  ts.flags.push_back(F("max_MAC_for_ER", Kind::Num, "X", "",
      "binary traits: exact test for MAC <= X", "MACCutoffforER"));
  ts.flags.push_back(F("is_noadjCov", Kind::Bool, "", "(R: TRUE)",
      "skip the covariate adjustment of the score", "isnoadjCov"));
  ts.flags.push_back(F("relatednessCutoff", Kind::Num, "X", "",
      "sparse GRM entries below this are dropped", "relatednessCutoff"));
  ts.flags.push_back(F("cateVarRatioMinMACVecExclude", Kind::NumList, "X,X,...", "",
      "MAC category lower bounds", "cateVarRatioMinMACVecExclude"));
  ts.flags.push_back(F("cateVarRatioMaxMACVecInclude", Kind::NumList, "X,...", "",
      "MAC category upper bounds", "cateVarRatioMaxMACVecInclude"));
  ts.flags.push_back(F("condition", Kind::Str, "ID[,ID...]", "",
      "conditioning marker IDs", "condition"));
  ts.flags.push_back(F("weights_for_condition", Kind::NumList, "X,X,...", "",
      "weights of the conditioning markers", "weights_for_condition"));
  S.push_back(ts);

  Section rg{"Region / gene tests", {}};
  rg.flags.push_back(F("groupFile", Kind::Path, "FILE", "",
      "group file; turns on region tests (one trait per run)", "groupFile"));
  rg.flags.push_back(F("annotation_in_groupTest", Kind::StrList, "A[,A;B...]",
      "(R: lof,missense;lof,missense;lof;synonymous)",
      "annotation sets, ',' between sets and ';' inside one", "annotationList"));
  rg.flags.push_back(F("maxMAF_in_groupTest", Kind::NumList, "X,X,...", "(R: 0.0001,0.001,0.01)",
      "max MAF cutoffs", "maxMAFList"));
  rg.flags.push_back(F("r.corr", Kind::Num, "0|1", "", "0: SKAT-O (+SKAT, Burden); 1: Burden only",
      "r_corr"));
  rg.flags.push_back(F("weights.beta", Kind::NumList, "A,B", "", "Beta(MAF; A, B) weights",
      "weights_beta"));
  rg.flags.push_back(F("MACCutoff_to_CollapseUltraRare", Kind::Num, "X", "",
      "collapse variants with MAC <= X", "MACCutoff_to_CollapseUltraRare"));
  rg.flags.push_back(F("markers_per_chunk_in_groupTest", Kind::Int, "N", "(R: 100)",
      "markers per chunk", "markers_per_chunk_in_groupTest"));
  rg.flags.push_back(F("groups_per_chunk", Kind::Int, "N", "", "regions per I/O chunk",
      "groups_per_chunk"));
  rg.flags.push_back(F("is_single_in_groupTest", Kind::Bool, "", "(R: FALSE)",
      "also write single-variant results", "isSingleInGroupTest"));
  rg.flags.push_back(F("is_output_markerList_in_groupTest", Kind::Bool, "", "",
      "write the marker list of each region", "isOutputMarkerList"));
  rg.flags.push_back(F("minGroupMAC_in_BurdenTest", Kind::Num, "X", "",
      "min group MAC for Burden-only tests", "min_gourpmac_for_burdenonly"));
  S.push_back(rg);

  common_flags(S, 2);
  annotate(S, {
    {"config", "(the file is loaded first)"},
    {"outDir", "models[].outputFile = DIR/<trait>.txt"},
    {"useGPU", "useGPU"},
    {"gpuDevice", "gpuDevice"},
    {"set", "KEY"},
    {"dryRun", "-"},
    {"plinkFile", "genoType: plink, plinkFile"},
    {"bedFile", "plinkFile (the shared prefix)"}, {"bimFile", "plinkFile"}, {"famFile", "plinkFile"},
    {"pgenPrefix", "genoType: pgen, pgenFile, pvarFile, psamFile"},
    {"pgenFile", "genoType: pgen, pgenFile"}, {"pvarFile", "pvarFile"}, {"psamFile", "psamFile"},
    {"bgenFile", "genoType: bgen, bgenFile"},
    {"bgenFileIndex", "- (checked)"},
    {"sampleFile", "bgenSampleFile"},
    {"vcfFile", "genoType: vcf, vcfFile"},
    {"vcfFileIndex", "- (not used)"},
    {"vcfField", "vcfField"},
    {"AlleleOrder", "AlleleOrder"},
    {"step1Dir", "models[]: traitName, modelFile, varianceRatioFile"},
    {"phenoCol", "models[] (subset)"},
    {"GMMATmodelFile", "modelFile"},
    {"varianceRatioFile", "varianceRatioFile"},
    {"SAIGEOutputFile", "outputFile"},
    {"outputFormat", "outputFormat"},
    {"sgsPrecision", "sgsPrecision"},
    {"is_overwrite_output", "- (checked)"},
    {"is_Firth_beta", "is_Firth_beta, isFirth"},
  }, {
    {"useGPU", "useGPU"}, {"gpuDevice", "gpuDevice"}, {"vcfField", "vcfField"},
    {"AlleleOrder", "AlleleOrder"}, {"outputFormat", "outputFormat"},
    {"sgsPrecision", "sgsPrecision"}, {"is_Firth_beta", "is_Firth_beta"},
  });
  return S;
}

std::vector<Rejected> step2_rejected() {
  const std::string s1 = std::string("it is stored in the model; give it to `") + PROG + " step1`";
  return {
    {"SPAcutoff", s1 + " --SPAcutoff", ""},
    {"is_fastTest", s1 + " --is_fastTest", ""},
    {"impute_method", s1 + " --impute_method", ""},
    {"sparseGRMFile", "the sparse GRM comes from the step-1 model (fit it with `" +
        std::string(PROG) + " step1 --sparseGRMFile ...`)", ""},
    {"sparseGRMSampleIDFile", "the sparse GRM comes from the step-1 model", ""},
    {"vcfFilters", "not implemented", ""},
    {"savFile", "SAV input is not supported; convert to VCF, BGEN or PGEN", ""},
    {"savFileIndex", "SAV input is not supported", ""},
    {"idstoIncludeFile", "not implemented; subset the genotype file first (plink2 --extract)", ""},
    {"rangestoIncludeFile", "not implemented; subset the genotype file first (plink2 --extract range)", ""},
    {"sampleFile_male", "X-chromosome male handling is not implemented", ""},
    {"X_PARregion", "X-chromosome male handling is not implemented", ""},
    {"is_rewrite_XnonPAR_forMales", "X-chromosome male handling is not implemented", "FALSE"},
    {"subSampleFile", "not implemented", ""},
    {"maxMAC_in_groupTest", "a MAC cap per group test is not implemented", "0"},
    {"is_no_weight_in_groupTest", "unweighted group tests are not implemented", "FALSE"},
  };
}

// ------------------------------------------------------------------ help
std::string self_dir();

// Defaults live in the engines only: --help asks them (`saige-null
// --print-defaults`, `saige-step2 --print-defaults`, each printing
// "default<TAB>key<TAB>value" from its own config parsing) and shows the value
// of the config key the flag sets. The front end writes no default into a config.
struct EngineDefaults {
  std::map<std::string, std::string> kv;
  std::string error;  // why there are none
};
EngineDefaults engine_defaults(int step) {
  EngineDefaults d;
  const std::string exe = self_dir() + "/" + (step == 1 ? ENGINE_STEP1 : ENGINE_STEP2);
  if (::access(exe.c_str(), X_OK) != 0) { d.error = exe + " not found"; return d; }
  FILE* p = ::popen((shell_quote(exe) + " --print-defaults 2>/dev/null").c_str(), "r");
  if (!p) { d.error = "cannot run " + exe; return d; }
  char buf[4096];
  while (std::fgets(buf, sizeof buf, p)) {
    std::string l(buf);
    while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back();
    if (!starts_with(l, "default\t")) continue;
    auto t = split(l, '\t');
    if (t.size() >= 3) d.kv[t[1]] = t[2];
  }
  ::pclose(p);
  if (d.kv.empty()) d.error = exe + " --print-defaults gave nothing";
  return d;
}
std::string shown_default(const Flag& f, const EngineDefaults& d) {
  if (f.dkey.empty()) return f.def;
  auto it = d.kv.find(f.dkey);
  std::string v;
  if (it == d.kv.end())
    v = "engine default (" + (d.error.empty() ? "key " + f.dkey + " not reported" : d.error) + ")";
  else {
    v = it->second;
    if (v == "true") v = "TRUE";
    else if (v == "false") v = "FALSE";
  }
  return f.def.empty() ? v : v + " " + f.def;
}

std::string flag_usage(const Flag& f) {
  std::string s = "--" + f.name;
  if (f.kind == Kind::Bool) s += "[=TRUE|FALSE]";
  else s += " " + f.meta;
  return s;
}

void print_help(int step) {
  const auto S = step == 1 ? step1_sections() : step2_sections();
  const EngineDefaults D = engine_defaults(step);
  const auto R = step == 1 ? step1_rejected() : step2_rejected();
  std::ostream& o = std::cout;
  if (step == 1) {
    o << "Usage: " << PROG << " step1 --plinkFile PREFIX --phenoFile FILE --phenoCol COL[,COL...] "
         "--outDir DIR [flags]\n\n"
         "Fits one null GLMM per phenotype column (step 1). Writes <outDir>/step1.yaml (the\n"
         "resolved config) and runs `" << ENGINE_STEP1 << " -c <outDir>/step1.yaml` on it.\n"
         "Outputs: <outDir>/models/<trait>/ and <outDir>/vr_<trait>.varianceRatio.txt.\n";
  } else {
    o << "Usage: " << PROG << " step2 --step1Dir DIR --outDir DIR (--plinkFile PREFIX | --pgenPrefix PREFIX |\n"
         "         --bgenFile FILE --sampleFile FILE | --vcfFile FILE) [flags]\n\n"
         "Tests every marker against the step-1 models (step 2). Writes <outDir>/step2.yaml\n"
         "(the resolved config) and runs `" << ENGINE_STEP2 << " <outDir>/step2.yaml` on it.\n"
         "Outputs: <outDir>/<trait>.txt per trait.\n";
  }
  o << "\nFlag names are R SAIGE's (1.5.2). [+] marks flags R does not have. A logical flag\n"
       "takes --flag=TRUE, --flag=FALSE, --flag TRUE or a bare --flag (= TRUE). Values: --flag=V\n"
       "or --flag V. Flags override the keys of --config. A flag that is not given is left out\n"
       "of the written config, so the engine's own default applies (shown below, as the engine\n"
       "reports it; where R's default differs it is noted).\n";
  for (const auto& sec : S) {
    o << "\n" << sec.title << ":\n";
    for (const auto& f : sec.flags) {
      std::string u = flag_usage(f);
      std::string line = "  " + u;
      const size_t col = 40;
      if (line.size() + 1 < col) line += std::string(col - line.size(), ' ');
      else line += "\n" + std::string(col, ' ');
      line += (f.ext ? "[+] " : "") + f.help + "  [default: " + shown_default(f, D) + "]";
      o << line << "\n";
    }
  }
  o << "\nR SAIGE flags that are refused (the run stops with the reason):\n";
  for (const auto& r : R) {
    o << "  --" << r.name << ": " << r.why;
    if (!r.accept_if.empty()) o << " (--" << r.name << "=" << r.accept_if << " is accepted)";
    o << "\n";
  }
  o << "\nOther:\n  -h, --help                            this text\n"
       "  --help-markdown                       the flag tables as Markdown\n";
}

std::string md_escape(std::string s) {
  std::string o;
  for (char c : s) { if (c == '|') o += "\\|"; else o += c; }
  return o;
}

// The flag tables of docs/step1.md and docs/step2.md.
void print_help_markdown(int step) {
  const auto S = step == 1 ? step1_sections() : step2_sections();
  const EngineDefaults D = engine_defaults(step);
  const auto R = step == 1 ? step1_rejected() : step2_rejected();
  std::ostream& o = std::cout;
  for (const auto& sec : S) {
    o << "\n**" << sec.title << "**\n\n| Flag | Config key | Default | Meaning |\n|---|---|---|---|\n";
    for (const auto& f : sec.flags)
      o << "| `" << md_escape(flag_usage(f)) << "`" << (f.ext ? " [+]" : "") << " | "
        << (f.key.empty() ? "" : (f.key[0] == '-' || f.key[0] == '(') ? md_escape(f.key)
                                                                       : "`" + md_escape(f.key) + "`")
        << " | " << md_escape(shown_default(f, D))
        << " | " << md_escape(f.help) << " |\n";
  }
  o << "\n**R SAIGE flags that are refused**\n\n| Flag | Why |\n|---|---|\n";
  for (const auto& r : R) {
    o << "| `--" << r.name << "` | " << md_escape(r.why);
    if (!r.accept_if.empty()) o << " (`--" << r.name << "=" << r.accept_if << "` is accepted)";
    o << " |\n";
  }
}

void print_main_help() {
  std::cout <<
    "Usage: " << PROG << " <command> [flags]\n\n"
    "Commands:\n"
    "  step1     fit the null models            (" << PROG << " step1 --help)\n"
    "  step2     test the markers               (" << PROG << " step2 --help)\n"
    "  sgs2txt   convert .sgs output to text    (" << PROG << " sgs2txt --help)\n\n"
    "  --version  print the version\n"
    "  --help     this text\n\n"
    "Every step1/step2 run writes its resolved YAML config into the output directory and\n"
    "runs the engine (" << ENGINE_STEP1 << ", " << ENGINE_STEP2 << ") on it; the engines are\n"
    "installed next to this program and can be run on a YAML config directly.\n";
}

// ------------------------------------------------------------------ parsing
struct Parsed {
  std::vector<std::pair<const Flag*, std::string>> items;
};

const Flag* find_flag(const std::vector<Section>& S, const std::string& name) {
  for (const auto& s : S)
    for (const auto& f : s.flags)
      if (f.name == name) return &f;
  return nullptr;
}

std::string suggestion(const std::vector<Section>& S, const std::string& name) {
  std::string best;
  size_t bd = 1000;
  for (const auto& s : S)
    for (const auto& f : s.flags) {
      if (lower(f.name) == lower(name)) return f.name;
      size_t d = edit_distance(lower(f.name), lower(name));
      if (d < bd) { bd = d; best = f.name; }
    }
  return bd <= 3 ? best : "";
}

Parsed parse_args(int step, const std::vector<std::string>& args, const std::vector<Section>& S,
                  const std::vector<Rejected>& R) {
  Parsed P;
  const std::string sub = step == 1 ? "step1" : "step2";
  for (size_t i = 0; i < args.size(); ++i) {
    const std::string& a = args[i];
    if (a == "-h" || a == "--help") { print_help(step); std::exit(0); }
    if (a == "--help-markdown") { print_help_markdown(step); std::exit(0); }
    if (!starts_with(a, "--") || a.size() == 2)
      throw UsageError("unexpected argument '" + a + "' (flags start with --)");
    std::string name = a.substr(2), val;
    bool hasEq = false;
    auto eq = name.find('=');
    if (eq != std::string::npos) { val = name.substr(eq + 1); name = name.substr(0, eq); hasEq = true; }

    const Flag* f = find_flag(S, name);
    if (!f) {
      for (const auto& r : R) {
        if (r.name != name) continue;
        // value of a refused flag (it may come as the next token)
        std::string rv = val;
        if (!hasEq && i + 1 < args.size() && !starts_with(args[i + 1], "--")) rv = args[++i];
        bool same = false;
        if (!r.accept_if.empty()) {
          bool b1, b2;
          if (parse_bool_lit(rv, b1) && parse_bool_lit(r.accept_if, b2)) same = (b1 == b2);
          else same = (trim(rv) == r.accept_if) ||
                      (is_number(rv) && is_number(r.accept_if) &&
                       std::stod(rv) == std::stod(r.accept_if));
        }
        if (same) goto next_arg;  // the R default: nothing to do
        throw UsageError("--" + name + " is an R SAIGE flag that " + PROG + " " + sub +
                         " does not support: " + r.why);
      }
      {
        std::string msg = "unknown flag --" + name + " for " + PROG + " " + sub;
        const std::string sg = suggestion(S, name);
        if (!sg.empty()) msg += " (did you mean --" + sg + "?)";
        msg += "; see " + std::string(PROG) + " " + sub + " --help";
        throw UsageError(msg);
      }
    }
    if (f->kind == Kind::Bool) {
      if (!hasEq) {
        bool b;
        if (i + 1 < args.size() && parse_bool_lit(args[i + 1], b)) val = args[++i];
        else val = "TRUE";
      }
      bool b;
      if (!parse_bool_lit(val, b))
        throw UsageError("--" + name + " expects TRUE or FALSE, got '" + val + "'");
    } else if (!hasEq) {
      if (i + 1 >= args.size() || starts_with(args[i + 1], "--"))
        throw UsageError("--" + name + " needs a value (" + f->meta + ")");
      val = args[++i];
    }
    P.items.push_back({f, val});
  next_arg:;
  }
  return P;
}

// ------------------------------------------------------------------ config load
const std::vector<std::string> S1_PATH_KEYS_CFGDIR = {
    "paths.bed", "paths.bim", "paths.fam", "paths.plinkFile", "paths.plinkfile",
    "paths.sparse_grm", "paths.sparse_grm_ids", "paths.out_prefix", "paths.out_prefix_vr"};
const std::vector<std::string> S1_PATH_KEYS_CWD = {"design.csv", "design.whitelist_ids"};
const std::vector<std::string> S2_PATH_KEYS = {
    "plinkFile", "vcfFile", "bgenFile", "bgenSampleFile", "pgenFile", "pvarFile", "psamFile",
    "groupFile", "modelFile", "varianceRatioFile", "outputFile", "checkpointDir"};

void absolutize(YAML::Node& cfg, const std::string& key, const fs::path& base) {
  const std::string v = get_str(cfg, key);
  if (!v.empty()) set_path(cfg, key, YAML::Node(abs_path(v, base)));
}

void load_config(Ctx& c) {
  if (c.configPath.empty()) return;
  if (!fs::exists(c.configPath)) throw UsageError("--config: file not found: " + c.configPath);
  try {
    c.cfg = YAML::LoadFile(c.configPath);
  } catch (const std::exception& e) {
    throw UsageError("--config " + c.configPath + ": " + e.what());
  }
  if (!c.cfg.IsMap()) throw UsageError("--config " + c.configPath + ": the top level must be a map");
  // Relative paths keep the meaning they have when the engine reads this file
  // in place: step 1 resolves paths.* against the config's directory and every
  // other path against the working directory; step 2 resolves all of them
  // against the working directory. The written config holds absolute paths.
  const fs::path cfgDir = fs::path(c.configPath).parent_path();
  const fs::path cwd = fs::current_path();
  if (c.step == 1) {
    for (const auto& k : S1_PATH_KEYS_CFGDIR) absolutize(c.cfg, k, cfgDir);
    for (const auto& k : S1_PATH_KEYS_CWD) absolutize(c.cfg, k, cwd);
    if (c.cfg["models"] && c.cfg["models"].IsSequence())
      for (auto m : c.cfg["models"]) {
        for (const char* k : {"out_prefix", "out_prefix_vr"})
          if (m[k] && m[k].IsScalar()) m[k] = abs_path(m[k].as<std::string>(), cwd);
      }
  } else {
    for (const auto& k : S2_PATH_KEYS) absolutize(c.cfg, k, cwd);
    if (c.cfg["models"] && c.cfg["models"].IsSequence())
      for (auto m : c.cfg["models"]) {
        for (const char* k : {"modelFile", "varianceRatioFile", "outputFile"})
          if (m[k] && m[k].IsScalar()) m[k] = abs_path(m[k].as<std::string>(), cwd);
      }
  }
}

// ------------------------------------------------------------------ resolve
void mkdirs(const std::string& d) {
  if (d.empty()) return;
  std::error_code ec;
  fs::create_directories(d, ec);
  if (ec) throw UsageError("cannot create directory " + d + ": " + ec.message());
}
std::string strip_ext(std::string p, const std::vector<std::string>& exts) {
  for (const auto& e : exts)
    if (ends_with(lower(p), e)) return p.substr(0, p.size() - e.size());
  return p;
}

std::string resolve_step1(Ctx& c) {
  YAML::Node& y = c.cfg;
  const YAML::Node& cy = c.cfg;  // existence checks: const access adds no placeholder keys
  std::string record;
  // ---- outputs
  if (!c.outDir.empty() && !c.outputPrefix.empty())
    throw UsageError("give either --outDir or --outputPrefix, not both");
  if (!c.outDir.empty()) {
    if (!c.outputPrefixVR.empty())
      throw UsageError("--outputPrefix_varRatio goes with --outputPrefix; with --outDir the "
                       "variance ratios are written to <outDir>/vr_<trait>.varianceRatio.txt");
    set_path(y, "paths.out_prefix", YAML::Node(c.outDir + "/models"));
    set_path(y, "paths.out_prefix_vr", YAML::Node(c.outDir + "/vr"));
    record = c.outDir + "/step1.yaml";
  } else if (!c.outputPrefix.empty()) {
    set_path(y, "paths.out_prefix", YAML::Node(c.outputPrefix));
    const std::string vr = c.outputPrefixVR.empty() ? c.outputPrefix : c.outputPrefixVR;
    set_path(y, "paths.out_prefix_vr", YAML::Node(vr));
    record = c.outputPrefix + ".step1.yaml";
  } else {
    throw UsageError("missing --outDir (where the models, variance ratios and step1.yaml go)");
  }

  // ---- phenotypes
  const bool hasModels = (bool)cy["models"];
  if (!c.phenoCols.empty()) {
    erase_path(y, "design.y_col");
    erase_path(y, "design.y_cols");
    if (hasModels) y.remove("models");
    if (!c.outDir.empty() || c.phenoCols.size() > 1) {
      YAML::Node seq(YAML::NodeType::Sequence);
      for (const auto& p : c.phenoCols) seq.push_back(p);
      set_path(y, "design.y_cols", seq);
    } else {
      set_path(y, "design.y_col", YAML::Node(c.phenoCols[0]));
    }
  } else if (has_path(y, "design.y_col") && !c.outDir.empty()) {
    // --outDir layout is models/<trait>/: the y_cols form gives it for one trait too
    YAML::Node seq(YAML::NodeType::Sequence);
    seq.push_back(get_str(y, "design.y_col"));
    erase_path(y, "design.y_col");
    set_path(y, "design.y_cols", seq);
  } else if (!has_path(y, "design.y_cols") && !hasModels && !has_path(y, "design.y_col")) {
    throw UsageError("missing --phenoCol (the phenotype column, or a comma list of them)");
  }

  // ---- genotypes
  if (!c.plinkFile.empty()) {
    if (!c.plinkTrioFromFlags)
      for (const char* k : {"paths.bed", "paths.bim", "paths.fam"}) erase_path(y, k);
    erase_path(y, "paths.plinkfile");
    set_path(y, "paths.plinkFile", YAML::Node(c.plinkFile));
  }
  const bool havePlink = has_path(y, "paths.plinkFile") || has_path(y, "paths.plinkfile") ||
      (has_path(y, "paths.bed") && has_path(y, "paths.bim") && has_path(y, "paths.fam"));
  if (!havePlink)
    throw UsageError("missing --plinkFile (or all three of --bedFile, --bimFile, --famFile)");
  if (!has_path(y, "design.csv"))
    throw UsageError("missing --phenoFile (the phenotype/covariate table)");

  // ---- sparse GRM build
  bool mk = false;
  { YAML::Node n = get_path(y, "fit.make_sparse_grm_only"); if (n.IsDefined() && n.IsScalar()) mk = n.as<bool>(); }
  if (mk && (!has_path(y, "paths.sparse_grm") || !has_path(y, "paths.sparse_grm_ids")))
    throw UsageError("--makeSparseGRMOnly needs --sparseGRMFile and --sparseGRMSampleIDFile "
                     "(the files to write)");
  bool useSp = false;
  { YAML::Node n = get_path(y, "fit.use_sparse_grm_to_fit"); if (n.IsDefined() && n.IsScalar()) useSp = n.as<bool>(); }
  if (useSp && !mk && (!has_path(y, "paths.sparse_grm") || !has_path(y, "paths.sparse_grm_ids")))
    throw UsageError("--useSparseGRMtoFitNULL=TRUE needs --sparseGRMFile and --sparseGRMSampleIDFile");

  // everything checked: create the output directories
  mkdirs(fs::path(record).parent_path().string());
  mkdirs(fs::path(get_str(y, "paths.out_prefix")).parent_path().string());
  mkdirs(fs::path(get_str(y, "paths.out_prefix_vr")).parent_path().string());
  return record;
}

struct TraitModel { std::string name, model, vr; };

std::vector<TraitModel> discover_step1(const std::string& dir) {
  if (!fs::is_directory(dir)) throw UsageError("--step1Dir: not a directory: " + dir);
  std::vector<TraitModel> out;
  const std::string rec = dir + "/step1.yaml";
  auto exists_model = [](const std::string& d) { return fs::exists(d + "/nullmodel.json"); };
  if (fs::exists(rec)) {
    YAML::Node y;
    try { y = YAML::LoadFile(rec); }
    catch (const std::exception& e) { throw UsageError("--step1Dir: cannot read " + rec + ": " + e.what()); }
    const std::string op = get_str(y, "paths.out_prefix");
    std::string opv = get_str(y, "paths.out_prefix_vr");
    if (opv.empty()) opv = op;
    std::vector<TraitModel> cand;
    YAML::Node ycols = get_path(y, "design.y_cols");
    if (ycols.IsDefined() && ycols.IsSequence()) {
      for (const auto& t : ycols) {
        const std::string n = t.as<std::string>();
        cand.push_back({n, op + "/" + n, opv + "_" + n + ".varianceRatio.txt"});
      }
    } else if (y["models"] && y["models"].IsSequence()) {
      for (const auto& m : y["models"]) {
        const std::string n = m["y_col"].as<std::string>();
        cand.push_back({n, m["out_prefix"] ? m["out_prefix"].as<std::string>() : op + "/" + n,
                        (m["out_prefix_vr"] ? m["out_prefix_vr"].as<std::string>() : opv + "_" + n) +
                            ".varianceRatio.txt"});
      }
    } else if (has_path(y, "design.y_col")) {
      const std::string n = get_str(y, "design.y_col");
      cand.push_back({n, op, opv + ".varianceRatio.txt"});
    }
    for (auto& t : cand) {
      // the directory may have been moved since: fall back to the standard layout
      if (!exists_model(t.model) && exists_model(dir + "/models/" + t.name)) {
        t.model = dir + "/models/" + t.name;
        t.vr = dir + "/vr_" + t.name + ".varianceRatio.txt";
      }
      out.push_back(t);
    }
  } else if (fs::is_directory(dir + "/models")) {
    std::vector<std::string> names;
    for (const auto& e : fs::directory_iterator(dir + "/models"))
      if (e.is_directory() && exists_model(e.path().string())) names.push_back(e.path().filename().string());
    std::sort(names.begin(), names.end());
    for (const auto& n : names)
      out.push_back({n, dir + "/models/" + n, dir + "/vr_" + n + ".varianceRatio.txt"});
  }
  if (out.empty())
    throw UsageError("--step1Dir " + dir + ": no step-1 models found (expected the step1.yaml and "
                     "models/<trait>/nullmodel.json that `" + std::string(PROG) +
                     " step1 --outDir " + dir + "` writes)");
  return out;
}

std::string resolve_step2(Ctx& c) {
  YAML::Node& y = c.cfg;
  const YAML::Node& cy = c.cfg;  // existence checks: const access adds no placeholder keys
  // ---- genotypes
  const bool fPlink = !c.s2plink.empty() || !c.s2bed.empty() || !c.s2bim.empty() || !c.s2fam.empty();
  const bool fPgen = !c.pgenPrefix.empty() || !c.pgenFile.empty() || !c.pvarFile.empty() || !c.psamFile.empty();
  const bool fBgen = !c.bgenFile.empty() || !c.sampleFile.empty() || !c.bgenIndex.empty();
  const bool fVcf = !c.vcfFile.empty() || !c.vcfIndex.empty();
  const int nFam = (int)fPlink + (int)fPgen + (int)fBgen + (int)fVcf;
  if (nFam > 1)
    throw UsageError("give the genotypes in one format: --plinkFile / --bedFile..., --pgenPrefix / "
                     "--pgenFile..., --bgenFile, or --vcfFile");
  auto clear_geno = [&]() {
    for (const char* k : {"plinkFile", "vcfFile", "bgenFile", "bgenSampleFile", "pgenFile", "pvarFile", "psamFile"})
      if (cy[k]) y.remove(k);
  };
  if (fPlink) {
    std::string prefix;
    if (!c.s2plink.empty()) {
      prefix = strip_ext(c.s2plink, {".bed", ".bim", ".fam"});
      if (!c.s2bed.empty() || !c.s2bim.empty() || !c.s2fam.empty())
        throw UsageError("give --plinkFile or --bedFile/--bimFile/--famFile, not both");
    } else {
      if (c.s2bed.empty() || c.s2bim.empty() || c.s2fam.empty())
        throw UsageError("--bedFile, --bimFile and --famFile go together (or use --plinkFile PREFIX)");
      const std::string a = strip_ext(c.s2bed, {".bed"}), b = strip_ext(c.s2bim, {".bim"}),
                        f = strip_ext(c.s2fam, {".fam"});
      if (!ends_with(lower(c.s2bed), ".bed") || !ends_with(lower(c.s2bim), ".bim") ||
          !ends_with(lower(c.s2fam), ".fam") || a != b || a != f)
        throw UsageError("this step 2 reads PREFIX.bed/.bim/.fam: --bedFile, --bimFile and --famFile "
                         "must share one prefix (or use --plinkFile PREFIX)");
      prefix = a;
    }
    clear_geno();
    y["genoType"] = "plink";
    y["plinkFile"] = prefix;
  } else if (fPgen) {
    std::string pre = c.pgenPrefix;
    if (!pre.empty() && !c.pgenFile.empty())
      throw UsageError("give --pgenPrefix or --pgenFile, not both");
    std::string pg = !c.pgenFile.empty() ? c.pgenFile : (pre.empty() ? "" : pre + ".pgen");
    if (pg.empty()) throw UsageError("--pvarFile/--psamFile need --pgenFile (or use --pgenPrefix)");
    if (pre.empty()) pre = strip_ext(pg, {".pgen"});
    clear_geno();
    y["genoType"] = "pgen";
    y["pgenFile"] = pg;
    y["pvarFile"] = c.pvarFile.empty() ? pre + ".pvar" : c.pvarFile;
    y["psamFile"] = c.psamFile.empty() ? pre + ".psam" : c.psamFile;
  } else if (fBgen) {
    if (c.bgenFile.empty()) throw UsageError("--sampleFile/--bgenFileIndex need --bgenFile");
    if (c.sampleFile.empty()) throw UsageError("--bgenFile needs --sampleFile (the BGEN .sample file)");
    if (!c.bgenIndex.empty() && c.bgenIndex != c.bgenFile + ".bgi")
      throw UsageError("--bgenFileIndex: this step 2 reads the index from <bgenFile>.bgi (" +
                       c.bgenFile + ".bgi); put or link the index there and drop the flag");
    clear_geno();
    y["genoType"] = "bgen";
    y["bgenFile"] = c.bgenFile;
    y["bgenSampleFile"] = c.sampleFile;
  } else if (fVcf) {
    if (c.vcfFile.empty()) throw UsageError("--vcfFileIndex needs --vcfFile");
    if (!c.vcfIndex.empty())
      c.notes.push_back("--vcfFileIndex is not used: the VCF is read from start to end");
    clear_geno();
    y["genoType"] = "vcf";
    y["vcfFile"] = c.vcfFile;
  } else {
    if (!cy["plinkFile"] && !cy["vcfFile"] && !cy["bgenFile"] && !cy["pgenFile"])
      throw UsageError("missing genotypes: --plinkFile PREFIX, --pgenPrefix PREFIX, "
                       "--bgenFile FILE --sampleFile FILE, or --vcfFile FILE");
  }

  // ---- models
  if (!c.step1Dir.empty() && (!c.gmmat.empty() || !c.vrFile.empty() || !c.outFile.empty()))
    throw UsageError("give --step1Dir or --GMMATmodelFile/--varianceRatioFile/--SAIGEOutputFile, not both");
  if (!c.phenoCols.empty() && c.step1Dir.empty())
    throw UsageError("--phenoCol in step 2 picks traits out of --step1Dir; it needs --step1Dir");
  std::vector<std::string> outputs;
  auto clear_models = [&]() {
    for (const char* k : {"models", "modelFile", "varianceRatioFile", "outputFile", "traitName"})
      if (cy[k]) y.remove(k);
  };
  if (!c.step1Dir.empty()) {
    if (c.outDir.empty()) throw UsageError("--step1Dir needs --outDir (where <trait>.txt go)");
    auto all = discover_step1(c.step1Dir);
    std::vector<TraitModel> sel;
    if (c.phenoCols.empty()) sel = all;
    else
      for (const auto& p : c.phenoCols) {
        auto it = std::find_if(all.begin(), all.end(), [&](const TraitModel& t) { return t.name == p; });
        if (it == all.end()) {
          std::string avail;
          for (const auto& t : all) avail += (avail.empty() ? "" : ", ") + t.name;
          throw UsageError("--phenoCol " + p + ": no such trait in " + c.step1Dir + " (it has: " + avail + ")");
        }
        sel.push_back(*it);
      }
    for (const auto& t : sel) {
      if (!fs::exists(t.model + "/nullmodel.json"))
        throw UsageError("--step1Dir: trait " + t.name + " has no model (" + t.model +
                         "/nullmodel.json missing); did its step 1 fail?");
      if (!fs::exists(t.vr))
        throw UsageError("--step1Dir: trait " + t.name + " has no variance-ratio file (" + t.vr + ")");
    }
    clear_models();
    YAML::Node seq(YAML::NodeType::Sequence);
    for (const auto& t : sel) {
      YAML::Node m(YAML::NodeType::Map);
      m["traitName"] = t.name;
      m["modelFile"] = t.model;
      m["varianceRatioFile"] = t.vr;
      m["outputFile"] = c.outDir + "/" + t.name + ".txt";
      outputs.push_back(c.outDir + "/" + t.name + ".txt");
      seq.push_back(m);
    }
    y["models"] = seq;
  } else if (!c.gmmat.empty() || !c.vrFile.empty()) {
    if (c.gmmat.empty()) throw UsageError("--varianceRatioFile needs --GMMATmodelFile");
    if (c.vrFile.empty()) throw UsageError("--GMMATmodelFile needs --varianceRatioFile");
    std::string out = c.outFile;
    if (out.empty()) {
      if (c.outDir.empty()) throw UsageError("give --SAIGEOutputFile or --outDir");
      out = c.outDir + "/" + fs::path(c.gmmat).lexically_normal().filename().string() + ".txt";
      if (ends_with(c.gmmat, "/"))
        out = c.outDir + "/" + fs::path(c.gmmat.substr(0, c.gmmat.size() - 1)).filename().string() + ".txt";
    }
    clear_models();
    y["modelFile"] = c.gmmat;
    y["varianceRatioFile"] = c.vrFile;
    y["outputFile"] = out;
    outputs.push_back(out);
  } else {
    if (!c.outFile.empty()) {
      if (!cy["modelFile"])
        throw UsageError("--SAIGEOutputFile goes with --GMMATmodelFile (or a config with modelFile:)");
      y["outputFile"] = c.outFile;
    }
    if (cy["models"] && cy["models"].IsSequence()) {
      for (const auto& m : cy["models"]) if (m["outputFile"]) outputs.push_back(m["outputFile"].as<std::string>());
    } else if (cy["modelFile"]) {
      if (cy["outputFile"]) outputs.push_back(cy["outputFile"].as<std::string>());
    } else {
      throw UsageError("missing models: --step1Dir DIR (every trait of a step-1 run), or "
                       "--GMMATmodelFile DIR --varianceRatioFile FILE");
    }
  }

  // ---- where the record goes
  std::string record;
  if (!c.outDir.empty()) record = c.outDir + "/step2.yaml";
  else if (!c.outFile.empty()) record = c.outFile + ".step2.yaml";
  else throw UsageError("missing --outDir (where the results and step2.yaml go)");
  if (!c.outDir.empty()) mkdirs(c.outDir);
  for (const auto& o : outputs) mkdirs(fs::path(o).parent_path().string());

  if (!c.overwriteOutput) {
    const bool sgs = cy["outputFormat"] && cy["outputFormat"].as<std::string>() == "sgs";
    for (const auto& o : outputs) {
      const std::string p = sgs ? o + ".sgs" : o;
      if (fs::exists(p))
        throw UsageError("--is_overwrite_output=FALSE and " + p + " exists");
    }
  }
  return record;
}

// ------------------------------------------------------------------ run
std::string self_dir() {
  std::error_code ec;
  fs::path p = fs::read_symlink("/proc/self/exe", ec);
  if (ec) return ".";
  return p.parent_path().string();
}

std::string engine_path(const char* name) {
  const std::string p = self_dir() + "/" + name;
  if (::access(p.c_str(), X_OK) != 0)
    throw std::runtime_error(std::string("engine not found: ") + p + " (" + PROG +
                             " runs the engines installed next to it; rebuild with `make`)");
  return p;
}

std::string command_line(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& a : v) s += (s.empty() ? "" : " ") + shell_quote(a);
  return s;
}

int run_step(int step, const std::vector<std::string>& args) {
  Ctx c;
  c.step = step;
  const auto S = step == 1 ? step1_sections() : step2_sections();
  const auto R = step == 1 ? step1_rejected() : step2_rejected();
  Parsed P = parse_args(step, args, S, R);
  if (P.items.empty()) { print_help(step); return 2; }

  // 1) --config first, wherever it appears; 2) every other flag, in order
  for (const auto& it : P.items)
    if (it.first->name == "config") c.configPath = abs_path(trim(it.second));
  load_config(c);
  if (step == 1)  // section order of the written file; empty ones are dropped below
    for (const char* k : {"paths", "design", "fit"})
      if (!c.cfg[k]) c.cfg[k] = YAML::Node(YAML::NodeType::Map);
  for (const auto& it : P.items) {
    if (it.first->name == "config") continue;
    c.given.insert(it.first->name);
    it.first->apply(c, it.second);
  }

  const std::string record = step == 1 ? resolve_step1(c) : resolve_step2(c);

  std::vector<std::string> cmd;
  if (step == 1) cmd = {engine_path(ENGINE_STEP1), "-c", record};
  else cmd = {engine_path(ENGINE_STEP2), record};
  std::string envPrefix;
  if (step == 1 && c.gpuDevice >= 0) envPrefix = "SAIGE_GPU_DEVICE=" + std::to_string(c.gpuDevice) + " ";

  // The record: the resolved config plus how it was made.
  if (step == 1)
    for (const char* k : {"paths", "design", "fit"})
      if (c.cfg[k] && c.cfg[k].IsMap() && c.cfg[k].size() == 0) c.cfg.remove(k);
  {
    std::vector<std::string> invoked = {PROG, step == 1 ? "step1" : "step2"};
    invoked.insert(invoked.end(), args.begin(), args.end());
    YAML::Emitter e;
    e << c.cfg;
    std::ofstream f(record);
    if (!f) throw std::runtime_error("cannot write " + record);
    f << "# Resolved step-" << step << " config, written by " << PROG << " " << SAIGE_CLI_VERSION << ".\n"
      << "# command: " << command_line(invoked) << "\n"
      << "# working directory: " << fs::current_path().string() << "\n"
      << "# same run by hand: " << envPrefix << command_line(cmd) << "\n"
      << e.c_str() << "\n";
    if (!f) throw std::runtime_error("cannot write " + record);
  }
  for (const auto& n : c.notes) std::cout << PROG << ": note: " << n << "\n";
  std::cout << PROG << ": config written to " << record << "\n";
  if (c.dryRun) {
    std::cout << envPrefix << command_line(cmd) << "\n";
    return 0;
  }
  std::cout << PROG << ": running " << envPrefix << command_line(cmd) << "\n" << std::flush;
  if (step == 1 && c.gpuDevice >= 0) setenv("SAIGE_GPU_DEVICE", std::to_string(c.gpuDevice).c_str(), 1);
  std::vector<char*> av;
  for (auto& s : cmd) av.push_back(const_cast<char*>(s.c_str()));
  av.push_back(nullptr);
  ::execv(av[0], av.data());
  throw std::runtime_error(std::string("cannot run ") + cmd[0] + ": " + std::strerror(errno));
}

int run_sgs2txt(const std::vector<std::string>& args) {
  if (args.empty() || args[0] == "--help" || args[0] == "-h") {
    std::cout <<
      "Usage: " << PROG << " sgs2txt [-m MARKERS.sgs] [-o OUT.txt] [-j N] TRAIT.sgs [TRAIT.sgs ...]\n\n"
      "Converts the .sgs files of `" << PROG << " step2 --outputFormat sgs` back to the text\n"
      "result files (byte-identical to --outputFormat text when written with sgsPrecision fp64).\n"
      "  (no -o)          write each trait's text to the path recorded in its .sgs\n"
      "                   (<outDir>/<trait>.txt of the step-2 run)\n"
      "  -o OUT           write OUT instead (one input file)\n"
      "  -m MARKERS.sgs   the shared markers file (default: the one recorded; use it after\n"
      "                   moving the files: <outDir>/<first trait>.txt.markers.sgs)\n"
      "  -j N             convert N traits at a time\n";
    return args.empty() ? 2 : 0;
  }
  std::vector<std::string> cmd = {engine_path(ENGINE_SGS2TXT)};
  cmd.insert(cmd.end(), args.begin(), args.end());
  std::vector<char*> av;
  for (auto& s : cmd) av.push_back(const_cast<char*>(s.c_str()));
  av.push_back(nullptr);
  ::execv(av[0], av.data());
  throw std::runtime_error(std::string("cannot run ") + cmd[0] + ": " + std::strerror(errno));
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  try {
    if (args.empty() || args[0] == "--help" || args[0] == "-h" || args[0] == "help") {
      if (args.size() >= 2 && args[0] == "help") {
        if (args[1] == "step1") { print_help(1); return 0; }
        if (args[1] == "step2") { print_help(2); return 0; }
        if (args[1] == "sgs2txt") return run_sgs2txt({"--help"});
      }
      print_main_help();
      return args.empty() ? 2 : 0;
    }
    if (args[0] == "--version" || args[0] == "-V" || args[0] == "version") {
      std::cout << PROG << " " << SAIGE_CLI_VERSION;
      if (std::strlen(SAIGE_CLI_BUILD)) std::cout << " (" << SAIGE_CLI_BUILD << ")";
      std::cout << "\n";
      return 0;
    }
    const std::string sub = args[0];
    std::vector<std::string> rest(args.begin() + 1, args.end());
    if (sub == "step1") return run_step(1, rest);
    if (sub == "step2") return run_step(2, rest);
    if (sub == "sgs2txt") return run_sgs2txt(rest);
    throw UsageError("unknown command '" + sub + "' (step1, step2, sgs2txt; see " + PROG + " --help)");
  } catch (const UsageError& e) {
    std::cerr << PROG << ": error: " << e.what() << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << PROG << ": error: " << e.what() << "\n";
    return 1;
  }
}
