// bed_parallel_test.cpp — PR-3 end-to-end + scaling harness.
//
// What it does:
//   1. Load IIDs from a saige-null cov_input.csv (or TSV pheno), build
//      FAM-ascending ptrsubSampleInGeno (same logic as bed_pipeline_test).
//   2. For each requested thread count, run parallel_decode_bed() over the
//      first `n_markers` of the BED and record wall time.
//   3. The nthreads=1 output is taken as ground truth; each N>1 run must
//      produce byte-identical PackedFlat + stats + orig_plink_idx.
//   4. First 20 markers' stats are additionally diffed against the saige-null
//      reference log (if supplied) for absolute correctness.
//
// Usage:
//   ./bed_parallel_test <bed> <fam> <iids.csv|tsv> <iid_col> <ref_log>
//                        [n_markers=20] [threads=1,4,16] [minMAF=0.01] [maxMiss=0.15]
//
// Example (matching Makefile target):
//   ./bed_parallel_test $UKB_BED $UKB_FAM \
//       .../ukb_ldl_dense_output_cov_input.csv iid \
//       .../benchmark_results/ukb_ldl/cpp_covT/stdout.log \
//       20 1,4,16
#include "bed_reader.hpp"
#include "marker_decoder.hpp"
#include "packed_store.hpp"
#include "parallel_decode.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// --- tiny helpers (duplicated with bed_pipeline_test; keep the two sandboxes
//     self-contained so either can be deleted without breaking the other) ---

std::vector<std::string> split_delim(const std::string& s, char d) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == d) { out.push_back(cur); cur.clear(); } else cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}
char detect_delim(const std::string& hdr) {
  return (hdr.find('\t') != std::string::npos) ? '\t' : ',';
}
std::vector<std::string> read_fam_iids(const std::string& p) {
  std::ifstream f(p);
  if (!f) throw std::runtime_error("open " + p);
  std::vector<std::string> out;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream is(line);
    std::string fid, iid;
    if (is >> fid >> iid) out.push_back(iid);
  }
  return out;
}
std::vector<std::string> read_iids(const std::string& p, const std::string& col) {
  std::ifstream f(p);
  if (!f) throw std::runtime_error("open " + p);
  std::string hdr; if (!std::getline(f, hdr)) throw std::runtime_error("empty " + p);
  const char d = detect_delim(hdr);
  auto cols = split_delim(hdr, d);
  int idx = -1;
  for (size_t i = 0; i < cols.size(); ++i) if (cols[i] == col) { idx = int(i); break; }
  if (idx < 0) throw std::runtime_error("input missing col '" + col + "'");
  std::vector<std::string> out;
  std::string line;
  while (std::getline(f, line)) {
    auto row = split_delim(line, d);
    if (int(row.size()) > idx) out.push_back(row[idx]);
  }
  return out;
}

struct RefLine {
  int   marker;
  float freq, maf, missRate;
  bool  passQC;
};
std::vector<RefLine> parse_ref_log(const std::string& p, int n) {
  std::ifstream f(p);
  if (!f) throw std::runtime_error("open " + p);
  std::regex re(R"(^Marker\s+(\d+):\s+freq=([-0-9.eE+]+),\s+maf=([-0-9.eE+]+),\s+missRate=([-0-9.eE+]+),\s+passQC=([01])\s*$)");
  std::string line;
  std::vector<RefLine> out;
  while (std::getline(f, line)) {
    std::smatch m;
    if (!std::regex_match(line, m, re)) continue;
    RefLine r;
    r.marker = std::stoi(m[1]);
    r.freq = std::stof(m[2]);
    r.maf = std::stof(m[3]);
    r.missRate = std::stof(m[4]);
    r.passQC = (m[5] == "1");
    if (r.marker < n) out.push_back(r);
  }
  std::sort(out.begin(), out.end(),
            [](const RefLine& a, const RefLine& b) { return a.marker < b.marker; });
  out.erase(std::unique(out.begin(), out.end(),
            [](const RefLine& a, const RefLine& b) { return a.marker == b.marker; }),
            out.end());
  return out;
}
bool approx_eq(float a, float b, float rel = 1e-4f, float abs_ = 1e-6f) {
  const float d = std::fabs(a - b);
  return d <= abs_ || d <= rel * std::max(std::fabs(a), std::fabs(b));
}

std::vector<int> parse_threads_list(const std::string& s) {
  std::vector<int> out;
  std::string cur;
  for (char c : s) {
    if (c == ',') { if (!cur.empty()) out.push_back(std::atoi(cur.c_str())); cur.clear(); }
    else cur.push_back(c);
  }
  if (!cur.empty()) out.push_back(std::atoi(cur.c_str()));
  return out;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr,
      "usage: %s <bed> <fam> <iids.csv|tsv> <iid_col> <ref_log> "
      "[n_markers=20] [threads=1,4,16] [minMAF=0.01] [maxMiss=0.15]\n",
      argv[0]);
    return 2;
  }
  const std::string bed_path = argv[1];
  const std::string fam_path = argv[2];
  const std::string iid_src  = argv[3];
  const std::string iid_col  = argv[4];
  const std::string log_path = argv[5];
  const std::size_t n_markers_arg = (argc > 6) ? std::strtoull(argv[6], nullptr, 10) : 20ULL;
  const std::string threads_str = (argc > 7) ? argv[7] : "1,4,16";
  const float min_maf   = (argc > 8) ? std::atof(argv[8]) : 0.01f;
  const float max_miss  = (argc > 9) ? std::atof(argv[9]) : 0.15f;

  try {
    // --- Inputs ------------------------------------------------------------
    std::fprintf(stderr, "[1/5] loading FAM + IID source...\n");
    auto fam_iids = read_fam_iids(fam_path);
    auto iids     = read_iids(iid_src, iid_col);
    const std::size_t N = fam_iids.size();
    std::unordered_map<std::string, int> fam_pos;
    fam_pos.reserve(N * 2);
    for (std::size_t i = 0; i < N; ++i) fam_pos.emplace(fam_iids[i], int(i) + 1);

    std::vector<bool> indicator(N, false);
    for (const auto& iid : iids) {
      auto it = fam_pos.find(iid);
      if (it != fam_pos.end()) indicator[it->second - 1] = true;
    }
    std::vector<int> ptrsub;
    ptrsub.reserve(iids.size());
    for (std::size_t i = 0; i < N; ++i)
      if (indicator[i]) ptrsub.push_back(int(i) + 1);
    const std::size_t Nnomissing = ptrsub.size();
    std::fprintf(stderr, "  N=%zu  Nnomissing=%zu\n", N, Nnomissing);

    // --- Ref log -----------------------------------------------------------
    std::fprintf(stderr, "[2/5] parsing reference log (first 20 markers)...\n");
    const auto ref = parse_ref_log(log_path, 20);
    std::fprintf(stderr, "  %zu ref lines\n", ref.size());

    // --- Reader geometry ---------------------------------------------------
    const std::vector<int> threads_list = parse_threads_list(threads_str);
    const int max_threads = *std::max_element(threads_list.begin(), threads_list.end());
    if (max_threads < 1) throw std::runtime_error("need at least one thread in list");

    // Determine M from the file.
    saige::BedReaderPool tmp(bed_path, N, 1);
    const std::size_t M_total = tmp.n_markers_from_file_size();
    const std::size_t n_markers = (n_markers_arg == 0) ? M_total
                                  : std::min(n_markers_arg, M_total);
    std::fprintf(stderr, "[3/5] BED geometry: M_total=%zu  running on %zu markers\n",
                 M_total, n_markers);

    // --- Reference run at nthreads=1 (ground truth) ------------------------
    std::fprintf(stderr, "[4/5] building serial reference (nthreads=1)...\n");
    saige::BedReaderPool reader1(bed_path, N, 1);
    auto t0 = std::chrono::steady_clock::now();
    auto ref_res = saige::parallel_decode_bed(
        reader1, ptrsub.data(), Nnomissing, n_markers, min_maf, max_miss, 1);
    auto t1 = std::chrono::steady_clock::now();
    const double sec1 = std::chrono::duration<double>(t1 - t0).count();
    std::fprintf(stderr, "  serial: %.3f s  (pass-QC=%zu, store.bytes=%zu)\n",
                 sec1, ref_res.store.n_stored(), ref_res.store.bytes());

    // Absolute check: first 20 markers vs saige-null log.
    int abs_mismatch = 0;
    for (int i = 0; i < 20 && static_cast<std::size_t>(i) < n_markers; ++i) {
      const saige::MarkerStats& s = ref_res.stats[i];
      const RefLine* r = nullptr;
      for (const auto& x : ref) if (x.marker == i) { r = &x; break; }
      if (!r) continue;
      const bool bad =
          !approx_eq(r->freq, s.altFreq) ||
          !approx_eq(r->maf, std::min(s.altFreq, 1.0f - s.altFreq)) ||
          !approx_eq(r->missRate, s.missingRate) ||
          (r->passQC != s.passQC);
      if (bad) ++abs_mismatch;
    }
    std::fprintf(stderr, "  absolute mismatches on first 20 markers: %d\n", abs_mismatch);
    if (abs_mismatch) {
      std::fprintf(stderr, "  (serial decode disagrees with saige-null log — bug upstream)\n");
      return 1;
    }

    // --- Multi-threaded runs vs reference ---------------------------------
    std::fprintf(stderr, "[5/5] scaling runs:\n");
    std::printf("\n%-6s %-12s %-10s %-16s %-16s %s\n",
                "nth", "wall(s)", "speedup", "store==ref?", "stats==ref?", "orig_idx==ref?");
    std::printf("%-6s %-12s %-10s %-16s %-16s %s\n",
                "---", "---", "---", "---", "---", "---");

    int fail_any = 0;
    for (int nt : threads_list) {
      if (nt < 1) continue;
      saige::BedReaderPool reader(bed_path, N, nt);
      auto ta = std::chrono::steady_clock::now();
      auto res = saige::parallel_decode_bed(
          reader, ptrsub.data(), Nnomissing, n_markers, min_maf, max_miss, nt);
      auto tb = std::chrono::steady_clock::now();
      const double sec = std::chrono::duration<double>(tb - ta).count();

      // Byte-level equality to serial reference.
      const bool eq_bytes = (res.store.n_stored() == ref_res.store.n_stored())
                         && (res.store.bytes()    == ref_res.store.bytes())
                         && (std::memcmp(res.store.raw(), ref_res.store.raw(),
                                         res.store.bytes()) == 0);
      const bool eq_orig = (res.orig_plink_idx == ref_res.orig_plink_idx);

      // Stats equality (field-by-field). Use strict equality on passQC + ints,
      // approx_eq on floats.
      bool eq_stats = (res.stats.size() == ref_res.stats.size());
      if (eq_stats) {
        for (std::size_t i = 0; i < res.stats.size(); ++i) {
          const auto& a = res.stats[i];
          const auto& b = ref_res.stats[i];
          if (a.passQC != b.passQC || a.alleleCount != b.alleleCount ||
              a.numMissing != b.numMissing || a.mac != b.mac ||
              !approx_eq(a.altFreq, b.altFreq, 1e-6f) ||
              !approx_eq(a.missingRate, b.missingRate, 1e-6f)) {
            eq_stats = false; break;
          }
        }
      }

      const double speedup = (sec > 0.0) ? (sec1 / sec) : 0.0;
      std::printf("%-6d %-12.3f %-10.2fx %-16s %-16s %s\n",
                  nt, sec, speedup,
                  eq_bytes ? "bit-identical" : "MISMATCH",
                  eq_stats ? "bit-identical" : "MISMATCH",
                  eq_orig  ? "bit-identical" : "MISMATCH");
      if (!(eq_bytes && eq_stats && eq_orig)) fail_any = 1;
    }

    std::printf("\n");
    if (fail_any) {
      std::fprintf(stderr, "bed_parallel_test: at least one nthreads run disagreed with the serial reference\n");
      return 1;
    }
    std::fprintf(stderr, "bed_parallel_test: all runs bit-identical to serial reference\n");
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "bed_parallel_test: %s\n", e.what());
    return 1;
  }
}
