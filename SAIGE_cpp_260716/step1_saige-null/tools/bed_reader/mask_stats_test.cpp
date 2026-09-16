// mask_stats_test.cpp — gates A1 and A2 of SCHEME_C_DESIGN.md §5, on real data.
//
// A1  For each trait t: decode the UNION of all traits' sample sets once with
//     exclusion sets, rebuild t's per-marker stats from the tallies, and compare
//     against decoding t's OWN sample set through the unmodified code path.
//     freq / invstd / fill / passQC / passVR / mac / numMissing / missingRate /
//     M_t must be identical BIT FOR BIT (fp32 bit patterns, not a tolerance).
//
// A2  For each trait: fill_U[j] + delta == fill_t[j] on every correction entry,
//     and the correction list is exactly the set of (union row, marker) cells
//     where the marker is missing, the row belongs to S_t, and the two fills
//     differ. Ground truth is brute-forced from a raw BED re-decode, not from
//     the decoder's own missing-cell output.
//
// Also checked along the way:
//   · the §2 union-pack rule (keep marker j iff some trait's passQC_t is true):
//     the decoder's keep flag == OR_t passQC_t from compute_trait_stats, and a
//     pack_if_keep decode stores exactly those markers, with byte-identical
//     rows to the ordinary decode.
//   · the decoder's per-marker missing-cell lists vs the brute-force decode.
//
// Sample sets come from a phenotype TSV the multi-trait tolerance gate already
// generated (tests/mt_tolerance/make_cases.py): S_t = FAM rows whose IID has a
// non-NA value in trait column t.
//
// usage:
//   ./mask_stats_test <bed> <fam> <pheno.tsv> [options]
//     --traits a,b,c      trait columns (default: all but IID and x*)
//     --markers N         first N markers only (default: all)
//     --threads T         decode threads (default 4)
//     --min-maf F         default 0.01
//     --max-miss F        default 0.15
//     --vr MIN[,MAX]      enable the variance-ratio rule (MAX -1 == open bin)
//     --vrdraw PATH       vrdraw binary (default /opt/saige/logs/mt_gate/data/vrdraw)
//     --synth P,PCT,SEED  ignore the pheno traits; build P traits, each dropping
//                         a random 1..PCT% of the FAM samples (the biobank-like
//                         regime scheme C actually targets)
// Exit 0 when every assertion holds.

#include "bed_reader.hpp"
#include "marker_decoder.hpp"
#include "mask_stats.hpp"
#include "packed_store.hpp"
#include "parallel_decode.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

std::vector<std::string> split_delim(const std::string& s, char d) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == d) { out.push_back(cur); cur.clear(); }
    else if (c != '\r') cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

char detect_delim(const std::string& hdr) {
  if (hdr.find('\t') != std::string::npos) return '\t';
  if (hdr.find(',')  != std::string::npos) return ',';
  return ' ';
}

std::vector<std::string> read_fam_iids(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("open " + path);
  std::vector<std::string> out;
  std::string line;
  while (std::getline(f, line)) {
    std::istringstream is(line);
    std::string fid, iid;
    if (!(is >> fid >> iid)) continue;
    out.push_back(iid);
  }
  return out;
}

bool is_na(const std::string& v) {
  return v.empty() || v == "NA" || v == "NaN" || v == "nan" || v == "-9" || v == ".";
}

// Bit-level float comparison — the point of gate A1.
inline uint32_t fbits(float x) {
  uint32_t u;
  std::memcpy(&u, &x, sizeof(u));
  return u;
}
inline bool fsame(float a, float b) { return fbits(a) == fbits(b); }

struct Fail {
  int  n = 0;
  int  shown = 0;
  void report(const char* what, const std::string& trait, std::size_t j,
              const std::string& got, const std::string& want) {
    ++n;
    if (shown < 12) {
      ++shown;
      std::printf("  FAIL [%s] trait=%s marker=%zu  scheme-C=%s  solo=%s\n",
                  what, trait.c_str(), j, got.c_str(), want.c_str());
    }
  }
};

std::string f2s(float x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.9g (0x%08x)", x, fbits(x));
  return buf;
}

std::vector<unsigned char> load_vr_draw(const std::string& bin, std::size_t M) {
  std::string cmd = bin + " " + std::to_string(M);
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) throw std::runtime_error("popen " + cmd);
  std::vector<unsigned char> drawn(M, 0);
  long idx;
  int n = 0;
  while (std::fscanf(p, "%ld", &idx) == 1) {
    if (idx >= 0 && idx < static_cast<long>(M)) { drawn[idx] = 1; ++n; }
  }
  pclose(p);
  if (n == 0) throw std::runtime_error("vrdraw produced nothing: " + cmd);
  std::fprintf(stderr, "  VR draw: %d of %zu markers\n", n, M);
  return drawn;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::fprintf(stderr,
      "usage: %s <bed> <fam> <pheno.tsv> [--traits a,b] [--markers N] "
      "[--threads T] [--min-maf F] [--max-miss F] [--vr MIN[,MAX]] "
      "[--vrdraw PATH] [--synth P,PCT,SEED]\n", argv[0]);
    return 2;
  }
  const std::string bed_path   = argv[1];
  const std::string fam_path   = argv[2];
  const std::string pheno_path = argv[3];

  std::string trait_arg, synth_arg, vr_arg;
  std::string vrdraw_bin = "/opt/saige/logs/mt_gate/data/vrdraw";
  long  n_markers_arg = -1;
  int   nthreads      = 4;
  float min_maf       = 0.01f;
  float max_miss      = 0.15f;
  for (int i = 4; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
      return argv[++i];
    };
    if      (a == "--traits")   trait_arg     = next();
    else if (a == "--markers")  n_markers_arg = std::atol(next().c_str());
    else if (a == "--threads")  nthreads      = std::atoi(next().c_str());
    else if (a == "--min-maf")  min_maf       = std::atof(next().c_str());
    else if (a == "--max-miss") max_miss      = std::atof(next().c_str());
    else if (a == "--vr")       vr_arg        = next();
    else if (a == "--vrdraw")   vrdraw_bin    = next();
    else if (a == "--synth")    synth_arg     = next();
    else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
  }

  try {
    // ---------------- sample sets --------------------------------------
    auto fam_iids = read_fam_iids(fam_path);
    const std::size_t N = fam_iids.size();
    std::unordered_map<std::string, int> fam_pos;
    fam_pos.reserve(N * 2);
    for (std::size_t i = 0; i < N; ++i) fam_pos.emplace(fam_iids[i], int(i));

    std::vector<std::string>       trait_names;
    std::vector<std::vector<char>> member;   // P × N, 1 iff FAM row ∈ S_t

    if (!synth_arg.empty()) {
      auto p = split_delim(synth_arg, ',');
      if (p.size() != 3) throw std::runtime_error("--synth wants P,PCT,SEED");
      const int P    = std::atoi(p[0].c_str());
      const int pct  = std::atoi(p[1].c_str());
      const unsigned sd = static_cast<unsigned>(std::atol(p[2].c_str()));
      std::mt19937 gen(sd);
      for (int t = 0; t < P; ++t) {
        // Each trait drops its own random slice, size cycling through 1..pct %.
        const double frac = 0.01 * static_cast<double>(1 + (t % pct));
        std::bernoulli_distribution drop(frac);
        std::vector<char> m(N, 1);
        std::size_t ndrop = 0;
        for (std::size_t i = 0; i < N; ++i)
          if (drop(gen)) { m[i] = 0; ++ndrop; }
        if (ndrop == N) throw std::runtime_error("synthetic trait dropped everyone");
        trait_names.push_back("s" + std::to_string(t));
        member.push_back(std::move(m));
      }
    } else {
      std::ifstream pf(pheno_path);
      if (!pf) throw std::runtime_error("open " + pheno_path);
      std::string hdr;
      if (!std::getline(pf, hdr)) throw std::runtime_error("empty " + pheno_path);
      const char delim = detect_delim(hdr);
      auto cols = split_delim(hdr, delim);
      int iid_col = -1;
      for (std::size_t c = 0; c < cols.size(); ++c)
        if (cols[c] == "IID" || cols[c] == "iid") iid_col = int(c);
      if (iid_col < 0) throw std::runtime_error("no IID column in " + pheno_path);

      std::vector<int> tcol;
      if (!trait_arg.empty()) {
        for (const auto& nm : split_delim(trait_arg, ',')) {
          int c = -1;
          for (std::size_t k = 0; k < cols.size(); ++k) if (cols[k] == nm) c = int(k);
          if (c < 0) throw std::runtime_error("no such column: " + nm);
          tcol.push_back(c);
          trait_names.push_back(nm);
        }
      } else {
        for (std::size_t c = 0; c < cols.size(); ++c) {
          if (int(c) == iid_col) continue;
          if (cols[c].empty() || cols[c][0] == 'x' || cols[c] == "FID") continue;
          tcol.push_back(int(c));
          trait_names.push_back(cols[c]);
        }
      }
      member.assign(tcol.size(), std::vector<char>(N, 0));
      std::string line;
      while (std::getline(pf, line)) {
        auto row = split_delim(line, delim);
        if (int(row.size()) <= iid_col) continue;
        auto it = fam_pos.find(row[iid_col]);
        if (it == fam_pos.end()) continue;
        for (std::size_t t = 0; t < tcol.size(); ++t)
          if (tcol[t] < int(row.size()) && !is_na(row[tcol[t]]))
            member[t][it->second] = 1;
      }
    }

    const int P = static_cast<int>(trait_names.size());
    if (P < 1) throw std::runtime_error("no traits");

    // Union, in ascending FAM order.
    std::vector<int> ptr_U;
    for (std::size_t i = 0; i < N; ++i) {
      for (int t = 0; t < P; ++t)
        if (member[t][i]) { ptr_U.push_back(int(i) + 1); break; }
    }
    const std::size_t nU = ptr_U.size();
    if (nU == 0) throw std::runtime_error("empty union");

    // Per-trait union-local membership, exclusion lists, own ptrsub.
    std::vector<std::vector<unsigned char>> in_trait(P, std::vector<unsigned char>(nU, 0));
    std::vector<std::vector<int>>           excl(P);
    std::vector<std::vector<int>>           ptr_t(P);
    std::vector<int>                        n_t(P);
    for (int t = 0; t < P; ++t) {
      for (std::size_t k = 0; k < nU; ++k) {
        if (member[t][ptr_U[k] - 1]) {
          in_trait[t][k] = 1;
          ptr_t[t].push_back(ptr_U[k]);
        } else {
          excl[t].push_back(int(k));
        }
      }
      n_t[t] = int(ptr_t[t].size());
      if (n_t[t] == 0) throw std::runtime_error("trait " + trait_names[t] + " has no samples");
    }

    saige::ExclusionSets es;
    es.P    = P;
    es.excl = excl.data();
    es.n_t  = n_t.data();

    // ---------------- decode --------------------------------------------
    saige::BedReaderPool reader(bed_path, N, nthreads);
    std::size_t M = reader.n_markers_from_file_size();
    if (n_markers_arg > 0) M = std::min<std::size_t>(M, std::size_t(n_markers_arg));

    saige::VarRatioRule vr;
    std::vector<unsigned char> vr_drawn;
    if (!vr_arg.empty()) {
      auto p = split_delim(vr_arg, ',');
      vr.enabled = true;
      vr.min_mac = std::atof(p[0].c_str());
      vr.max_mac = (p.size() > 1) ? std::atof(p[1].c_str()) : -1.0f;
      vr_drawn   = load_vr_draw(vrdraw_bin, M);
    }

    std::printf("data   : %s\n", bed_path.c_str());
    std::printf("N(fam) = %zu   M = %zu   P = %d   |U| = %zu\n", N, M, P, nU);
    for (int t = 0; t < P; ++t)
      std::printf("  trait %-6s n_t=%7d  excluded from U = %zu (%.2f%%)\n",
                  trait_names[t].c_str(), n_t[t], excl[t].size(),
                  100.0 * double(excl[t].size()) / double(nU));
    std::printf("QC     : minMAF=%g maxMiss=%g  VR=%s\n", min_maf, max_miss,
                vr.enabled ? "on" : "off");

    std::fprintf(stderr, "[1/4] union decode (%d threads)...\n", nthreads);
    saige::ParallelDecodeAux pa;
    pa.excl                  = &es;
    pa.collect_tally         = true;
    pa.collect_missing_cells = true;
    pa.collect_keep          = true;
    auto uni = saige::parallel_decode_bed(
        reader, ptr_U.data(), nU, M, min_maf, max_miss, nthreads,
        vr, vr_drawn.empty() ? nullptr : vr_drawn.data(), &pa);

    std::vector<int> fill_U(M);
    for (std::size_t j = 0; j < M; ++j) fill_U[j] = uni.stats[j].fillin;

    // ---------------- A1 -------------------------------------------------
    std::fprintf(stderr, "[2/4] A1: per-trait stats vs solo decode...\n");
    Fail a1;
    std::vector<saige::TraitMarkerStats> ts(P);
    long a1_compared = 0;
    for (int t = 0; t < P; ++t) {
      saige::compute_trait_stats(uni.stats.data(), M, uni.tally.data(), P, t,
                                 n_t[t], min_maf, max_miss,
                                 vr, vr_drawn.empty() ? nullptr : vr_drawn.data(),
                                 ts[t]);

      // Reference: the trait's own sample set, unmodified code path.
      auto solo = saige::parallel_decode_bed(
          reader, ptr_t[t].data(), std::size_t(n_t[t]), M, min_maf, max_miss,
          nthreads, vr, vr_drawn.empty() ? nullptr : vr_drawn.data(), nullptr);

      int M_solo = 0;
      for (std::size_t j = 0; j < M; ++j) {
        const saige::MarkerStats& r = solo.stats[j];
        const float inv_ref = r.passQC
            ? ((std::sqrt(2.0f * r.altFreq * (1.0f - r.altFreq)) == 0.0f)
                 ? 0.0f : 1.0f / std::sqrt(2.0f * r.altFreq * (1.0f - r.altFreq)))
            : 0.0f;
        if (r.passQC) ++M_solo;
        ++a1_compared;

        if (!fsame(ts[t].freq[j], r.altFreq))
          a1.report("freq", trait_names[t], j, f2s(ts[t].freq[j]), f2s(r.altFreq));
        if (!fsame(ts[t].invstd[j], inv_ref))
          a1.report("invstd", trait_names[t], j, f2s(ts[t].invstd[j]), f2s(inv_ref));
        if (!fsame(ts[t].missingRate[j], r.missingRate))
          a1.report("missingRate", trait_names[t], j,
                    f2s(ts[t].missingRate[j]), f2s(r.missingRate));
        if (ts[t].fill[j] != r.fillin)
          a1.report("fill", trait_names[t], j, std::to_string(ts[t].fill[j]),
                    std::to_string(r.fillin));
        if ((ts[t].passQC[j] != 0) != r.passQC)
          a1.report("passQC", trait_names[t], j, std::to_string(int(ts[t].passQC[j])),
                    std::to_string(int(r.passQC)));
        if ((ts[t].passVR[j] != 0) != bool(solo.passVR[j]))
          a1.report("passVR", trait_names[t], j, std::to_string(int(ts[t].passVR[j])),
                    std::to_string(int(solo.passVR[j])));
        if (ts[t].mac[j] != r.mac)
          a1.report("mac", trait_names[t], j, std::to_string(ts[t].mac[j]),
                    std::to_string(r.mac));
        if (ts[t].numMissing[j] != r.numMissing)
          a1.report("numMissing", trait_names[t], j, std::to_string(ts[t].numMissing[j]),
                    std::to_string(r.numMissing));
        if (ts[t].alleleCount[j] != r.alleleCount)
          a1.report("alleleCount", trait_names[t], j,
                    std::to_string(ts[t].alleleCount[j]), std::to_string(r.alleleCount));
      }
      if (ts[t].M_t != M_solo)
        a1.report("M_t", trait_names[t], 0, std::to_string(ts[t].M_t),
                  std::to_string(M_solo));
      std::printf("  A1 %-6s  M_t=%d (solo %d)  markers=%zu  %s\n",
                  trait_names[t].c_str(), ts[t].M_t, M_solo, M,
                  a1.n == 0 ? "bit-identical" : "MISMATCH");
    }

    // ---------------- §2 union-keep rule ---------------------------------
    std::fprintf(stderr, "[3/4] union-pack keep rule + pack_if_keep store...\n");
    int keep_mismatch = 0;
    std::size_t n_keep = 0, n_union_passQC = 0;
    std::vector<std::size_t> keep_idx;
    for (std::size_t j = 0; j < M; ++j) {
      bool any = false;
      for (int t = 0; t < P; ++t) if (ts[t].passQC[j]) { any = true; break; }
      if (any) { ++n_keep; keep_idx.push_back(j); }
      if (uni.stats[j].passQC) ++n_union_passQC;
      if ((uni.keep_union[j] != 0) != any) {
        if (keep_mismatch < 8)
          std::printf("  FAIL [keep] marker=%zu decoder=%d  OR_t passQC_t=%d\n",
                      j, int(uni.keep_union[j]), int(any));
        ++keep_mismatch;
      }
    }
    std::printf("  keep (any trait passQC) = %zu / %zu   union's own passQC = %zu\n",
                n_keep, M, n_union_passQC);

    // pack_if_keep decode: same rows, byte-identical to the ordinary pack.
    saige::ParallelDecodeAux pk = pa;
    pk.pack_if_keep = true;
    auto uni_keep = saige::parallel_decode_bed(
        reader, ptr_U.data(), nU, M, min_maf, max_miss, nthreads,
        vr, vr_drawn.empty() ? nullptr : vr_drawn.data(), &pk);
    int pack_fail = 0;
    if (uni_keep.store.n_stored() != n_keep) {
      std::printf("  FAIL [pack] store holds %zu rows, keep rule says %zu\n",
                  uni_keep.store.n_stored(), n_keep);
      ++pack_fail;
    } else {
      // index list matches, and every row also present in the ordinary pack is
      // byte-identical there (the packed content does not depend on the rule).
      std::unordered_map<std::size_t, std::size_t> ord;
      for (std::size_t r = 0; r < uni.orig_plink_idx.size(); ++r)
        ord.emplace(uni.orig_plink_idx[r], r);
      const std::size_t nb = uni_keep.store.nbyte();
      for (std::size_t r = 0; r < n_keep; ++r) {
        if (uni_keep.orig_plink_idx[r] != keep_idx[r]) {
          if (pack_fail < 8)
            std::printf("  FAIL [pack] row %zu holds marker %zu, expected %zu\n",
                        r, uni_keep.orig_plink_idx[r], keep_idx[r]);
          ++pack_fail;
          continue;
        }
        auto it = ord.find(keep_idx[r]);
        if (it == ord.end()) continue;   // union's own QC dropped it; nothing to compare
        if (std::memcmp(uni_keep.store.span(r).data,
                        uni.store.span(it->second).data, nb) != 0) {
          if (pack_fail < 8)
            std::printf("  FAIL [pack] marker %zu packed bytes differ\n", keep_idx[r]);
          ++pack_fail;
        }
      }
    }

    // ---------------- A2 + missing-cell ground truth ---------------------
    std::fprintf(stderr, "[4/4] A2: fill corrections vs brute-force decode...\n");
    std::vector<saige::FillCorrection> corr(P);
    for (int t = 0; t < P; ++t)
      saige::build_fill_corrections(fill_U.data(), ts[t].fill.data(), M,
                                    uni.missing_cells, in_trait[t].data(), corr[t]);

    // Brute force: re-decode the BED from scratch, build the expected
    // correction lists and the expected missing-cell lists.
    std::vector<std::vector<int>>   exp_row(P), exp_col(P);
    std::vector<std::vector<float>> exp_delta(P);
    saige::BedLut lut;
    saige::build_bed_lookup(lut);
    std::vector<unsigned char> raw(reader.n_bytes_per_marker());
    std::vector<int>           geno(N, 0);
    std::vector<int>           miss_rows;
    int mcell_fail = 0;
    long n_missing_cells = 0;
    for (std::size_t j = 0; j < M; ++j) {
      reader.read_marker(0, j, raw.data());
      for (std::size_t b = 0; b < raw.size(); ++b) {
        const int* l = lut[raw[b]];
        for (int q = 0; q < 4; ++q) {
          const std::size_t fi = b * 4 + q;
          if (fi >= N) break;
          geno[fi] = l[q];
        }
      }
      miss_rows.clear();
      for (std::size_t k = 0; k < nU; ++k)
        if (geno[ptr_U[k] - 1] == 3) miss_rows.push_back(int(k));
      n_missing_cells += long(miss_rows.size());

      if (uni.missing_cells[j] != miss_rows) {
        if (mcell_fail < 8)
          std::printf("  FAIL [missing_cells] marker=%zu decoder has %zu rows, "
                      "brute force %zu\n", j, uni.missing_cells[j].size(),
                      miss_rows.size());
        ++mcell_fail;
      }
      for (int t = 0; t < P; ++t) {
        const int d = ts[t].fill[j] - fill_U[j];
        if (d == 0) continue;
        for (int k : miss_rows) {
          if (!in_trait[t][k]) continue;
          exp_row[t].push_back(k);
          exp_col[t].push_back(int(j));
          exp_delta[t].push_back(float(d));
        }
      }
    }

    Fail a2;
    for (int t = 0; t < P; ++t) {
      if (corr[t].size() != exp_col[t].size()) {
        a2.report("corr size", trait_names[t], 0,
                  std::to_string(corr[t].size()), std::to_string(exp_col[t].size()));
      } else {
        for (std::size_t e = 0; e < corr[t].size(); ++e) {
          if (corr[t].row[e] != exp_row[t][e] || corr[t].col[e] != exp_col[t][e] ||
              corr[t].delta[e] != exp_delta[t][e]) {
            a2.report("corr entry", trait_names[t], std::size_t(corr[t].col[e]),
                      "(" + std::to_string(corr[t].row[e]) + "," +
                      std::to_string(corr[t].col[e]) + "," +
                      std::to_string(corr[t].delta[e]) + ")",
                      "(" + std::to_string(exp_row[t][e]) + "," +
                      std::to_string(exp_col[t][e]) + "," +
                      std::to_string(exp_delta[t][e]) + ")");
            continue;
          }
          // the contract's own identity
          const int j = corr[t].col[e];
          if (fill_U[j] + int(corr[t].delta[e]) != ts[t].fill[j])
            a2.report("fill_U+delta", trait_names[t], std::size_t(j),
                      std::to_string(fill_U[j] + int(corr[t].delta[e])),
                      std::to_string(ts[t].fill[j]));
          if (!in_trait[t][corr[t].row[e]])
            a2.report("row not in S_t", trait_names[t], std::size_t(j),
                      std::to_string(corr[t].row[e]), "in S_t");
        }
      }
      std::size_t ncols = 0;
      for (std::size_t j = 0; j < M; ++j) if (ts[t].fill[j] != fill_U[j]) ++ncols;
      std::printf("  A2 %-6s  corrections=%zu over %zu markers with fill_t != fill_U  %s\n",
                  trait_names[t].c_str(), corr[t].size(), ncols,
                  a2.n == 0 ? "exact" : "MISMATCH");
    }

    std::printf("\nmissing cells in the union pack: %ld  (over %zu markers)\n",
                n_missing_cells, M);
    std::printf("A1 mismatches                  : %d  (%ld marker-trait comparisons, "
                "9 fields each)\n", a1.n, a1_compared);
    std::printf("A2 mismatches                  : %d\n", a2.n);
    std::printf("union keep-rule mismatches     : %d\n", keep_mismatch);
    std::printf("pack_if_keep store mismatches  : %d\n", pack_fail);
    std::printf("missing-cell list mismatches   : %d\n", mcell_fail);

    const bool ok = (a1.n == 0 && a2.n == 0 && keep_mismatch == 0 &&
                     pack_fail == 0 && mcell_fail == 0);
    std::printf("%s\n", ok ? "MASK_STATS_TEST PASS" : "MASK_STATS_TEST FAIL");
    return ok ? 0 : 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "mask_stats_test: %s\n", e.what());
    return 2;
  }
}
