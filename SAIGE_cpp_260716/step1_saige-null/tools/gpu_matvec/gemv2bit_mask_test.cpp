// gemv2bit_mask_test.cpp — acceptance B1/B2/B3 for scheme C stage B
// (optimization/missing_mt/SCHEME_C_DESIGN.md §5).
//
// This test talks to the tier-4 kernels directly (gemv2bit.hpp), not through
// the gpu_matvec facade: bind_trait's mask and correction lists have no facade
// entry point yet, and the point of the exercise is the kernel identity.
//
//   B1  bind_trait with no mask and no corrections must be BIT-IDENTICAL to
//       the pre-scheme-C kernel. The comparison is against float bits captured
//       from the old code (git ce9bacd6) into /opt/saige/logs/scheme_c/b1_ref/,
//       not against anything this binary computes.
//   B2  masked + corrected K·b against a float64 CPU reference written
//       straight from §1's identity.
//   B3  masked K·b on the union matrix vs the same phenotype's K·b on a matrix
//       built from only its own rows — the two differ ONLY by fp32 reduction
//       order (§4), so this measures that order effect and nothing else.
//   T   wall-time cost of the mask + correction scatter-adds at a realistic
//       size, since the design budgets ~5% for masking.
//
// Usage: ./gemv2bit_mask_test [--ref DIR] [--no-timing]
#include "gemv2bit.hpp"
#include "g2b_testdata.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace g2btest;
namespace g2b = saige::gpu::g2b;

namespace {

int g_fail = 0;

const char* kDefaultRef = "/opt/saige/logs/scheme_c/b1_ref";

// ------------------------------------------------------------------ helpers
bool load_ref(const std::string& dir, const std::string& name,
              std::size_t n, std::vector<float>& out) {
  const std::string path = dir + "/" + name + ".bin";
  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) return false;
  out.resize(n);
  const std::size_t got = std::fread(out.data(), sizeof(float), n, f);
  std::fclose(f);
  return got == n;
}

// Bit-level comparison: memcmp per word so a NaN or a −0/+0 flip is caught too.
int bitdiff(const std::vector<float>& a, const std::vector<float>& b) {
  int n = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::memcmp(&a[i], &b[i], sizeof(float)) != 0) ++n;
  return n;
}

struct Rel { double max_abs, rel_max, rel_l2; };

// ref is float64, got is float32 — compare in double.
Rel relerr(const std::vector<double>& ref, const std::vector<float>& got) {
  double max_abs = 0, max_ref = 0, num = 0, den = 0;
  for (std::size_t i = 0; i < ref.size(); ++i) {
    const double d = static_cast<double>(got[i]) - ref[i];
    max_abs = std::max(max_abs, std::fabs(d));
    max_ref = std::max(max_ref, std::fabs(ref[i]));
    num += d * d;
    den += ref[i] * ref[i];
  }
  return {max_abs, max_abs / std::max(1e-300, max_ref),
          std::sqrt(num) / std::max(1e-300, std::sqrt(den))};
}

// ------------------------------------------------- the phenotype's own view
struct Trait {
  std::vector<int>   mask;                      // union rows NOT owned, ascending
  std::vector<unsigned char> is_mask;           // N flags
  std::vector<int>   crow, ccol;                // corrections, sorted by (col,row)
  std::vector<float> cdel;
  std::vector<float> freq, invstd;              // length M; invstd 0 = QC dropped
  float              inv_M = 1.f;
  int                M_t = 0;
};

// ~20% of rows masked, `n_corr` correction cells, ~5% of markers QC-dropped.
// Δ is drawn only from the values that keep g+Δ inside {0,1,2}: that is exactly
// the real constraint (Δ = fill_t − fill_U, both in {0,1,2}) and it is what
// makes B3's own-rows matrix well defined.
Trait make_trait(const std::vector<unsigned char>& packed, std::size_t stride,
                 int N, int M, const std::vector<float>& freq0,
                 const std::vector<float>& invstd0,
                 double mask_frac, int n_corr, double drop_frac,
                 std::uint64_t seed) {
  Trait t;
  Rng rng(seed);
  t.is_mask.assign(N, 0);
  const int want = static_cast<int>(mask_frac * N);
  int have = 0;
  while (have < want) {
    const int i = rng.below(N);
    if (!t.is_mask[i]) { t.is_mask[i] = 1; ++have; }
  }
  for (int i = 0; i < N; ++i) if (t.is_mask[i]) t.mask.push_back(i);

  t.freq = freq0;
  t.invstd = invstd0;
  const int ndrop = static_cast<int>(drop_frac * M);
  for (int k = 0; k < ndrop; ++k) t.invstd[rng.below(M)] = 0.0f;
  t.M_t = 0;
  for (int j = 0; j < M; ++j) if (t.invstd[j] != 0.0f) ++t.M_t;
  t.inv_M = 1.0f / static_cast<float>(t.M_t);

  // Distinct (col,row) cells on unmasked rows, then sorted by (col,row) as
  // §3.2 promises the real builder will deliver.
  std::vector<std::pair<int,int>> cells;   // (col, row)
  cells.reserve(n_corr);
  while (static_cast<int>(cells.size()) < n_corr) {
    const int i = rng.below(N);
    if (t.is_mask[i]) continue;
    const int j = rng.below(M);
    cells.emplace_back(j, i);
  }
  std::sort(cells.begin(), cells.end());
  cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
  for (auto& cj : cells) {
    const int g = geno_at(packed.data(), stride, cj.first, cj.second);
    // g=0 → Δ∈{1,2}; g=1 → Δ∈{−1,1}; g=2 → Δ∈{−2,−1}
    const int opt[3][2] = {{1,2},{-1,1},{-2,-1}};
    const int d = opt[g][rng.below(2)];
    t.ccol.push_back(cj.first);
    t.crow.push_back(cj.second);
    t.cdel.push_back(static_cast<float>(d));
  }
  return t;
}

// ------------------------------------- float64 reference, straight from §1
// raw[j] = Σ_i g[i,j] x[i] + Σ_{(i,j)} Δ·x[i]
// Y[j]   = invstd[j]·(raw[j] − 2 f[j] S),  S = Σ_i x[i]
// w[j]   = invstd[j]·Y[j]
// Z[i]   = Σ_j g[i,j] w[j] + Σ_{(i,j)} Δ·w[j] − C,  C = Σ_j 2 f[j] w[j]
// ret    = Z · inv_M,  and 0 on the masked rows (see gemv2bit.hpp on why).
void ref_matvec_range(const std::vector<unsigned char>& packed, std::size_t stride,
                      int N, int M, int j0, int jn,
                      const Trait& t, const std::vector<float>& b,
                      std::vector<double>& ret);

void ref_matvec(const std::vector<unsigned char>& packed, std::size_t stride,
                int N, int M, const Trait& t, const std::vector<float>& b,
                std::vector<double>& ret) {
  ref_matvec_range(packed, stride, N, M, 0, M, t, b, ret);
}

// [j0, j0+jn) only. Corrections whose marker falls outside the range are
// skipped — matching what the kernel does, which is what makes a LOCO slice
// consistent with the full call.
void ref_matvec_range(const std::vector<unsigned char>& packed, std::size_t stride,
                      int N, int M, int j0, int jn,
                      const Trait& t, const std::vector<float>& b,
                      std::vector<double>& ret) {
  (void)M;
  std::vector<double> x(N);
  double S = 0;
  for (int i = 0; i < N; ++i) {
    x[i] = t.is_mask.empty() || !t.is_mask[i] ? static_cast<double>(b[i]) : 0.0;
    S += x[i];
  }

  std::vector<double> w(j0 + jn, 0.0);
  for (int j = j0; j < j0 + jn; ++j) {
    double raw = 0;
    for (int i = 0; i < N; ++i)
      raw += static_cast<double>(geno_at(packed.data(), stride, j, i)) * x[i];
    for (std::size_t k = 0; k < t.cdel.size(); ++k)
      if (t.ccol[k] == j) raw += static_cast<double>(t.cdel[k]) * x[t.crow[k]];
    const double y = static_cast<double>(t.invstd[j]) *
                     (raw - 2.0 * static_cast<double>(t.freq[j]) * S);
    w[j] = static_cast<double>(t.invstd[j]) * y;
  }

  double C = 0;
  for (int j = j0; j < j0 + jn; ++j) C += 2.0 * static_cast<double>(t.freq[j]) * w[j];

  ret.assign(N, 0.0);
  for (int j = j0; j < j0 + jn; ++j) {
    if (w[j] == 0.0) continue;
    for (int i = 0; i < N; ++i)
      ret[i] += static_cast<double>(geno_at(packed.data(), stride, j, i)) * w[j];
  }
  for (std::size_t k = 0; k < t.cdel.size(); ++k)
    if (t.ccol[k] >= j0 && t.ccol[k] < j0 + jn)
      ret[t.crow[k]] += static_cast<double>(t.cdel[k]) * w[t.ccol[k]];
  for (int i = 0; i < N; ++i) {
    ret[i] = (ret[i] - C) * static_cast<double>(t.inv_M);
    if (!t.is_mask.empty() && t.is_mask[i]) ret[i] = 0.0;
  }
}

}  // namespace

// ============================================================== B1
static void b1(const std::string& refdir) {
  std::printf("\n================ B1: degenerate bind == pre-scheme-C bits ==========\n");
  int ncase = 0;
  const Case* cs = cases(&ncase);
  const int kNcol[] = {1, 2, 3, 5, 8, 9};
  int total_words = 0, total_diff = 0, missing = 0;

  for (int k = 0; k < ncase; ++k) {
    const int N = cs[k].N, M = cs[k].M;
    std::vector<unsigned char> packed; std::size_t stride = 0;
    std::vector<float> freq, invstd;
    gen_matrix(N, M, cs[k].seed, packed, stride, freq, invstd);

    g2b::Ctx* c = g2b::create(packed.data(), stride, N, M);
    if (!c) { std::fprintf(stderr, "B1 create failed (case %d)\n", k); ++g_fail; return; }
    g2b::TraitBind* t = g2b::bind_trait(c, freq.data(), invstd.data(), nullptr, 0,
                                        nullptr, nullptr, nullptr, 0);
    if (!t) { std::fprintf(stderr, "B1 bind failed (case %d)\n", k); ++g_fail; g2b::destroy(c); return; }

    const std::string tag = "case" + std::to_string(k);
    std::vector<float> x; gen_vec(N, cs[k].seed ^ 0xB1ull, x);
    std::vector<float> got(N, 0.f), ref;

    auto check = [&](const std::string& name, const std::vector<float>& v) {
      if (!load_ref(refdir, name, v.size(), ref)) {
        std::printf("  %-20s REFERENCE MISSING (%s/%s.bin)\n", name.c_str(),
                    refdir.c_str(), name.c_str());
        ++missing; return;
      }
      const int d = bitdiff(ref, v);
      total_words += (int)v.size(); total_diff += d;
      std::printf("  %-20s n=%-7zu words differing = %d  %s\n", name.c_str(),
                  v.size(), d, d == 0 ? "bit-identical" : "*** MISMATCH ***");
      if (d) ++g_fail;
    };

    if (!g2b::matvec_range(c, t, 0, M, 1.0f/(float)M, x.data(), got.data())) { ++g_fail; }
    check(tag + "_range_full", got);

    // Same bind, different range, different 1/M_t — which is the whole reason
    // inv_M is a call argument and not part of the bind (a LOCO run does this
    // 23 times per phenotype).
    const int j0 = M/3, jn = M/2;
    std::fill(got.begin(), got.end(), 0.f);
    if (!g2b::matvec_range(c, t, j0, jn, 1.0f/(float)jn, x.data(), got.data())) { ++g_fail; }
    check(tag + "_range_sub", got);

    for (int nc : kNcol) {
      std::vector<float> X; gen_vec((std::size_t)N*nc, cs[k].seed ^ (0xC0ull + nc), X);
      std::vector<float> KU((std::size_t)N*nc, 0.f);
      if (!g2b::matvec_mat(c, t, nc, 1.0f/(float)M, X.data(), KU.data())) { ++g_fail; }
      check(tag + "_mat" + std::to_string(nc), KU);
    }

    g2b::unbind_trait(t);
    g2b::destroy(c);
  }
  if (missing)
    std::printf("  %d reference file(s) missing — B1 INCONCLUSIVE, rerun the capture\n", missing);
  std::printf("B1: %d/%d words differ from the pre-scheme-C baseline  %s\n",
              total_diff, total_words,
              (total_diff == 0 && missing == 0) ? "PASS" : "FAIL");
  if (missing) ++g_fail;
}

// ============================================================== B2
static void b2() {
  std::printf("\n================ B2: masked + corrected vs float64 reference =======\n");
  const int N = 20000, M = 3000;
  std::vector<unsigned char> packed; std::size_t stride = 0;
  std::vector<float> freq, invstd;
  gen_matrix(N, M, 0x51A6E2ull, packed, stride, freq, invstd);

  Trait tr = make_trait(packed, stride, N, M, freq, invstd,
                        0.20, 400, 0.05, 0xB2005EEDull);
  std::printf("  N=%d M=%d  masked rows=%zu (%.1f%%)  corrections=%zu  "
              "QC-dropped markers=%d  M_t=%d\n",
              N, M, tr.mask.size(), 100.0*tr.mask.size()/N, tr.cdel.size(),
              M - tr.M_t, tr.M_t);

  g2b::Ctx* c = g2b::create(packed.data(), stride, N, M);
  g2b::TraitBind* t = g2b::bind_trait(
      c, tr.freq.data(), tr.invstd.data(),
      tr.mask.data(), (int)tr.mask.size(),
      tr.crow.data(), tr.ccol.data(), tr.cdel.data(), (int)tr.cdel.size());
  if (!c || !t) { std::fprintf(stderr, "B2 setup failed\n"); ++g_fail; return; }

  std::vector<float> b; gen_vec(N, 0xB2B2ull, b);
  std::vector<double> ref;
  ref_matvec(packed, stride, N, M, tr, b, ref);

  std::vector<float> got(N, 0.f);
  if (!g2b::matvec_range(c, t, 0, M, tr.inv_M, b.data(), got.data())) { ++g_fail; }
  Rel r = relerr(ref, got);
  // GATED ON rel_L2 ONLY. rel_max is printed because it is the number that
  // says how bad a single element can get, but it is not a threshold: fp32
  // reduction error grows like √M, so a max-norm gate calibrated at M=3000
  // would fail at UKB's M=1.1e5 for no reason other than the problem being
  // bigger. Same reason §4 refuses to gate anything bit-exact here.
  std::printf("  matvec_range : max|Δ|=%.4g  rel_max=%.4g (ungated)  rel_L2=%.4g  %s\n",
              r.max_abs, r.rel_max, r.rel_l2, r.rel_l2 <= 1e-6 ? "PASS" : "FAIL(rel_L2>1e-6)");
  if (r.rel_l2 > 1e-6) ++g_fail;

  // Sub-range under mask + corrections: the marker CSR is offset to j0 and the
  // sample CSR is filtered by the range test, so both need their own check.
  {
    const int j0 = M/3, jn = M/2;
    std::vector<double> rs;
    Trait trng = tr;
    {
      int mt = 0;
      for (int j = j0; j < j0+jn; ++j) if (tr.invstd[j] != 0.0f) ++mt;
      trng.inv_M = 1.0f/(float)mt;
    }
    ref_matvec_range(packed, stride, N, M, j0, jn, trng, b, rs);
    std::vector<float> gs(N, 0.f);
    // Same bind, a different range and a different inv_M.
    int Mt_rng = 0;
    for (int j = j0; j < j0+jn; ++j) if (tr.invstd[j] != 0.0f) ++Mt_rng;
    if (!g2b::matvec_range(c, t, j0, jn, 1.0f/(float)Mt_rng, b.data(), gs.data())) ++g_fail;
    Rel r2 = relerr(rs, gs);
    std::printf("  range[%d,%d) : max|Δ|=%.4g  rel_max=%.4g (ungated)  rel_L2=%.4g  %s\n",
                j0, j0+jn, r2.max_abs, r2.rel_max, r2.rel_l2,
                r2.rel_l2 <= 1e-6 ? "PASS" : "FAIL(rel_L2>1e-6)");
    if (r2.rel_l2 > 1e-6) ++g_fail;
  }

  // The caller must not have to pre-mask, and must get zeros back on the rows
  // it does not own. Both are contract, so both get asserted.
  int nonzero_masked = 0;
  for (int i : tr.mask) if (got[i] != 0.0f) ++nonzero_masked;
  std::printf("  masked rows nonzero in ret: %d  %s\n", nonzero_masked,
              nonzero_masked ? "*** FAIL ***" : "PASS");
  if (nonzero_masked) ++g_fail;

  // Determinism with the scatter-adds live.
  std::vector<float> got2(N, 0.f);
  g2b::matvec_range(c, t, 0, M, tr.inv_M, b.data(), got2.data());
  const int bd = bitdiff(got, got2);
  std::printf("  inter-run: %d/%d words differ  %s\n", bd, N,
              bd ? "*** NONDETERMINISTIC ***" : "bit-identical");
  if (bd) ++g_fail;

  // Negative control: a gate that cannot fail is not a gate. Re-run the
  // reference with the corrections deleted and with nothing masked, and
  // require the GPU result NOT to match either — otherwise B2 would pass just
  // as happily on a kernel that silently ignored both lists.
  {
    Trait nocorr = tr; nocorr.crow.clear(); nocorr.ccol.clear(); nocorr.cdel.clear();
    std::vector<double> r2;
    ref_matvec(packed, stride, N, M, nocorr, b, r2);
    Rel rn = relerr(r2, got);
    std::printf("  neg-control (no corrections in ref): rel_L2=%.4g  %s\n",
                rn.rel_l2, rn.rel_l2 > 1e-5 ? "PASS (corrections do bite)"
                                            : "*** FAIL: corrections are a no-op ***");
    if (rn.rel_l2 <= 1e-5) ++g_fail;

    Trait nomask = tr; nomask.is_mask.assign(N, 0);
    ref_matvec(packed, stride, N, M, nomask, b, r2);
    Rel rm = relerr(r2, got);
    std::printf("  neg-control (no mask in ref):        rel_L2=%.4g  %s\n",
                rm.rel_l2, rm.rel_l2 > 1e-5 ? "PASS (mask does bite)"
                                            : "*** FAIL: mask is a no-op ***");
    if (rm.rel_l2 <= 1e-5) ++g_fail;
  }

  // Same check through the batch path, so the _mc finish kernels' scatter-adds
  // are covered too. ncol=3 pads to NC=4 — the zero columns must stay zero even
  // with corrections live.
  const int nc = 3;
  std::vector<float> X((std::size_t)N*nc), KU((std::size_t)N*nc, 0.f);
  for (int cc = 0; cc < nc; ++cc) {
    std::vector<float> col; gen_vec(N, 0xB2C0ull + cc, col);
    std::copy(col.begin(), col.end(), X.begin() + (std::size_t)cc*N);
  }
  if (!g2b::matvec_mat(c, t, nc, tr.inv_M, X.data(), KU.data())) { ++g_fail; }
  double worst_l2 = 0, worst_max = 0;
  for (int cc = 0; cc < nc; ++cc) {
    std::vector<float> col(X.begin() + (std::size_t)cc*N,
                           X.begin() + (std::size_t)(cc+1)*N);
    ref_matvec(packed, stride, N, M, tr, col, ref);
    std::vector<float> gc(KU.begin() + (std::size_t)cc*N,
                          KU.begin() + (std::size_t)(cc+1)*N);
    Rel rc = relerr(ref, gc);
    worst_l2 = std::max(worst_l2, rc.rel_l2);
    worst_max = std::max(worst_max, rc.rel_max);
  }
  std::printf("  matvec_mat(3): worst rel_max=%.4g (ungated)  worst rel_L2=%.4g  %s\n",
              worst_max, worst_l2, worst_l2 <= 1e-6 ? "PASS" : "FAIL(rel_L2>1e-6)");
  if (worst_l2 > 1e-6) ++g_fail;

  g2b::unbind_trait(t);
  g2b::destroy(c);
}

// ============================================================== B3
static void b3() {
  std::printf("\n================ B3: union+mask vs own-rows matrix =================\n");
  const int N = 20000, M = 3000;
  std::vector<unsigned char> packed; std::size_t stride = 0;
  std::vector<float> freq, invstd;
  gen_matrix(N, M, 0x51A6E2ull, packed, stride, freq, invstd);

  Trait tr = make_trait(packed, stride, N, M, freq, invstd,
                        0.20, 400, 0.05, 0xB3005EEDull);

  // Own-rows matrix: this phenotype's rows only, with its OWN fill baked in
  // (g + Δ at every correction cell), so it needs neither mask nor corrections.
  std::vector<int> local(N, -1);
  int Nt = 0;
  for (int i = 0; i < N; ++i) if (!tr.is_mask[i]) local[i] = Nt++;
  const std::size_t stride_t = nbyte_of(Nt);
  std::vector<unsigned char> own((std::size_t)M * stride_t, 0xFF);
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < N; ++i)
      if (local[i] >= 0)
        set_geno(own.data(), stride_t, j, local[i],
                 geno_at(packed.data(), stride, j, i));
  for (std::size_t k = 0; k < tr.cdel.size(); ++k) {
    const int i = tr.crow[k], j = tr.ccol[k];
    const int g = geno_at(packed.data(), stride, j, i) + (int)tr.cdel[k];
    set_geno(own.data(), stride_t, j, local[i], g);
  }
  std::printf("  union N=%d -> own N_t=%d  M=%d  corrections=%zu  M_t=%d\n",
              N, Nt, M, tr.cdel.size(), tr.M_t);

  std::vector<float> b; gen_vec(N, 0xB3B3ull, b);
  std::vector<float> b_own(Nt);
  for (int i = 0; i < N; ++i) if (local[i] >= 0) b_own[local[i]] = b[i];

  g2b::Ctx* cu = g2b::create(packed.data(), stride, N, M);
  g2b::TraitBind* tu = g2b::bind_trait(
      cu, tr.freq.data(), tr.invstd.data(),
      tr.mask.data(), (int)tr.mask.size(),
      tr.crow.data(), tr.ccol.data(), tr.cdel.data(), (int)tr.cdel.size());
  g2b::Ctx* co = g2b::create(own.data(), stride_t, Nt, M);
  g2b::TraitBind* to = g2b::bind_trait(co, tr.freq.data(), tr.invstd.data(),
                                       nullptr, 0,
                                       nullptr, nullptr, nullptr, 0);
  if (!cu || !tu || !co || !to) { std::fprintf(stderr, "B3 setup failed\n"); ++g_fail; return; }

  std::vector<float> ru(N, 0.f), ro(Nt, 0.f);
  if (!g2b::matvec_range(cu, tu, 0, M, tr.inv_M, b.data(), ru.data())) ++g_fail;
  if (!g2b::matvec_range(co, to, 0, M, tr.inv_M, b_own.data(), ro.data())) ++g_fail;

  double max_abs = 0, max_ref = 0, num = 0, den = 0;
  int arg = -1;
  for (int i = 0; i < N; ++i) {
    if (local[i] < 0) continue;
    const double d = (double)ru[i] - (double)ro[local[i]];
    if (std::fabs(d) > max_abs) { max_abs = std::fabs(d); arg = i; }
    max_ref = std::max(max_ref, std::fabs((double)ro[local[i]]));
    num += d*d; den += (double)ro[local[i]]*(double)ro[local[i]];
  }
  const double rel_max = max_abs / std::max(1e-300, max_ref);
  const double rel_l2  = std::sqrt(num) / std::max(1e-300, std::sqrt(den));
  std::printf("  max|Δ|=%.6g at union row %d (value %.6g)  rel_max=%.4g  rel_L2=%.4g\n",
              max_abs, arg, arg >= 0 ? ro[local[arg]] : 0.f, rel_max, rel_l2);
  std::printf("  B3 %s (threshold rel ≤ 1e-5)\n",
              (rel_max <= 1e-5 && rel_l2 <= 1e-5) ? "PASS" : "FAIL");
  if (rel_max > 1e-5 || rel_l2 > 1e-5) ++g_fail;

  g2b::unbind_trait(tu); g2b::destroy(cu);
  g2b::unbind_trait(to); g2b::destroy(co);
}

// ============================================================== timing
// Fast packed generator — this matrix is only ever fed to the GPU, so it does
// not need to match the B1 baseline and can be built a byte at a time.
static void gen_packed_fast(int N, int M, std::uint64_t seed,
                            std::vector<unsigned char>& packed, std::size_t& stride) {
  stride = nbyte_of(N);
  packed.assign((std::size_t)M * stride, 0xFF);
  Rng rng(seed);
  static const unsigned char kCode[3] = {0x0, 0x2, 0x3};
  for (std::size_t p = 0; p < packed.size(); ++p) {
    const std::uint32_t r = rng.u32();
    unsigned char by = 0;
    for (int k = 0; k < 4; ++k) by |= (unsigned char)(kCode[(r >> (8*k)) % 3] << (2*k));
    packed[p] = by;
  }
}

static void timing() {
  std::printf("\n================ T: mask + correction overhead =====================\n");
  const int N = 50000, M = 40000;
  std::vector<unsigned char> packed; std::size_t stride = 0;
  gen_packed_fast(N, M, 0x7717ull, packed, stride);
  std::vector<float> freq(M), invstd(M);
  for (int j = 0; j < M; ++j) {
    freq[j] = 0.03f + (j % 100) * 0.0037f;
    invstd[j] = 1.0f/std::sqrt(2.0f*freq[j]*(1.0f-freq[j]));
  }
  Trait tr = make_trait(packed, stride, N, M, freq, invstd, 0.20, 400, 0.0, 0x7718ull);

  g2b::Ctx* c = g2b::create(packed.data(), stride, N, M);
  if (!c) { std::fprintf(stderr, "timing create failed (VRAM?)\n"); ++g_fail; return; }
  g2b::TraitBind* plain = g2b::bind_trait(c, freq.data(), invstd.data(), nullptr, 0,
                                          nullptr, nullptr, nullptr, 0);
  g2b::TraitBind* full = g2b::bind_trait(
      c, tr.freq.data(), tr.invstd.data(),
      tr.mask.data(), (int)tr.mask.size(),
      tr.crow.data(), tr.ccol.data(), tr.cdel.data(), (int)tr.cdel.size());
  g2b::TraitBind* maskonly = g2b::bind_trait(
      c, tr.freq.data(), tr.invstd.data(),
      tr.mask.data(), (int)tr.mask.size(), nullptr, nullptr, nullptr, 0);
  if (!plain || !full || !maskonly) { std::fprintf(stderr, "timing bind failed\n"); ++g_fail; return; }

  std::vector<float> x; gen_vec(N, 0x7719ull, x);
  std::vector<float> out(N, 0.f);

  std::printf("  N=%d M=%d  packed=%zu MB  masked=%zu rows  corrections=%zu\n",
              N, M, ((std::size_t)M*stride) >> 20, tr.mask.size(), tr.cdel.size());

  // Variants are INTERLEAVED, one matvec each per round. Measured back to back
  // in three separate loops the numbers moved by ±10% — enough to make the
  // corrected bind look faster than the plain one — because the V100's SM
  // clock drifts over a multi-second run. Round-robin puts every variant on
  // the same clock trajectory, and the per-variant min is reported alongside
  // the mean as the drift-free figure.
  g2b::TraitBind* v[3]   = {plain, maskonly, full};
  const char*     nm[3]  = {"no mask, no corrections", "mask only",
                            "mask + 400 corrections"};
  double tot[3] = {0,0,0}, best[3] = {1e30,1e30,1e30};
  const int reps = 60;
  for (int r = 0; r < 5; ++r)
    for (int k = 0; k < 3; ++k) g2b::matvec_range(c, v[k], 0, M, 1.0f/(float)M, x.data(), out.data());
  for (int r = 0; r < reps; ++r) {
    for (int k = 0; k < 3; ++k) {
      auto t0 = std::chrono::steady_clock::now();
      g2b::matvec_range(c, v[k], 0, M, 1.0f/(float)M, x.data(), out.data());
      auto t1 = std::chrono::steady_clock::now();
      const double ms = std::chrono::duration<double, std::milli>(t1-t0).count();
      tot[k] += ms;
      best[k] = std::min(best[k], ms);
    }
  }
  for (int k = 0; k < 3; ++k)
    std::printf("  %-28s mean %7.3f ms   min %7.3f ms\n", nm[k], tot[k]/reps, best[k]);
  std::printf("  mask overhead            mean %+.2f%%   min %+.2f%%\n",
              100.0*(tot[1]-tot[0])/tot[0], 100.0*(best[1]-best[0])/best[0]);
  std::printf("  mask+correction overhead mean %+.2f%%   min %+.2f%%\n",
              100.0*(tot[2]-tot[0])/tot[0], 100.0*(best[2]-best[0])/best[0]);

  // How the correction lists scale. The design says "a few hundred" cells, but
  // that is a guess about real missingness — the pass-2 CSR walk is serial
  // within a sample, so a list concentrated on few rows is the shape that
  // would hurt. Random cells are the benign case; report it as the floor.
  std::printf("  -- correction count sweep (mask fixed at 20%%, cells random) --\n");
  for (int nc : {2000, 20000, 200000}) {
    Trait tv = make_trait(packed, stride, N, M, freq, invstd, 0.20, nc, 0.0,
                          0x7720ull + nc);
    g2b::TraitBind* tb2 = g2b::bind_trait(
        c, tv.freq.data(), tv.invstd.data(),
        tv.mask.data(), (int)tv.mask.size(),
        tv.crow.data(), tv.ccol.data(), tv.cdel.data(), (int)tv.cdel.size());
    if (!tb2) { std::fprintf(stderr, "  sweep bind failed at %d\n", nc); ++g_fail; continue; }
    double bmin = 1e30, bref = 1e30;
    for (int r = 0; r < 5; ++r) {
      g2b::matvec_range(c, plain, 0, M, 1.0f/(float)M, x.data(), out.data());
      g2b::matvec_range(c, tb2,   0, M, tv.inv_M, x.data(), out.data());
    }
    for (int r = 0; r < 25; ++r) {
      auto a0 = std::chrono::steady_clock::now();
      g2b::matvec_range(c, plain, 0, M, 1.0f/(float)M, x.data(), out.data());
      auto a1 = std::chrono::steady_clock::now();
      g2b::matvec_range(c, tb2, 0, M, tv.inv_M, x.data(), out.data());
      auto a2 = std::chrono::steady_clock::now();
      bref = std::min(bref, std::chrono::duration<double,std::milli>(a1-a0).count());
      bmin = std::min(bmin, std::chrono::duration<double,std::milli>(a2-a1).count());
    }
    std::printf("  %7zu corrections   min %7.3f ms  vs plain %7.3f ms  %+.2f%%  (bind %zu MB)\n",
                tv.cdel.size(), bmin, bref, 100.0*(bmin-bref)/bref,
                g2b::bind_bytes(N, M, (int)tv.mask.size(), (int)tv.cdel.size()) >> 20);
    g2b::unbind_trait(tb2);
  }

  g2b::unbind_trait(plain); g2b::unbind_trait(full); g2b::unbind_trait(maskonly);
  g2b::destroy(c);
}

// ====================================== S: several binds on one matrix
// The whole point of the split is that switching phenotypes does not touch the
// matrix. Bind two phenotypes to one Ctx, interleave their matvecs, and require
// each to return exactly what it returned before the other existed.
static void coexist() {
  std::printf("\n================ S: two binds on one Ctx ===========================\n");
  const int N = 20000, M = 3000;
  std::vector<unsigned char> packed; std::size_t stride = 0;
  std::vector<float> freq, invstd;
  gen_matrix(N, M, 0x51A6E2ull, packed, stride, freq, invstd);

  Trait a = make_trait(packed, stride, N, M, freq, invstd, 0.20, 400, 0.05, 0x5A11ull);
  Trait bb = make_trait(packed, stride, N, M, freq, invstd, 0.07, 150, 0.02, 0x5A22ull);

  g2b::Ctx* c = g2b::create(packed.data(), stride, N, M);
  g2b::TraitBind* ta = g2b::bind_trait(c, a.freq.data(), a.invstd.data(),
                                       a.mask.data(), (int)a.mask.size(),
                                       a.crow.data(), a.ccol.data(), a.cdel.data(),
                                       (int)a.cdel.size());
  std::vector<float> x; gen_vec(N, 0x5A33ull, x);
  std::vector<float> ra0(N, 0.f), ra1(N, 0.f), rb(N, 0.f);
  if (!c || !ta || !g2b::matvec_range(c, ta, 0, M, a.inv_M, x.data(), ra0.data())) {
    std::fprintf(stderr, "S setup failed\n"); ++g_fail; return;
  }

  // Second bind arrives AFTER the first has already run.
  g2b::TraitBind* tb = g2b::bind_trait(c, bb.freq.data(), bb.invstd.data(),
                                       bb.mask.data(), (int)bb.mask.size(),
                                       bb.crow.data(), bb.ccol.data(), bb.cdel.data(),
                                       (int)bb.cdel.size());
  if (!tb) { std::fprintf(stderr, "S second bind failed\n"); ++g_fail; return; }
  g2b::matvec_range(c, tb, 0, M, bb.inv_M, x.data(), rb.data());
  g2b::matvec_range(c, ta, 0, M, a.inv_M,  x.data(), ra1.data());

  const int bd = bitdiff(ra0, ra1);
  std::printf("  phenotype A before/after B ran: %d/%d words differ  %s\n",
              bd, N, bd ? "*** FAIL ***" : "bit-identical");
  if (bd) ++g_fail;

  int same = 0;
  for (int i = 0; i < N; ++i) if (std::memcmp(&ra1[i], &rb[i], sizeof(float)) == 0) ++same;
  std::printf("  A vs B outputs identical on %d/%d rows  %s\n", same, N,
              same < N/2 ? "PASS (binds really differ)"
                         : "*** FAIL: the bind is being ignored ***");
  if (same >= N/2) ++g_fail;

  g2b::unbind_trait(ta);
  // B must still work after A is gone — unbind frees only its own buffers.
  std::vector<float> rb2(N, 0.f);
  g2b::matvec_range(c, tb, 0, M, bb.inv_M, x.data(), rb2.data());
  const int bd2 = bitdiff(rb, rb2);
  std::printf("  phenotype B after A unbound:   %d/%d words differ  %s\n",
              bd2, N, bd2 ? "*** FAIL ***" : "bit-identical");
  if (bd2) ++g_fail;

  g2b::unbind_trait(tb);
  g2b::destroy(c);
}

int main(int argc, char** argv) {
  std::string refdir = kDefaultRef;
  bool do_timing = true;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "--ref") && i+1 < argc) refdir = argv[++i];
    else if (!std::strcmp(argv[i], "--no-timing")) do_timing = false;
  }
  b1(refdir);
  b2();
  b3();
  coexist();
  if (do_timing) timing();
  std::printf("\n%s (%d failing check(s))\n", g_fail ? "FAILED" : "ALL CHECKS PASSED", g_fail);
  return g_fail ? 1 : 0;
}
