// spa_gpu_test -- a device SPA against the shipped CPU SPA, pair by pair, on
// identical inputs. Simulated data (no individual-level data goes near this
// box): N samples, p covariates, T binary traits with chosen case fractions,
// M markers spanning MAC 1 to MAF 0.5 with missing calls, allele flips and the
// MAC-gated .clean() zeroing, all expressed through the PLINK 2-bit codes +
// 4-entry dosage table the step-2 reducer stages. For every (marker, trait)
// pair a score statistic z is drawn so that the full-N and the carriers-only
// variants, the linear and the log domain, out-of-support roots and the K2
// guard are all exercised.
//
//   CPU  getadjGFast's arithmetic -> g~, m1, the variant rule, NAmu / NAsigma,
//        q / qinv, then spa.cpp's SPA / SPA_fast (and the inner getroot /
//        Get_Saddle_Prob calls replayed to expose roots, iteration counts and
//        the saddle flags), then getMarkerPval's host post-rules (quantile
//        step, p == 0, the printed string)
//   GPU  Tstat / var1 / var2 / pno handed to the device SPA, then the same
//        post-rules. Two implementations can be driven:
//          --impl mine    gpu/spa_gpu (this library), fast = -1 so the device
//                         applies the variant rule itself
//          --impl integ   the integrator's gpu/gpu_spa.cu, compiled against
//                         integ_shim/ (Makefile target spa_gpu_test_integ);
//                         it takes the variant from the caller
//
// Gate: variant, convergence, saddle flags, out-of-support flags and
// iteration counts must be identical for every pair; p-values are reported
// as max relative difference and max |delta log10 p|; printed strings must
// agree. Non-zero exit on any routing mismatch.
//
//   spa_gpu_test [--n N] [--p P] [--traits r1,r2,...] [--markers M] [--pairs K]
//                [--seed S] [--threads T] [--blocks B] [--erfc-mode 0|1]
//                [--device D] [--tsv FILE] [--bench] [--erfc-sweep] [--no-gate]
//                [--impl mine|integ] [--zmix a,b,c,d]
//                (zmix: weights of the |z| bands [2,5) [5,15) [15,37) [37,50))
#include <armadillo>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "spa.hpp"
#include "spa_binary.hpp"
#include "UTIL.hpp"
#include "score_format.hpp"
#include "spa_gpu.hpp"
#ifdef WITH_INTEGRATOR
#include "gpu_step2.hpp"   // integ_shim/
#include "gpu_spa.hpp"     // the integrator's header (snapshot dir)
#endif

using namespace saige::spa_gpu;

namespace {

struct Opts {
    int N = 50000, p = 3, M = 1500, maxPairs = 30000, seed = 7, threads = 4, blocks = 256;
    int erfcMode = 1, device = 0;
    std::vector<double> rates = {0.5, 0.10, 0.05, 0.01};
    std::vector<double> zmix = {0.45, 0.25, 0.15, 0.15};
    std::string tsv, impl = "mine";
    bool bench = false, erfcSweep = false, gate = true;
};

std::vector<double> parseList(const char* s)
{
    std::vector<double> v; std::string cur;
    for (const char* c = s;; ++c) {
        if (*c == ',' || *c == 0) { if (!cur.empty()) v.push_back(std::atof(cur.c_str())); cur.clear(); if (!*c) break; }
        else cur += *c;
    }
    return v;
}

// ---------------------------------------------------------------------------
// Simulated null models and markers
// ---------------------------------------------------------------------------
struct Trait {
    double rate;
    arma::vec mu;          // N
    arma::mat XV;          // p x N
    arma::mat XXVX_inv;    // N x p
};

struct Marker {
    std::vector<unsigned char> packed;   // bpv bytes
    double fd[4];                        // code -> dosage
    arma::vec G;                         // dense, from fd
    arma::uvec nz, zero;
    double mac;                          // post-impute minor allele count
    bool flip, cleaned; int nMissing;
};

// PLINK: 00 hom A1 (alt here) = 2, 01 missing, 10 het = 1, 11 hom A2 = 0.
constexpr int DMAP[4] = {2, -1, 1, 0};
constexpr unsigned CODE_OF_DOSAGE[3] = {3u, 2u, 0u};
constexpr unsigned CODE_MISSING = 1u;

Marker makeMarker(int N, std::mt19937_64& rng, int kind, double maf, int mac, double missRate, bool flipAlt, bool clean)
{
    Marker mk;
    std::vector<unsigned> code(N, 3u);
    std::uniform_real_distribution<double> U(0, 1);
    if (kind == 0) {
        // exactly `mac` alt alleles: hets, with a hom every third carrier when mac >= 4
        std::vector<int> idx(N); for (int i = 0; i < N; ++i) idx[i] = i;
        std::shuffle(idx.begin(), idx.end(), rng);
        int left = mac, k = 0;
        while (left > 0) {
            if (left >= 4 && (k % 3 == 2)) { code[idx[k]] = CODE_OF_DOSAGE[2]; left -= 2; }
            else { code[idx[k]] = CODE_OF_DOSAGE[1]; left -= 1; }
            ++k;
        }
    } else {
        const double f = flipAlt ? 1.0 - maf : maf;
        std::binomial_distribution<int> B(2, f);
        for (int i = 0; i < N; ++i) code[i] = CODE_OF_DOSAGE[B(rng)];
    }
    mk.nMissing = 0;
    if (missRate > 0) for (int i = 0; i < N; ++i) if (U(rng) < missRate) { code[i] = CODE_MISSING; ++mk.nMissing; }
    // finalizeFusedStats semantics
    unsigned long long counts[4] = {0, 0, 0, 0};
    for (int i = 0; i < N; ++i) counts[code[i]]++;
    double altCount = 2.0 * counts[0] + 1.0 * counts[2];
    const int nCalled = N - mk.nMissing;
    double altFreq = nCalled > 0 ? altCount / (2.0 * nCalled) : 0.0;
    mk.flip = (altFreq > 0.5);
    double af = mk.flip ? 1 - altFreq : altFreq;
    double imputeG = 0, MAC = mk.flip ? (2.0 * nCalled - altCount) : altCount;
    if (mk.nMissing > 0) { imputeG = 2 * af; MAC = MAC + imputeG * mk.nMissing; }
    mk.mac = MAC;
    const double cutoff = clean ? 0.2 : 0.0, macCutoff = 10.0;
    const bool doClean = (cutoff > 0) && (MAC <= macCutoff);
    mk.cleaned = doClean;
    for (int c = 0; c < 4; ++c) {
        double d;
        if ((unsigned)c == CODE_MISSING) d = imputeG;
        else { const double d0 = DMAP[c]; d = mk.flip ? (2 - d0) : d0; }
        if (doClean && std::abs(d) <= cutoff) d = 0;
        mk.fd[c] = d;
    }
    const std::size_t bpv = (std::size_t)(N + 3) / 4;
    mk.packed.assign(bpv, 0);
    for (int i = 0; i < N; ++i) mk.packed[i >> 2] |= (unsigned char)(code[i] << (2 * (i & 3)));
    mk.G.set_size(N);
    std::vector<arma::uword> nz, zr;
    for (int i = 0; i < N; ++i) {
        const double d = mk.fd[code[i]];
        mk.G[i] = d;
        if (d != 0.0) nz.push_back(i); else zr.push_back(i);
    }
    mk.nz = arma::uvec(nz); mk.zero = arma::uvec(zr);
    return mk;
}

// ---------------------------------------------------------------------------
// The CPU side of one pair, as getMarkerPval does it
// ---------------------------------------------------------------------------
struct CpuRes {
    bool fast, logp;
    double Tstat, var1, var2, pno;
    double m1, q, qinv;
    RootResult r1, r2;
    SaddleResult s1, s2;   // valid when r1 & r2 converged
    double pval; bool conv;        // from SPA / SPA_fast
    bool convF; std::string str;   // after the host post-rules
    double z; int nnz;
};

// getMarkerPval's post-SPA rules: the quantile step (which can withdraw
// convergence), the p == 0 rule, the printed string.
bool postRules(double p, bool islog, bool conv, std::string& str, const std::string& noSpaStr)
{
    if (conv) {
        boost::math::normal ns;
        try {
            if (islog) {
                const double half = std::exp(p / 2.0);
                if (half > 0 && half < 1) (void)boost::math::quantile(complement(ns, half));
                else conv = false;
            } else {
                (void)boost::math::quantile(complement(ns, p / 2));
            }
        } catch (const std::overflow_error&) { conv = false; }
    }
    if (!islog && p == 0) conv = false;
    if (conv) {
        char buf[100];
        if (!islog) std::snprintf(buf, sizeof(buf), "%.6E", p);
        else {
            const double l10 = p / std::log(10.0);
            int ex = (int)std::floor(l10);
            double fr = std::pow(10.0, l10 - ex);
            if (fr >= 9.95) { fr = 1; ex++; }
            std::snprintf(buf, sizeof(buf), "%.1fE%d", fr, ex);
        }
        str = buf;
    } else str = noSpaStr;
    return conv;
}

void cpuPair(const Trait& T, const Marker& M, double z, double vr, CpuRes& r, std::string& noSpaStr)
{
    const int N = (int)T.mu.n_elem;
    // getadjGFast
    arma::vec XVG(T.XV.n_rows, arma::fill::zeros);
    for (arma::uword i = 0; i < M.nz.n_elem; ++i) XVG += T.XV.col(M.nz(i)) * M.G(M.nz(i));
    arma::vec gt = M.G - T.XXVX_inv * XVG;
    const double m1 = arma::dot(T.mu, gt);
    arma::vec mu2 = T.mu % (1 - T.mu);
    const double var2 = arma::dot(mu2, arma::square(gt));
    const double var1 = var2 * vr;
    const double S = z * std::sqrt(var1);
    double Beta, seBeta, pval, Tstat, v1, v2; bool islogp;
    SAIGE::format_score_result(S, var1, var2, Beta, seBeta, noSpaStr, pval, islogp, Tstat, v1, v2);
    r.Tstat = Tstat; r.var1 = v1; r.var2 = v2; r.pno = pval; r.logp = islogp; r.z = z;
    r.nnz = (int)M.nz.n_elem; r.m1 = m1;
    const double pz = double(M.zero.n_elem) / N;
    r.fast = (pz >= 0.5);
    arma::vec gNB, gNA, muNB, muNA; double NAmu = 0, NAsigma = 0;
    if (r.fast) {
        gNB = gt(M.nz); gNA = gt(M.zero); muNB = T.mu(M.nz); muNA = T.mu(M.zero);
        NAmu = m1 - arma::dot(gNB, muNB);
        NAsigma = v2 - arma::sum(muNB % (1 - muNB) % arma::pow(gNB, 2));
    }
    const double tol1 = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
    double q = Tstat / std::sqrt(v1 / v2) + m1, qinv;
    if ((q - m1) > 0) qinv = -1 * std::abs(q - m1) + m1;
    else if ((q - m1) == 0) qinv = m1;
    else qinv = std::abs(q - m1) + m1;
    r.q = q; r.qinv = qinv;
    arma::vec mu = T.mu;   // the shipped functions take non-const refs
    if (r.fast) {
        r.r1 = getroot_K1_fast_Binom(0, mu, gt, q,    gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1);
        r.r2 = getroot_K1_fast_Binom(0, mu, gt, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1);
        if (r.r1.Isconverge && r.r2.Isconverge) {
            r.s1 = Get_Saddle_Prob_fast_Binom(r.r1.root, mu, gt, q,    gNA, gNB, muNA, muNB, NAmu, NAsigma, islogp);
            r.s2 = Get_Saddle_Prob_fast_Binom(r.r2.root, mu, gt, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, islogp);
        }
        SPA_fast(mu, gt, q, qinv, pval, islogp, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, "binary", r.pval, r.conv);
    } else {
        r.r1 = getroot_K1_Binom(0, mu, gt, q, tol1);
        r.r2 = getroot_K1_Binom(0, mu, gt, qinv, tol1);
        if (r.r1.Isconverge && r.r2.Isconverge) {
            r.s1 = Get_Saddle_Prob_Binom(r.r1.root, mu, gt, q, islogp);
            r.s2 = Get_Saddle_Prob_Binom(r.r2.root, mu, gt, qinv, islogp);
        }
        SPA(mu, gt, q, qinv, pval, tol1, islogp, "binary", r.pval, r.conv);
    }
    r.convF = postRules(r.pval, islogp, r.conv, r.str, noSpaStr);
}

// ---------------------------------------------------------------------------
// Device implementations behind one interface
// ---------------------------------------------------------------------------
struct GRes {
    double pval; int conv; bool fast; int s1, s2; int niter1, niter2; double root1, root2; double m1; int nnz;
    unsigned status; int reason1, reason2; bool hasBits;   // status / reasons only from spa_gpu
};

struct Impl {
    virtual ~Impl() {}
    virtual const char* name() const = 0;
    // Solve the pairs idx[] (fastMode -1: device decides when supported, else the CPU's choice).
    virtual bool solve(const std::vector<int>& idx, int fastMode, std::vector<GRes>& res,
                       double* kernelSec, double* h2dSec, double* d2hSec) = 0;
    virtual bool decidesVariant() const = 0;
};

struct PairSrc {
    const std::vector<CpuRes>* cres; const std::vector<int>* pm; const std::vector<int>* pt;
};

struct MineImpl : Impl {
    Spa* S = nullptr; Geno geno; PairSrc src; int maxPairs = 0;
    const char* name() const override { return "spa_gpu (this library)"; }
    bool decidesVariant() const override { return true; }
    bool solve(const std::vector<int>& idx, int fastMode, std::vector<GRes>& res, double* ks, double* hs, double* ds) override
    {
        res.resize(idx.size());
        double k0, h0, d0, k1, h1, d1; long long n0, n1;
        timings(S, &k0, &h0, &d0, &n0);
        for (std::size_t s0 = 0; s0 < idx.size(); s0 += maxPairs) {
            const int kc = (int)std::min<std::size_t>(maxPairs, idx.size() - s0);
            PairIn* pi = in(S);
            for (int k = 0; k < kc; ++k) {
                const int j = idx[s0 + k]; const CpuRes& c = (*src.cres)[j];
                pi[k].slot = (*src.pm)[j]; pi[k].trait = (*src.pt)[j];
                pi[k].fast = (fastMode < 0) ? -1 : (c.fast ? 1 : 0);
                pi[k].logp = c.logp ? 1 : 0;
                pi[k].Tstat = c.Tstat; pi[k].var1 = c.var1; pi[k].var2 = c.var2; pi[k].pno = c.pno;
            }
            if (!run(S, geno, kc)) { std::fprintf(stderr, "spa_gpu::run failed: %s\n", lastError()); return false; }
            const PairOut* po = out(S);
            for (int k = 0; k < kc; ++k) {
                GRes& g = res[s0 + k];
                g.pval = po[k].pval; g.conv = po[k].conv; g.fast = (po[k].status & ST_FAST) != 0;
                g.s1 = po[k].s1; g.s2 = po[k].s2; g.niter1 = po[k].niter1; g.niter2 = po[k].niter2;
                g.root1 = po[k].root1; g.root2 = po[k].root2; g.m1 = po[k].m1; g.nnz = po[k].nnz;
                g.status = po[k].status; g.reason1 = po[k].reason1; g.reason2 = po[k].reason2; g.hasBits = true;
            }
        }
        timings(S, &k1, &h1, &d1, &n1);
        if (ks) *ks = k1 - k0; if (hs) *hs = h1 - h0; if (ds) *ds = d1 - d0;
        return true;
    }
};

#ifdef WITH_INTEGRATOR
struct IntegImpl : Impl {
    saige::gpu2::Spa* S = nullptr; saige::gpu2::Reducer R; PairSrc src; int maxPairs = 0;
    const char* name() const override { return "integrator gpu/gpu_spa.cu (snapshot)"; }
    bool decidesVariant() const override { return false; }
    bool solve(const std::vector<int>& idx, int, std::vector<GRes>& res, double* ks, double* hs, double* ds) override
    {
        res.resize(idx.size());
        double k0, k1; long long n0, n1;
        saige::gpu2::spaTimings(S, &k0, &n0);
        auto w0 = std::chrono::steady_clock::now();
        for (std::size_t s0 = 0; s0 < idx.size(); s0 += maxPairs) {
            const int kc = (int)std::min<std::size_t>(maxPairs, idx.size() - s0);
            saige::gpu2::SpaPairIn* pi = saige::gpu2::spaIn(S);
            for (int k = 0; k < kc; ++k) {
                const int j = idx[s0 + k]; const CpuRes& c = (*src.cres)[j];
                pi[k].slot = (*src.pm)[j]; pi[k].trait = (*src.pt)[j];
                pi[k].fast = c.fast ? 1 : 0; pi[k].logp = c.logp ? 1 : 0;
                pi[k].Tstat = c.Tstat; pi[k].var1 = c.var1; pi[k].var2 = c.var2; pi[k].pno = c.pno;
            }
            if (!saige::gpu2::spaRun(S, &R, kc)) { std::fprintf(stderr, "gpu2::spaRun failed\n"); return false; }
            const saige::gpu2::SpaPairOut* po = saige::gpu2::spaOut(S);
            for (int k = 0; k < kc; ++k) {
                GRes& g = res[s0 + k]; const CpuRes& c = (*src.cres)[idx[s0 + k]];
                g.pval = po[k].pval; g.conv = po[k].conv; g.fast = c.fast;
                g.s1 = po[k].s1; g.s2 = po[k].s2; g.niter1 = po[k].niter1; g.niter2 = po[k].niter2;
                g.root1 = po[k].root1; g.root2 = po[k].root2; g.m1 = po[k].m1; g.nnz = -1;
                g.status = 0; g.reason1 = g.reason2 = 0; g.hasBits = false;
            }
        }
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count();
        saige::gpu2::spaTimings(S, &k1, &n1);
        if (ks) *ks = k1 - k0; if (hs) *hs = wall - (k1 - k0); if (ds) *ds = 0;   // its timer covers the kernel only
        return true;
    }
};
#endif

// ---------------------------------------------------------------------------
// Comparison bookkeeping
// ---------------------------------------------------------------------------
struct Tally {
    long n = 0, mmVariant = 0, mmConv = 0, mmSaddle = 0, mmInf = 0, mmIter = 0, mmConvF = 0, mmStr = 0, mmZero = 0;
    double maxRel = 0, maxDlog10 = 0, maxDroot = 0;
    long nConv = 0, nInf = 0, nSaddleFail = 0, nRootFail = 0;
};

struct Cmp {
    bool variantOk, convOk, saddleOk, infOk, iterOk, convFOk, strOk, zeroOk;
    double rel, dlog10, droot;
};

Cmp compare(const CpuRes& c, const GRes& g, bool gConvF, const std::string& gStr)
{
    Cmp o;
    o.variantOk = (g.fast == c.fast);
    o.convOk = ((g.conv != 0) == c.conv);
    const bool bothRoots = c.r1.Isconverge && c.r2.Isconverge;
    if (bothRoots) o.saddleOk = (g.s1 == (c.s1.isSaddle ? 1 : 0)) && (g.s2 == (c.s2.isSaddle ? 1 : 0));
    else o.saddleOk = (g.s1 == -1 && g.s2 == -1);
    const bool cInf1 = std::isinf(c.r1.root) && c.r1.niter == 0, cInf2 = std::isinf(c.r2.root) && c.r2.niter == 0;
    const bool gInf1 = g.hasBits ? ((g.status & ST_ROOT1_INF) != 0) : (std::isinf(g.root1) && g.niter1 == 0);
    const bool gInf2 = g.hasBits ? ((g.status & ST_ROOT2_INF) != 0) : (std::isinf(g.root2) && g.niter2 == 0);
    o.infOk = (gInf1 == cInf1) && (gInf2 == cInf2);
    if (g.hasBits) o.infOk = o.infOk && (((g.status & ST_ROOT1_FAIL) != 0) == !c.r1.Isconverge) && (((g.status & ST_ROOT2_FAIL) != 0) == !c.r2.Isconverge);
    o.iterOk = (g.niter1 == c.r1.niter) && (g.niter2 == c.r2.niter);
    o.convFOk = (gConvF == c.convF);
    o.strOk = (gStr == c.str);
    o.droot = 0;
    if (std::isfinite(c.r1.root) && std::isfinite(g.root1)) o.droot = std::max(o.droot, std::fabs(c.r1.root - g.root1));
    if (std::isfinite(c.r2.root) && std::isfinite(g.root2)) o.droot = std::max(o.droot, std::fabs(c.r2.root - g.root2));
    // p-value: in the log domain pval is log p, so |dp| / ln 10 is |delta log10 p|
    const double pc = c.pval, pg = g.pval;
    o.zeroOk = true; o.rel = 0; o.dlog10 = 0;
    if (c.logp) {
        if (std::isinf(pc) || std::isinf(pg) || std::isnan(pc) || std::isnan(pg)) {
            o.zeroOk = (std::isinf(pc) == std::isinf(pg)) && (std::isnan(pc) == std::isnan(pg)) && (!std::isinf(pc) || (pc < 0) == (pg < 0));
        } else {
            o.dlog10 = std::fabs(pc - pg) / std::log(10.0);
            o.rel = std::fabs(std::expm1(pg - pc));
        }
    } else {
        if (pc == 0 || pg == 0) o.zeroOk = (pc == 0 && pg == 0);
        else {
            o.rel = std::fabs(pc - pg) / std::fabs(pc);
            o.dlog10 = std::fabs(std::log10(pc) - std::log10(pg));
        }
    }
    return o;
}

void tally(Tally& t, const CpuRes& c, const Cmp& o)
{
    t.n++;
    if (!o.variantOk) t.mmVariant++;
    if (!o.convOk) t.mmConv++;
    if (!o.saddleOk) t.mmSaddle++;
    if (!o.infOk) t.mmInf++;
    if (!o.iterOk) t.mmIter++;
    if (!o.convFOk) t.mmConvF++;
    if (!o.strOk) t.mmStr++;
    if (!o.zeroOk) t.mmZero++;
    t.maxRel = std::max(t.maxRel, o.rel);
    t.maxDlog10 = std::max(t.maxDlog10, o.dlog10);
    t.maxDroot = std::max(t.maxDroot, o.droot);
    if (c.conv) t.nConv++;
    if ((std::isinf(c.r1.root) && c.r1.niter == 0) || (std::isinf(c.r2.root) && c.r2.niter == 0)) t.nInf++;
    if (!c.r1.Isconverge || !c.r2.Isconverge) t.nRootFail++;
    else if (!c.s1.isSaddle || !c.s2.isSaddle) t.nSaddleFail++;
}

const char* zBand(double z)
{
    z = std::fabs(z);
    if (z < 5) return "|z| [2,5)";
    if (z < 15) return "|z| [5,15)";
    if (z < 37) return "|z| [15,37)";
    return "|z| >= 37";
}
const char* macBand(double mac)
{
    if (mac <= 10) return "MAC <= 10";
    if (mac <= 100) return "MAC (10,100]";
    return "MAC > 100";
}

void printTally(const std::string& name, const Tally& t)
{
    std::printf("  %-34s n=%7ld  conv=%7ld inf=%5ld rootFail=%4ld saddleFail=%4ld | mismatch variant=%ld conv=%ld saddle=%ld inf/fail=%ld iter=%ld convF=%ld str=%ld zero=%ld | maxRel=%.2e maxDlog10=%.2e maxDroot=%.2e\n",
        name.c_str(), t.n, t.nConv, t.nInf, t.nRootFail, t.nSaddleFail,
        t.mmVariant, t.mmConv, t.mmSaddle, t.mmInf, t.mmIter, t.mmConvF, t.mmStr, t.mmZero,
        t.maxRel, t.maxDlog10, t.maxDroot);
}

// ---------------------------------------------------------------------------
// The tail study: Boost's erfc on the host, the ported one and CUDA's on the
// device, over a grid of arguments, as the normal upper tail p = erfc(x)/2.
// ---------------------------------------------------------------------------
int erfcSweep(const Opts& o)
{
    std::vector<double> z;
    for (double x = -3.0; x < 26.5; x += 1e-3) z.push_back(x);
    for (double x = 26.5; x < 28.2; x += 1e-6) z.push_back(x);          // the subnormal zone, finely
    std::mt19937_64 rng(3); std::uniform_real_distribution<double> U(-3, 28);
    for (int i = 0; i < 300000; ++i) z.push_back(U(rng));
    z.push_back(28.0); z.push_back(27.999999); z.push_back(27.2); z.push_back(0.0); z.push_back(0.5); z.push_back(1.5); z.push_back(2.5); z.push_back(4.5);
    z.push_back(std::numeric_limits<double>::infinity()); z.push_back(-std::numeric_limits<double>::infinity());
    const int n = (int)z.size();
    std::vector<double> dev[3], host(n);
    for (auto& v : dev) v.resize(n);
    if (!debugErfc(o.device, z.data(), n, dev[0].data(), dev[1].data(), dev[2].data())) { std::fprintf(stderr, "debugErfc: %s\n", lastError()); return 2; }
    for (int i = 0; i < n; ++i) host[i] = std::isnan(z[i]) ? NAN : boost::math::erfc(z[i]);
    const char* names[3] = {"scaled port (mode 1, default)", "literal port (mode 2)", "CUDA libm erfc (mode 0)"};
    auto ulpDiff = [](double a, double b) -> double {
        if (a == b) return 0;
        if (!std::isfinite(a) || !std::isfinite(b) || a == 0 || b == 0) return std::numeric_limits<double>::infinity();
        return std::fabs(a - b) / std::fabs(std::nexttoward(a, b) - a);
    };
    const double dmin = std::numeric_limits<double>::denorm_min(), nmin = std::numeric_limits<double>::min();
    struct Band { const char* name; double lo, hi; long n = 0, exact[3] = {0, 0, 0}, zeroMis[3] = {0, 0, 0}; double maxUlp[3] = {0, 0, 0}, maxGrid[3] = {0, 0, 0}; };
    std::vector<Band> bands = { {"z < 0.5", -4, 0.5}, {"[0.5,1.5)", 0.5, 1.5}, {"[1.5,2.5)", 1.5, 2.5}, {"[2.5,4.5)", 2.5, 4.5},
                                {"[4.5,26.55) normal p", 4.5, 26.55}, {"[26.55,27.3) subnormal p", 26.55, 27.3}, {"[27.3,28)", 27.3, 28.0}, {"z >= 28", 28.0, 1e300} };
    double zZero[4] = {NAN, NAN, NAN, NAN};   // host, modes
    double logMaxDiff[3] = {0, 0, 0};          // |log(dev) - log(host)| where both non-zero (the log-domain p)
    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(z[i])) continue;
        const double h = host[i];
        if (std::isnan(zZero[0]) && h == 0 && z[i] > 0) zZero[0] = z[i];
        for (int m = 0; m < 3; ++m) {
            const double d = dev[m][i];
            if (std::isnan(zZero[m + 1]) && d == 0 && z[i] > 0) zZero[m + 1] = z[i];
            if (h != 0 && d != 0) logMaxDiff[m] = std::max(logMaxDiff[m], std::fabs(std::log(d) - std::log(h)));
        }
        for (Band& b : bands) if (z[i] >= b.lo && z[i] < b.hi) {
            b.n++;
            for (int m = 0; m < 3; ++m) {
                const double d = dev[m][i];
                if (d == h) b.exact[m]++;
                if ((h == 0) != (d == 0)) b.zeroMis[m]++;
                const bool sub = (h != 0 && std::fabs(h) < nmin) || (d != 0 && std::fabs(d) < nmin);
                if (sub) b.maxGrid[m] = std::max(b.maxGrid[m], std::fabs(d - h) / dmin);
                else if (h != 0 && d != 0) { const double u = ulpDiff(d, h); if (std::isfinite(u)) b.maxUlp[m] = std::max(b.maxUlp[m], u); }
            }
        }
    }
    std::printf("erfc tail study on the device: %d arguments; host = boost::math::erfc (Boost 1.85, evaluated in long double, i.e. correctly rounded, 0 from z = %.7f)\n", n, zZero[0]);
    std::printf("  first positive z with erfc == 0: host %.7f | %s %.7f | %s %.7f | %s %.7f\n", zZero[0], names[0], zZero[1], names[1], zZero[2], names[2], zZero[3]);
    for (int m = 0; m < 3; ++m) std::printf("  max |log(dev) - log(host)| over arguments where both are non-zero: %s %.3e\n", names[m], logMaxDiff[m]);
    std::printf("  per band: n, then for each mode [bit-identical, max ulp (normal), max subnormal grid units, zero/non-zero disagreements]\n");
    std::printf("  %-26s %8s | %-32s | %-32s | %-32s\n", "band", "n", names[0], names[1], names[2]);
    for (const Band& b : bands) {
        std::printf("  %-26s %8ld |", b.name, b.n);
        for (int m = 0; m < 3; ++m) std::printf(" %7ld %6.2f %8.0f %7ld |", b.exact[m], b.maxUlp[m], b.maxGrid[m], b.zeroMis[m]);
        std::printf("\n");
    }
    std::printf("  p-value terms: p = erfc(z)/2, Z = z*sqrt2: z = 26.55 -> p %.3e (Z %.2f); z = 27.226 -> host 0 (Z %.3f); z = 28 (Z %.3f)\n",
                boost::math::erfc(26.55) / 2, 26.55 * std::sqrt(2.0), 27.226 * std::sqrt(2.0), 28 * std::sqrt(2.0));
    for (int i = n - 10; i < n; ++i) std::printf("  z=%-12.6g host=%-24.17g scaled=%-24.17g literal=%-24.17g cuda=%-24.17g\n", z[i], host[i], dev[0][i], dev[1][i], dev[2][i]);
    return 0;
}

}  // namespace

int main(int argc, char** argv)
{
    Opts o;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> const char* { if (i + 1 >= argc) { std::fprintf(stderr, "%s needs a value\n", a.c_str()); std::exit(2); } return argv[++i]; };
        if (a == "--n") o.N = std::atoi(next());
        else if (a == "--p") o.p = std::atoi(next());
        else if (a == "--markers") o.M = std::atoi(next());
        else if (a == "--pairs") o.maxPairs = std::atoi(next());
        else if (a == "--seed") o.seed = std::atoi(next());
        else if (a == "--threads") o.threads = std::atoi(next());
        else if (a == "--blocks") o.blocks = std::atoi(next());
        else if (a == "--erfc-mode") o.erfcMode = std::atoi(next());
        else if (a == "--device") o.device = std::atoi(next());
        else if (a == "--traits") o.rates = parseList(next());
        else if (a == "--zmix") o.zmix = parseList(next());
        else if (a == "--tsv") o.tsv = next();
        else if (a == "--impl") o.impl = next();
        else if (a == "--bench") o.bench = true;
        else if (a == "--erfc-sweep") o.erfcSweep = true;
        else if (a == "--no-gate") o.gate = false;
        else { std::fprintf(stderr, "unknown option %s\n", a.c_str()); return 2; }
    }
    if (o.erfcSweep) return erfcSweep(o);
    if (o.p < 1 || o.p > PMAX) { std::fprintf(stderr, "--p must be 1..%d\n", PMAX); return 2; }
#ifndef WITH_INTEGRATOR
    if (o.impl == "integ") { std::fprintf(stderr, "--impl integ needs the spa_gpu_test_integ build\n"); return 2; }
#endif
    omp_set_num_threads(o.threads);
    const int N = o.N, p = o.p, T = (int)o.rates.size();
    std::mt19937_64 rng(o.seed);
    std::normal_distribution<double> Nrm(0, 1);
    std::uniform_real_distribution<double> U(0, 1);

    // ---- covariates and traits
    arma::mat X(N, p, arma::fill::ones);
    for (int i = 0; i < N; ++i) {
        if (p > 1) X(i, 1) = Nrm(rng);
        if (p > 2) X(i, 2) = (U(rng) < 0.5) ? 1.0 : 0.0;
        for (int j = 3; j < p; ++j) X(i, j) = Nrm(rng);
    }
    std::vector<Trait> traits(T);
    for (int t = 0; t < T; ++t) {
        Trait& tr = traits[t];
        tr.rate = o.rates[t];
        arma::vec a(p, arma::fill::zeros);
        for (int j = 1; j < p; ++j) a[j] = 0.3 * Nrm(rng);
        a[0] = std::log(tr.rate / (1 - tr.rate));
        arma::vec eta = X * a;
        tr.mu = 1.0 / (1.0 + arma::exp(-eta));
        arma::vec W = tr.mu % (1 - tr.mu);
        arma::mat XW = X.each_col() % W;
        tr.XV = XW.t();                                   // p x N
        arma::mat XVX = X.t() * XW;
        tr.XXVX_inv = X * arma::inv_sympd(XVX);           // N x p
    }

    // ---- markers
    const int fixedMac[] = {1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500};
    const int nFixed = (int)(sizeof(fixedMac) / sizeof(int));
    std::vector<Marker> markers; markers.reserve(o.M);
    for (int m = 0; m < o.M; ++m) {
        const bool missing = (m % 10 == 3);
        const bool flip = (m % 10 == 7);
        const bool clean = (m % 10 == 5);
        if (m < 2 * nFixed) {
            markers.push_back(makeMarker(N, rng, 0, 0, fixedMac[m % nFixed], (m >= nFixed && missing) ? 0.02 : 0.0, false, clean));
        } else {
            const double maf = std::exp(std::log(2.0 / N) + U(rng) * (std::log(0.5) - std::log(2.0 / N)));
            markers.push_back(makeMarker(N, rng, 1, maf, 0, missing ? 0.02 : 0.0, flip, clean));
        }
    }
    const std::size_t bpv = (std::size_t)(N + 3) / 4;
    std::vector<unsigned char> packedAll((std::size_t)o.M * bpv);
    std::vector<double> lutAll((std::size_t)o.M * 4);
    for (int m = 0; m < o.M; ++m) {
        std::memcpy(packedAll.data() + (std::size_t)m * bpv, markers[m].packed.data(), bpv);
        for (int c = 0; c < 4; ++c) lutAll[4 * m + c] = markers[m].fd[c];
    }

    // ---- pairs
    std::vector<int> pm, pt; std::vector<double> pz, pvr;
    {
        std::vector<std::pair<int, int>> all;
        for (int m = 0; m < o.M; ++m) for (int t = 0; t < T; ++t) all.push_back({m, t});
        std::shuffle(all.begin(), all.end(), rng);
        if ((int)all.size() > o.maxPairs) all.resize(o.maxPairs);
        // every fixed-MAC marker gets every trait, so low-MAC coverage does not depend on the shuffle
        std::vector<char> have((std::size_t)o.M * T, 0);
        for (auto& pr : all) have[(std::size_t)pr.first * T + pr.second] = 1;
        for (int m = 0; m < std::min(o.M, 2 * nFixed); ++m) for (int t = 0; t < T; ++t) if (!have[(std::size_t)m * T + t]) all.push_back({m, t});
        const double wsum = o.zmix[0] + o.zmix[1] + o.zmix[2] + o.zmix[3];
        for (auto& pr : all) {
            pm.push_back(pr.first); pt.push_back(pr.second);
            const double u = U(rng) * wsum; double zlo, zhi;
            if (u < o.zmix[0]) { zlo = 2; zhi = 5; }
            else if (u < o.zmix[0] + o.zmix[1]) { zlo = 5; zhi = 15; }
            else if (u < o.zmix[0] + o.zmix[1] + o.zmix[2]) { zlo = 15; zhi = 37; }
            else { zlo = 37; zhi = 50; }
            pz.push_back((U(rng) < 0.5 ? -1 : 1) * (zlo + U(rng) * (zhi - zlo)));
            pvr.push_back(0.85 + 0.3 * U(rng));
        }
    }
    const int K = (int)pm.size();
    std::printf("spa_gpu_test: impl=%s N=%d p=%d traits=%d (case rates", o.impl.c_str(), N, p, T);
    for (double r : o.rates) std::printf(" %g", r);
    std::printf(") markers=%d pairs=%d seed=%d threads=%d blocks=%d erfcMode=%d\n", o.M, K, o.seed, o.threads, o.blocks, o.erfcMode);

    // ---- CPU reference
    std::vector<CpuRes> cres(K);
    std::vector<std::string> noSpa(K);
    auto t0 = std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic, 16)
    for (int k = 0; k < K; ++k) cpuPair(traits[pt[k]], markers[pm[k]], pz[k], pvr[k], cres[k], noSpa[k]);
    const double tCpu = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    long nFastCpu = 0, nLogCpu = 0;
    for (auto& c : cres) { nFastCpu += c.fast; nLogCpu += c.logp; }
    std::printf("CPU reference: %.1f s on %d threads (%.3f ms/pair); fast variant %ld, log domain %ld\n", tCpu, o.threads, tCpu / K * 1e3, nFastCpu, nLogCpu);

    // ---- device: this library always provides the genotype upload
    std::vector<TraitArgs> ta(T);
    for (int t = 0; t < T; ++t) { ta[t].mu = traits[t].mu.memptr(); ta[t].XV = traits[t].XV.memptr(); ta[t].XXVX_inv = traits[t].XXVX_inv.memptr(); ta[t].p = p; }
    CreateArgs ca;
    ca.device = o.device; ca.N = N; ca.nTraits = T; ca.traits = ta.data();
    ca.maxPairs = std::min(K, 65536); ca.tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
    ca.maxiter = 1000; ca.blocks = o.blocks; ca.erfcMode = o.erfcMode;
    Spa* S = create(ca);
    if (!S) { std::fprintf(stderr, "spa_gpu::create failed: %s\n", lastError()); return 2; }
    Geno geno;
    if (!uploadPacked(S, packedAll.data(), bpv, lutAll.data(), o.M, &geno)) { std::fprintf(stderr, "uploadPacked: %s\n", lastError()); return 2; }
    std::printf("GPU: spa_gpu device bytes %.1f MB (scratch %d blocks x 16N)\n", deviceBytes(S) / 1048576.0, o.blocks);

    PairSrc src{&cres, &pm, &pt};
    std::unique_ptr<Impl> impl;
    if (o.impl == "mine") {
        auto* m = new MineImpl(); m->S = S; m->geno = geno; m->src = src; m->maxPairs = ca.maxPairs; impl.reset(m);
    }
#ifdef WITH_INTEGRATOR
    else if (o.impl == "integ") {
        auto* m = new IntegImpl();
        std::vector<saige::gpu2::SpaTraitArgs> ita(T);
        for (int t = 0; t < T; ++t) { ita[t].mu = ta[t].mu; ita[t].XV = ta[t].XV; ita[t].XXVX_inv = ta[t].XXVX_inv; ita[t].p = p; }
        saige::gpu2::SpaCreateArgs ia;
        ia.device = o.device; ia.N = N; ia.nTraits = T; ia.traits = ita.data(); ia.maxPairs = ca.maxPairs;
        ia.tol = ca.tol; ia.maxiter = 1000; ia.blocks = o.blocks;
        m->S = saige::gpu2::spaCreate(ia);
        if (!m->S) { std::fprintf(stderr, "gpu2::spaCreate failed\n"); return 2; }
        m->R.pk = geno.packed; m->R.lut = geno.lut; m->R.bpv = geno.bpv;
        m->src = src; m->maxPairs = ca.maxPairs; impl.reset(m);
        std::printf("GPU: integrator kernel created, device bytes %.1f MB\n", saige::gpu2::spaDeviceBytes(m->S) / 1048576.0);
    }
#endif
    else { std::fprintf(stderr, "unknown --impl %s\n", o.impl.c_str()); return 2; }

    std::vector<int> allIdx(K); for (int k = 0; k < K; ++k) allIdx[k] = k;
    std::vector<GRes> gres, gres2;
    double k1, h1, d1, k2, h2, d2;
    if (!impl->solve(allIdx, -1, gres, &k1, &h1, &d1)) return 2;    // device decides the variant (where supported)
    if (!impl->solve(allIdx, 0, gres2, &k2, &h2, &d2)) return 2;    // the caller decides, as the integrator does
    long diffDecide = 0;
    for (int k = 0; k < K; ++k) if (std::memcmp(&gres[k], &gres2[k], sizeof(GRes)) != 0) diffDecide++;
    std::printf("GPU [%s]: %d pairs twice: kernel %.3f + %.3f s (%.2f us/pair), H2D %.3f s, D2H %.3f s; results with the variant %s vs given by the caller: %ld pairs differ\n",
                impl->name(), K, k1, k2, (k1 + k2) / (2.0 * K) * 1e6, h1 + h2, d1 + d2, impl->decidesVariant() ? "decided on the device" : "(no device decision)", diffDecide);

    // ---- compare
    int rc = 0;
    if (o.gate) {
        Tally all;
        std::map<std::string, Tally> by;
        std::vector<int> bad;
        FILE* ft = o.tsv.empty() ? nullptr : std::fopen(o.tsv.c_str(), "w");
        if (ft) std::fprintf(ft, "m\tt\trate\tmac\tz\tfast_c\tfast_g\tlogp\tTstat\tvar1\tvar2\tpno\tnnz_c\tnnz_g\tm1_c\tm1_g\tr1c\tr1g\tr2c\tr2g\tn1c\tn1g\tn2c\tn2g\ts1c\ts1g\ts2c\ts2g\tconvc\tconvg\tstatus\treason1\treason2\tpc\tpg\tconvFc\tconvFg\tstrc\tstrg\n");
        for (int k = 0; k < K; ++k) {
            const CpuRes& c = cres[k]; const GRes& g = gres[k];
            std::string gStr; const bool gConvF = postRules(g.pval, c.logp, g.conv != 0, gStr, noSpa[k]);
            const Cmp cm = compare(c, g, gConvF, gStr);
            tally(all, c, cm);
            tally(by[std::string("variant: ") + (c.fast ? "fast (carriers)" : "full N")], c, cm);
            tally(by[std::string("domain: ") + (c.logp ? "log-p" : "linear")], c, cm);
            char buf[64]; std::snprintf(buf, sizeof(buf), "case rate %g", traits[pt[k]].rate); tally(by[buf], c, cm);
            tally(by[zBand(c.z)], c, cm);
            tally(by[macBand(markers[pm[k]].mac)], c, cm);
            if (markers[pm[k]].nMissing > 0) tally(by["markers with missing calls"], c, cm);
            if (markers[pm[k]].flip) tally(by["markers with allele flip"], c, cm);
            if (markers[pm[k]].cleaned) tally(by["markers with .clean() zeroing"], c, cm);
            const bool routingOk = cm.variantOk && cm.convOk && cm.saddleOk && cm.infOk && cm.iterOk && cm.convFOk && cm.strOk && cm.zeroOk;
            if (!routingOk) bad.push_back(k);
            if (ft) std::fprintf(ft, "%d\t%d\t%g\t%g\t%.6g\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%u\t%d\t%d\t%.17g\t%.17g\t%d\t%d\t%s\t%s\n",
                pm[k], pt[k], traits[pt[k]].rate, markers[pm[k]].mac, c.z, (int)c.fast, (int)g.fast, (int)c.logp,
                c.Tstat, c.var1, c.var2, c.pno, c.nnz, g.nnz, c.m1, g.m1, c.r1.root, g.root1, c.r2.root, g.root2, c.r1.niter, g.niter1, c.r2.niter, g.niter2,
                (c.r1.Isconverge && c.r2.Isconverge) ? (int)c.s1.isSaddle : -1, g.s1, (c.r1.Isconverge && c.r2.Isconverge) ? (int)c.s2.isSaddle : -1, g.s2,
                (int)c.conv, g.conv, g.status, g.reason1, g.reason2, c.pval, g.pval, (int)c.convF, (int)gConvF, c.str.c_str(), gStr.c_str());
        }
        if (ft) std::fclose(ft);
        std::printf("\nACCURACY GATE [%s] (CPU = shipped spa.cpp / spa_binary.cpp; identical inputs)\n", impl->name());
        printTally("all pairs", all);
        for (auto& kv : by) printTally(kv.first, kv.second);
        const long mm = all.mmVariant + all.mmConv + all.mmSaddle + all.mmInf + all.mmIter + all.mmConvF + all.mmStr + all.mmZero;
        std::printf("\nRESULT [%s]: %s -- routing mismatches %ld of %ld pairs; max rel |dp| = %.3e, max |dlog10 p| = %.3e, max |droot| = %.3e; printed strings differ on %ld pairs\n",
                    impl->name(), mm == 0 ? "PASS" : "FAIL", mm, all.n, all.maxRel, all.maxDlog10, all.maxDroot, all.mmStr);
        for (std::size_t i = 0; i < bad.size() && i < 20; ++i) {
            const int k = bad[i]; const CpuRes& c = cres[k]; const GRes& g = gres[k];
            std::string gStr; postRules(g.pval, c.logp, g.conv != 0, gStr, noSpa[k]);
            std::printf("  mismatch pair %d: m=%d t=%d mac=%g z=%.4g fast c/g=%d/%d logp=%d | roots c (%.17g,%d,%d) (%.17g,%d,%d) g (%.17g,%d) (%.17g,%d) status=%u reasons=%d,%d | saddle c=%d,%d g=%d,%d | conv c/g=%d/%d p c/g=%.17g/%.17g str %s / %s\n",
                k, pm[k], pt[k], markers[pm[k]].mac, c.z, (int)c.fast, (int)g.fast, (int)c.logp,
                c.r1.root, c.r1.niter, (int)c.r1.Isconverge, c.r2.root, c.r2.niter, (int)c.r2.Isconverge, g.root1, g.niter1, g.root2, g.niter2, g.status, g.reason1, g.reason2,
                (c.r1.Isconverge && c.r2.Isconverge) ? (int)c.s1.isSaddle : -1, (c.r1.Isconverge && c.r2.Isconverge) ? (int)c.s2.isSaddle : -1, g.s1, g.s2,
                (int)c.conv, g.conv, c.pval, g.pval, c.str.c_str(), gStr.c_str());
        }
        if (mm != 0) rc = 1;
    }

    // ---- throughput
    if (o.bench) {
        std::printf("\nTHROUGHPUT [%s] (device kernel time from events; pair table %zu B/pair in, results %zu B/pair out; best of 3)\n", impl->name(), sizeof(PairIn), sizeof(PairOut));
        auto bench = [&](const char* name, const std::vector<int>& idx) {
            if (idx.empty()) return;
            std::vector<int> sub(idx.begin(), idx.begin() + std::min((int)idx.size(), ca.maxPairs));
            double best = 1e30, bh = 0, bd = 0;
            for (int rep = 0; rep < 3; ++rep) {
                std::vector<GRes> r; double ks, hs, ds;
                if (!impl->solve(sub, 0, r, &ks, &hs, &ds)) std::exit(2);
                if (ks < best) { best = ks; bh = hs; bd = ds; }
            }
            double itSum = 0; for (int k : sub) itSum += cres[k].r1.niter + cres[k].r2.niter;
            const int n = (int)sub.size();
            std::printf("  %-36s pairs=%6d  kernel %.2f us/pair  (%.0f pairs/s)   H2D %.3f us/pair  D2H %.3f us/pair   mean Newton steps per pair %.2f\n",
                        name, n, best / n * 1e6, n / best, bh / n * 1e6, bd / n * 1e6, itSum / n);
        };
        std::vector<int> fullIdx, fastIdx, linIdx, logIdx;
        std::map<int, std::vector<int>> byTrait;
        for (int k = 0; k < K; ++k) {
            (cres[k].fast ? fastIdx : fullIdx).push_back(k);
            (cres[k].logp ? logIdx : linIdx).push_back(k);
            byTrait[pt[k]].push_back(k);
        }
        bench("all pairs (mix)", allIdx);
        bench("full-N variant only", fullIdx);
        bench("fast (carriers) variant only", fastIdx);
        bench("linear domain only", linIdx);
        bench("log domain only", logIdx);
        for (auto& kv : byTrait) { char b[64]; std::snprintf(b, sizeof(b), "case rate %g (both variants)", traits[kv.first].rate); bench(b, kv.second); }
        for (auto& kv : byTrait) {
            std::vector<int> f, g;
            for (int k : kv.second) (cres[k].fast ? f : g).push_back(k);
            char b[64];
            std::snprintf(b, sizeof(b), "case rate %g, full N", traits[kv.first].rate); bench(b, g);
            std::snprintf(b, sizeof(b), "case rate %g, fast", traits[kv.first].rate); bench(b, f);
        }
    }
#ifdef WITH_INTEGRATOR
    if (auto* m = dynamic_cast<IntegImpl*>(impl.get())) saige::gpu2::spaDestroy(m->S);
#endif
    impl.reset();
    destroy(S);
    return rc;
}
