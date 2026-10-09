// stats_selfcheck.cpp -- see stats_selfcheck.hpp. The Makefile compiles this
// translation unit with -ffp-contract=off: the "no fma" candidates below must
// stay unfused, and the fma candidates spell std::fma out.

#include "stats_selfcheck.hpp"
#include "score_format.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace SAIGE {
namespace devstats {

double dotPattern(int t_pat, const double* a, const double* b, int p)
{
    if (p == 1) return (t_pat == 1) ? std::fma(a[0], b[0], 0.0) : a[0] * b[0];
    double acc;
    int k = 2;
    switch (t_pat) {
    case 1:
        acc = std::fma(a[0], b[0], 0.0);
        acc = std::fma(a[1], b[1], acc);
        for (; k < p; ++k) acc = std::fma(a[k], b[k], acc);
        break;
    case 2:
        acc = std::fma(a[0], b[0], a[1] * b[1]);
        for (; k < p; ++k) acc = std::fma(a[k], b[k], acc);
        break;
    case 3:
        acc = a[0] * b[0] + a[1] * b[1];
        for (; k < p; ++k) acc = std::fma(a[k], b[k], acc);
        break;
    case 4:
        acc = std::fma(a[0], b[0], a[1] * b[1]);
        for (; k < p; ++k) acc = acc + a[k] * b[k];
        break;
    default:
        acc = a[0] * b[0] + a[1] * b[1];
        for (; k < p; ++k) acc = acc + a[k] * b[k];
        break;
    }
    return acc;
}

double hadSumPattern(int t_pat, const double* x, const double* y, int p)
{
    double v1 = 0.0, v2 = 0.0;
    int i = 0;
    if (t_pat == 1) {
        for (; i + 1 < p; i += 2) { v1 = std::fma(x[i], y[i], v1); v2 = std::fma(x[i + 1], y[i + 1], v2); }
        if (i < p) v1 = std::fma(x[i], y[i], v1);
    } else {
        for (; i + 1 < p; i += 2) { v1 = v1 + x[i] * y[i]; v2 = v2 + x[i + 1] * y[i + 1]; }
        if (i < p) v1 = v1 + x[i] * y[i];
    }
    return v1 + v2;
}

void refPairStats(const Patterns& t_pat, int p, const double* XVX, const double* Sa,
                  double tau0, const double* z, const double* w, double gr, double g2,
                  double& S, double& var2)
{
    double xz[saige::gpu2::STATS_PMAX];
    for (int i = 0; i < p; ++i) xz[i] = dotPattern(t_pat.xz, XVX + i * saige::gpu2::STATS_PMAX, z, p);
    const double zxz = hadSumPattern(t_pat.zxz, z, xz, p);
    const double saz = dotPattern(t_pat.saz, Sa, z, p);
    const double gwz = hadSumPattern(t_pat.gwz, w, z, p);
    S    = (gr - saz) / tau0;
    var2 = (zxz + g2) - 2.0 * gwz;
}

namespace {

inline bool sameBits(double a, double b)
{
    return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof(double)) == 0;
}

// The trait's constants as the device holds them.
struct TraitConst {
    int p = 0;
    double XVX[saige::gpu2::STATS_PMAX * saige::gpu2::STATS_PMAX];
    double Sa[saige::gpu2::STATS_PMAX];
    double tau0 = 1.0;
};
bool traitConst(const MTContext& ctx, int t, TraitConst& c, std::string& why)
{
    const TraitMeta& M = ctx.meta[t];
    c.p = M.p; c.tau0 = M.tau0;
    if (M.p < 1 || M.p > saige::gpu2::STATS_PMAX) { why = "trait '" + M.name + "': p outside 1.." + std::to_string(saige::gpu2::STATS_PMAX); return false; }
    if ((int)ctx.XVX[t].n_rows != M.p || (int)ctx.XVX[t].n_cols != M.p || (int)ctx.S_a[t].n_elem != M.p) {
        why = "trait '" + M.name + "': XVX / S_a are not p x p / p"; return false;
    }
    std::memset(c.XVX, 0, sizeof(c.XVX));
    std::memset(c.Sa, 0, sizeof(c.Sa));
    for (int i = 0; i < M.p; ++i) {
        c.Sa[i] = ctx.S_a[t](i);
        for (int k = 0; k < M.p; ++k) c.XVX[i * saige::gpu2::STATS_PMAX + k] = ctx.XVX[t](i, k);
    }
    return true;
}

// Random columns of the shape the host tail reads. In real data the three
// terms of var2 = Z'XVX Z + (g^2)'mu2 - 2 W'Z are of one size (W'Z and Z'XVX Z
// are the same quadratic form up to rounding, the covariate adjustment), so
// W is drawn around XVX Z, (g^2)'mu2 around 1.2..2.2 x Z'XVX Z, and g'res so
// that stat = S^2 / var2 spans 0 .. ~80: every term's last bit then reaches
// var2, which is what tells the contraction orders apart.
struct RandCols {
    std::mt19937_64 rng{20261009ull};
    std::normal_distribution<double> nd{0.0, 1.0};
    std::uniform_real_distribution<double> ud{0.0, 1.0};
};
inline double zScale(const TraitConst& c)
{
    double s = 0.0;
    for (int i = 0; i < c.p; ++i) s += c.Sa[i] * c.Sa[i];
    s = std::sqrt(s);
    return (s > 1e-300 && std::isfinite(s)) ? 5.0 / s : 1.0;
}
// One pair's columns: z, w (p each), gr, g2. sMode 0: g'res such that S is of
// the size of S_a'Z itself, so the last bit of S_a'Z reaches S (that is what
// tells its contraction order apart; with the model's small S_a a
// statistic-sized S would swallow it); 1: g'res such that stat = S^2 / var2
// spans 0 .. ~80 (realistic gates). The probes alternate the two.
inline void drawPair(RandCols& rc, const TraitConst& c, double* z, double* w, double& gr, double& g2, int sMode)
{
    const double sZ = zScale(c);
    double xz[saige::gpu2::STATS_PMAX];
    double zxz = 0.0, saz = 0.0;
    for (int i = 0; i < c.p; ++i) z[i] = rc.nd(rc.rng) * sZ;
    for (int i = 0; i < c.p; ++i) {
        double a = 0.0;
        for (int k = 0; k < c.p; ++k) a += c.XVX[i * saige::gpu2::STATS_PMAX + k] * z[k];
        xz[i] = a;
        zxz += z[i] * a;
        saz += c.Sa[i] * z[i];
    }
    for (int i = 0; i < c.p; ++i) w[i] = xz[i] * (1.0 + 0.1 * rc.nd(rc.rng));
    const double zabs = std::fabs(zxz) + 1e-300;
    g2 = zabs * (1.2 + rc.ud(rc.rng));
    const double var2 = g2 - zxz;
    const double Starget = (sMode == 0)
        ? rc.nd(rc.rng) * (std::fabs(saz) + 1e-300) * (0.3 + 3.0 * rc.ud(rc.rng))
        : rc.nd(rc.rng) * std::sqrt(std::fabs(var2)) * (0.5 + 2.5 * rc.ud(rc.rng));
    gr = saz + Starget * c.tau0;
}

}  // namespace

std::string identifyHostPatterns(const MTContext& ctx, const std::vector<int>& traits,
                                 Patterns& out, std::string& detail)
{
    out = Patterns();
    if (traits.empty()) return "no binary trait";
    if (ctx.sampleSetsDiffer) return "the models do not share one sample list";
    std::vector<TraitConst> tc(ctx.P);
    for (int t : traits) { std::string why; if (!traitConst(ctx, t, tc[t], why)) return why; }

    const int B = 128;
    RandCols rc;
    MTScratch scr;
    scr.Zall.set_size(ctx.sumP, B);
    scr.GWbin.set_size(std::max(ctx.sumPbin, 1), B);
    scr.GR.set_size(B, ctx.P);
    scr.G2Mu2.set_size(B, std::max(ctx.nBin, 1));
    scr.Gsq.set_size(B);
    if (ctx.sumPbin != (int)scr.GWbin.n_rows || ctx.nBin != (int)scr.G2Mu2.n_cols) return "no binary stack block";
    MTBlockResult res;
    res.resize(B, ctx.P);
    arma::mat VR(B, ctx.P);
    scr.Zall.zeros(); scr.GWbin.zeros(); scr.GR.zeros(); scr.G2Mu2.zeros(); scr.Gsq.ones(); VR.ones();
    for (int t : traits) {
        const TraitMeta& M = ctx.meta[t];
        for (int j = 0; j < B; ++j) {
            double z[saige::gpu2::STATS_PMAX], w[saige::gpu2::STATS_PMAX], gr, g2;
            drawPair(rc, tc[t], z, w, gr, g2, j & 1);
            for (int i = 0; i < M.p; ++i) {
                scr.Zall(M.colOff + i, j)  = z[i];
                scr.GWbin(M.binOff + i, j) = w[i];
            }
            scr.GR(j, t) = gr;
            scr.G2Mu2(j, M.binIdx) = g2;
            VR(j, t) = 0.9 + 0.2 * rc.ud(rc.rng);
        }
    }
    try {
        scoreTestBatchMTBinPre(ctx, traits, 0, B, VR, scr, res);
    } catch (const std::exception& e) {
        return std::string("scoreTestBatchMTBinPre refused the probe block: ") + e.what();
    }

    long mismS[5] = {0, 0, 0, 0, 0};
    long mismV[5][2][2];
    std::memset(mismV, 0, sizeof(mismV));
    long nPairs = 0;
    for (int t : traits) {
        const TraitMeta& M = ctx.meta[t];
        const TraitConst& c = tc[t];
        for (int j = 0; j < B; ++j) {
            double z[saige::gpu2::STATS_PMAX], w[saige::gpu2::STATS_PMAX], xz[saige::gpu2::STATS_PMAX];
            for (int i = 0; i < M.p; ++i) { z[i] = scr.Zall(M.colOff + i, j); w[i] = scr.GWbin(M.binOff + i, j); }
            const double gr = scr.GR(j, t), g2 = scr.G2Mu2(j, M.binIdx);
            const double Sh = res.Tstat(j, t), v2h = res.var2(j, t);
            ++nPairs;
            for (int ps = 0; ps < 5; ++ps) {
                const double S = (gr - dotPattern(ps, c.Sa, z, M.p)) / M.tau0;
                if (!sameBits(S, Sh)) ++mismS[ps];
            }
            for (int px = 0; px < 5; ++px) {
                for (int i = 0; i < M.p; ++i) xz[i] = dotPattern(px, c.XVX + i * saige::gpu2::STATS_PMAX, z, M.p);
                for (int pz = 0; pz < 2; ++pz) {
                    const double zxz = hadSumPattern(pz, z, xz, M.p);
                    for (int pg = 0; pg < 2; ++pg) {
                        const double gwz = hadSumPattern(pg, w, z, M.p);
                        const double v2 = (zxz + g2) - 2.0 * gwz;
                        if (!sameBits(v2, v2h)) ++mismV[px][pz][pg];
                    }
                }
            }
        }
    }
    int saz = -1, nSaz = 0;
    for (int ps = 0; ps < 5; ++ps) if (mismS[ps] == 0) { if (saz < 0) saz = ps; ++nSaz; }
    int xz = -1, zxz = -1, gwz = -1;
    for (int px = 0; px < 5 && xz < 0; ++px)
        for (int pz = 0; pz < 2 && xz < 0; ++pz)
            for (int pg = 0; pg < 2 && xz < 0; ++pg)
                if (mismV[px][pz][pg] == 0) { xz = px; zxz = pz; gwz = pg; }
    std::ostringstream d;
    d << nPairs << " probe pairs; S_a'Z pattern mismatches";
    for (int ps = 0; ps < 5; ++ps) d << " " << mismS[ps];
    d << "; var2 pattern (XVX Z, Z%XZ sum, W%Z sum) mismatches";
    for (int px = 0; px < 5; ++px) for (int pz = 0; pz < 2; ++pz) for (int pg = 0; pg < 2; ++pg)
        d << " " << px << pz << pg << ":" << mismV[px][pz][pg];
    detail = d.str();
    if (saz < 0) return "no candidate order reproduces the host's S_a'Z (" + detail + ")";
    if (xz < 0) return "no candidate order reproduces the host's var2 (" + detail + ")";
    // p >= 3 separates every candidate on these probes; several matches mean
    // the probe could not see the last bit, and a guess is not good enough.
    int nV = 0;
    for (int px = 0; px < 5; ++px) for (int pz = 0; pz < 2; ++pz) for (int pg = 0; pg < 2; ++pg) if (mismV[px][pz][pg] == 0) ++nV;
    int pmax = 0;
    for (int t : traits) pmax = std::max(pmax, ctx.meta[t].p);
    if (pmax >= 3 && (nSaz > 1 || nV > 1))
        return "the probe does not separate the candidate orders (" + detail + ")";
    out.saz = saz; out.xz = xz; out.zxz = zxz; out.gwz = gwz;
    return std::string();
}

void fillStatsTraits(const MTContext& ctx, const std::vector<int>& binTraits,
                     const std::vector<char>& enabled, int oW, int oR,
                     std::vector<saige::gpu2::StatsTrait>& out)
{
    out.assign((std::size_t)ctx.nBin, saige::gpu2::StatsTrait());
    for (int t : binTraits) {
        const TraitMeta& M = ctx.meta[t];
        saige::gpu2::StatsTrait& T = out[(std::size_t)M.binIdx];
        T.p = M.p; T.rowZ = M.colOff; T.rowW = oW + M.binOff; T.rowGR = oR + t; T.colG2 = M.binIdx;
        T.enabled = enabled[(std::size_t)M.binIdx] ? 1 : 0;
        T.isFirth = M.is_Firth_beta ? 1 : 0;
        T.isFast = M.isFastTest ? 1 : 0;
        T.tau0 = M.tau0; T.spaCut = M.SPA_Cutoff; T.firthCut = M.pCutoffforFirth; T.fastCut = M.pval_cutoff_for_fastTest;
        std::memset(T.XVX, 0, sizeof(T.XVX));
        std::memset(T.Sa, 0, sizeof(T.Sa));
        if (M.p >= 1 && M.p <= saige::gpu2::STATS_PMAX) {
            for (int i = 0; i < M.p; ++i) {
                T.Sa[i] = ctx.S_a[t](i);
                for (int k = 0; k < M.p; ++k) T.XVX[i * saige::gpu2::STATS_PMAX + k] = ctx.XVX[t](i, k);
            }
        } else {
            T.enabled = 0;
        }
    }
}

std::string deviceSelfTest(saige::gpu2::Reducer* R, const MTContext& ctx,
                           const std::vector<int>& binTraits, const std::vector<char>& enabled,
                           int K1, int K2, int oW, int oR, int maxSlots, int Bblk,
                           double statCutoff, double pRelTol, std::string& detail)
{
    using namespace saige::gpu2;
    std::vector<int> traits;
    for (int t : binTraits) if (enabled[(std::size_t)ctx.meta[t].binIdx]) traits.push_back(t);
    if (traits.empty()) return "no binary trait enabled";
    std::vector<TraitConst> tc(ctx.P);
    for (int t : traits) { std::string why; if (!traitConst(ctx, t, tc[t], why)) return why; }

    // ---- (a) the device's erfc against boost's tail ----
    {
        const int n = (1 << 17) + 2048;
        std::vector<double> st((std::size_t)n), pd((std::size_t)n);
        std::mt19937_64 rng(20261010ull);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < (1 << 17); ++i) {
            const double r = u(rng);
            st[(std::size_t)i] = (i & 1) ? r * statCutoff : std::exp(std::log(1e-12) * r) * statCutoff;
        }
        for (int i = 0; i < 1024; ++i) st[(std::size_t)(1 << 17) + i] = statCutoff * (1.0 - 1e-9 * (i + 1));
        for (int i = 0; i < 1024; ++i) st[(std::size_t)(1 << 17) + 1024 + i] = 1e-300 * std::pow(10.0, i * 0.29);
        if (!statsErfcTest(R, st.data(), n, pd.data())) return "the device erfc probe did not run (" + statsLastError() + ")";
        double maxRel = 0.0;
        long nBad = 0;
        for (int i = 0; i < n; ++i) {
            const double s = st[(std::size_t)i];
            boost::math::chi_squared chisq(1);
            const double pb = boost::math::cdf(complement(chisq, s));
            const double rel = std::fabs(pd[(std::size_t)i] - pb) / pb;
            if (!(rel <= pRelTol * 0.5)) ++nBad;
            if (rel > maxRel) maxRel = rel;
        }
        std::ostringstream d;
        d << "device erfc vs boost on " << n << " stats in [0, " << statCutoff << "): max rel diff " << maxRel;
        detail = d.str();
        if (nBad > 0) {
            std::ostringstream w;
            w << nBad << " of " << n << " device p-values differ from boost by more than " << pRelTol * 0.5 << " relative (max " << maxRel << ")";
            return w.str();
        }
    }

    // ---- (b) the kernel against scoreTestBatchMTBinPre, bit for bit ----
    const int nS = std::min(maxSlots, 4 * Bblk);
    std::vector<double> C1((std::size_t)maxSlots * K1, 0.0), C2((std::size_t)maxSlots * std::max(K2, 1), 0.0);
    std::vector<double> vr((std::size_t)nS * ctx.nBin, 1.0);
    RandCols rc;
    rc.rng.seed(20261011ull);
    for (int t : traits) {
        const TraitMeta& M = ctx.meta[t];
        for (int s = 0; s < nS; ++s) {
            double z[saige::gpu2::STATS_PMAX], w[saige::gpu2::STATS_PMAX], gr, g2;
            drawPair(rc, tc[t], z, w, gr, g2, s & 1);
            // a few exactly-degenerate columns too (zero columns, as an unused slot)
            const bool zero = (s % 97 == 5);
            for (int i = 0; i < M.p; ++i) {
                C1[(std::size_t)(M.colOff + i) * maxSlots + s]      = zero ? 0.0 : z[i];
                C1[(std::size_t)(oW + M.binOff + i) * maxSlots + s] = zero ? 0.0 : w[i];
            }
            C1[(std::size_t)(oR + t) * maxSlots + s] = zero ? 0.0 : gr;
            C2[(std::size_t)M.binIdx * maxSlots + s] = zero ? 0.0 : g2;
            vr[(std::size_t)s * ctx.nBin + M.binIdx] = 0.9 + 0.2 * rc.ud(rc.rng);
        }
    }
    if (!statsSelfTest(R, C1.data(), C2.data(), vr.data(), nS)) return "the device stats probe did not run (" + statsLastError() + ")";
    const double* dS = statsS(R, 0);
    const double* dV = statsVar2(R, 0);
    const double* dP = statsP(R, 0);
    const unsigned char* dF = statsFlags(R, 0);
    if (!dS || !dV || !dP || !dF) return "no device stats buffers";

    MTScratch scr;
    scr.Zall.set_size(ctx.sumP, Bblk);
    scr.GWbin.set_size(ctx.sumPbin, Bblk);
    scr.GR.set_size(Bblk, ctx.P);
    scr.G2Mu2.set_size(Bblk, ctx.nBin);
    scr.Gsq.set_size(Bblk);
    MTBlockResult res;
    res.resize(Bblk, ctx.P);
    arma::mat VR(Bblk, ctx.P);
    long nPairs = 0, mismS = 0, mismV = 0, nHost = 0, gateBad = 0, nGate = 0;
    long reason[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double maxRel = 0.0;
    for (int j0 = 0; j0 < nS; j0 += Bblk) {
        const int nb = std::min(Bblk, nS - j0);
        scr.Zall.zeros(); scr.GWbin.zeros(); scr.GR.zeros(); scr.G2Mu2.zeros(); scr.Gsq.ones(); VR.ones();
        for (int r = 0; r < ctx.sumP; ++r)
            for (int j = 0; j < nb; ++j) scr.Zall(r, j) = C1[(std::size_t)r * maxSlots + j0 + j];
        for (int r = 0; r < ctx.sumPbin; ++r)
            for (int j = 0; j < nb; ++j) scr.GWbin(r, j) = C1[(std::size_t)(oW + r) * maxSlots + j0 + j];
        for (int t = 0; t < ctx.P; ++t)
            for (int j = 0; j < nb; ++j) scr.GR(j, t) = C1[(std::size_t)(oR + t) * maxSlots + j0 + j];
        for (int b = 0; b < ctx.nBin; ++b)
            for (int j = 0; j < nb; ++j) scr.G2Mu2(j, b) = C2[(std::size_t)b * maxSlots + j0 + j];
        for (int t : traits)
            for (int j = 0; j < nb; ++j) VR(j, t) = vr[(std::size_t)(j0 + j) * ctx.nBin + ctx.meta[t].binIdx];
        try {
            scoreTestBatchMTBinPre(ctx, traits, 0, nb, VR, scr, res);
        } catch (const std::exception& e) {
            return std::string("scoreTestBatchMTBinPre refused the probe block: ") + e.what();
        }
        for (int t : traits) {
            const TraitMeta& M = ctx.meta[t];
            for (int j = 0; j < nb; ++j) {
                const std::size_t o = (std::size_t)(j0 + j) * ctx.nBin + M.binIdx;
                ++nPairs;
                const bool badS = !sameBits(dS[o], res.Tstat(j, t));
                const bool badV = !sameBits(dV[o], res.var2(j, t));
                if (badS) ++mismS;
                if (badV) ++mismV;
                if ((badS || badV) && std::getenv("SAIGE_DEVSTATS_DEBUG")) {
                    // the pair's inputs and both results, for a look at the order
                    std::printf("  devstats mismatch trait %s slot %d:%s%s\n", M.name.c_str(), j0 + j, badS ? " S" : "", badV ? " var2" : "");
                    std::printf("    tau0 %.17g gr %.17g g2 %.17g vr %.17g\n", M.tau0, C1[(std::size_t)(oR + t) * maxSlots + j0 + j],
                                C2[(std::size_t)M.binIdx * maxSlots + j0 + j], vr[(std::size_t)(j0 + j) * ctx.nBin + M.binIdx]);
                    for (int i = 0; i < M.p; ++i)
                        std::printf("    i %d z %.17g w %.17g Sa %.17g XVX row %.17g %.17g %.17g\n", i,
                                    C1[(std::size_t)(M.colOff + i) * maxSlots + j0 + j], C1[(std::size_t)(oW + M.binOff + i) * maxSlots + j0 + j],
                                    tc[t].Sa[i], tc[t].XVX[i * saige::gpu2::STATS_PMAX], tc[t].XVX[i * saige::gpu2::STATS_PMAX + 1],
                                    M.p > 2 ? tc[t].XVX[i * saige::gpu2::STATS_PMAX + 2] : 0.0);
                    std::printf("    host S %.17g var2 %.17g | device S %.17g var2 %.17g\n", res.Tstat(j, t), res.var2(j, t), dS[o], dV[o]);
                    double zz[saige::gpu2::STATS_PMAX];
                    for (int i = 0; i < M.p; ++i) zz[i] = C1[(std::size_t)(M.colOff + i) * maxSlots + j0 + j];
                    for (int ps = 0; ps < 5; ++ps)
                        std::printf("    saz pattern %d -> S %.17g\n", ps,
                                    (C1[(std::size_t)(oR + t) * maxSlots + j0 + j] - dotPattern(ps, tc[t].Sa, zz, M.p)) / M.tau0);
                }
                const unsigned f = dF[o];
                if (f & STATS_HOST) { ++nHost; ++reason[(f >> STATS_REASON_SHIFT) & 7u]; continue; }
                ++nGate;
                const double ph = res.pvalRaw(j, t);
                const bool islog = res.pvalIsLog[(std::size_t)t][(std::size_t)j] != 0;
                const double rel = islog ? 1.0 : std::fabs(dP[o] - ph) / ph;
                if (rel > maxRel) maxRel = rel;
                const double sd = res.StdStat(j, t);
                const bool spaH = !std::isnan(sd) && sd > M.SPA_Cutoff;
                const bool firthH = M.is_Firth_beta && (islog ? ph <= std::log(M.pCutoffforFirth) : ph <= M.pCutoffforFirth);
                bool fastH = false;
                if (M.isFastTest) {
                    double pnum;
                    try { pnum = std::stod(res.pvalStr[(std::size_t)t][(std::size_t)j]); }
                    catch (...) { pnum = 0; }
                    fastH = pnum < M.pval_cutoff_for_fastTest;
                }
                if (spaH != ((f & STATS_SPA) != 0) || firthH != ((f & STATS_FIRTH) != 0) ||
                    fastH != ((f & STATS_FAST) != 0) || islog || !(rel <= pRelTol))
                    ++gateBad;
            }
        }
    }
    std::ostringstream d;
    d << detail << "; kernel vs host on " << nPairs << " random pairs: S bit-mismatches " << mismS
      << ", var2 " << mismV << "; " << nHost << " flagged for the host (tail " << reason[STATS_R_TAIL]
      << ", degenerate " << reason[STATS_R_DEGEN] << ", print " << reason[STATS_R_PRINT]
      << ", Firth cutoff " << reason[STATS_R_FIRTH] << ", fast cutoff " << reason[STATS_R_FASTC] << "), "
      << nGate << " device-decided: gate / p disagreements " << gateBad << ", max p rel diff " << maxRel;
    detail = d.str();
    if (mismS > 0 || mismV > 0) return "device S / var2 differ from scoreTestBatchMTBinPre in some bit (" + detail + ")";
    if (gateBad > 0) return "device gate bits or p disagree with the host's (" + detail + ")";
    return std::string();
}

}  // namespace devstats
}  // namespace SAIGE
