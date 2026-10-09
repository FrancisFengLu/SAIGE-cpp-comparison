// spa_fp32_check -- gpuPrecisionSPA fp32 against fp64, pair by pair, for both
// device SPA implementations (gpuSpaImpl: lib = spa_gpu.hpp, own =
// gpu_spa.hpp), with the production CPU SPA alongside. Same inputs and
// selection as tools/spa_gpu_check.cpp (which it is derived from): the shipped
// PLINK decode, scoreTestFast's Tstat / var1 / var2 / pval_noadj, |z| > zcut,
// or --sweep z1,z2,... to force +-z on every (marker, trait) -- which is how
// p far below 1e-100 and the log-p domain (|z| > ~37) are reached.
//
// For each implementation it runs the fp64 and the fp32 kernel on the same
// pair table and reports, overall and by |z| bin / variant / domain:
// convergence, saddle-flag and final-convergence (after the host post-rules)
// mismatches, iteration-count mismatches, the root difference, the relative
// p difference (|dln p| in the log domain), p == 0 / non-finite results, and
// the worst pairs with the library's status words.
//
//   spa_fp32_check BED MODELDIR:VRFILE[,...] [--markers M] [--zcut 2]
//                  [--maxpairs K] [--sweep z1,z2,...] [--worst 10] [--flipmu]
//
// --flipmu hands every SPA (CPU and device) mu' = 1 - mu instead of the
// model's mu (Tstat / var come from the model as usual), so the fp32 path's
// mu > 1/2 branch is exercised on data whose case rates are all <= 1/2.
// --full forces the full-N variant on every pair (as a sparse-GRM model does).
#include <armadillo>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "../../genotype_reader.hpp"
#include "../../null_model_loader.hpp"
#include "../../saige_test.hpp"
#include "../../score_format.hpp"
#include "../../spa.hpp"
#include "../../spa_binary.hpp"
#include "../../UTIL.hpp"
#include "../gpu_step2.hpp"
#include "../gpu_spa.hpp"
#include "spa_gpu.hpp"

namespace {

struct Model {
    NullModelData nm;
    SAIGE::SAIGEClass* obj = nullptr;
    arma::vec mu;
    bool singleVR = false;
};

std::vector<std::string> split(const std::string& s, char c)
{
    std::vector<std::string> v; std::string cur;
    for (char ch : s) { if (ch == c) { v.push_back(cur); cur.clear(); } else cur += ch; }
    v.push_back(cur); return v;
}

Model* loadModel(const std::string& spec)
{
    auto a = split(spec, ':');
    if (a.size() != 2) { std::fprintf(stderr, "model spec must be DIR:VRFILE\n"); std::exit(2); }
    Model* M = new Model();
    M->nm = loadNullModel(a[0], a[1]);
    NullModelData& nm = M->nm;
    M->obj = new SAIGE::SAIGEClass(nm.XVX, nm.XXVX_inv, nm.XV, nm.XVX_inv_XV, nm.Sigma_iXXSigma_iX, nm.X,
        nm.S_a, nm.res, nm.mu2, nm.mu, nm.varRatio_sparse, nm.varRatio_null, nm.varRatio_null_noXadj,
        nm.cateVarRatioMinMACVecExclude, nm.cateVarRatioMaxMACVecInclude, nm.SPA_Cutoff, nm.tauvec,
        nm.traitType, nm.y, nm.impute_method, nm.flagSparseGRM, nm.isFastTest, nm.isnoadjCov,
        nm.pval_cutoff_for_fastTest, nm.locationMat, nm.valueVec, nm.dimNum, nm.isCondition,
        nm.condition_genoIndex, nm.is_Firth_beta, nm.pCutoffforFirth, nm.offset, nm.resout);
    // As mainMarkerMT: one category -> the single ratio; several -> per MAC.
    M->singleVR = (M->obj->m_varRatio_null.n_elem == 1);
    if (M->singleVR) M->obj->assignSingleVarianceRatio(false, false);
    M->obj->get_mu(M->mu);
    return M;
}

double traitVR(Model* M, double MAC)
{
    if (M->singleVR) return M->obj->computeSingleVarianceRatio(false, false);
    bool has;
    return M->obj->computeVarianceRatio(MAC, false, false, has);
}

// The host post-SPA rules of getMarkerPval: the quantile step (which can
// withdraw convergence), the p == 0 rule, and the printed string.
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
    } else {
        str = noSpaStr;
    }
    return conv;
}

// One pair: inputs, the CPU result, and per (impl, precision) the device result.
struct Res { double r1, r2, p; int n1, n2, s1, s2, conv, convF; unsigned status; int rs1, rs2; std::string str; };
struct Rec {
    int m, t; double z; int fast, logp;
    double Tstat, var1, var2, pno;
    std::string noSpa;
    double pc; int convc, convFc; std::string strc;
    Res g[2][2];   // [impl: 0 lib, 1 own][prec: 0 fp64, 1 fp32]
};

struct Cmp {
    long n = 0, conv = 0, sadd = 0, convF = 0, iter = 0, strd = 0, pzero = 0, pnonfin = 0, prel1e6 = 0;
    double maxrel = 0, maxrelSmall = 0, maxroot = 0, maxrelConv = 0;
};

// Relative p difference: |a - b| / |a| in the linear domain, |a - b| (= |d ln p|) in the log domain.
double prel(double a, double b, bool logp)
{
    if (a == b) return 0;
    if (std::isnan(a) && std::isnan(b)) return 0;
    if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::infinity();
    if (logp) return std::fabs(a - b);
    return std::fabs(a - b) / std::max(std::fabs(a), 5e-324);
}

void tally(Cmp& c, const Rec& r, const Res& A, const Res& B)
{
    c.n++;
    if (A.conv != B.conv) c.conv++;
    if (A.s1 != B.s1 || A.s2 != B.s2) c.sadd++;
    if (A.convF != B.convF) c.convF++;
    if (A.n1 != B.n1 || A.n2 != B.n2) c.iter++;
    if (A.str != B.str) c.strd++;
    if (!r.logp && B.p == 0 && A.p != 0) c.pzero++;
    if (!std::isfinite(B.p) && std::isfinite(A.p)) c.pnonfin++;
    auto rd = [](double a, double b) { if (std::isinf(a) && std::isinf(b) && (a > 0) == (b > 0)) return 0.0;
                                       if (std::isinf(a) != std::isinf(b)) return std::numeric_limits<double>::infinity();
                                       return std::fabs(a - b) / std::max(1.0, std::fabs(a)); };
    const double dr = std::max(rd(A.r1, B.r1), rd(A.r2, B.r2));
    if (A.conv && B.conv && dr > c.maxroot) c.maxroot = dr;
    const double rel = prel(A.p, B.p, r.logp);
    if (rel > c.maxrel) c.maxrel = rel;
    if (rel > 1e-6) c.prel1e6++;
    if (A.convF && B.convF && rel > c.maxrelConv) c.maxrelConv = rel;
    const bool small = r.logp || A.p < 1e-5;
    if (small && rel > c.maxrelSmall) c.maxrelSmall = rel;
}

void printCmp(const char* name, const Cmp& c)
{
    if (c.n == 0) return;
    std::printf("  %-24s n=%7ld conv %4ld saddle %4ld convFinal %4ld iter %6ld string %6ld p==0 %ld nonfinite %ld | "
                "max rel p %.2e (both SPA %.2e, p<1e-5/logp %.2e) rel>1e-6 %ld | max rel root %.1e\n",
                name, c.n, c.conv, c.sadd, c.convF, c.iter, c.strd, c.pzero, c.pnonfin, c.maxrel, c.maxrelConv,
                c.maxrelSmall, c.prel1e6, c.maxroot);
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: spa_fp32_check BED MODELDIR:VRFILE[,...] [--markers M] [--zcut Z] [--maxpairs K] [--sweep z,..] [--worst W] [--flipmu]\n"); return 2; }
    const std::string bed = argv[1];
    long nMarkers = -1; double zcut = 2.0; long maxPairs = 20000; std::vector<double> sweep; int nWorst = 10; bool flipmu = false, forceFull = false;
    for (int i = 3; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--markers" && i + 1 < argc) nMarkers = std::atol(argv[++i]);
        else if (a == "--zcut" && i + 1 < argc) zcut = std::atof(argv[++i]);
        else if (a == "--maxpairs" && i + 1 < argc) maxPairs = std::atol(argv[++i]);
        else if (a == "--sweep" && i + 1 < argc) for (auto& s : split(argv[++i], ',')) sweep.push_back(std::atof(s.c_str()));
        else if (a == "--worst" && i + 1 < argc) nWorst = std::atoi(argv[++i]);
        else if (a == "--flipmu") flipmu = true;
        else if (a == "--full") forceFull = true;
        else { std::fprintf(stderr, "bad arg %s\n", a.c_str()); return 2; }
    }
    std::vector<Model*> models;
    for (auto& s : split(argv[2], ',')) models.push_back(loadModel(s));
    const int T = (int)models.size();
    if (flipmu) for (auto* Mo : models) Mo->mu = 1.0 - Mo->mu;
    std::vector<std::string> ids = models[0]->nm.sampleIDs;
    setPLINKobjInCPP(bed + ".bim", bed + ".fam", bed + ".bed", ids, "alt-first");
    const int n = (int)ptr_gPLINKobj->getN();
    const long M = (long)ptr_gPLINKobj->getM();
    if (nMarkers < 0 || nMarkers > M) nMarkers = M;
    for (auto* Mo : models) if (Mo->obj->m_n != n) { std::fprintf(stderr, "model n != reader n\n"); return 2; }
    const int imputeCase = string_to_case.at(models[0]->nm.impute_method);
    std::fprintf(stderr, "N=%d markers=%ld (using %ld) traits=%d sweep=%zu\n", n, M, nMarkers, T, sweep.size());

    const int slotsCap = 4096;
    std::vector<double> B1((std::size_t)n, 0.0);
    saige::gpu2::CreateArgs ca; ca.device = 0; ca.N = n; ca.K1 = 1; ca.B1 = B1.data(); ca.maxSlots = slotsCap; ca.precision = saige::gpu2::Prec::FP64;
    saige::gpu2::Reducer* R = saige::gpu2::create(ca);
    if (!R) { std::fprintf(stderr, "reducer create failed\n"); return 3; }
    const double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
    const int maxP = (int)std::min<long>(maxPairs, (long)slotsCap * T * (long)std::max<size_t>(1, 2 * sweep.size()));
    // own
    std::vector<saige::gpu2::SpaTraitArgs> ta(T);
    std::vector<saige::spa_gpu::TraitArgs> tl(T);
    for (int t = 0; t < T; t++) {
        ta[t].mu = models[t]->mu.memptr(); ta[t].XV = models[t]->obj->m_XV.memptr(); ta[t].XXVX_inv = models[t]->obj->m_XXVX_inv.memptr(); ta[t].p = models[t]->obj->m_p;
        tl[t].mu = ta[t].mu; tl[t].XV = ta[t].XV; tl[t].XXVX_inv = ta[t].XXVX_inv; tl[t].p = ta[t].p;
    }
    saige::gpu2::Spa* SO[2]; saige::spa_gpu::Spa* SL[2];
    for (int pr = 0; pr < 2; pr++) {
        const auto prec = pr ? saige::gpu2::Prec::FP32 : saige::gpu2::Prec::FP64;
        saige::gpu2::SpaCreateArgs sa; sa.device = 0; sa.N = n; sa.nTraits = T; sa.traits = ta.data();
        sa.maxPairs = maxP; sa.tol = tol; sa.maxiter = 1000; sa.blocks = 256; sa.precision = prec;
        SO[pr] = saige::gpu2::spaCreate(sa);
        saige::spa_gpu::CreateArgs la; la.device = 0; la.N = n; la.nTraits = T; la.traits = tl.data();
        la.maxPairs = maxP; la.tol = tol; la.maxiter = 1000; la.blocks = 256; la.precision = prec;
        la.minBlocksPerSM = 3;   // main.cpp's default with useGPU
        SL[pr] = saige::spa_gpu::create(la);
        if (!SO[pr] || !SL[pr]) { std::fprintf(stderr, "spa create failed: %s\n", saige::spa_gpu::lastError()); return 3; }
    }
    const std::size_t bpv = saige::gpu2::bytesPerSlot(R);
    unsigned char* hPk = saige::gpu2::packed(R);
    double* hLu = saige::gpu2::lut(R);

    std::vector<Rec> recs;
    for (long m0 = 0; m0 < nMarkers && (long)recs.size() < maxPairs; m0 += slotsCap) {
        const long m1 = std::min(nMarkers, m0 + (long)slotsCap);
        const int nS = (int)(m1 - m0);
        std::vector<arma::vec> gs(nS); std::vector<arma::uvec> izs(nS), inzs(nS);
        std::vector<double> af(nS), mac(nS); std::vector<char> okm(nS, 0);
        #pragma omp parallel for schedule(dynamic, 8)
        for (int s = 0; s < nS; s++) {
            PLINK::PlinkClass::FusedMarkerStats fs;
            if (!ptr_gPLINKobj->getOneMarkerFusedStats_ts((uint64_t)(m0 + s), fs)) continue;
            double MAF = std::min(fs.altFreq, 1 - fs.altFreq);
            double MAC = MAF * n * (1 - fs.missingRate) * 2;
            if (MAC < 1 || fs.missingRate > 0.15) continue;
            PLINK::finalizeFusedStats(fs, imputeCase, 0.2, 10.0, MAC);
            ptr_gPLINKobj->fillOneMarkerFusedDense_ts(fs, gs[s], izs[s], inzs[s]);
            ptr_gPLINKobj->copyFusedPacked_ts(fs, hPk + (std::size_t)s * bpv);
            for (int k = 0; k < 4; k++) hLu[(std::size_t)s * 4 + k] = fs.fd[k];
            af[s] = fs.altFreq_post; okm[s] = 1;
            mac[s] = std::min(fs.altCounts_post, 2.0 * n - fs.altCounts_post);
        }
        for (int s = 0; s < nS; s++) if (!okm[s]) { std::memset(hPk + (std::size_t)s * bpv, 0, (n + 3) / 4); for (int k = 0; k < 4; k++) hLu[(std::size_t)s * 4 + k] = 0; }
        if (!saige::gpu2::reduce(R, nS)) { std::fprintf(stderr, "reduce failed\n"); return 3; }

        struct Job { int s, t; double z; };
        std::vector<Job> jobs;
        for (int s = 0; s < nS; s++) if (okm[s]) for (int t = 0; t < T; t++) {
            if (sweep.empty()) jobs.push_back({s, t, 0.0});
            else for (double z : sweep) { jobs.push_back({s, t, z}); jobs.push_back({s, t, -z}); }
        }
        std::vector<Rec> batch(jobs.size()); std::vector<char> keep(jobs.size(), 0);
        #pragma omp parallel for schedule(dynamic, 16)
        for (std::size_t j = 0; j < jobs.size(); j++) {
            const Job& J = jobs[j]; Model* Mo = models[J.t]; SAIGE::SAIGEClass* obj = Mo->obj;
            arma::vec& g = gs[J.s]; arma::uvec& inz = inzs[J.s]; arma::uvec& iz = izs[J.s];
            SAIGE::PerMarkerCtx ctx; ctx.flagSparseGRM_cur = false; ctx.isnoadjCov_cur = false; ctx.varRatioVal = traitVR(Mo, mac[J.s]);
            double Beta, seBeta, pval, Tstat, var1, var2; bool islogp; std::string pstr;
            obj->scoreTestFast(g, inz, Beta, seBeta, pstr, pval, islogp, af[J.s], Tstat, var1, var2, ctx);
            if (J.z != 0.0) {
                Tstat = J.z * std::sqrt(var1);
                const double stat = Tstat * Tstat / var1;
                boost::math::chi_squared chisq(1);
                pval = boost::math::cdf(complement(chisq, stat)); islogp = false;
                char buf[100];
                if (pval != 0) std::snprintf(buf, sizeof(buf), "%.6E", pval);
                else { pval = SAIGE::log_chisq1_uppertail(stat); islogp = true;
                       const double l10 = pval / std::log(10.0); int ex = (int)std::floor(l10); double fr = std::pow(10.0, l10 - ex);
                       if (fr >= 9.95) { fr = 1; ex++; } std::snprintf(buf, sizeof(buf), "%.1fE%d", fr, ex); }
                pstr = buf;
            }
            const double StdStat = std::fabs(Tstat) / std::sqrt(var1);
            if (!(StdStat > zcut) && J.z == 0.0) continue;
            Rec& r = batch[j]; keep[j] = 1;
            r.m = (int)(m0 + J.s); r.t = J.t; r.z = (J.z != 0.0) ? J.z : StdStat; r.logp = islogp;
            r.Tstat = Tstat; r.var1 = var1; r.var2 = var2; r.pno = pval; r.noSpa = pstr;
            arma::vec gt(n); obj->getadjGFast(g, gt, inz);
            const double m1 = arma::dot(Mo->mu, gt);
            const bool fast = !forceFull && (double(iz.n_elem) / n >= 0.5);
            r.fast = fast;
            double q = Tstat / std::sqrt(var1 / var2) + m1, qinv;
            if ((q - m1) > 0) qinv = -1 * std::fabs(q - m1) + m1; else if ((q - m1) == 0) qinv = m1; else qinv = std::fabs(q - m1) + m1;
            double pw; bool cw;
            if (fast) {
                arma::vec gNB = gt(inz), gNA = gt(iz), muNB = Mo->mu(inz), muNA = Mo->mu(iz);
                const double NAmu = m1 - arma::dot(gNB, muNB), NAsigma = var2 - arma::sum(muNB % (1 - muNB) % arma::pow(gNB, 2));
                SPA_fast(Mo->mu, gt, q, qinv, pval, islogp, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol, "binary", pw, cw);
            } else {
                SPA(Mo->mu, gt, q, qinv, pval, tol, islogp, "binary", pw, cw);
            }
            r.pc = pw; r.convc = cw;
            r.convFc = postRules(pw, islogp, cw, r.strc, pstr);
        }
        std::vector<int> idx;
        for (std::size_t j = 0; j < jobs.size(); j++) if (keep[j]) idx.push_back((int)j);
        saige::spa_gpu::Geno geno;
        geno.packed = (const unsigned char*)saige::gpu2::devicePacked(R); geno.bpv = bpv;
        geno.lut = (const double*)saige::gpu2::deviceLut(R);
        for (std::size_t k0 = 0; k0 < idx.size(); k0 += (std::size_t)maxP) {
            const int nk = (int)std::min<std::size_t>(idx.size() - k0, (std::size_t)maxP);
            for (int pr = 0; pr < 2; pr++) {
                saige::gpu2::SpaPairIn* in = saige::gpu2::spaIn(SO[pr]);
                saige::spa_gpu::PairIn* li = saige::spa_gpu::in(SL[pr]);
                for (int k = 0; k < nk; k++) {
                    const Rec& r = batch[idx[k0 + k]];
                    in[k].slot = r.m - (int)m0; in[k].trait = r.t; in[k].fast = r.fast; in[k].logp = r.logp;
                    in[k].Tstat = r.Tstat; in[k].var1 = r.var1; in[k].var2 = r.var2; in[k].pno = r.pno;
                    li[k].slot = in[k].slot; li[k].trait = r.t; li[k].fast = r.fast; li[k].logp = r.logp;
                    li[k].Tstat = r.Tstat; li[k].var1 = r.var1; li[k].var2 = r.var2; li[k].pno = r.pno;
                }
                if (!saige::gpu2::spaRun(SO[pr], R, nk)) { std::fprintf(stderr, "spaRun failed\n"); return 3; }
                if (!saige::spa_gpu::run(SL[pr], geno, nk)) { std::fprintf(stderr, "lib run failed: %s\n", saige::spa_gpu::lastError()); return 3; }
                const saige::gpu2::SpaPairOut* o = saige::gpu2::spaOut(SO[pr]);
                const saige::spa_gpu::PairOut* lo = saige::spa_gpu::out(SL[pr]);
                for (int k = 0; k < nk; k++) {
                    Rec& r = batch[idx[k0 + k]];
                    Res& L = r.g[0][pr];
                    L.r1 = lo[k].root1; L.r2 = lo[k].root2; L.n1 = lo[k].niter1; L.n2 = lo[k].niter2;
                    L.s1 = lo[k].s1; L.s2 = lo[k].s2; L.conv = lo[k].conv; L.p = lo[k].pval;
                    L.status = lo[k].status; L.rs1 = lo[k].reason1; L.rs2 = lo[k].reason2;
                    L.convF = postRules(L.p, r.logp != 0, L.conv != 0, L.str, r.noSpa);
                    Res& O = r.g[1][pr];
                    O.r1 = o[k].root1; O.r2 = o[k].root2; O.n1 = o[k].niter1; O.n2 = o[k].niter2;
                    O.s1 = o[k].s1; O.s2 = o[k].s2; O.conv = o[k].conv; O.p = o[k].pval;
                    O.status = 0; O.rs1 = O.rs2 = 0;
                    O.convF = postRules(O.p, r.logp != 0, O.conv != 0, O.str, r.noSpa);
                }
            }
        }
        for (int j : idx) { recs.push_back(batch[j]); if ((long)recs.size() >= maxPairs) break; }
        std::fprintf(stderr, "  markers %ld-%ld: %zu pairs (total kept %zu)\n", m0, m1, idx.size(), recs.size());
    }

    const char* implName[2] = {"lib", "own"};
    long bad = 0;
    for (int im = 0; im < 2; im++) {
        std::map<std::string, Cmp> bins; Cmp all, cpu64;
        Res cpuRes;
        for (const Rec& r : recs) {
            const Res& A = r.g[im][0]; const Res& B = r.g[im][1];
            tally(all, r, A, B);
            tally(bins[r.fast ? "variant: fast" : "variant: full N"], r, A, B);
            tally(bins[r.logp ? "domain: log-p" : "domain: linear"], r, A, B);
            const double z = std::fabs(r.z);
            const char* zb = z < 3 ? "|z| [2,3)" : z < 6 ? "|z| [3,6)" : z < 10 ? "|z| [6,10)" : z < 20 ? "|z| [10,20)" :
                             z < 37 ? "|z| [20,37)" : "|z| >= 37";
            tally(bins[std::string("z: ") + zb], r, A, B);
            if (r.m >= 0) { Res C; C.p = r.pc; C.conv = r.convc; C.convF = r.convFc; C.str = r.strc; C.r1 = A.r1; C.r2 = A.r2;
                            C.n1 = A.n1; C.n2 = A.n2; C.s1 = A.s1; C.s2 = A.s2; tally(cpu64, r, C, B); }
        }
        std::printf("== %s: fp32 vs fp64 (%zu pairs, N=%d, traits=%d%s)\n", implName[im], recs.size(), n, T, sweep.empty() ? "" : ", sweep");
        printCmp("all", all);
        for (auto& kv : bins) printCmp(kv.first.c_str(), kv.second);
        printCmp("fp32 vs CPU fp64", cpu64);
        bad += all.convF + all.pzero + all.pnonfin;
        double kt[2]; long long np[2];
        for (int pr = 0; pr < 2; pr++) {
            if (im == 0) saige::spa_gpu::timings(SL[pr], &kt[pr], nullptr, nullptr, &np[pr]);
            else saige::gpu2::spaTimings(SO[pr], &kt[pr], &np[pr]);
        }
        std::printf("  kernel: fp64 %.2f us/pair, fp32 %.2f us/pair (%lld pairs)\n", kt[0] / std::max(1LL, np[0]) * 1e6,
                    kt[1] / std::max(1LL, np[1]) * 1e6, np[0]);
        // worst pairs
        std::vector<std::pair<double, int>> w;
        for (std::size_t i = 0; i < recs.size(); i++) {
            const Rec& r = recs[i];
            double rel = prel(r.g[im][0].p, r.g[im][1].p, r.logp);
            if (r.g[im][0].convF != r.g[im][1].convF) rel = 1e300;
            w.push_back({rel, (int)i});
        }
        std::sort(w.begin(), w.end(), [](auto& a, auto& b) { return a.first > b.first; });
        for (int i = 0; i < nWorst && i < (int)w.size() && w[i].first > 0; i++) {
            const Rec& r = recs[w[i].second];
            std::printf("   worst m=%d t=%d z=%.3f fast=%d logp=%d pno=%.6g CPU p=%.10g conv=%d\n", r.m, r.t, r.z, r.fast, r.logp, r.pno, r.pc, r.convc);
            for (int pr = 0; pr < 2; pr++) {
                const Res& X = r.g[im][pr];
                std::printf("     %s p=%.10g conv=%d/%d s=%d,%d n=%d,%d root=%.10g,%.10g status=0x%x reason=%d,%d\n", pr ? "fp32" : "fp64",
                            X.p, X.conv, X.convF, X.s1, X.s2, X.n1, X.n2, X.r1, X.r2, X.status, X.rs1, X.rs2);
            }
        }
    }
    if (flipmu) std::printf("(--flipmu: SPA used mu' = 1 - mu)\n");
    std::printf("RESULT: %s (final convergence identical, no p == 0 or non-finite p introduced by fp32)\n", bad ? "FAIL" : "PASS");
    for (int pr = 0; pr < 2; pr++) { saige::gpu2::spaDestroy(SO[pr]); saige::spa_gpu::destroy(SL[pr]); }
    saige::gpu2::destroy(R);
    return bad ? 1 : 0;
}
