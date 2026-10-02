// spa_gpu_check -- the device saddlepoint kernel (gpu/gpu_spa.hpp) against the
// production CPU SPA, pair by pair, on identical inputs.
//
// Links the step-2 objects (everything but main.o) and the CUDA objects, so
// getadjGFast / scoreTestFast / getroot_K1_* / Get_Saddle_Prob_* / SPA /
// SPA_fast are the shipped functions, and the genotype column comes from the
// shipped PLINK fused decode (the same packed bytes and dosage table the GPU
// loop stages). For every selected (marker, trait) pair:
//
//   CPU   scoreTestFast -> Tstat / var1 / var2 / pval_noadj / islogp, then the
//         SPA block of getMarkerPval replayed with the roots exposed (the two
//         getroot calls, the two saddle probabilities, the SPA / SPA_fast
//         wrapper) and the host post-rules (quantile step, p == 0, the string)
//   GPU   the same Tstat / var1 / var2 / pval_noadj handed to spaRun(); the
//         kernel recomputes g~, m1, q, qinv, the bounds, the fast-variant
//         terms, the roots and the tails; then the same host post-rules
//
// and reports, by |z| bin and by variant: iteration-count mismatches, root
// differences, convergence / saddle-flag / final-convergence mismatches, the
// relative p-value difference, and whether the printed "%.6E" / "%.1fE%d"
// strings agree.
//
//   spa_gpu_check BED MODELDIR:VRFILE[,...] [--markers M] [--zcut 2] [--maxpairs K]
//                 [--sweep z1,z2,...] [--tsv OUT]
//
// --sweep replaces Tstat by +-z sqrt(var1) for every listed z on every
// (marker, trait) of the first --markers markers, recomputing pval_noadj from
// the chi-square tail exactly as format_score_result does -- which is how the
// log-p path (p below ~1e-300, |z| above ~37) and the out-of-support roots
// are reached on data that has no such pair.
#include <armadillo>
#include <omp.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "../genotype_reader.hpp"
#include "../null_model_loader.hpp"
#include "../saige_test.hpp"
#include "../score_format.hpp"
#include "../spa.hpp"
#include "../spa_binary.hpp"
#include "../UTIL.hpp"
#include "../gpu/gpu_step2.hpp"
#include "../gpu/gpu_spa.hpp"

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

struct Rec {
    int m, t; double z; int fast, logp;
    double Tstat, var1, var2, pno;
    // CPU
    double root1c, root2c; int n1c, n2c, s1c, s2c, convc; double pc; int convFc; std::string strc;
    // GPU
    double root1g, root2g; int n1g, n2g, s1g, s2g, convg; double pg; int convFg; std::string strg;
};

struct Bin { long n = 0, iter = 0, root = 0, sadd = 0, conv = 0, convF = 0, str = 0, prel = 0; double maxrel = 0, maxroot = 0; };

void tally(Bin& b, const Rec& r)
{
    b.n++;
    if (r.n1c != r.n1g || r.n2c != r.n2g) b.iter++;
    auto rd = [](double a, double c) { if (std::isinf(a) && std::isinf(c) && (a > 0) == (c > 0)) return 0.0;
                                       if (std::isinf(a) != std::isinf(c)) return std::numeric_limits<double>::infinity(); return std::fabs(a - c); };
    const double dr = std::max(rd(r.root1c, r.root1g), rd(r.root2c, r.root2g));
    if (dr > 1e-12) b.root++;
    if (dr > b.maxroot) b.maxroot = dr;
    if (r.s1c != r.s1g || r.s2c != r.s2g) b.sadd++;
    if (r.convc != r.convg) b.conv++;
    if (r.convFc != r.convFg) b.convF++;
    if (r.strc != r.strg) b.str++;
    double rel;
    if (std::isnan(r.pc) && std::isnan(r.pg)) rel = 0;
    else if (std::isinf(r.pc) && std::isinf(r.pg) && (r.pc > 0) == (r.pg > 0)) rel = 0;
    else rel = std::fabs(r.pc - r.pg) / std::max(std::fabs(r.pc), r.logp ? 1e-300 : 5e-324);
    if (rel > 1e-10) b.prel++;
    if (rel > b.maxrel) b.maxrel = rel;
}

void printBin(const char* name, const Bin& b)
{
    if (b.n == 0) { std::printf("  %-26s n=%8ld\n", name, b.n); return; }
    std::printf("  %-26s n=%8ld  iter-mismatch %6ld  |droot|>1e-12 %6ld (max %.1e)  saddle-flag %5ld  conv %5ld  conv_final %5ld  p rel>1e-10 %6ld (max %.1e)  string-differs %6ld\n",
                name, b.n, b.iter, b.root, b.maxroot, b.sadd, b.conv, b.convF, b.prel, b.maxrel, b.str);
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: spa_gpu_check BED MODELDIR:VRFILE[,...] [--markers M] [--zcut Z] [--maxpairs K] [--sweep z,..] [--tsv OUT]\n"); return 2; }
    const std::string bed = argv[1];
    long nMarkers = -1; double zcut = 2.0; long maxPairs = 20000; std::vector<double> sweep; std::string tsv;
    for (int i = 3; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--markers" && i + 1 < argc) nMarkers = std::atol(argv[++i]);
        else if (a == "--zcut" && i + 1 < argc) zcut = std::atof(argv[++i]);
        else if (a == "--maxpairs" && i + 1 < argc) maxPairs = std::atol(argv[++i]);
        else if (a == "--sweep" && i + 1 < argc) for (auto& s : split(argv[++i], ',')) sweep.push_back(std::atof(s.c_str()));
        else if (a == "--tsv" && i + 1 < argc) tsv = argv[++i];
        else { std::fprintf(stderr, "bad arg %s\n", a.c_str()); return 2; }
    }
    std::vector<Model*> models;
    for (auto& s : split(argv[2], ',')) models.push_back(loadModel(s));
    const int T = (int)models.size();
    std::vector<std::string> ids = models[0]->nm.sampleIDs;
    setPLINKobjInCPP(bed + ".bim", bed + ".fam", bed + ".bed", ids, "alt-first");
    const int n = (int)ptr_gPLINKobj->getN();
    const long M = (long)ptr_gPLINKobj->getM();
    if (nMarkers < 0 || nMarkers > M) nMarkers = M;
    for (auto* Mo : models) if (Mo->obj->m_n != n) { std::fprintf(stderr, "model n != reader n\n"); return 2; }
    const int imputeCase = string_to_case.at(models[0]->nm.impute_method);
    std::fprintf(stderr, "N=%d markers=%ld (using %ld) traits=%d sweep=%zu\n", n, M, nMarkers, T, sweep.size());

    // ---- device: reducer with a dummy 1-column right operand, SPA with the traits ----
    const int slotsCap = 4096;
    std::vector<double> B1((std::size_t)n, 0.0);
    saige::gpu2::CreateArgs ca; ca.device = 0; ca.N = n; ca.K1 = 1; ca.B1 = B1.data(); ca.maxSlots = slotsCap; ca.fp64 = true;
    saige::gpu2::Reducer* R = saige::gpu2::create(ca);
    if (!R) { std::fprintf(stderr, "reducer create failed\n"); return 3; }
    std::vector<saige::gpu2::SpaTraitArgs> ta(T);
    for (int t = 0; t < T; t++) { ta[t].mu = models[t]->mu.memptr(); ta[t].XV = models[t]->obj->m_XV.memptr(); ta[t].XXVX_inv = models[t]->obj->m_XXVX_inv.memptr(); ta[t].p = models[t]->obj->m_p; }
    const double tol = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
    saige::gpu2::SpaCreateArgs sa; sa.device = 0; sa.N = n; sa.nTraits = T; sa.traits = ta.data();
    sa.maxPairs = (int)std::min<long>(maxPairs, (long)slotsCap * T * (long)std::max<size_t>(1, 2 * sweep.size())); sa.tol = tol; sa.maxiter = 1000; sa.blocks = 256;
    saige::gpu2::Spa* SP = saige::gpu2::spaCreate(sa);
    if (!SP) { std::fprintf(stderr, "spa create failed\n"); return 3; }
    const std::size_t bpv = saige::gpu2::bytesPerSlot(R);
    unsigned char* hPk = saige::gpu2::packed(R);
    double* hLu = saige::gpu2::lut(R);

    std::vector<Rec> recs;
    std::vector<double> gpuUs;
    long pairsTotal = 0;
    std::FILE* ft = tsv.empty() ? nullptr : std::fopen(tsv.c_str(), "w");
    if (ft) std::fprintf(ft, "m\tt\tz\tfast\tlogp\tTstat\tvar1\tvar2\tpno\tr1c\tr1g\tr2c\tr2g\tn1c\tn1g\tn2c\tn2g\tconvc\tconvg\tpc\tpg\tconvFc\tconvFg\tstrc\tstrg\n");

    // ---- batches of markers ----
    for (long m0 = 0; m0 < nMarkers && (long)recs.size() < maxPairs; m0 += slotsCap) {
        const long m1 = std::min(nMarkers, m0 + (long)slotsCap);
        const int nS = (int)(m1 - m0);
        // stage + CPU per marker
        std::vector<arma::vec> gs(nS); std::vector<arma::uvec> izs(nS), inzs(nS);
        std::vector<double> af(nS), mac(nS), fds((std::size_t)nS * 4); std::vector<char> okm(nS, 0);
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
            for (int k = 0; k < 4; k++) { fds[(std::size_t)s * 4 + k] = fs.fd[k]; hLu[(std::size_t)s * 4 + k] = fs.fd[k]; }
            af[s] = fs.altFreq_post; okm[s] = 1;
            mac[s] = std::min(fs.altCounts_post, 2.0 * n - fs.altCounts_post);
        }
        for (int s = 0; s < nS; s++) if (!okm[s]) { std::memset(hPk + (std::size_t)s * bpv, 0, (n + 3) / 4); for (int k = 0; k < 4; k++) hLu[(std::size_t)s * 4 + k] = 0; }
        if (!saige::gpu2::reduce(R, nS)) { std::fprintf(stderr, "reduce failed\n"); return 3; }

        // CPU: score test and SPA replica per (marker, trait[, z])
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
            r.Tstat = Tstat; r.var1 = var1; r.var2 = var2; r.pno = pval;
            // --- the SPA block of getMarkerPval ---
            arma::vec gt(n); obj->getadjGFast(g, gt, inz);
            const double m1 = arma::dot(Mo->mu, gt);
            const double pz = double(iz.n_elem) / n;
            const bool fast = (pz >= 0.5);
            r.fast = fast;
            double q = Tstat / std::sqrt(var1 / var2) + m1, qinv;
            if ((q - m1) > 0) qinv = -1 * std::fabs(q - m1) + m1; else if ((q - m1) == 0) qinv = m1; else qinv = std::fabs(q - m1) + m1;
            arma::vec gNB, gNA, muNB, muNA; double NAmu = 0, NAsigma = 0;
            if (fast) { gNB = gt(inz); gNA = gt(iz); muNB = Mo->mu(inz); muNA = Mo->mu(iz); NAmu = m1 - arma::dot(gNB, muNB); NAsigma = var2 - arma::sum(muNB % (1 - muNB) % arma::pow(gNB, 2)); }
            RootResult r1, r2;
            if (fast) { r1 = getroot_K1_fast_Binom(0, Mo->mu, gt, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol); r2 = getroot_K1_fast_Binom(0, Mo->mu, gt, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol); }
            else { r1 = getroot_K1_Binom(0, Mo->mu, gt, q, tol); r2 = getroot_K1_Binom(0, Mo->mu, gt, qinv, tol); }
            r.root1c = r1.root; r.root2c = r2.root; r.n1c = r1.niter; r.n2c = r2.niter; r.s1c = -1; r.s2c = -1;
            if (r1.Isconverge && r2.Isconverge) {
                SaddleResult a, b;
                if (fast) { a = Get_Saddle_Prob_fast_Binom(r1.root, Mo->mu, gt, q, gNA, gNB, muNA, muNB, NAmu, NAsigma, islogp); b = Get_Saddle_Prob_fast_Binom(r2.root, Mo->mu, gt, qinv, gNA, gNB, muNA, muNB, NAmu, NAsigma, islogp); }
                else { a = Get_Saddle_Prob_Binom(r1.root, Mo->mu, gt, q, islogp); b = Get_Saddle_Prob_Binom(r2.root, Mo->mu, gt, qinv, islogp); }
                r.s1c = a.isSaddle; r.s2c = b.isSaddle;
            }
            double pw; bool cw;
            if (fast) SPA_fast(Mo->mu, gt, q, qinv, pval, islogp, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol, "binary", pw, cw);
            else SPA(Mo->mu, gt, q, qinv, pval, tol, islogp, "binary", pw, cw);
            r.pc = pw; r.convc = cw;
            r.convFc = postRules(pw, islogp, cw, r.strc, pstr);
            r.strg = pstr;   // placeholder, replaced below
        }
        // GPU
        std::vector<int> idx;
        for (std::size_t j = 0; j < jobs.size(); j++) if (keep[j]) idx.push_back((int)j);
        saige::gpu2::SpaPairIn* in = saige::gpu2::spaIn(SP);
        for (std::size_t k0 = 0; k0 < idx.size(); k0 += (std::size_t)sa.maxPairs) {
            const int nk = (int)std::min<std::size_t>(idx.size() - k0, (std::size_t)sa.maxPairs);
            for (int k = 0; k < nk; k++) {
                const Rec& r = batch[idx[k0 + k]];
                in[k].slot = r.m - (int)m0; in[k].trait = r.t; in[k].fast = r.fast; in[k].logp = r.logp;
                in[k].Tstat = r.Tstat; in[k].var1 = r.var1; in[k].var2 = r.var2; in[k].pno = r.pno;
            }
            const double t0 = omp_get_wtime();
            if (!saige::gpu2::spaRun(SP, R, nk)) { std::fprintf(stderr, "spaRun failed\n"); return 3; }
            gpuUs.push_back((omp_get_wtime() - t0) / nk * 1e6);
            const saige::gpu2::SpaPairOut* o = saige::gpu2::spaOut(SP);
            for (int k = 0; k < nk; k++) {
                Rec& r = batch[idx[k0 + k]];
                r.root1g = o[k].root1; r.root2g = o[k].root2; r.n1g = o[k].niter1; r.n2g = o[k].niter2;
                r.s1g = o[k].s1; r.s2g = o[k].s2; r.convg = o[k].conv; r.pg = o[k].pval;
                std::string noSpa = r.strg;
                r.convFg = postRules(r.pg, r.logp != 0, r.convg != 0, r.strg, noSpa);
            }
        }
        for (int j : idx) { recs.push_back(batch[j]); if ((long)recs.size() >= maxPairs) break; }
        pairsTotal += (long)idx.size();
        std::fprintf(stderr, "  markers %ld-%ld: %zu pairs (total kept %zu)\n", m0, m1, idx.size(), recs.size());
    }

    // ---- report ----
    std::map<std::string, Bin> bins;
    Bin all;
    for (const Rec& r : recs) {
        tally(all, r);
        tally(bins[r.fast ? "variant: fast (carriers)" : "variant: full N"], r);
        tally(bins[r.logp ? "domain: log-p" : "domain: linear"], r);
        const double z = std::fabs(r.z);
        const char* zb = z < 3 ? "|z| in [2,3)" : z < 4 ? "|z| in [3,4)" : z < 6 ? "|z| in [4,6)" : z < 10 ? "|z| in [6,10)" : z < 20 ? "|z| in [10,20)" :
                         z < 30 ? "|z| in [20,30)" : z < 37 ? "|z| in [30,37)" : z < 40 ? "|z| in [37,40)" : "|z| >= 40";
        tally(bins[std::string("z: ") + zb], r);
        if (std::isinf(r.root1c) || std::isinf(r.root2c)) tally(bins["root: out of support (inf)"], r);
        if (ft) std::fprintf(ft, "%d\t%d\t%.3f\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%d\t%d\t%s\t%s\n",
            r.m, r.t, r.z, r.fast, r.logp, r.Tstat, r.var1, r.var2, r.pno, r.root1c, r.root1g, r.root2c, r.root2g, r.n1c, r.n1g, r.n2c, r.n2g, r.convc, r.convg, r.pc, r.pg, r.convFc, r.convFg, r.strc.c_str(), r.strg.c_str());
    }
    if (ft) std::fclose(ft);
    double us = 0; for (double v : gpuUs) us += v; if (!gpuUs.empty()) us /= gpuUs.size();
    std::printf("spa_gpu_check: %zu pairs (N=%d, traits=%d, zcut %.2f%s), GPU %.1f us/pair incl. H2D/D2H\n", recs.size(), n, T, zcut, sweep.empty() ? "" : ", sweep", us);
    printBin("all", all);
    for (auto& kv : bins) printBin(kv.first.c_str(), kv.second);
    const bool pass = (all.iter == 0 && all.sadd == 0 && all.conv == 0 && all.convF == 0 && all.prel == 0);
    std::printf("RESULT: %s (iteration counts, saddle flags, convergence, final convergence all identical; p rel <= 1e-10)\n", pass ? "PASS" : "FAIL");
    saige::gpu2::spaDestroy(SP);
    saige::gpu2::destroy(R);
    return pass ? 0 : 1;
}
