// main.cpp -- step1bench: a SAIGE step-1 null-model benchmark whose only
// variable is how Sigma is inverted.
//
// Usage is in README.md. The order of operations is deliberate:
//   1. read the GRM and the design
//   2. run the Sigma self-check (all solvers must form the SAME Sigma)
//   3. build the chosen solver, run AI-REML, report
// The self-check comes before any timer starts, and uses throwaway solver
// instances, so it cannot flatter or penalise the measured run.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "aireml.hpp"
#include "bench_io.hpp"
#include "sigma_solver.hpp"

using namespace step1bench;

namespace {

struct Args {
    std::string grm, pheno, trait = "quant", phenoCol, covarCols, phenoIdCol = "IID";
    std::string solver = "block", trace = "hutchinson", report;
    int probes = 30, maxiter = 20, threads = 1;
    double tol = 0.02;
    unsigned long long seed = 1;
    bool selfCheckOnly = false;
};

[[noreturn]] void usage(int code) {
    std::fprintf(code ? stderr : stdout,
"step1bench -- SAIGE step-1 null model, three interchangeable Sigma solvers\n"
"\n"
"  --grm <prefix>        sparse GRM: <prefix>.mtx + <prefix>.ids   (required)\n"
"  --pheno <file>        phenotype table, header line              (required)\n"
"  --pheno-col <name>    the trait column                          (required)\n"
"  --trait quant|binary  default quant\n"
"  --covar-cols a,b,c    covariates; an intercept is always added\n"
"  --pheno-id-col <name> id column in --pheno (default IID)\n"
"  --solver block|cholmod|cholmod_reuse|superlu   default block\n"
"  --trace exact|hutchinson                       default hutchinson\n"
"  --probes N            Hutchinson probes per iteration (default 30)\n"
"  --maxiter N           AI-REML iterations (default 20)\n"
"  --tol X               relative change on tau (default 0.02, SAIGE's)\n"
"  --seed N              probe RNG seed (default 1)\n"
"  --threads N           OpenMP threads (default 1)\n"
"  --report <file.json>  write the machine-readable report too\n"
"  --self-check-only     run the Sigma self-check and exit\n");
    std::exit(code);
}

std::string need(int& i, int argc, char** argv) {
    if (i + 1 >= argc) { std::fprintf(stderr, "step1bench: %s needs a value\n", argv[i]); usage(2); }
    return std::string(argv[++i]);
}

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        const std::string k = argv[i];
        if      (k == "--grm")            a.grm = need(i, argc, argv);
        else if (k == "--pheno")          a.pheno = need(i, argc, argv);
        else if (k == "--pheno-col")      a.phenoCol = need(i, argc, argv);
        else if (k == "--pheno-id-col")   a.phenoIdCol = need(i, argc, argv);
        else if (k == "--trait")          a.trait = need(i, argc, argv);
        else if (k == "--covar-cols")     a.covarCols = need(i, argc, argv);
        else if (k == "--solver")         a.solver = need(i, argc, argv);
        else if (k == "--trace")          a.trace = need(i, argc, argv);
        else if (k == "--probes")         a.probes = std::stoi(need(i, argc, argv));
        else if (k == "--maxiter")        a.maxiter = std::stoi(need(i, argc, argv));
        else if (k == "--tol")            a.tol = std::stod(need(i, argc, argv));
        else if (k == "--seed")           a.seed = std::stoull(need(i, argc, argv));
        else if (k == "--threads")        a.threads = std::stoi(need(i, argc, argv));
        else if (k == "--report")         a.report = need(i, argc, argv);
        else if (k == "--self-check-only") a.selfCheckOnly = true;
        else if (k == "-h" || k == "--help") usage(0);
        else { std::fprintf(stderr, "step1bench: unknown option %s\n", argv[i]); usage(2); }
    }
    if (a.grm.empty() || a.pheno.empty() || a.phenoCol.empty()) {
        std::fprintf(stderr, "step1bench: --grm, --pheno and --pheno-col are required\n");
        usage(2);
    }
    if (a.trait != "quant" && a.trait != "binary") {
        std::fprintf(stderr, "step1bench: --trait must be quant or binary\n"); std::exit(2);
    }
    if (a.trace != "exact" && a.trace != "hutchinson") {
        std::fprintf(stderr, "step1bench: --trace must be exact or hutchinson\n"); std::exit(2);
    }
    if (a.probes < 1)  { std::fprintf(stderr, "step1bench: --probes must be >= 1\n");  std::exit(2); }
    if (a.maxiter < 1) { std::fprintf(stderr, "step1bench: --maxiter must be >= 1\n"); std::exit(2); }
    return a;
}

double max_abs_diff(const arma::sp_mat& A, const arma::sp_mat& B) {
    arma::sp_mat D = A - B;
    double m = 0.0;
    for (arma::sp_mat::const_iterator it = D.begin(); it != D.end(); ++it)
        m = std::max(m, std::fabs(*it));
    return m;
}

// The fairness gate. All solvers must form the same Sigma; if any of them does
// not, every number downstream is meaningless, so this runs on every
// invocation and refuses to continue.
struct SelfCheck {
    bool ok = false;
    double maxDiff = 0.0;
    std::map<std::string, double> residual;    // solver -> max|Sigma x - b| / max|b|
};

SelfCheck run_self_check(const SparseGRM& g, unsigned long long seed) {
    SelfCheck sc;
    const arma::uword n = (arma::uword)g.n;

    std::vector<std::string> names = {"block", "cholmod", "superlu"};
    std::vector<std::unique_ptr<SigmaSolver>> sv;
    for (const auto& nm : names) {
        auto s = make_solver(nm);
        if (!s->build(g.loc, g.val, g.n)) {
            std::fprintf(stderr, "step1bench: self-check: %s failed to build\n", nm.c_str());
            return sc;
        }
        sv.push_back(std::move(s));
    }

    // Two probe points: W == 1 (the quantitative case) and a small, uneven W
    // (the binary case, where tau0/w_i is large and the fp32 rounding of that
    // product is what the arithmetic has to agree on).
    std::vector<std::pair<arma::fvec, arma::fvec>> pts;
    {
        arma::fvec w1(n, arma::fill::ones), t1(2);
        t1(0) = 0.712932f; t1(1) = 0.295917f;
        pts.emplace_back(w1, t1);

        arma::arma_rng::set_seed(seed);
        arma::fvec w2 = 0.02f + 0.23f * arma::randu<arma::fvec>(n);
        arma::fvec t2(2); t2(0) = 1.0f; t2(1) = 0.210022f;
        pts.emplace_back(w2, t2);
    }

    double worst = 0.0;
    for (const auto& p : pts) {
        std::vector<arma::sp_mat> S;
        for (auto& s : sv) S.push_back(s->debugSigma(p.first, p.second));
        for (size_t i = 0; i < S.size(); i++) {
            if (S[i].n_nonzero != S[0].n_nonzero) {
                std::fprintf(stderr,
                    "step1bench: self-check: %s formed Sigma with %llu nonzeros, "
                    "%s with %llu\n", names[i].c_str(),
                    (unsigned long long)S[i].n_nonzero, names[0].c_str(),
                    (unsigned long long)S[0].n_nonzero);
            }
            worst = std::max(worst, max_abs_diff(S[i], S[0]));
        }
    }
    sc.maxDiff = worst;
    sc.ok = (worst == 0.0);

    // Numeric half: the reconstruction above proves the three form the same
    // matrix; this proves each solver's machinery actually inverts it.
    {
        const arma::fvec& w = pts[1].first;
        const arma::fvec& t = pts[1].second;
        arma::sp_mat Sref = sv[2]->debugSigma(w, t);
        arma::arma_rng::set_seed(seed + 7);
        arma::fvec b = arma::randn<arma::fvec>(n);
        const double bmax = arma::max(arma::abs(b));
        for (size_t i = 0; i < sv.size(); i++) {
            sv[i]->refresh(w, t);
            arma::fvec x = sv[i]->solve(b);
            arma::vec r = Sref * arma::conv_to<arma::vec>::from(x)
                        - arma::conv_to<arma::vec>::from(b);
            sc.residual[names[i]] = arma::max(arma::abs(r)) / bmax;
        }
    }
    return sc;
}

std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) { if (c == '"' || c == '\\') o.push_back('\\'); o.push_back(c); }
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    const double t_start = now_s();
    Args a = parse(argc, argv);

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, a.threads));
#endif
    setenv("OMP_NUM_THREADS", std::to_string(std::max(1, a.threads)).c_str(), 1);

    try {
        std::printf("step1bench\n");
        std::printf("  grm       %s\n", a.grm.c_str());
        std::printf("  pheno     %s  [%s]\n", a.pheno.c_str(), a.phenoCol.c_str());
        std::printf("  trait     %s\n", a.trait.c_str());
        std::printf("  solver    %s\n", a.solver.c_str());
        std::printf("  trace     %s%s\n", a.trace.c_str(),
                    a.trace == "hutchinson" ? "" : " (no probes)");
        if (a.trace == "hutchinson") std::printf("  probes    %d\n", a.probes);
        std::printf("  threads   %d\n", std::max(1, a.threads));
        std::fflush(stdout);

        const double t_read0 = now_s();
        SparseGRM g = read_sparse_grm(a.grm);
        PhenoTable ph = read_pheno(a.pheno, a.phenoIdCol);
        arma::fvec y; arma::fmat X;
        assemble_design(g, ph, a.phenoCol, split_commas(a.covarCols), y, X);
        const double t_read = now_s() - t_read0;

        const BlockSummary bsum = summarize_blocks(g.loc, g.n);
        std::printf("\n  N = %d, nnz(Psi) = %llu (both triangles), X = %llu x %llu\n",
                    g.n, (unsigned long long)g.val.n_elem,
                    (unsigned long long)X.n_rows, (unsigned long long)X.n_cols);
        std::printf("  connected components: %d blocks, max %d\n",
                    bsum.nblocks, bsum.maxBlock);
        std::printf("  block-size histogram:");
        for (auto& kv : bsum.hist) std::printf(" %d:%lld", kv.first, kv.second);
        std::printf("\n  read %.3fs\n", t_read);

        // --- gate 1 -------------------------------------------------------
        std::printf("\n[self-check] forming Sigma three ways at two (w,tau) points...\n");
        std::fflush(stdout);
        SelfCheck sc = run_self_check(g, a.seed);
        std::printf("[self-check] max|Sigma_block - Sigma_ref| over block/cholmod/superlu "
                    "= %.17g  %s\n", sc.maxDiff, sc.ok ? "PASS" : "FAIL");
        for (auto& kv : sc.residual)
            std::printf("[self-check] %-14s max|Sigma x - b| / max|b| = %.3e\n",
                        kv.first.c_str(), kv.second);
        std::fflush(stdout);
        if (!sc.ok) {
            std::fprintf(stderr, "step1bench: solvers do not form the same Sigma; "
                                 "refusing to benchmark.\n");
            return 3;
        }
        if (a.selfCheckOnly) return 0;

        if (a.trace == "exact" && a.solver == "superlu") {
            std::fprintf(stderr,
                "step1bench: --trace exact is not available with --solver superlu.\n"
                "  arma::spsolve never materialises Sigma^-1, so there is no entry to\n"
                "  contract against Psi. Use --trace hutchinson (that combination IS the\n"
                "  SAIGE baseline).\n");
            return 2;
        }

        // --- the measured run ---------------------------------------------
        auto raw = make_solver(a.solver);
        if (!raw) { std::fprintf(stderr, "step1bench: unknown solver %s\n", a.solver.c_str()); return 2; }
        TimedSolver sol(std::move(raw));

        const double t_fit0 = now_s();
        if (!sol.build(g.loc, g.val, g.n)) {
            std::fprintf(stderr, "step1bench: solver %s refused to build on this GRM\n",
                         a.solver.c_str());
            return 4;
        }
        if (a.trace == "exact" && !sol.hasExplicitInverse()) {
            std::fprintf(stderr, "step1bench: solver %s has no explicit inverse\n",
                         a.solver.c_str());
            return 2;
        }

        arma::sp_mat psi(g.loc, g.val, (arma::uword)g.n, (arma::uword)g.n);

        FitConfig cfg;
        cfg.binary     = (a.trait == "binary");
        cfg.exactTrace = (a.trace == "exact");
        cfg.probes     = a.probes;
        cfg.maxiter    = a.maxiter;
        cfg.tol        = a.tol;
        cfg.seed       = a.seed;

        std::printf("\n[fit] %s / %s / %s\n", a.solver.c_str(), a.trait.c_str(), a.trace.c_str());
        std::fflush(stdout);
        FitResult fr = run_aireml(sol, psi, y, X, cfg);
        const double t_fit = now_s() - t_fit0;

        // --- report --------------------------------------------------------
        const StageAcc& acc = sol.acc();
        static const char* stageName[ST_N] = {"build", "refresh", "solve", "trace"};

        std::printf("\n  iter  tau0        tau1        score0        score1        "
                    "trPA0         trPA1         relchg     inner\n");
        for (const auto& r : fr.traj)
            std::printf("  %4d  %-11.6f %-11.6f %-13.5g %-13.5g %-13.5g %-13.5g %-10.3e %d\n",
                        (int)r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], (int)r[8]);

        std::printf("\n  tau        = [%.6f, %.6f]\n", (double)fr.tau(0), (double)fr.tau(1));
        std::printf("  alpha      = [");
        for (arma::uword k = 0; k < fr.alpha.n_elem; k++)
            std::printf("%s%.6f", k ? ", " : "", (double)fr.alpha(k));
        std::printf("]\n");
        std::printf("  iterations = %d%s, inner IRLS solves of (Y,X) = %lld\n",
                    fr.iters, fr.converged ? " (converged)" : " (NOT converged)",
                    fr.innerIRLS);

        std::printf("\n  stage      calls        seconds\n");
        double tot = 0.0;
        for (int s = 0; s < ST_N; s++) {
            std::printf("  %-9s  %-11lld  %.4f\n", stageName[s], acc.calls[s], acc.secs[s]);
            tot += acc.secs[s];
        }
        std::printf("  %-9s  %-11s  %.4f\n", "stages", "", tot);
        std::printf("  %-9s  %-11s  %.4f\n", "fit wall", "", t_fit);
        std::printf("  %-9s  %-11s  %.4f\n", "total", "", now_s() - t_start);

        const long long nref = sol.refreshCount(), nreuse = sol.reuseCount();
        if (nref + nreuse > 0)
            std::printf("\n  (w,tau) cache: %lld rebuilds, %lld reuses (%.1f%% hit)\n",
                        nref, nreuse, 100.0 * (double)nreuse / (double)(nref + nreuse));
        else
            std::printf("\n  (w,tau) cache: none (this solver does not cache)\n");
        const std::string extra = sol.extraReport();
        if (!extra.empty()) std::printf("  %s\n", extra.c_str());

        if (!a.report.empty()) {
            std::ofstream f(a.report);
            f.setf(std::ios::fixed);
            f << "{\n";
            f << "  \"grm\": \"" << json_escape(a.grm) << "\",\n";
            f << "  \"pheno\": \"" << json_escape(a.pheno) << "\",\n";
            f << "  \"pheno_col\": \"" << json_escape(a.phenoCol) << "\",\n";
            f << "  \"trait\": \"" << a.trait << "\",\n";
            f << "  \"solver\": \"" << a.solver << "\",\n";
            f << "  \"trace\": \"" << a.trace << "\",\n";
            f << "  \"probes\": " << a.probes << ",\n";
            f << "  \"threads\": " << std::max(1, a.threads) << ",\n";
            f << "  \"seed\": " << a.seed << ",\n";
            f << "  \"tol\": " << std::scientific << a.tol << std::fixed << ",\n";
            f << "  \"n\": " << g.n << ",\n";
            f << "  \"nnz\": " << g.val.n_elem << ",\n";
            f << "  \"nblocks\": " << bsum.nblocks << ",\n";
            f << "  \"max_block\": " << bsum.maxBlock << ",\n";
            f << "  \"block_hist\": {";
            for (size_t k = 0; k < bsum.hist.size(); k++)
                f << (k ? ", " : "") << "\"" << bsum.hist[k].first << "\": " << bsum.hist[k].second;
            f << "},\n";
            f << "  \"selfcheck_max_sigma_diff\": " << std::scientific << sc.maxDiff << ",\n";
            f << "  \"selfcheck_residual\": {";
            { bool first = true;
              for (auto& kv : sc.residual) { f << (first ? "" : ", ") << "\"" << kv.first
                                               << "\": " << kv.second; first = false; } }
            f << "},\n" << std::fixed;
            f << "  \"tau\": [" << std::setprecision(9) << (double)fr.tau(0) << ", "
              << (double)fr.tau(1) << "],\n";
            f << "  \"alpha\": [";
            for (arma::uword k = 0; k < fr.alpha.n_elem; k++)
                f << (k ? ", " : "") << (double)fr.alpha(k);
            f << "],\n" << std::setprecision(6);
            f << "  \"iterations\": " << fr.iters << ",\n";
            f << "  \"converged\": " << (fr.converged ? "true" : "false") << ",\n";
            f << "  \"stages\": {\n";
            for (int s = 0; s < ST_N; s++)
                f << "    \"" << stageName[s] << "\": {\"calls\": " << acc.calls[s]
                  << ", \"seconds\": " << acc.secs[s] << "}"
                  << (s + 1 < ST_N ? ",\n" : "\n");
            f << "  },\n";
            f << "  \"fit_seconds\": " << t_fit << ",\n";
            f << "  \"total_seconds\": " << (now_s() - t_start) << ",\n";
            f << "  \"cache_rebuilds\": " << nref << ",\n";
            f << "  \"cache_reuses\": " << nreuse << ",\n";
            f << "  \"extra\": \"" << json_escape(extra) << "\"\n";
            f << "}\n";
            std::printf("  report -> %s\n", a.report.c_str());
        }
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "step1bench: %s\n", e.what());
        return 1;
    }
}
