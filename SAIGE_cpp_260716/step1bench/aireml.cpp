#include "aireml.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>

namespace step1bench {

namespace {

// ---------------------------------------------------------------------------
// Every Sigma^-1 v in this file goes through here, so all three solvers see
// the same call pattern: refresh before EVERY solve, exactly as SAIGE's
// gen_spsolve_v4 calls blocksigma::refresh on each of its 199-261 calls.
// Solvers that cache on (w, tau) return from refresh immediately and report
// the hit; SuperLU has no cache and pays the rebuild inside solve().
arma::fvec sigma_solve(SigmaSolver& sol, const arma::fvec& w,
                       const arma::fvec& tau, const arma::fvec& b) {
    sol.refresh(w, tau);
    return sol.solve(b);
}

// Psi * v, in fp64, rounded once on return. Identical for every solver.
arma::fvec psi_mul(const arma::sp_mat& psi, const arma::fvec& v) {
    arma::vec vd = arma::conv_to<arma::vec>::from(v);
    arma::vec r  = psi * vd;
    return arma::conv_to<arma::fvec>::from(r);
}

arma::fmat psi_mul(const arma::sp_mat& psi, const arma::fmat& V) {
    arma::mat Vd = arma::conv_to<arma::mat>::from(V);
    arma::mat R  = psi * Vd;
    return arma::conv_to<arma::fmat>::from(R);
}

arma::fmat inv_psd_or_pinv(const arma::fmat& A) {
    try { return arma::inv_sympd(arma::symmatu(A)); }
    catch (const std::exception&) { return arma::pinv(arma::symmatu(A)); }
}

// ---------------------------------------------------------------------------
// The covariate correction, written out rather than folded into a helper that
// hides it.
//
//   P = Sigma^-1 - Sigma^-1 X (X' Sigma^-1 X)^-1 X' Sigma^-1
//
// so for any v, with Sigma_iX = Sigma^-1 X and cov = (X' Sigma_iX)^-1,
//
//   P v = Sigma^-1 v  -  Sigma_iX * ( cov * ( Sigma_iX' * v ) )
//
// The second term is the whole of the correction. Nothing else in this file
// touches X.
arma::fvec apply_P(const arma::fvec& Sigma_iv, const arma::fmat& Sigma_iX,
                   const arma::fmat& cov, const arma::fvec& v) {
    const arma::fvec Xtv     = Sigma_iX.t() * v;     // p
    const arma::fvec covXtv  = cov * Xtv;            // p
    const arma::fvec correct = Sigma_iX * covXtv;    // n
    return Sigma_iv - correct;
}

struct Coef {
    arma::fvec Sigma_iY;
    arma::fmat Sigma_iX;
    arma::fmat cov;
    arma::fvec alpha;
    arma::fvec eta;
};

Coef get_coefficients(SigmaSolver& sol, const arma::fvec& Y, const arma::fmat& X,
                      const arma::fvec& w, const arma::fvec& tau) {
    Coef c;
    c.Sigma_iY = sigma_solve(sol, w, tau, Y);
    c.Sigma_iX.set_size(X.n_rows, X.n_cols);
    for (arma::uword j = 0; j < X.n_cols; j++)
        c.Sigma_iX.col(j) = sigma_solve(sol, w, tau, arma::fvec(X.col(j)));
    c.cov   = inv_psd_or_pinv(X.t() * c.Sigma_iX);
    c.alpha = c.cov * (c.Sigma_iX.t() * Y);
    // eta = Y - tau0 * (Sigma^-1 Y - Sigma^-1 X alpha) / w, the BLUP identity
    // SAIGE uses (saige_ai.cpp::finishCoefficients_cpp).
    c.eta   = Y - tau(0) * (c.Sigma_iY - c.Sigma_iX * c.alpha) / w;
    return c;
}

// ---------------------------------------------------------------------------
// tr(P A0) and tr(P A1) with A0 = dSigma/dtau0 = diag(1/w) and A1 = Psi.
struct Traces { double trPA0 = 0.0, trPA1 = 0.0; };

// Exact, from an explicit Sigma^-1:
//   tr(P Psi) = tr(Sigma^-1 Psi) - tr( cov * Sigma_iX' * (Psi Sigma_iX) )
//   tr(P)     = tr(Sigma^-1)     - tr( cov * Sigma_iX' * Sigma_iX )
// The second identity is tr(P A0) only when w == 1 (quantitative here); for a
// general w it would need sum_i (Sigma^-1)_ii / w_i, which SigmaSolver::traces
// does not expose. Binary never asks for it (tau0 is fixed), and the caller
// checks.
Traces exact_traces(SigmaSolver& sol, const arma::sp_mat& psi,
                    const arma::fmat& Sigma_iX, const arma::fmat& cov,
                    bool wantTrPA0) {
    double trSiPsi = 0.0, trSi = 0.0;
    sol.traces(&trSiPsi, &trSi);

    const arma::mat SX  = arma::conv_to<arma::mat>::from(Sigma_iX);
    const arma::mat C   = arma::conv_to<arma::mat>::from(cov);
    const arma::mat PsiX = psi * SX;                        // n x p, fp64
    Traces t;
    t.trPA1 = trSiPsi - arma::trace(C * (SX.t() * PsiX));
    if (wantTrPA0) t.trPA0 = trSi - arma::trace(C * (SX.t() * SX));
    return t;
}

// Hutchinson, the SAIGE baseline. `probes` fixed -- no CV growth loop, so the
// work is identical for every solver.
Traces hutchinson_traces(SigmaSolver& sol, const arma::sp_mat& psi,
                         const arma::fvec& w, const arma::fvec& tau,
                         const arma::fmat& Sigma_iX, const arma::fmat& cov,
                         int probes, unsigned long long seed, bool wantTrPA0) {
    const arma::uword n = w.n_elem;
    std::mt19937_64 rng(seed);
    double s0 = 0.0, s1 = 0.0;
    arma::fvec u((arma::uword)n);
    for (int k = 0; k < probes; k++) {
        for (arma::uword i = 0; i < n; i++) u(i) = (rng() & 1ULL) ? 1.0f : -1.0f;
        const arma::fvec Siu = sigma_solve(sol, w, tau, u);
        const arma::fvec Pu  = apply_P(Siu, Sigma_iX, cov, u);
        const arma::fvec A1u = psi_mul(psi, u);
        s1 += (double)arma::dot(A1u, Pu);
        if (wantTrPA0) {
            const arma::fvec A0u = u / w;                  // diag(1/w) u
            s0 += (double)arma::dot(A0u, Pu);
        }
    }
    Traces t;
    t.trPA1 = s1 / probes;
    t.trPA0 = wantTrPA0 ? s0 / probes : 0.0;
    return t;
}

Traces get_traces(SigmaSolver& sol, const arma::sp_mat& psi,
                  const arma::fvec& w, const arma::fvec& tau,
                  const arma::fmat& Sigma_iX, const arma::fmat& cov,
                  const FitConfig& cfg, int it, bool wantTrPA0) {
    StageScope _(ST_TRACE);
    if (cfg.exactTrace) {
        // The exact path needs the inverse for THIS (w, tau); refresh here so
        // the requirement is explicit rather than inherited from whatever solve
        // happened last.
        sol.refresh(w, tau);
        return exact_traces(sol, psi, Sigma_iX, cov, wantTrPA0);
    }
    // Seeded per iteration so the probe stream at iteration k is the same for
    // every solver even if they take different numbers of iterations.
    return hutchinson_traces(sol, psi, w, tau, Sigma_iX, cov, cfg.probes,
                             cfg.seed + 1000ULL * (unsigned long long)it, wantTrPA0);
}

double rel_change(const arma::fvec& a, const arma::fvec& b, double tol) {
    double m = 0.0;
    for (arma::uword i = 0; i < a.n_elem; i++) {
        const double d = std::fabs((double)a(i) - (double)b(i));
        const double s = std::fabs((double)a(i)) + std::fabs((double)b(i)) + tol;
        m = std::max(m, d / s);
    }
    return m;
}

// Plain logistic regression, no random effect, to start eta somewhere sane.
// Uses no Sigma solve, so it costs the benchmark nothing and is identical for
// every solver.
arma::fvec logistic_init_eta(const arma::fvec& y, const arma::fmat& X) {
    arma::vec  yd = arma::conv_to<arma::vec>::from(y);
    arma::mat  Xd = arma::conv_to<arma::mat>::from(X);
    arma::vec  beta(Xd.n_cols, arma::fill::zeros);
    const double ybar = arma::mean(yd);
    beta(0) = std::log(ybar / (1.0 - ybar));
    for (int it = 0; it < 25; it++) {
        arma::vec eta = Xd * beta;
        arma::vec mu  = 1.0 / (1.0 + arma::exp(-eta));
        arma::vec wv  = arma::clamp(mu % (1.0 - mu), 1e-10, 1.0);
        arma::vec z   = eta + (yd - mu) / wv;
        arma::mat XtW = Xd.t() * arma::diagmat(wv);
        arma::vec nb  = arma::solve(XtW * Xd, XtW * z, arma::solve_opts::likely_sympd);
        const double mv = arma::max(arma::abs(nb - beta));
        beta = nb;
        if (mv < 1e-10) break;
    }
    return arma::conv_to<arma::fvec>::from(Xd * beta);
}

}  // namespace

// ===========================================================================
FitResult run_aireml(SigmaSolver& sol, const arma::sp_mat& psi,
                     const arma::fvec& y, const arma::fmat& X,
                     const FitConfig& cfg) {
    const arma::uword n = y.n_elem;
    FitResult R;
    const double t_begin = now_s();

    arma::fvec tau(2);
    arma::fvec w((arma::uword)n, arma::fill::ones);
    arma::fvec Y = y;
    arma::fvec eta((arma::uword)n, arma::fill::zeros);

    if (!cfg.binary) {
        // --- quantitative: W === 1, Y === y, both variance components free.
        // Initial tau: split the OLS residual variance evenly, the standard
        // GMMAT/AI-REML start. SAIGE starts at [1, 0] and spends its first
        // iteration on a "conservative" tau^2*score/n step to leave the
        // boundary; starting off the boundary is simpler and reaches the same
        // optimum, which is what the gates check.
        const arma::mat  Xd = arma::conv_to<arma::mat>::from(X);
        const arma::vec  yd = arma::conv_to<arma::vec>::from(y);
        const arma::vec  b  = arma::solve(Xd.t() * Xd, Xd.t() * yd);
        const arma::vec  r  = yd - Xd * b;
        const double s2 = arma::dot(r, r) / (double)(n - X.n_cols);
        tau(0) = (float)(0.5 * s2);
        tau(1) = (float)(0.5 * s2);
    } else {
        // --- binary: tau0 fixed at 1, tau1 free.
        tau(0) = 1.0f;
        tau(1) = 0.5f;
        eta = logistic_init_eta(y, X);
    }

    Coef coef;
    arma::fvec tau_prev = tau;

    for (int it = 0; it < cfg.maxiter; it++) {
        // ---------------------------------------------------------------
        // Inner loop: IRLS with tau HELD FIXED.
        //
        // Direction confirmed against ../step1_saige-null/saige_ai.cpp:52
        // ("the inner IRLS loop calls it on every iteration of every outer
        // AI-REML iteration") and loco_engine.cpp:243 ("PQL Newton loop with
        // tau held FIXED"). Quantitative degenerates: for a gaussian link
        // W === 1 and Y === y independently of eta, so one pass is exact and
        // the loop below runs once.
        // ---------------------------------------------------------------
        StageScope _solveStage(ST_SOLVE);
        arma::fvec alpha_prev;
        int inner = 0;
        for (; inner < (cfg.binary ? cfg.maxiter : 1); inner++) {
            if (cfg.binary) {
                const arma::fvec mu    = 1.0f / (1.0f + arma::exp(-eta));
                const arma::fvec muEta = mu % (1.0f - mu);
                w = arma::clamp(muEta, 1e-6f, 1.0f);
                Y = eta + (y - mu) / w;
            } else {
                w.ones();
                Y = y;
            }
            coef = get_coefficients(sol, Y, X, w, tau);
            R.innerIRLS++;
            if (!cfg.binary) { eta = coef.eta; break; }
            eta = coef.eta;
            if (inner > 0 && rel_change(coef.alpha, alpha_prev, cfg.tolCoef) < cfg.tolCoef) {
                inner++;
                break;
            }
            alpha_prev = coef.alpha;
        }

        // Rebuild W, Y at the converged eta and re-solve, so Sigma_iX and cov
        // belong to the SAME W the AI step uses. (Upstream SAIGE leaves them
        // on the previous W; that mixed operator is a faithfulness detail this
        // harness does not need and does not want, because it would make the
        // AI step depend on iteration history.)
        if (cfg.binary) {
            const arma::fvec mu    = 1.0f / (1.0f + arma::exp(-eta));
            const arma::fvec muEta = mu % (1.0f - mu);
            w = arma::clamp(muEta, 1e-6f, 1.0f);
            Y = eta + (y - mu) / w;
            coef = get_coefficients(sol, Y, X, w, tau);
            R.innerIRLS++;
        }

        // ---------------------------------------------------------------
        // AI-REML step.
        // ---------------------------------------------------------------
        const arma::fvec PY   = apply_P(coef.Sigma_iY, coef.Sigma_iX, coef.cov, Y);
        const arma::fvec A1PY = psi_mul(psi, PY);
        const arma::fvec A0PY = PY / w;                    // diag(1/w) PY

        const bool wantTrPA0 = !cfg.binary;
        const Traces tr = get_traces(sol, psi, w, tau, coef.Sigma_iX, coef.cov,
                                     cfg, it, wantTrPA0);

        const double yPA1Py = (double)arma::dot(PY, A1PY);
        const double score1 = yPA1Py - tr.trPA1;

        arma::fvec tau_new = tau;
        double score0 = 0.0;

        if (!cfg.binary) {
            const double yPA0Py = (double)arma::dot(PY, A0PY);
            score0 = yPA0Py - tr.trPA0;

            const arma::fvec PA0PY = apply_P(sigma_solve(sol, w, tau, A0PY),
                                             coef.Sigma_iX, coef.cov, A0PY);
            const arma::fvec PA1PY = apply_P(sigma_solve(sol, w, tau, A1PY),
                                             coef.Sigma_iX, coef.cov, A1PY);
            arma::mat AI(2, 2);
            AI(0, 0) = (double)arma::dot(A0PY, PA0PY);
            AI(1, 1) = (double)arma::dot(A1PY, PA1PY);
            AI(0, 1) = (double)arma::dot(A0PY, PA1PY);
            AI(1, 0) = AI(0, 1);

            arma::vec s(2); s(0) = score0; s(1) = score1;
            arma::vec delta = arma::solve(arma::symmatu(AI), s);

            double step = 1.0;
            tau_new(0) = (float)((double)tau(0) + delta(0));
            tau_new(1) = (float)((double)tau(1) + delta(1));
            while ((tau_new(0) < 0.0f || tau_new(1) < 0.0f) && step > 1e-10) {
                step *= 0.5;
                tau_new(0) = (float)((double)tau(0) + step * delta(0));
                tau_new(1) = (float)((double)tau(1) + step * delta(1));
            }
            tau_new(0) = std::max(0.0f, tau_new(0));
            tau_new(1) = std::max(0.0f, tau_new(1));
        } else {
            const arma::fvec PA1PY = apply_P(sigma_solve(sol, w, tau, A1PY),
                                             coef.Sigma_iX, coef.cov, A1PY);
            const double AI = std::max(1e-12, (double)arma::dot(A1PY, PA1PY));
            double t1 = (double)tau(1) + score1 / AI;
            double step = 1.0;
            while (t1 < 0.0 && step > 1e-10) {
                step *= 0.5;
                t1 = (double)tau(1) + step * score1 / AI;
            }
            tau_new(0) = 1.0f;                 // fixed
            tau_new(1) = (float)std::max(0.0, t1);
        }

        const double rc = rel_change(tau_new, tau, cfg.tol);
        R.traj.push_back({(double)it, (double)tau_new(0), (double)tau_new(1),
                          score0, score1, tr.trPA0, tr.trPA1, rc, (double)inner});
        tau_prev = tau;
        tau = tau_new;
        R.iters = it + 1;
        if (rc < cfg.tol) { R.converged = true; break; }
        if (!cfg.binary && (tau(0) <= 0.0f || tau(1) <= 0.0f)) break;
        if (cfg.binary && tau(1) <= 0.0f) break;
    }

    // Final fixed effects at the converged tau.
    {
        StageScope _(ST_SOLVE);
        coef = get_coefficients(sol, Y, X, w, tau);
        R.innerIRLS++;
    }
    R.tau   = tau;
    R.alpha = coef.alpha;
    R.fitSeconds = now_s() - t_begin;
    (void)tau_prev;
    return R;
}

}  // namespace step1bench
