// Standalone port of SKAT/Burden/SKAT-O tests for SAIGE Step 2
//
// Implements:
//   1. Davies method (characteristic function inversion) for mixture of chi-squared
//   2. Liu moment-matching method (fallback when Davies fails)
//   3. SKAT-O optimal.adj method with Gauss-Kronrod adaptive quadrature
//   4. get_SKAT_pvalue() matching SAIGE's R wrapper
//
// Original R sources:
//   - SAIGE/R/SAIGE_SPATest_Region_Func.R :: get_SKAT_pvalue()
//   - SKAT:::Met_SKAT_Get_Pvalue (R/SKAT_Optimal_Adj.R)
//   - SKAT:::SKAT_Optimal_PValue_Davies
//   - SKAT:::SKAT_liu.MOD.Lambda
//
// Conversions from R/Rcpp:
//   - R::pchisq  -> boost::math::cdf(chi_squared)
//   - R::qnorm   -> boost::math::quantile(normal)
//   - R::qchisq  -> boost::math::quantile(chi_squared)
//   - Rcpp::List  -> C++ structs
//   - Rcpp::stop  -> throw std::runtime_error

#include "skat.hpp"

#include <armadillo>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cfloat>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>


// ============================================================
// Davies qfc algorithm (AS 155) for P(Q > q)
// ============================================================
// Computes P(Q > q) where Q = sum(lambda_j * chi2(1)) using
// numerical inversion of the characteristic function.
//
// This is the classic "qfc" algorithm by Robert B. Davies (1980),
// as distributed in the CompQuadForm R package by P. Duchesne
// and P. Lafaye de Micheaux, and bundled in the SKAT R package.
//
// The code below is adapted from SKAT/src/qfc.cpp with:
//   - Global statics moved into a struct for reentrancy
//   - R.h / Rmath.h dependencies removed
//   - setjmp/longjmp replaced with a flag-based early exit
//   - Wrapped in the qfc_impl namespace
//
// Falls back to Liu moment-matching if qfc fails.
// ============================================================

namespace qfc_impl {

#define QFC_PI 3.14159265358979
#define QFC_LOG28 0.0866  /* log(2.0) / 8.0 */

struct QfcState {
    double sigsq, lmax, lmin, mean, c;
    double intl, ersm;
    int count, r, lim;
    bool ndtsrt, fail, exceeded;
    int *n, *th;
    double *lb, *nc;
};

static double exp1(double x) {
    return x < -50.0 ? 0.0 : exp(x);
}

static double square(double x) { return x * x; }
static double cube(double x) { return x * x * x; }

static void counter(QfcState& S) {
    S.count++;
    if (S.count > S.lim) S.exceeded = true;
}

static double log1(double x, bool first) {
    if (fabs(x) > 0.1) {
        return first ? log(1.0 + x) : (log(1.0 + x) - x);
    } else {
        double s, s1, term, y, k;
        y = x / (2.0 + x);  term = 2.0 * cube(y);  k = 3.0;
        s = (first ? 2.0 : -x) * y;
        y = square(y);
        for (s1 = s + term / k; s1 != s; s1 = s + term / k)
        { k = k + 2.0; term = term * y; s = s1; }
        return s;
    }
}

static void order(QfcState& S) {
    for (int j = 0; j < S.r; j++) {
        double lj = fabs(S.lb[j]);
        int k;
        for (k = j - 1; k >= 0; k--) {
            if (lj > fabs(S.lb[S.th[k]])) S.th[k + 1] = S.th[k];
            else break;
        }
        S.th[k + 1] = j;
    }
    S.ndtsrt = false;
}

static double errbd(QfcState& S, double u, double* cx) {
    double sum1, lj, ncj, x, y, xconst;
    int j, nj;
    counter(S);
    if (S.exceeded) return 0.0;
    xconst = u * S.sigsq;  sum1 = u * xconst;  u = 2.0 * u;
    for (j = S.r - 1; j >= 0; j--) {
        nj = S.n[j]; lj = S.lb[j]; ncj = S.nc[j];
        x = u * lj; y = 1.0 - x;
        xconst = xconst + lj * (ncj / y + nj) / y;
        sum1 = sum1 + ncj * square(x / y)
            + nj * (square(x) / y + log1(-x, false));
    }
    *cx = xconst;
    return exp1(-0.5 * sum1);
}

static double ctff(QfcState& S, double accx, double* upn) {
    double u1, u2, u, rb, xconst, c1, c2;
    u2 = *upn;   u1 = 0.0;  c1 = S.mean;
    rb = 2.0 * ((u2 > 0.0) ? S.lmax : S.lmin);
    for (u = u2 / (1.0 + u2 * rb); errbd(S, u, &c2) > accx;
         u = u2 / (1.0 + u2 * rb)) {
        if (S.exceeded) return c1;
        u1 = u2;  c1 = c2;  u2 = 2.0 * u2;
    }
    if (S.exceeded) return c1;
    for (u = (c1 - S.mean) / (c2 - S.mean); u < 0.9;
         u = (c1 - S.mean) / (c2 - S.mean)) {
        if (S.exceeded) return c2;
        u = (u1 + u2) / 2.0;
        if (errbd(S, u / (1.0 + u * rb), &xconst) > accx)
            { u1 = u; c1 = xconst; }
        else
            { u2 = u; c2 = xconst; }
    }
    *upn = u2;
    return c2;
}

static double truncation(QfcState& S, double u, double tausq) {
    double sum1, sum2, prod1, prod2, prod3, lj, ncj, x, y, err1, err2;
    int j, nj, s;
    counter(S);
    if (S.exceeded) return 1.0;
    sum1 = 0.0; prod2 = 0.0; prod3 = 0.0; s = 0;
    sum2 = (S.sigsq + tausq) * square(u); prod1 = 2.0 * sum2;
    u = 2.0 * u;
    for (j = 0; j < S.r; j++) {
        lj = S.lb[j]; ncj = S.nc[j]; nj = S.n[j];
        x = square(u * lj);
        sum1 = sum1 + ncj * x / (1.0 + x);
        if (x > 1.0) {
            prod2 = prod2 + nj * log(x);
            prod3 = prod3 + nj * log1(x, true);
            s = s + nj;
        } else {
            prod1 = prod1 + nj * log1(x, true);
        }
    }
    sum1 = 0.5 * sum1;
    prod2 = prod1 + prod2; prod3 = prod1 + prod3;
    x = exp1(-sum1 - 0.25 * prod2) / QFC_PI;
    y = exp1(-sum1 - 0.25 * prod3) / QFC_PI;
    err1 = (s == 0) ? 1.0 : x * 2.0 / s;
    err2 = (prod3 > 1.0) ? 2.5 * y : 1.0;
    if (err2 < err1) err1 = err2;
    x = 0.5 * sum2;
    err2 = (x <= y) ? 1.0 : y / x;
    return (err1 < err2) ? err1 : err2;
}

static void findu(QfcState& S, double* utx, double accx) {
    double u, ut;
    static double divis[] = {2.0, 1.4, 1.2, 1.1};
    ut = *utx; u = ut / 4.0;
    if (truncation(S, u, 0.0) > accx) {
        for (u = ut; truncation(S, u, 0.0) > accx; u = ut) {
            if (S.exceeded) { *utx = ut; return; }
            ut = ut * 4.0;
        }
    } else {
        ut = u;
        for (u = u / 4.0; truncation(S, u, 0.0) <= accx; u = u / 4.0) {
            if (S.exceeded) { *utx = ut; return; }
            ut = u;
        }
    }
    for (int i = 0; i < 4; i++) {
        if (S.exceeded) { *utx = ut; return; }
        u = ut / divis[i];
        if (truncation(S, u, 0.0) <= accx) ut = u;
    }
    *utx = ut;
}

static void do_integrate(QfcState& S, int nterm, double interv, double tausq, bool mainx) {
    double inpi, u, sum1, sum2, sum3, x, y, z;
    int k, j, nj;
    inpi = interv / QFC_PI;
    for (k = nterm; k >= 0; k--) {
        u = (k + 0.5) * interv;
        sum1 = -2.0 * u * S.c; sum2 = fabs(sum1);
        sum3 = -0.5 * S.sigsq * square(u);
        for (j = S.r - 1; j >= 0; j--) {
            nj = S.n[j]; x = 2.0 * S.lb[j] * u; y = square(x);
            sum3 = sum3 - 0.25 * nj * log1(y, true);
            y = S.nc[j] * x / (1.0 + y);
            z = nj * atan(x) + y;
            sum1 = sum1 + z; sum2 = sum2 + fabs(z);
            sum3 = sum3 - 0.5 * x * y;
        }
        x = inpi * exp1(sum3) / u;
        if (!mainx)
            x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
        sum1 = sin(0.5 * sum1) * x; sum2 = 0.5 * sum2 * x;
        S.intl = S.intl + sum1; S.ersm = S.ersm + sum2;
    }
}

static double cfe(QfcState& S, double x) {
    double axl, axl1, axl2, sxl, sum1, lj;
    int j, k, t;
    counter(S);
    if (S.exceeded) return 1.0;
    if (S.ndtsrt) order(S);
    axl = fabs(x); sxl = (x > 0.0) ? 1.0 : -1.0; sum1 = 0.0;
    for (j = S.r - 1; j >= 0; j--) {
        t = S.th[j];
        if (S.lb[t] * sxl > 0.0) {
            lj = fabs(S.lb[t]);
            axl1 = axl - lj * (S.n[t] + S.nc[t]); axl2 = lj / QFC_LOG28;
            if (axl1 > axl2) { axl = axl1; }
            else {
                if (axl > axl2) axl = axl2;
                sum1 = (axl - axl1) / lj;
                for (k = j - 1; k >= 0; k--)
                    sum1 = sum1 + (S.n[S.th[k]] + S.nc[S.th[k]]);
                goto done;
            }
        }
    }
done:
    if (sum1 > 100.0) { S.fail = true; return 1.0; }
    return pow(2.0, (sum1 / 4.0)) / (QFC_PI * square(axl));
}

// Main qfc entry point.
// Computes P(Q < c) where Q = sum(lb[j] * chi2(n[j], nc[j])) + sigma*N(0,1).
// For our use: n[j]=1, nc[j]=0, sigma=0.
// Returns qfval = P(Q < c). On failure, ifault != 0.
static void qfc(double* lb1, double* nc1, int* n1, int r1,
                double sigma, double c1, int lim1, double acc_in,
                double* trace, int* ifault, double* res) {
    QfcState S;
    S.exceeded = false;
    S.r = r1; S.lim = lim1; S.c = c1;
    S.n = n1; S.lb = lb1; S.nc = nc1;
    for (int j = 0; j < 7; j++) trace[j] = 0.0;
    *ifault = 0; S.count = 0;
    S.intl = 0.0; S.ersm = 0.0;
    double qfval = -1.0;
    double acc1 = acc_in;
    S.ndtsrt = true; S.fail = false;
    double xlim = (double)S.lim;

    S.th = (int*)malloc(S.r * sizeof(int));
    if (!S.th) { *ifault = 5; res[0] = -1.0; return; }

    // find mean, sd, max and min of lb; validate parameters
    S.sigsq = square(sigma); double sd = S.sigsq;
    S.lmax = 0.0; S.lmin = 0.0; S.mean = 0.0;
    for (int j = 0; j < S.r; j++) {
        int nj = S.n[j]; double lj = S.lb[j]; double ncj = S.nc[j];
        if (nj < 0 || ncj < 0.0) { *ifault = 3; goto endofproc; }
        sd = sd + square(lj) * (2 * nj + 4.0 * ncj);
        S.mean = S.mean + lj * (nj + ncj);
        if (S.lmax < lj) S.lmax = lj;
        else if (S.lmin > lj) S.lmin = lj;
    }
    if (sd == 0.0) { qfval = (S.c > 0.0) ? 1.0 : 0.0; goto endofproc; }
    if (S.lmin == 0.0 && S.lmax == 0.0 && sigma == 0.0) { *ifault = 3; goto endofproc; }
    sd = sqrt(sd);
    {
        double almx = (S.lmax < -S.lmin) ? -S.lmin : S.lmax;
        double utx = 16.0 / sd;
        double up = 4.5 / sd;
        double un = -up;

        // truncation point with no convergence factor
        findu(S, &utx, 0.5 * acc1);
        if (S.exceeded) { *ifault = 4; goto endofproc; }

        // does convergence factor help?
        if (S.c != 0.0 && (almx > 0.07 * sd)) {
            double tausq = 0.25 * acc1 / cfe(S, S.c);
            if (S.exceeded) { *ifault = 4; goto endofproc; }
            if (S.fail) S.fail = false;
            else if (truncation(S, utx, tausq) < 0.2 * acc1) {
                if (S.exceeded) { *ifault = 4; goto endofproc; }
                S.sigsq = S.sigsq + tausq;
                findu(S, &utx, 0.25 * acc1);
                if (S.exceeded) { *ifault = 4; goto endofproc; }
                trace[5] = sqrt(tausq);
            }
        }
        trace[4] = utx; acc1 = 0.5 * acc1;

        // find range of distribution
        double d1, d2, intv, xnt, xntm;
    l1:
        d1 = ctff(S, acc1, &up) - S.c;
        if (S.exceeded) { *ifault = 4; goto endofproc; }
        if (d1 < 0.0) { qfval = 1.0; goto endofproc; }
        d2 = S.c - ctff(S, acc1, &un);
        if (S.exceeded) { *ifault = 4; goto endofproc; }
        if (d2 < 0.0) { qfval = 0.0; goto endofproc; }

        // find integration interval
        intv = 2.0 * QFC_PI / ((d1 > d2) ? d1 : d2);
        xnt = utx / intv; xntm = 3.0 / sqrt(acc1);
        if (xnt > xntm * 1.5) {
            // parameters for auxiliary integration
            if (xntm > xlim) { *ifault = 1; goto endofproc; }
            int ntm = (int)floor(xntm + 0.5);
            double intv1 = utx / ntm;
            double x2 = 2.0 * QFC_PI / intv1;
            if (x2 <= fabs(S.c)) goto l2;
            double tausq = 0.33 * acc1 / (1.1 * (cfe(S, S.c - x2) + cfe(S, S.c + x2)));
            if (S.exceeded) { *ifault = 4; goto endofproc; }
            if (S.fail) goto l2;
            acc1 = 0.67 * acc1;
            do_integrate(S, ntm, intv1, tausq, false);
            xlim = xlim - xntm; S.sigsq = S.sigsq + tausq;
            trace[2] = trace[2] + 1; trace[1] = trace[1] + ntm + 1;
            findu(S, &utx, 0.25 * acc1);
            if (S.exceeded) { *ifault = 4; goto endofproc; }
            acc1 = 0.75 * acc1;
            goto l1;
        }

    l2:
        trace[3] = intv;
        if (xnt > xlim) { *ifault = 1; goto endofproc; }
        {
            int nt = (int)floor(xnt + 0.5);
            do_integrate(S, nt, intv, 0.0, true);
            trace[2] = trace[2] + 1; trace[1] = trace[1] + nt + 1;
        }
        qfval = 0.5 - S.intl;
        trace[0] = S.ersm;

        // test whether round-off error could be significant
        {
            static int rats[] = {1, 2, 4, 8};
            double up2 = S.ersm;
            double x = up2 + acc_in / 10.0;
            for (int j = 0; j < 4; j++) {
                if (rats[j] * x == rats[j] * up2) *ifault = 2;
            }
        }
    }

endofproc:
    free(S.th);
    trace[6] = (double)S.count;
    res[0] = qfval;
}

#undef QFC_PI
#undef QFC_LOG28

} // namespace qfc_impl


double davies_pvalue(double q, const arma::vec& lambda) {
    // Filter out zero/near-zero eigenvalues
    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    int n = (int)lam.n_elem;

    if (n == 0) {
        return (q > 0.0) ? 0.0 : 1.0;
    }

    // Special case: single eigenvalue -> simple scaled chi-squared(1)
    if (n == 1) {
        double l1 = lam(0);
        if (l1 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        double q_scaled = q / l1;
        if (q_scaled <= 0.0) return 1.0;
        try {
            boost::math::chi_squared chi2_1(1.0);
            return boost::math::cdf(boost::math::complement(chi2_1, q_scaled));
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Call qfc (Davies 1980, AS 155) via the SKAT/CompQuadForm implementation.
    // Interface matches R's:
    //   SKAT_davies(Q, lambda, h=rep(1,r), delta=rep(0,r), sigma=0, lim=10000, acc=1e-6)
    //   out$Qq = 1 - out$res   (qfc returns P(Q<q), we need P(Q>q))
    std::vector<double> lb(lam.memptr(), lam.memptr() + n);
    std::vector<double> nc(n, 0.0);   // non-centrality parameters (all zero)
    std::vector<int> df(n, 1);        // degrees of freedom (all 1)
    double sigma = 0.0;
    int lim = 10000;
    double acc = 1e-6;
    double trace[7];
    int ifault = 0;
    double res = 0.0;

    qfc_impl::qfc(lb.data(), nc.data(), df.data(), n,
                   sigma, q, lim, acc, trace, &ifault, &res);

    // qfc returns P(Q < q) in res; we need P(Q > q)
    double pval = 1.0 - res;

    // Check convergence (matches R's Get_Davies_PVal logic):
    //   if ifault != 0 -> fall back to Liu
    //   if p > 1 or p <= 0 -> fall back to Liu
    if (ifault == 0 && std::isfinite(pval) && pval > 0.0 && pval <= 1.0) {
        return pval;
    }

    // Fall back to Liu moment-matching
    return liu_pvalue(q, lam);
}


// Strict variant of davies_pvalue for use inside the SKAT-O integrand.
// Returns NaN when Davies fails (ifault != 0) instead of silently falling
// back to Liu.  This matches R's behavior where the integrand calls stop()
// on ifault != 0, causing try(integrate(...)) to fail and trigger the Liu
// integration fallback path.
double davies_pvalue_strict(double q, const arma::vec& lambda) {
    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    int n = (int)lam.n_elem;

    if (n == 0) {
        return (q > 0.0) ? 0.0 : 1.0;
    }

    if (n == 1) {
        double l1 = lam(0);
        if (l1 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        double q_scaled = q / l1;
        if (q_scaled <= 0.0) return 1.0;
        try {
            boost::math::chi_squared chi2_1(1.0);
            return boost::math::cdf(boost::math::complement(chi2_1, q_scaled));
        } catch (...) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    std::vector<double> lb(lam.memptr(), lam.memptr() + n);
    std::vector<double> nc(n, 0.0);
    std::vector<int> df(n, 1);
    double sigma = 0.0;
    int lim = 10000;
    double acc = 1e-6;
    double trace[7];
    int ifault = 0;
    double res = 0.0;

    qfc_impl::qfc(lb.data(), nc.data(), df.data(), n,
                   sigma, q, lim, acc, trace, &ifault, &res);

    double pval = 1.0 - res;

    // Strict check: if ifault != 0 or pval out of range, return NaN
    // (do NOT fall back to Liu — caller will detect and handle)
    if (ifault != 0 || !std::isfinite(pval) || pval <= 0.0 || pval > 1.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return pval;
}


// ============================================================
// Liu moment-matching method
// ============================================================
// Approximates P(Q > q) where Q = sum(lambda_i * chi2(1))
// by matching the first 4 cumulants to a scaled noncentral chi-squared.
//
// From SKAT:::SKAT_liu.MOD.Lambda:
//   c_k = sum(lambda^k) for k=1,2,3,4
//   s1 = c3/c2^1.5, s2 = c4/c2^2
//   if s1^2 > s2: use noncentral chi-squared
//   else: use central chi-squared
//
// The normalized statistic is mapped to the chi-squared distribution.

double liu_pvalue(double q, const arma::vec& lambda) {
    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    int n = (int)lam.n_elem;

    if (n == 0) return (q > 0.0) ? 0.0 : 1.0;

    // Cumulant sums: c_k = sum(lambda^k)
    double c1 = arma::sum(lam);
    double c2 = arma::sum(lam % lam);
    double c3 = arma::sum(lam % lam % lam);
    double c4 = arma::sum(lam % lam % lam % lam);

    if (c2 <= 0.0) return (q > c1) ? 0.0 : 1.0;

    // Skewness and kurtosis ratios
    double s1 = c3 / std::pow(c2, 1.5);
    double s2 = c4 / (c2 * c2);

    double muQ = c1;
    double sigmaQ = std::sqrt(2.0 * c2);

    // Fit chi-squared approximation parameters
    double l, delta;  // df and noncentrality

    if (s1 * s1 > s2) {
        // Noncentral chi-squared approximation
        double a = 1.0 / (s1 - std::sqrt(s1 * s1 - s2));
        delta = s1 * a * a * a - a * a;
        l = a * a - 2.0 * delta;
    } else {
        // Central chi-squared approximation
        delta = 0.0;
        l = 1.0 / s2;
    }

    // Ensure valid parameters
    if (l <= 0.0) l = 1.0;
    if (delta < 0.0) delta = 0.0;

    // Normalize Q to the chi-squared scale
    double Qnorm = (q - muQ) / sigmaQ * std::sqrt(2.0 * l + 4.0 * delta) + l + delta;

    if (Qnorm <= 0.0) return 1.0;

    // Compute p-value
    try {
        if (delta > 1e-6) {
            boost::math::non_central_chi_squared ncchi2(l, delta);
            return boost::math::cdf(boost::math::complement(ncchi2, Qnorm));
        } else {
            boost::math::chi_squared chi2(l);
            return boost::math::cdf(boost::math::complement(chi2, Qnorm));
        }
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}


// ============================================================
// Internal: Liu method returning all parameters (for SKAT-O)
// ============================================================

struct LiuParams {
    double pvalue;
    double muQ;
    double sigmaQ;
    double l;       // df
    double delta;   // noncentrality
};

static LiuParams liu_params(const arma::vec& lambda) {
    LiuParams res;
    res.pvalue = 1.0;
    res.muQ = 0.0;
    res.sigmaQ = 1.0;
    res.l = 1.0;
    res.delta = 0.0;

    arma::vec lam = lambda(arma::find(arma::abs(lambda) > 1e-10));
    if (lam.n_elem == 0) return res;

    double c1 = arma::sum(lam);
    double c2 = arma::sum(lam % lam);
    double c3 = arma::sum(lam % lam % lam);
    double c4 = arma::sum(lam % lam % lam % lam);

    if (c2 <= 0.0) return res;

    double s1 = c3 / std::pow(c2, 1.5);
    double s2 = c4 / (c2 * c2);

    res.muQ = c1;
    res.sigmaQ = std::sqrt(2.0 * c2);

    if (s1 * s1 > s2) {
        double a = 1.0 / (s1 - std::sqrt(s1 * s1 - s2));
        res.delta = s1 * a * a * a - a * a;
        res.l = a * a - 2.0 * res.delta;
    } else {
        res.delta = 0.0;
        res.l = 1.0 / s2;
    }

    if (res.l <= 0.0) res.l = 1.0;
    if (res.delta < 0.0) res.delta = 0.0;

    return res;
}


// ============================================================
// Helper: compute eigenvalues of Phi_rho for a given rho
// ============================================================
// Q(rho) = Score' * R_rho * Score where R_rho = (1-rho)*I + rho*11'
// Under H0, Q(rho) ~ sum(lambda_j * chi2(1)) where lambda_j = eigenvalues
// of R_rho^{1/2} * Phi * R_rho^{1/2}

// R's Get_Lambda eigenvalue filter:
//   lambda1 = eigen(K)$values
//   IDX1 = which(lambda1 >= 0)
//   IDX2 = which(lambda1 > mean(lambda1[IDX1]) / 1e5)
//   lambda = lambda1[IDX2]
static arma::vec filter_eigenvalues_R(const arma::vec& ev_all) {
    arma::vec ev_pos = ev_all(arma::find(ev_all >= 0));
    double mean_pos = 0.0;
    if (ev_pos.n_elem > 0) {
        mean_pos = arma::mean(ev_pos);
    }
    double threshold = mean_pos / 1e5;
    return ev_all(arma::find(ev_all > threshold));
}

static arma::vec compute_phi_rho_eigenvalues(const arma::mat& Phi, double rho, int m) {
    if (std::abs(rho) < 1e-10) {
        // SKAT: eigenvalues of Phi
        arma::vec ev;
        arma::eig_sym(ev, Phi);
        return filter_eigenvalues_R(ev);
    }

    if (std::abs(rho - 1.0) < 1e-10) {
        // Burden: single eigenvalue = sum(Phi)
        double total = arma::accu(Phi);
        if (total > 1e-10) {
            arma::vec result(1);
            result(0) = total;
            return result;
        }
        return arma::vec();
    }

    // General case: R_rho = (1-rho)*I + rho*11'
    // Eigenvalues of R_rho: (1-rho) [multiplicity m-1] and (1-rho+m*rho) [multiplicity 1]
    // sqrt(R_rho) has same eigenvectors with sqrt eigenvalues
    double ev_small = 1.0 - rho;
    double ev_large = 1.0 - rho + m * rho;

    arma::vec ones_v = arma::ones<arma::vec>(m);
    arma::mat P_proj = ones_v * ones_v.t() / (double)m;  // projection onto 1
    arma::mat sqrtR = std::sqrt(ev_small) * (arma::eye<arma::mat>(m, m) - P_proj)
                    + std::sqrt(ev_large) * P_proj;

    arma::mat Phi_rho = sqrtR * Phi * sqrtR;
    Phi_rho = 0.5 * (Phi_rho + Phi_rho.t());  // ensure numerical symmetry

    arma::vec ev;
    arma::eig_sym(ev, Phi_rho);
    return filter_eigenvalues_R(ev);
}


// ============================================================
// Gauss-Kronrod adaptive quadrature (15-point / 7-point)
// ============================================================
// This replicates R's integrate() function which uses adaptive
// Gauss-Kronrod quadrature. R's integrate() uses QUADPACK's dqags
// which employs the 21-point Gauss-Kronrod rule (G10/K21) with
// embedded 10-point Gauss rule for automatic error estimation
// and adaptive subdivision.

namespace gauss_kronrod {

// ============================================================
// Faithful C++ port of R's QUADPACK dqags algorithm from
// R source: src/appl/integrate.c (GPL-2+, QUADPACK is public domain)
//
// Ported functions:
//   rdqk21  -> qk21   (21-point Gauss-Kronrod rule)
//   rdqpsrt -> qpsrt  (interval error sorting)
//   rdqelg  -> qelg   (Wynn epsilon extrapolation)
//   rdqagse -> integrate (main adaptive integrator with extrapolation)
//
// The code preserves R's exact variable names and control flow
// to ensure bit-identical numerical results.
// ============================================================

// Adaptive quadrature result (public interface)
struct QuadResult {
    double value;
    double error;
    int neval;
    bool converged;
};

// 21-point Gauss-Kronrod nodes and weights (on [-1, 1])
// Identical to R's rdqk21 constants (evaluated with 80 decimal digit arithmetic
// by L. W. Fullerton, Bell Labs, Nov. 1981)
static const double xgk[11] = {
    .995657163025808080735527280689003,
    .973906528517171720077964012084452,
    .930157491355708226001207180059508,
    .865063366688984510732096688423493,
    .780817726586416897063717578345042,
    .679409568299024406234327365114874,
    .562757134668604683339000099272694,
    .433395394129247190799265943165784,
    .294392862701460198131126603103866,
    .14887433898163121088482600112972,
    0.
};

static const double wgk[11] = {
    .011694638867371874278064396062192,
    .03255816230796472747881897245939,
    .05475589657435199603138130024458,
    .07503967481091995276704314091619,
    .093125454583697605535065465083366,
    .109387158802297641899210590325805,
    .123491976262065851077958109831074,
    .134709217311473325928054001771707,
    .142775938577060080797094273138717,
    .147739104901338491374841515972068,
    .149445554002916905664936468389821
};

static const double wg[5] = {
    .066671344308688137593568809893332,
    .149451349150580593145776339657697,
    .219086362515982043995534934228163,
    .269266719309996355091226921569469,
    .295524224714752870173892994651338
};

// ----------------------------------------------------------------
// rdqk21: 21-point Gauss-Kronrod rule on [a, b]
// Faithful port of R's rdqk21 from integrate.c
//
// Returns: result, abserr, resabs (integral of |f|), resasc (integral of |f - mean|)
// ----------------------------------------------------------------
static void qk21(const std::function<double(double)>& f,
                  double a, double b,
                  double& result, double& abserr,
                  double& resabs, double& resasc) {
    double fv1[10], fv2[10];
    double absc, resg, resk, fsum, fval1, fval2;
    double hlgth, centr, reskh, uflow;
    double fc, epmach, dhlgth;
    int j, jtw, jtwm1;

    epmach = DBL_EPSILON;
    uflow = DBL_MIN;

    centr = (a + b) * .5;
    hlgth = (b - a) * .5;
    dhlgth = fabs(hlgth);

    // Compute the 21-point Kronrod approximation to the integral,
    // and estimate the absolute error.
    // R's version vectorizes f calls, but we call f one point at a time.
    resg = 0.;
    fc = f(centr);
    resk = wgk[10] * fc;
    resabs = fabs(resk);

    for (j = 1; j <= 5; ++j) {
        jtw = j << 1;   // 2, 4, 6, 8, 10 -> xgk index jtw-1 = 1,3,5,7,9
        absc = hlgth * xgk[jtw - 1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtw - 1] = fval1;
        fv2[jtw - 1] = fval2;
        fsum = fval1 + fval2;
        resg += wg[j - 1] * fsum;
        resk += wgk[jtw - 1] * fsum;
        resabs += wgk[jtw - 1] * (fabs(fval1) + fabs(fval2));
    }

    for (j = 1; j <= 5; ++j) {
        jtwm1 = (j << 1) - 1;  // 1, 3, 5, 7, 9 -> xgk index jtwm1-1 = 0,2,4,6,8
        absc = hlgth * xgk[jtwm1 - 1];
        fval1 = f(centr - absc);
        fval2 = f(centr + absc);
        fv1[jtwm1 - 1] = fval1;
        fv2[jtwm1 - 1] = fval2;
        fsum = fval1 + fval2;
        resk += wgk[jtwm1 - 1] * fsum;
        resabs += wgk[jtwm1 - 1] * (fabs(fval1) + fabs(fval2));
    }

    reskh = resk * .5;
    resasc = wgk[10] * fabs(fc - reskh);
    for (j = 1; j <= 10; ++j) {
        resasc += wgk[j - 1] * (fabs(fv1[j - 1] - reskh) +
                                 fabs(fv2[j - 1] - reskh));
    }

    result = resk * hlgth;
    resabs *= dhlgth;
    resasc *= dhlgth;
    abserr = fabs((resk - resg) * hlgth);
    if (resasc != 0. && abserr != 0.) {
        abserr = resasc * std::min(1., pow(abserr * 200. / resasc, 1.5));
    }
    if (resabs > uflow / (epmach * 50.)) {
        abserr = std::max(epmach * 50. * resabs, abserr);
    }
}

// ----------------------------------------------------------------
// rdqpsrt: maintain descending ordering of error estimates
// Faithful port of R's rdqpsrt from integrate.c
//
// Uses 1-based indexing internally (Fortran heritage).
// Caller passes 1-based arrays (elist, iord).
// ----------------------------------------------------------------
static void qpsrt(int limit, int last, int& maxerr,
                   double& ermax, double* elist, int* iord, int& nrmax) {
    int i, j, k, ido, jbnd, isucc, jupbn;
    double errmin, errmax_local;

    // elist and iord are 1-based (caller passes &array[1])
    // so elist[1..last] and iord[1..last] are valid

    if (last <= 2) {
        iord[1] = 1;
        iord[2] = 2;
        goto Last;
    }

    errmax_local = elist[maxerr];
    if (nrmax > 1) {
        ido = nrmax - 1;
        for (i = 1; i <= ido; ++i) {
            isucc = iord[nrmax - 1];
            if (errmax_local <= elist[isucc])
                break;
            iord[nrmax] = isucc;
            --nrmax;
        }
    }

    if (last > limit / 2 + 2)
        jupbn = limit + 3 - last;
    else
        jupbn = last;

    errmin = elist[last];

    jbnd = jupbn - 1;
    for (i = nrmax + 1; i <= jbnd; ++i) {
        isucc = iord[i];
        if (errmax_local >= elist[isucc]) {
            // Insert errmax here, then insert errmin bottom-up
            iord[i - 1] = maxerr;
            for (j = i, k = jbnd; j <= jbnd; j++, k--) {
                isucc = iord[k];
                if (errmin < elist[isucc]) {
                    iord[k + 1] = last;
                    goto Last;
                }
                iord[k + 1] = isucc;
            }
            iord[i] = last;
            goto Last;
        }
        iord[i - 1] = isucc;
    }

    iord[jbnd] = maxerr;
    iord[jupbn] = last;

Last:
    maxerr = iord[nrmax];
    ermax = elist[maxerr];
}

// ----------------------------------------------------------------
// rdqelg: Wynn epsilon extrapolation algorithm
// Faithful port of R's rdqelg from integrate.c
//
// Determines the limit of a given sequence of approximations using
// the epsilon algorithm of P. Wynn. The condensed epsilon table is
// computed; only elements needed for the next diagonal are preserved.
//
// Parameters (1-based Fortran convention internally):
//   n      - epstab[n] contains the new element (modified on output)
//   epstab - double[52], elements of the two lower diagonals (0-based in C)
//   result - resulting approximation
//   abserr - estimate of absolute error
//   res3la - double[3], last 3 results (0-based in C, but used 1-based internally)
//   nres   - number of calls (should be 0 at first call)
// ----------------------------------------------------------------
static void qelg(int& n, double* epstab, double& result,
                  double& abserr, double* res3la, int& nres) {
    int i__, indx, ib, ib2, ie, k1, k2, k3, num, newelm, limexp;
    double delta1, delta2, delta3, e0, e1, e1abs, e2, e3, epmach, epsinf;
    double oflow, ss, res;
    double errA, err1, err2, err3, tol1, tol2, tol3;

    // Note: R's code uses 1-based indexing via "--epstab; --res3la;"
    // adjustments. We use 1-based access by subtracting 1 from pointer.
    // To match R exactly, we replicate the pointer adjustment.
    double* epstab1 = epstab - 1;  // epstab1[1] == epstab[0]
    double* res3la1 = res3la - 1;  // res3la1[1] == res3la[0]

    epmach = DBL_EPSILON;
    oflow = DBL_MAX;
    ++nres;
    abserr = oflow;
    result = epstab1[n];
    if (n < 3) {
        goto L100;
    }
    limexp = 50;
    epstab1[n + 2] = epstab1[n];
    newelm = (n - 1) / 2;
    epstab1[n] = oflow;
    num = n;
    k1 = n;
    for (i__ = 1; i__ <= newelm; ++i__) {
        k2 = k1 - 1;
        k3 = k1 - 2;
        res = epstab1[k1 + 2];
        e0 = epstab1[k3];
        e1 = epstab1[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2 - e1;
        err2 = fabs(delta2);
        tol2 = std::max(fabs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3 = fabs(delta3);
        tol3 = std::max(e1abs, fabs(e0)) * epmach;
        if (err2 <= tol2 && err3 <= tol3) {
            // e0, e1 and e2 are equal to within machine accuracy,
            // convergence is assumed.
            result = res;
            abserr = err2 + err3;
            goto L100;
        }

        e3 = epstab1[k1];
        epstab1[k1] = e1;
        delta1 = e1 - e3;
        err1 = fabs(delta1);
        tol1 = std::max(e1abs, fabs(e3)) * epmach;

        // If two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        if (err1 > tol1 && err2 > tol2 && err3 > tol3) {
            ss = 1. / delta1 + 1. / delta2 - 1. / delta3;
            epsinf = fabs(ss * e1);

            // Test to detect irregular behaviour in the table, and
            // eventually omit a part of the table adjusting the value of n.
            if (epsinf > 1e-4) {
                goto L30;
            }
        }

        n = i__ + i__ - 1;
        goto L50;

    L30:
        // Compute a new element and eventually adjust the value of result.
        res = e1 + 1. / ss;
        epstab1[k1] = res;
        k1 += -2;
        errA = err2 + fabs(res - e2) + err3;
        if (errA <= abserr) {
            abserr = errA;
            result = res;
        }
    }

    // Shift the table.
L50:
    if (n == limexp) {
        n = (limexp / 2 << 1) - 1;
    }

    if (num / 2 << 1 == num) ib = 2; else ib = 1;
    ie = newelm + 1;
    for (i__ = 1; i__ <= ie; ++i__) {
        ib2 = ib + 2;
        epstab1[ib] = epstab1[ib2];
        ib = ib2;
    }
    if (num != n) {
        indx = num - n + 1;
        for (i__ = 1; i__ <= n; ++i__) {
            epstab1[i__] = epstab1[indx];
            ++indx;
        }
    }
    if (nres >= 4) {
        abserr = fabs(result - res3la1[3]) +
                 fabs(result - res3la1[2]) +
                 fabs(result - res3la1[1]);
        res3la1[1] = res3la1[2];
        res3la1[2] = res3la1[3];
        res3la1[3] = result;
    } else {
        res3la1[nres] = result;
        abserr = oflow;
    }

L100:
    abserr = std::max(abserr, epmach * 5. * fabs(result));
}

// ----------------------------------------------------------------
// rdqagse: main adaptive integrator with Wynn epsilon extrapolation
// Faithful port of R's rdqagse from integrate.c
//
// This is the core QUADPACK dqags algorithm that R's integrate() calls.
// It uses adaptive subdivision with 21-point Gauss-Kronrod rule and
// Wynn epsilon extrapolation for convergence acceleration.
// ----------------------------------------------------------------
static QuadResult integrate(const std::function<double(double)>& f,
                             double a, double b,
                             double rel_tol = 1.220703e-4,
                             double abs_tol = 1e-25,
                             int limit = 1000) {
    // Allocate work arrays (1-based indexing, so size limit+1)
    std::vector<double> alist(limit + 1);
    std::vector<double> blist(limit + 1);
    std::vector<double> rlist(limit + 1);
    std::vector<double> elist(limit + 1);
    std::vector<int>    iord(limit + 1);

    // Local variables (matching R's rdqagse exactly)
    bool noext, extrap;
    int k, ksgn, nres;
    int ierro;
    int ktmin, nrmax;
    int iroff1, iroff2, iroff3;
    int id;
    int numrl2;
    int jupbnd;
    int maxerr;
    int last;
    int ier;
    double res3la[3];
    double rlist2[52];
    double abseps, area, area1, area2, area12, dres, epmach;
    double a1, a2, b1, b2, defabs, defab1, defab2, oflow, uflow, resabs, reseps;
    double error1, error2, erro12, errbnd, erlast, errmax, errsum;
    double result_val, abserr_val;
    int neval;

    double correc = 0.0, erlarg = 0.0, ertest = 0.0, small = 0.0;

    epmach = DBL_EPSILON;

    // Test on validity of parameters
    ier = 0;
    neval = 0;
    last = 0;
    result_val = 0.;
    abserr_val = 0.;
    alist[1] = a;
    blist[1] = b;
    rlist[1] = 0.;
    elist[1] = 0.;
    if (abs_tol <= 0. && rel_tol < std::max(epmach * 50., 5e-29)) {
        ier = 6;
        // Return with error
        QuadResult res;
        res.value = result_val;
        res.error = abserr_val;
        res.neval = neval;
        res.converged = false;
        return res;
    }

    // First approximation to the integral
    uflow = DBL_MIN;
    oflow = DBL_MAX;
    ierro = 0;
    {
        double defabs_init, resabs_init;
        qk21(f, a, b, result_val, abserr_val, defabs_init, resabs_init);
        defabs = defabs_init;
        resabs = resabs_init;
    }

    // Test on accuracy
    dres = fabs(result_val);
    errbnd = std::max(abs_tol, rel_tol * dres);
    last = 1;
    rlist[1] = result_val;
    elist[1] = abserr_val;
    iord[1] = 1;
    if (abserr_val <= epmach * 100. * defabs && abserr_val > errbnd)
        ier = 2;
    if (limit == 1)
        ier = 1;
    if (ier != 0 || (abserr_val <= errbnd && abserr_val != resabs)
        || abserr_val == 0.) {
        // Jump to final result
        neval = last * 42 - 21;
        QuadResult res;
        res.value = result_val;
        res.error = abserr_val;
        res.neval = neval;
        res.converged = (ier == 0);
        return res;
    }

    // Initialization
    rlist2[0] = result_val;
    errmax = abserr_val;
    maxerr = 1;
    area = result_val;
    errsum = abserr_val;
    abserr_val = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = false;
    noext = false;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1. - epmach * 50.) * defabs) {
        ksgn = 1;
    }

    // Main do-loop
    for (last = 2; last <= limit; ++last) {

        // Bisect the subinterval with the nrmax-th largest error estimate.
        a1 = alist[maxerr];
        b1 = (alist[maxerr] + blist[maxerr]) * .5;
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        qk21(f, a1, b1, area1, error1, resabs, defab1);
        qk21(f, a2, b2, area2, error2, resabs, defab2);

        // Improve previous approximations to integral
        // and error and test for accuracy.
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area = area + area12 - rlist[maxerr];
        if (!(defab1 == error1 || defab2 == error2)) {

            if (fabs(rlist[maxerr] - area12) <= fabs(area12) * 1e-5 &&
                erro12 >= errmax * .99) {
                if (extrap)
                    ++iroff2;
                else
                    ++iroff1;
            }
            if (last > 10 && erro12 > errmax)
                ++iroff3;
        }
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = std::max(abs_tol, rel_tol * fabs(area));

        // Test for roundoff error and eventually set error flag.
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20)
            ier = 2;
        if (iroff2 >= 5)
            ierro = 3;

        // Set error flag in the case that the number of subintervals equals limit.
        if (last == limit)
            ier = 1;

        // Set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        if (std::max(fabs(a1), fabs(b2)) <=
            (epmach * 100. + 1.) * (fabs(a2) + uflow * 1e3)) {
            ier = 4;
        }

        // Append the newly-created intervals to the list.
        if (error2 > error1) {
            alist[maxerr] = a2;
            alist[last] = a1;
            blist[last] = b1;
            rlist[maxerr] = area2;
            rlist[last] = area1;
            elist[maxerr] = error2;
            elist[last] = error1;
        } else {
            alist[last] = a2;
            blist[maxerr] = b1;
            blist[last] = b2;
            elist[maxerr] = error1;
            elist[last] = error2;
        }

        // Call qpsrt to maintain the descending ordering in the list of
        // error estimates and select the subinterval with nrmax-th largest
        // error estimate (to be bisected next).
        // qpsrt uses 1-based indexing on elist[] and iord[].
        // Our vectors are sized (limit+1) with valid indices [1..last].
        qpsrt(limit, last, maxerr, errmax, elist.data(), iord.data(), nrmax);

        if (errsum <= errbnd) goto L115;
        if (ier != 0) break;
        if (last == 2) {
            small = fabs(b - a) * .375;
            erlarg = errsum;
            ertest = errbnd;
            rlist2[1] = area;
            continue;
        }
        if (noext) continue;

        erlarg -= erlast;
        if (fabs(b1 - a1) > small) {
            erlarg += erro12;
        }
        if (!extrap) {
            // Test whether the interval to be bisected next is the
            // smallest interval.
            if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                continue;
            }
            extrap = true;
            nrmax = 2;
        }

        if (ierro != 3 && erlarg > ertest) {
            // The smallest interval has the largest error.
            // Before bisecting, decrease the sum of the errors over the
            // larger intervals (erlarg) and perform extrapolation.
            id = nrmax;
            jupbnd = last;
            if (last > limit / 2 + 2) {
                jupbnd = limit + 3 - last;
            }
            for (k = id; k <= jupbnd; ++k) {
                maxerr = iord[nrmax];
                errmax = elist[maxerr];
                if (fabs(blist[maxerr] - alist[maxerr]) > small) {
                    goto L90;
                }
                ++nrmax;
            }
        }

        // Perform extrapolation.
        ++numrl2;
        rlist2[numrl2 - 1] = area;
        qelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ++ktmin;
        if (ktmin > 5 && abserr_val < errsum * .001) {
            ier = 5;
        }
        if (abseps < abserr_val) {
            ktmin = 0;
            abserr_val = abseps;
            result_val = reseps;
            correc = erlarg;
            ertest = std::max(abs_tol, rel_tol * fabs(reseps));
            if (abserr_val <= ertest) {
                break;
            }
        }

        // Prepare bisection of the smallest interval.
        if (numrl2 == 1) {
            noext = true;
        }
        if (ier == 5) {
            break;
        }
        maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = false;
        small *= .5;
        erlarg = errsum;
    L90:
        ;
    }

    // Set final result and error estimate.
    if (abserr_val == oflow) goto L115;
    if (ier + ierro == 0) goto L110;
    if (ierro == 3)
        abserr_val += correc;
    if (ier == 0)
        ier = 3;
    if (result_val == 0. || area == 0.) {
        if (abserr_val > errsum) goto L115;
        if (area == 0.) goto L130;
    } else {
        if (abserr_val / fabs(result_val) > errsum / fabs(area))
            goto L115;
    }

L110:
    // Test on divergence.
    if (ksgn == -1 && std::max(fabs(result_val), fabs(area)) <= defabs * .01) {
        goto L130;
    }
    if (.01 > result_val / area || result_val / area > 100. || errsum > fabs(area)) {
        ier = 5;
    }
    goto L130;

L115:
    // Compute global integral sum.
    result_val = 0.;
    for (k = 1; k <= last; ++k)
        result_val += rlist[k];
    abserr_val = errsum;

L130:
    if (ier > 2)
        /* L140: */ ;
    neval = last * 42 - 21;

    QuadResult res;
    res.value = result_val;
    res.error = abserr_val;
    res.neval = neval;
    res.converged = (ier == 0);
    return res;
}

} // namespace gauss_kronrod


// ============================================================
// SKAT-O: optimal.adj p-value
// ============================================================
// Implements SKAT:::SKAT_Optimal_PValue_Davies / optimal.adj
//
// The R SKAT package algorithm:
// 1. For each rho in grid, compute per-rho p-value via Davies/Liu
// 2. Find T_min = min(p-values across rho)
// 3. For each rho, compute tau_rho = threshold Q value such that
//    P(Q_rho > tau_rho) = T_min (i.e., tau = quantile at T_min)
// 4. Compute P(SKAT-O) = integral_0^inf f(x) * P(any rho exceeds tau | Z=x) dx
//    where f(x) is the chi2(1) density and Z is the leading eigenvalue component
// 5. The integration uses adaptive Gauss-Kronrod quadrature (like R's integrate())
//
// This implementation closely follows the R SKAT source.

double SKATO_optimal_pvalue(const arma::vec& Score,
                             const arma::mat& Phi,
                             const arma::vec& rho_vec,
                             const arma::vec& pval_each,
                             const arma::vec& q_each) {
    // ================================================================
    // Exact port of R's SKAT package SKAT-O optimal.adj algorithm:
    //   SKAT_META_Optimal -> SKAT_META_Optimal_Get_Pvalue
    //     -> SKAT_META_Optimal_Param (shared decomposition)
    //     -> SKAT_Optimal_Each_Q (per-rho p-values + tau_rho)
    //     -> SKAT_Optimal_PValue_Davies (integration)
    // ================================================================

    int n_rho = (int)rho_vec.n_elem;
    int m = (int)Score.n_elem;

    // Cap rho at 0.999 (R: r.corr[r.corr >= 0.999] <- 0.999)
    arma::vec rho_capped = rho_vec;
    for (int k = 0; k < n_rho; k++) {
        if (rho_capped(k) >= 0.999) rho_capped(k) = 0.999;
    }

    // Find minimum p-value across the rho grid
    double min_pval = arma::min(pval_each);
    if (!std::isfinite(min_pval) || min_pval <= 0.0) {
        return min_pval;
    }

    // Special case: single variant (m <= 1)
    if (m <= 1) {
        return min_pval;
    }

    // ================================================================
    // R passes Phi/2 internally: SKAT_META_Optimal calls
    //   SKAT_META_Optimal_Get_Pvalue(Q.all, Phi/2, r.all, ...)
    // Q.all = Q/2 where Q = Score' * R_rho * Score
    // Since p-values are scale-invariant, we work with Phi/2.
    // ================================================================
    arma::mat Phi_half = Phi * 0.5;
    int p_m = m;  // number of markers

    // ================================================================
    // SKAT_META_Optimal_Param(Phi_half, r.all)
    // Shared decomposition of Phi/2 into orthogonal components
    // ================================================================

    // Z.item1.1 = Phi_half %*% rep(1, p.m)  -- column sums vector
    arma::vec ones_m = arma::ones<arma::vec>(p_m);
    arma::vec Z_item1_1 = Phi_half * ones_m;

    // sum(Phi_half) = total of all elements
    double sum_Phi_half = arma::accu(Phi_half);

    // ZMZ = Z.item1.1 %*% t(Z.item1.1) / sum(Phi_half)
    arma::mat ZMZ;
    if (std::abs(sum_Phi_half) > 1e-300) {
        ZMZ = (Z_item1_1 * Z_item1_1.t()) / sum_Phi_half;
    } else {
        ZMZ = arma::zeros<arma::mat>(p_m, p_m);
    }

    // W3.2.t = Phi_half - ZMZ  (orthogonal component)
    arma::mat W3_2_t = Phi_half - ZMZ;

    // Get_Lambda(W3.2.t) with R's eigenvalue filter:
    //   lambda1 = eigen(K)$values
    //   IDX1 = which(lambda1 >= 0)
    //   IDX2 = which(lambda1 > mean(lambda1[IDX1]) / 1e5)
    //   lambda = lambda1[IDX2]
    arma::vec lambda_all;
    arma::eig_sym(lambda_all, W3_2_t);

    // R's Get_Lambda filter
    arma::vec lambda_pos = lambda_all(arma::find(lambda_all >= 0));
    double mean_pos = 0.0;
    if (lambda_pos.n_elem > 0) {
        mean_pos = arma::mean(lambda_pos);
    }
    double threshold = mean_pos / 1e5;
    arma::vec lambda_shared = lambda_all(arma::find(lambda_all > threshold));

    // W3.3.item = sum(ZMZ * (Phi_half - ZMZ)) * 4
    double W3_3_item = arma::accu(ZMZ % (Phi_half - ZMZ)) * 4.0;

    // z_mean_2 = sum(Phi_half) / p.m^2
    double z_mean_2 = sum_Phi_half / ((double)p_m * (double)p_m);

    // tau1 = sum(Phi_half %*% Phi_half) / p.m^2 / z_mean_2
    // Note: sum(A %*% B) = accu(A %*% B) in R = sum of all elements of the product
    // More efficiently: sum(Phi_half %*% Phi_half) = trace(Phi_half' * Phi_half * ones)
    // Actually R's sum(A) for matrix = accu(A), so sum(Phi_half %*% Phi_half) = accu(Phi_half * Phi_half)
    double sum_PhiPhi = arma::accu(Phi_half * Phi_half);
    double tau1 = 0.0;
    if (std::abs(z_mean_2) > 1e-300) {
        tau1 = sum_PhiPhi / ((double)p_m * (double)p_m) / z_mean_2;
    }

    // MuQ = sum(lambda_shared)
    double MuQ = arma::sum(lambda_shared);

    // VarQ = sum(lambda^2) * 2 + W3.3.item
    double VarQ = arma::sum(lambda_shared % lambda_shared) * 2.0 + W3_3_item;

    // VarRemain = W3.3.item  (the cross-term variance)
    // Actually from R code: param.m$VarRemain = W3.3.item
    double VarRemain = W3_3_item;

    // tau[i] = p.m^2 * r.corr[i] * z_mean_2 + tau1 * (1 - r.corr[i])
    arma::vec tau(n_rho);
    for (int k = 0; k < n_rho; k++) {
        tau(k) = (double)p_m * (double)p_m * rho_capped(k) * z_mean_2 + tau1 * (1.0 - rho_capped(k));
    }

    // ================================================================
    // SKAT_Optimal_Each_Q: per-rho p-values and tau_rho
    // ================================================================
    // For each rho, compute per-rho eigenvalues via Cholesky:
    //   R.M = (1-rho)*diag(m) + rho*matrix(1,m,m)
    //   L = chol(R.M)  (upper Cholesky)
    //   Phi_rho = L %*% Phi_half %*% t(L)
    //   lambda_temp = Get_Lambda(Phi_rho)
    //   c1[k] = sum(lambda_temp^k)
    //   param.temp = Get_Liu_Params_Mod(c1)
    //   pval = Get_PValue.Lambda(lambda_temp, Q)$p.value  (Davies with Liu fallback)
    //
    // After finding pmin = min(pval across rhos):
    //   q.org = qchisq(1 - pmin, df = l)  *** CENTRAL chi-squared, no ncp! ***
    //   q.q = (q.org - l)/sqrt(2*l) * sqrt(VarQ) + MuQ  (back to Q scale)

    // Q.all in R = Q/2 where Q = Score' * R_rho * Score
    // Q(rho) = (1-rho)*Score'*Score + rho*(sum(Score))^2
    // Q.all(rho) = Q(rho)/2
    double sum_score = arma::sum(Score);
    arma::vec Q_all(n_rho);
    for (int k = 0; k < n_rho; k++) {
        double Q_rho = (1.0 - rho_capped(k)) * arma::dot(Score, Score) + rho_capped(k) * sum_score * sum_score;
        Q_all(k) = Q_rho / 2.0;
    }

    // Per-rho eigenvalues and p-values
    std::vector<arma::vec> lambda_rho_vec(n_rho);
    arma::vec pval_rho(n_rho);
    std::vector<LiuParams> liu_rho_vec(n_rho);

    for (int k = 0; k < n_rho; k++) {
        // Build R.M = (1-rho)*I + rho*11'
        double rho = rho_capped(k);
        arma::mat R_M = (1.0 - rho) * arma::eye<arma::mat>(p_m, p_m) + rho * arma::ones<arma::mat>(p_m, p_m);

        // L = chol(R.M) -- upper Cholesky: R.M = L' * L (R convention: chol returns upper)
        arma::mat L;
        bool chol_ok = arma::chol(L, R_M);  // L is upper triangular
        if (!chol_ok) {
            // Fallback: use identity
            L = arma::eye<arma::mat>(p_m, p_m);
        }

        // Phi_rho = L %*% Phi_half %*% t(L)
        arma::mat Phi_rho = L * Phi_half * L.t();
        Phi_rho = 0.5 * (Phi_rho + Phi_rho.t());  // ensure symmetry

        // Get_Lambda with R's eigenvalue filter
        arma::vec ev_all;
        arma::eig_sym(ev_all, Phi_rho);

        arma::vec ev_pos = ev_all(arma::find(ev_all >= 0));
        double mean_pos_k = 0.0;
        if (ev_pos.n_elem > 0) {
            mean_pos_k = arma::mean(ev_pos);
        }
        double thresh_k = mean_pos_k / 1e5;
        arma::vec lambda_temp = ev_all(arma::find(ev_all > thresh_k));

        lambda_rho_vec[k] = lambda_temp;

        if (lambda_temp.n_elem == 0) {
            pval_rho(k) = 1.0;
            liu_rho_vec[k] = LiuParams{1.0, 0.0, 1.0, 1.0, 0.0};
            continue;
        }

        // Liu params from per-rho eigenvalues (needed for pmin_q mapping)
        liu_rho_vec[k] = liu_params(lambda_temp);

        // Per-rho p-value via Davies (with Liu fallback).
        // SAIGE uses method="optimal.adj", which triggers:
        //   pval[, i] <- Get_PValue.Lambda(lambda.temp, Q)$p.value
        // Get_PValue.Lambda uses Davies (acc=1e-6) with Liu fallback when
        // ifault != 0 or p-value out of (0,1].
        {
            double pv = davies_pvalue(Q_all(k), lambda_temp);
            if (std::isnan(pv) || pv <= 0.0 || pv > 1.0) {
                pv = liu_pvalue(Q_all(k), lambda_temp);
            }
            if (std::isnan(pv) || pv < 0.0) pv = 1.0;
            pval_rho(k) = pv;
        }
    }

    // pmin = min(pval_rho)
    double pmin = arma::min(pval_rho);

    // ================================================================
    // Compute pmin.q: the Q threshold corresponding to pmin for each rho
    // R code:
    //   q.org = qchisq(1 - pmin, df = param.mat[,6])  *** CENTRAL, no ncp ***
    //   q.q = (q.org - param.mat[,6])/sqrt(2*param.mat[,6]) * sqrt(VarQ) + MuQ
    //
    // param.mat[,6] = l (df) from Liu params for each rho
    // ================================================================
    arma::vec pmin_q(n_rho);
    for (int k = 0; k < n_rho; k++) {
        double l_k = liu_rho_vec[k].l;
        if (l_k <= 0.0) l_k = 1.0;

        try {
            // CENTRAL chi-squared quantile (no noncentrality parameter!)
            boost::math::chi_squared chi2_df(l_k);
            double q_org = boost::math::quantile(boost::math::complement(chi2_df, pmin));

            // Map back to Q scale using PER-RHO varQ and muQ
            // R: q.q = (q.org - df)/sqrt(2*df) * sqrt(varQ) + muQ
            // where varQ = param.mat[i,2], muQ = param.mat[i,1], df = param.mat[i,3]
            // These are the per-rho Liu parameters
            double denom = std::sqrt(2.0 * l_k);
            if (denom > 0.0) {
                double varQ_k = liu_rho_vec[k].sigmaQ * liu_rho_vec[k].sigmaQ;
                pmin_q(k) = (q_org - l_k) / denom * std::sqrt(varQ_k) + liu_rho_vec[k].muQ;
            } else {
                pmin_q(k) = Q_all(k);
            }
        } catch (...) {
            pmin_q(k) = Q_all(k);
        }
    }

    // ================================================================
    // SKAT_Optimal_PValue_Davies: integration
    // ================================================================
    // R's integrand (SKAT_Optimal_Integrate_Func_Davies):
    //   temp1 = tau %x% t(x)        # tau[k] * x for each rho and quadrature point
    //   temp = (pmin.q - temp1)/(1-r.all)   # adjusted threshold per rho
    //   temp.min = apply(temp, 2, min)       # min across rhos
    //   For each x:
    //     min1 = temp.min[i]
    //     min1.temp = min1 - MuQ
    //     sd1 = sqrt((VarQ - VarRemain)/VarQ)
    //     min1.st = min1.temp * sd1 + MuQ
    //     dav.re = SKAT_davies(min1.st, lambda_shared, acc=10^-6)
    //     temp = dav.re$Qq   # P(Q > q) = tail probability
    //     re[i] = (1 - temp) * dchisq(x[i], df=1)   # CDF * density

    double sd1 = 0.0;
    if (VarQ > 1e-300) {
        double ratio = (VarQ - VarRemain) / VarQ;
        if (ratio > 0.0) {
            sd1 = std::sqrt(ratio);
        }
    }

    boost::math::chi_squared chi2_1(1.0);

    // R early-out threshold: if (min1 > sum(param.m$lambda) * 10^4) temp <- 0
    double sum_lambda_shared = arma::sum(lambda_shared);

    // Track whether any Davies call fails inside the integrand
    bool davies_failed_in_integrand = false;

    auto integrand_func = [&](double x) -> double {
        if (x <= 0.0) return 0.0;

        // chi2(1) density at x
        double fx;
        try {
            fx = boost::math::pdf(chi2_1, x);
        } catch (...) {
            return 0.0;
        }
        if (fx < 1e-300) return 0.0;

        // Compute min threshold across rhos:
        //   temp = (pmin.q[k] - tau[k]*x) / (1 - rho[k])  for each k
        //   temp.min = min(temp)
        double temp_min = std::numeric_limits<double>::infinity();
        for (int k = 0; k < n_rho; k++) {
            double denom = 1.0 - rho_capped(k);
            if (std::abs(denom) < 1e-15) continue;  // skip rho~1
            double val = (pmin_q(k) - tau(k) * x) / denom;
            if (val < temp_min) temp_min = val;
        }

        if (!std::isfinite(temp_min)) return 0.0;

        // R early-out: if (min1 > sum(param.m$lambda) * 10^4) temp <- 0
        // tail_prob = 0 means CDF = 1, so return 1.0 * fx
        if (temp_min > sum_lambda_shared * 1e4) {
            return fx;  // CDF = 1, tail_prob = 0
        }

        // Standardize: min1.st = (min1 - MuQ) * sd1 + MuQ
        double min1_temp = temp_min - MuQ;
        double min1_st = min1_temp * sd1 + MuQ;

        // ONE Davies call with shared lambda (strict: NaN on ifault != 0)
        // In R, ifault != 0 calls stop() inside the integrand, which causes
        // try(integrate(...)) to fail entirely, falling back to Liu integration.
        // We replicate this by using davies_pvalue_strict() which returns NaN
        // on failure, then set the flag so the caller switches to Liu integration.
        double tail_prob;
        if (lambda_shared.n_elem > 0) {
            tail_prob = davies_pvalue_strict(min1_st, lambda_shared);
            if (std::isnan(tail_prob) || tail_prob < 0.0 || tail_prob > 1.0) {
                // Davies failed — flag for Liu fallback and abort integration
                davies_failed_in_integrand = true;
                return 0.0;  // Return 0 to not corrupt the integral; it will be discarded
            }
        } else {
            tail_prob = 0.0;
        }

        // re[i] = (1 - tail_prob) * dchisq(x, 1)  =  CDF * density
        double cdf_val = 1.0 - tail_prob;
        if (cdf_val < 0.0) cdf_val = 0.0;
        if (cdf_val > 1.0) cdf_val = 1.0;

        return cdf_val * fx;
    };

    // R: integrate(integrand, 0, 40, subdivisions=1000, abs.tol=10^-25)
    // R uses rel.tol = .Machine$double.eps^0.25 = 1.220703e-4
    gauss_kronrod::QuadResult qr = gauss_kronrod::integrate(
        integrand_func,
        0.0, 40.0,
        1.220703e-4,  // rel_tol (R's .Machine$double.eps^0.25)
        1e-25,        // abs_tol (R uses 10^-25)
        1000          // subdivisions (R's subdivisions=1000)
    );

    // ================================================================
    // Davies fail -> Liu integration fallback
    // ================================================================
    // R: if Davies fails (ifault != 0) inside the integrand, the try()
    // wrapper catches the stop() and falls back to SKAT_Optimal_PValue_Liu.
    // R's Liu fallback uses param.m$Df = 12/KerQ = sum(lambda^2)^2/sum(lambda^4)
    // with CENTRAL chi-squared (no noncentrality parameter):
    //   temp.q = (temp.min - MuQ)/sqrt(VarQ) * sqrt(2*Df) + Df
    //   re = pchisq(temp.q, df=Df) * dchisq(x, 1)
    // Note: R uses temp.min directly (NOT the sd1-standardized min1.st).
    if (davies_failed_in_integrand || !std::isfinite(qr.value) ||
        qr.value < 0.0 || qr.value > 1.0) {
        // Compute Df = sum(lambda^2)^2 / sum(lambda^4) for the shared lambda
        double sum_lam2 = arma::sum(lambda_shared % lambda_shared);
        double sum_lam4 = arma::sum(lambda_shared % lambda_shared % lambda_shared % lambda_shared);
        double Df_shared = (sum_lam4 > 0.0) ? (sum_lam2 * sum_lam2 / sum_lam4) : 1.0;

        auto liu_integrand_func = [&](double x) -> double {
            if (x <= 0.0) return 0.0;

            double fx;
            try {
                fx = boost::math::pdf(chi2_1, x);
            } catch (...) {
                return 0.0;
            }
            if (fx < 1e-300) return 0.0;

            double temp_min_l = std::numeric_limits<double>::infinity();
            for (int k = 0; k < n_rho; k++) {
                double denom = 1.0 - rho_capped(k);
                if (std::abs(denom) < 1e-15) continue;
                double val = (pmin_q(k) - tau(k) * x) / denom;
                if (val < temp_min_l) temp_min_l = val;
            }
            if (!std::isfinite(temp_min_l)) return 0.0;

            if (temp_min_l > sum_lambda_shared * 1e4) {
                return fx;
            }

            // R's Liu fallback: use temp.min directly (no sd1 standardization)
            // temp.q = (temp.min - MuQ) / sqrt(VarQ) * sqrt(2*Df) + Df
            double temp_q = (temp_min_l - MuQ) / std::sqrt(VarQ) * std::sqrt(2.0 * Df_shared) + Df_shared;

            double cdf_val_l;
            if (temp_q <= 0.0) {
                cdf_val_l = 0.0;
            } else {
                try {
                    boost::math::chi_squared chi2_df(Df_shared);
                    cdf_val_l = boost::math::cdf(chi2_df, temp_q);
                } catch (...) {
                    cdf_val_l = 0.0;
                }
            }

            if (cdf_val_l < 0.0) cdf_val_l = 0.0;
            if (cdf_val_l > 1.0) cdf_val_l = 1.0;

            return cdf_val_l * fx;
        };

        qr = gauss_kronrod::integrate(
            liu_integrand_func,
            0.0, 40.0,
            1.220703e-4,
            1e-25,
            2000          // R's SKAT_Optimal_PValue_Liu uses subdivisions=2000
        );
    }

    // R: pvalue = 1 - re[[1]]  (convert CDF integral to tail probability)
    double p_skato = 1.0 - qr.value;

    // R: if (pmin * length(r.all) < pvalue) pvalue = pmin * length(r.all)
    // Bonferroni bound
    double bonf = pmin * (double)n_rho;
    if (bonf < p_skato) {
        p_skato = bonf;
    }

    // ================================================================
    // Multi correction from SKAT_META_Optimal_Get_Pvalue:
    //   multi = 3; if (length(r.all) < 3) multi = 2
    //   pval1 = min(pval.each) * multi
    //   if (pval[i] <= 0 || any(pval.each <= 0)) pval[i] = pval1
    //   if (pval[i] == 0) pval[i] = min(pval.each[pval.each > 0])
    // ================================================================
    int multi = (n_rho < 3) ? 2 : 3;
    double pval1 = pmin * (double)multi;

    bool any_pval_le_zero = false;
    for (int k = 0; k < n_rho; k++) {
        if (pval_rho(k) <= 0.0) { any_pval_le_zero = true; break; }
    }

    if (p_skato <= 0.0 || any_pval_le_zero) {
        p_skato = pval1;
    }

    if (p_skato == 0.0) {
        // Find min of positive p-values
        double min_pos = std::numeric_limits<double>::infinity();
        for (int k = 0; k < n_rho; k++) {
            if (pval_rho(k) > 0.0 && pval_rho(k) < min_pos) {
                min_pos = pval_rho(k);
            }
        }
        if (std::isfinite(min_pos)) {
            p_skato = min_pos;
        }
    }

    // Clamp to valid range
    if (p_skato < 0.0) p_skato = 0.0;
    if (p_skato > 1.0) p_skato = 1.0;

    return p_skato;
}


// ============================================================
// get_SKAT_pvalue: Main entry point
// ============================================================
// Direct port of SAIGE's get_SKAT_pvalue R function.
// Calls SKAT/Burden/SKAT-O tests depending on regionTestType.
//
// Arguments:
//   Score: weighted score vector (m x 1), weights already applied
//   Phi: weighted variance-covariance matrix (m x m), weights already applied
//   r_corr: correlation parameter grid for SKAT-O
//   regionTestType: "SKAT", "BURDEN", or "SKATO"

SKATResult get_SKAT_pvalue(const arma::vec& Score,
                            const arma::mat& Phi,
                            const arma::vec& r_corr,
                            const std::string& regionTestType) {
    SKATResult result;
    result.pvalue_SKATO = std::numeric_limits<double>::quiet_NaN();
    result.pvalue_Burden = std::numeric_limits<double>::quiet_NaN();
    result.pvalue_SKAT = std::numeric_limits<double>::quiet_NaN();
    result.beta_Burden = std::numeric_limits<double>::quiet_NaN();
    result.se_Burden = std::numeric_limits<double>::quiet_NaN();

    int m = (int)Score.n_elem;
    if (m == 0) return result;

    // ----------------------------------------------------------------
    // Burden effect size: BETA_Burden = sum(Score) / sum(diag(Phi))
    // This matches: BETA_Burden = sum(Score)/(sum(diag(Phi))) in R
    // ----------------------------------------------------------------
    double sum_score = arma::sum(Score);
    double sum_diag_phi = arma::trace(Phi);
    if (sum_diag_phi > 0.0) {
        result.beta_Burden = sum_score / sum_diag_phi;
    }

    // ----------------------------------------------------------------
    // SKAT test (rho = 0): Q_SKAT = Score' * Score
    // Under H0, Q_SKAT ~ sum(lambda_j * chi2(1))
    // where lambda_j = eigenvalues of Phi
    // ----------------------------------------------------------------
    {
        double Q_SKAT = arma::dot(Score, Score);
        arma::vec eigvals;
        bool ok = arma::eig_sym(eigvals, Phi);
        if (ok) {
            arma::vec lam = filter_eigenvalues_R(eigvals);

            if (lam.n_elem > 0) {
                result.pvalue_SKAT = davies_pvalue(Q_SKAT, lam);
                if (std::isnan(result.pvalue_SKAT) || result.pvalue_SKAT < 0.0) {
                    result.pvalue_SKAT = liu_pvalue(Q_SKAT, lam);
                }
            } else {
                result.pvalue_SKAT = 1.0;
            }
        }
    }

    // ----------------------------------------------------------------
    // Burden test (rho = 1): Q_Burden = (sum(Score))^2
    // Under H0, Q_Burden / sum(Phi) ~ chi2(1)
    // ----------------------------------------------------------------
    {
        double Q_Burden = sum_score * sum_score;
        double Var_Burden = arma::accu(Phi);  // sum of all elements

        if (Var_Burden > 0.0) {
            try {
                boost::math::chi_squared chi2_1(1.0);
                result.pvalue_Burden = boost::math::cdf(
                    boost::math::complement(chi2_1, Q_Burden / Var_Burden));
            } catch (...) {
                result.pvalue_Burden = std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            result.pvalue_Burden = 1.0;
        }
    }

    // ----------------------------------------------------------------
    // Burden SE: SE = |BETA| / |qnorm(p/2)|
    // Matches R: SE_Burden = abs(BETA_Burden/qnorm(p_burden/2))
    // ----------------------------------------------------------------
    if (std::isfinite(result.pvalue_Burden) && result.pvalue_Burden > 0.0 &&
        std::isfinite(result.beta_Burden)) {
        try {
            boost::math::normal norm(0, 1);
            double z = boost::math::quantile(norm, result.pvalue_Burden / 2.0);
            if (z != 0.0) {
                result.se_Burden = std::abs(result.beta_Burden / z);
            }
        } catch (...) {
            result.se_Burden = std::numeric_limits<double>::quiet_NaN();
        }
    }

    // ----------------------------------------------------------------
    // SKAT-O or single test dispatch
    // ----------------------------------------------------------------
    if (regionTestType == "BURDEN") {
        result.pvalue_SKATO = result.pvalue_Burden;
    } else if (regionTestType == "SKAT") {
        result.pvalue_SKATO = result.pvalue_SKAT;
    } else {
        // SKAT-O: optimal rho search

        // R rank check: if (ncol(Phi) <= 10) { if (qr(Phi)$rank <= 1) r.corr = 0 }
        // Forces pure SKAT for rank-1 Phi (skip SKAT-O integration)
        bool force_skat_only = false;
        if (m <= 10) {
            arma::uword phi_rank = arma::rank(Phi);
            if (phi_rank <= 1) {
                force_skat_only = true;
            }
        }

        if (force_skat_only) {
            // r.corr = 0 means pure SKAT
            result.pvalue_SKATO = result.pvalue_SKAT;
        } else {

        // Default rho grid: matches SKAT:::SKAT_Check_Method("optimal.adj", 0)
        // R returns: {0, 0.01, 0.04, 0.09, 0.25, 0.5, 1}
        arma::vec rho_vec;
        if (r_corr.n_elem > 0) {
            rho_vec = r_corr;
        } else {
            rho_vec = arma::vec({0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0});
        }

        int n_rho = (int)rho_vec.n_elem;
        arma::vec pval_each(n_rho);
        arma::vec q_each(n_rho);

        // R caps rho >= 0.999 to 0.999 inside SKAT_META_Optimal, then
        // computes per-rho p-values via Cholesky + eigenvalues + Davies.
        // The Burden p-value (rho=1 slot) is computed at rho=0.999, NOT
        // via the exact chi2(1) formula. We replicate this here.
        arma::mat Phi_half = Phi * 0.5;
        for (int k = 0; k < n_rho; k++) {
            double rho = rho_vec(k);
            // Cap rho >= 0.999 to 0.999, matching R's SKAT_META_Optimal
            double rho_capped = (rho >= 0.999) ? 0.999 : rho;

            // Q(rho) = (1-rho_capped) * Score'*Score + rho_capped * (sum(Score))^2
            double Q_rho = (1.0 - rho_capped) * arma::dot(Score, Score) + rho_capped * sum_score * sum_score;
            q_each(k) = Q_rho;

            if (std::abs(rho_capped) < 1e-10) {
                pval_each(k) = result.pvalue_SKAT;
            } else {
                // For all rho > 0 (including rho=1 capped to 0.999):
                // R.M = (1-rho)*I + rho*11'
                // L = chol(R.M)
                // Phi_rho = L * (Phi/2) * L'
                // eigenvalues -> Davies/Liu p-value
                arma::mat R_M = (1.0 - rho_capped) * arma::eye<arma::mat>(m, m)
                              + rho_capped * arma::ones<arma::mat>(m, m);
                arma::mat L;
                bool chol_ok = arma::chol(L, R_M);  // upper Cholesky
                if (!chol_ok) {
                    L = arma::eye<arma::mat>(m, m);
                }
                arma::mat Phi_rho = L * Phi_half * L.t();
                Phi_rho = 0.5 * (Phi_rho + Phi_rho.t());  // ensure symmetry

                arma::vec ev_all;
                arma::eig_sym(ev_all, Phi_rho);
                arma::vec ev_rho = filter_eigenvalues_R(ev_all);

                if (ev_rho.n_elem > 0) {
                    // Q.all in R = Q/2
                    double Q_half = Q_rho / 2.0;
                    pval_each(k) = davies_pvalue(Q_half, ev_rho);
                    if (std::isnan(pval_each(k)) || pval_each(k) < 0.0) {
                        pval_each(k) = liu_pvalue(Q_half, ev_rho);
                    }
                } else {
                    pval_each(k) = 1.0;
                }
            }

            if (std::isnan(pval_each(k))) {
                pval_each(k) = 1.0;
            }
        }

        // Check if both rho=0 and rho=1 are in the grid
        bool has_rho0 = false, has_rho1 = false;
        for (int k = 0; k < n_rho; k++) {
            if (std::abs(rho_vec(k)) < 1e-10) has_rho0 = true;
            if (std::abs(rho_vec(k) - 1.0) < 1e-10) has_rho1 = true;
        }

        if (!has_rho0 || !has_rho1) {
            // Edge case: grid does not contain both endpoints
            double min_pval = arma::min(pval_each);
            result.pvalue_SKATO = min_pval;
            result.pvalue_SKAT = min_pval;
            result.pvalue_Burden = min_pval;
        } else {
            // Compute optimal SKAT-O p-value
            result.pvalue_SKATO = SKATO_optimal_pvalue(Score, Phi, rho_vec, pval_each, q_each);

            // Extract SKAT and Burden p-values from the grid
            for (int k = 0; k < n_rho; k++) {
                if (std::abs(rho_vec(k)) < 1e-10) {
                    result.pvalue_SKAT = pval_each(k);
                }
                if (std::abs(rho_vec(k) - 1.0) < 1e-10) {
                    result.pvalue_Burden = pval_each(k);
                }
            }
        }

        } // end force_skat_only else

        // Recompute SE_Burden with final burden p-value
        if (std::isfinite(result.pvalue_Burden) && result.pvalue_Burden > 0.0 &&
            std::isfinite(result.beta_Burden)) {
            try {
                boost::math::normal norm(0, 1);
                double z = boost::math::quantile(norm, result.pvalue_Burden / 2.0);
                if (z != 0.0) {
                    result.se_Burden = std::abs(result.beta_Burden / z);
                }
            } catch (...) {
                result.se_Burden = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    return result;
}
