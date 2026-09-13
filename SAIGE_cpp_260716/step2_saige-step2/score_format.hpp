// score_format.hpp — shared chi-square(1) p-value formatting for the score test.
//
// Pure extraction (2026-09-13, multi-trait Phase 0): the ~30 lines below were
// duplicated verbatim in SAIGEClass::scoreTest, ::scoreTestFast and
// ::scoreTestFast_noadjCov. Not one operator is changed; the three call sites
// now share one definition so the scalar path and the future multi-trait batch
// kernel (scoreTestBatchMT) cannot drift apart in p-value formatting.
//
// Deliberately NOT applied to scoreTestFast_block / *_fused: those two set
// stat = 0.0 on the degenerate-variance branch (var1 <= DBL_MIN) whereas the
// three scalar functions leave stat untouched there, which changes seBeta.
// Unifying that difference would be a behaviour change, not an extraction, so
// it is left alone (see MULTITRAIT_DESIGN.md section 5).

#ifndef SAIGE_SCORE_FORMAT_HPP
#define SAIGE_SCORE_FORMAT_HPP

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

namespace SAIGE {

// P1 fix (2026-05-09): direct log P(X > stat) for chi-square(df=1) without
// underflow. The previous code computed log(boost::math::cdf(...)) which
// underflows to -inf when the cdf itself underflows below ~1e-300 (i.e. at
// genuinely significant markers). This matches R's pchisq(q, 1, lower.tail=FALSE,
// log.p=TRUE).
//
// chi-square(df=1) upper tail: P(X > stat) = erfc(sqrt(stat/2))
// Boost provides erfc; for very large argument it underflows too, so use
// the asymptotic expansion log(erfc(z)) ≈ -z² - log(z) - 0.5*log(π) for z > 6.
inline double log_chisq1_uppertail(double stat) {
    if (!std::isfinite(stat) || stat <= 0.0) return 0.0;  // log(p=1) = 0
    const double z = std::sqrt(stat / 2.0);
    if (z > 6.0) {
        // first-order asymptotic for erfc(z), accurate to ~1e-3 already at z=4
        // and machine precision at z=10+.
        return -z * z - std::log(z) - 0.5 * std::log(boost::math::constants::pi<double>());
    }
    const double e = boost::math::erfc(z);
    if (e <= 0.0 || !std::isfinite(e)) {
        // fallback: same asymptotic if boost's erfc itself returned 0
        return -z * z - std::log(z) - 0.5 * std::log(boost::math::constants::pi<double>());
    }
    return std::log(e);
}

// One (marker, trait) score-test result: turn (S, var1, var2) into the
// p-value string / Beta / seBeta the output writer expects.
//
// Byte-for-byte the code that used to live inline in scoreTestFast:
//   - var1 <= DBL_MIN          -> pval = 1, stat left as computed
//   - stat not finite          -> pval = 1, stat = 0
//   - pval != 0                -> "%.6E"        , islogp = false
//   - pval == 0 (underflow)    -> "%.1fE%d" of the log10 upper tail,
//                                 pval replaced by the natural log, islogp = true
inline void format_score_result(double S, double var1, double var2,
                                double& t_Beta,
                                double& t_seBeta,
                                std::string& t_pval_str,
                                double& t_pval,
                                bool& t_islogp,
                                double& t_Tstat,
                                double& t_var1,
                                double& t_var2)
{
    double stat = S * S / var1;
    if (var1 <= std::numeric_limits<double>::min()) {
        t_pval = 1;
    } else {
        if (!std::isnan(stat) && std::isfinite(stat)) {
            boost::math::chi_squared chisq_dist(1);
            t_pval = boost::math::cdf(complement(chisq_dist, stat));
        } else {
            t_pval = 1;
            stat = 0.0;
        }
    }
    char pValueBuf[100];
    if (t_pval != 0) {
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    } else {
        // P1 fix: direct log(upper-tail) avoiding cdf underflow at p<1e-300
        double logp = log_chisq1_uppertail(stat);
        double log10p = logp / (std::log(10.0));
        int exponent = std::floor(log10p);
        double fraction = std::pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
            fraction = 1;
            exponent++;
        }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp;
        t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S / var1;
    t_seBeta = std::fabs(t_Beta) / std::sqrt(std::fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}

}  // namespace SAIGE

#endif  // SAIGE_SCORE_FORMAT_HPP
