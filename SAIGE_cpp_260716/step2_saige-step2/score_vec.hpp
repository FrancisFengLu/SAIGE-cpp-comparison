// score_vec.hpp -- block-at-a-time chi-square(1) upper tail and "%.6E" for the
// quantitative multi-trait batch path (config key mtVecQuantStats).
//
// WHAT THIS REPLACES
//
// format_score_result (score_format.hpp) turns one (marker, trait) pair into
// its printed row.  At M2 = 10^6 markers x P = 128 traits that is 1.28e8 calls,
// each of which runs boost::math::cdf(complement(chi_squared(1), stat)) and one
// sprintf("%.6E").  Measured on this machine (micro/erfc_study.cpp, 5e6 draws
// of stat ~ chi-square(1), single thread, -O3 -march=native):
//
//   boost cdf(complement(chi_squared(1), stat))   538.5 ns
//   std::erfc(sqrt(stat/2))                        18.4 ns
//   boost::math::erfc(sqrt(stat/2))                69.9 ns
//   libmvec _ZGVdN4v_erfc, 4 wide                   5.1 ns
//   snprintf("%.6E")                              258.9 / 280.2 ns
//   std::to_chars(scientific, 6)                   83.3 ns
//
// The chi-square(1) upper tail IS erfc: P(X > stat) = erfc(sqrt(stat/2)).  So
// the 538 ns is not the special function, it is boost's incomplete-gamma
// machinery plus its default promote_double<true> policy, which evaluates the
// whole thing in long double.
//
// WHY THERE IS STILL A BOOST FALLBACK
//
// That long-double promotion is not free accuracy theatre: erfc(z) ~
// exp(-z^2)/(z sqrt(pi)), so a double erfc loses relative precision roughly in
// proportion to z^2 = stat/2, while boost's does not.  Measured against boost
// over stat values chosen to land in each decade of p (micro/erfc_study.cpp):
//
//   p in            max rel err of std::erfc / libmvec vs boost
//   [1e+0 , 1e-10]  4.2e-15
//   [1e-10, 1e-20]  5.9e-15
//   [1e-100,1e-110] 2.8e-14
//   [1e-290,1e-300] 9.3e-14
//
// and the two disagree about WHERE the tail underflows to exactly 0 (boost at
// stat = 1482.5168, double erfc at stat = 1482.5537), which is the test
// format_score_result uses to switch to the "%.1fE%d" log form.
//
// So the block path keeps the vectorised erfc only where it is provably
// indistinguishable and hands everything else back to format_score_result:
// every pair with stat >= qchisq(1e-5, 1, lower=FALSE) = 19.5114 -- i.e. every
// p-value below 1e-5, which is the entire range a GWAS reader cares about --
// is computed by boost itself and is bit-identical to a switch-off run.  Above
// that the two agree to <= 4.2e-15 relative, i.e. |d(-log10 p)| <= 1.8e-15.
// Under the null 1e-5 of pairs take the fallback; at 538 ns each that is
// 5.4 us per million pairs, which does not show up.
//
// THE STRING IS NOT VECTORISABLE
//
// "%.6E" is per element by construction.  std::to_chars(chars_format::scientific,
// 6) is specified to produce printf's "e" format and libstdc++ implements it
// with Ryu instead of the printf machinery: 3.4x faster and byte-identical
// (verified over 2e7 values spanning p = 1 down to 1e-320, subnormals
// included, tests/mtvec_fmt_check.cpp).  Where <charconv> has no floating
// point support (libc++ before 14) this falls back to snprintf.

#ifndef SAIGE_SCORE_VEC_HPP
#define SAIGE_SCORE_VEC_HPP

#include <cmath>
#include <cstdio>

#include <boost/math/distributions/chi_squared.hpp>

#if defined(__has_include)
#  if __has_include(<charconv>)
#    include <charconv>
#  endif
#endif
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#  define SAIGE_VEC_TO_CHARS 1
#endif

// libmvec's 4-wide double erfc.  glibc 2.35+ on x86-64; the Makefile probes
// for it (tools/mvec_probe.cpp) and only then defines SAIGE_HAVE_MVEC_ERFC and
// links -lmvec.  Without it the loop below is the scalar libm erfc, which is
// 18.4 ns instead of 5.1 ns -- 13 ns out of a ~1500 ns per-pair tail, so the
// probe failing costs under 1%.
#if defined(SAIGE_HAVE_MVEC_ERFC) && defined(__x86_64__) && defined(__AVX2__)
#  include <immintrin.h>
extern "C" __m256d _ZGVdN4v_erfc(__m256d);
#  define SAIGE_USE_MVEC_ERFC 1
#endif

namespace SAIGE {

// Pairs whose p-value would land below this are recomputed with boost, so the
// tail stays bit-identical to a switch-off run.  See the header comment.
constexpr double MT_VEC_EXACT_BELOW_P = 1e-5;

// stat >= this goes to boost.  qchisq(1e-5, df=1, lower.tail=FALSE) = 19.51142.
inline double mtVecStatCutoff()
{
    static const double c = boost::math::quantile(
        complement(boost::math::chi_squared(1), MT_VEC_EXACT_BELOW_P));
    return c;
}

// out[i] = erfc(z[i]).  Contiguous, no aliasing.
inline void mtVecErfc(const double* t_z, double* t_out, int t_n)
{
    int i = 0;
#ifdef SAIGE_USE_MVEC_ERFC
    for (; i + 4 <= t_n; i += 4)
        _mm256_storeu_pd(t_out + i, _ZGVdN4v_erfc(_mm256_loadu_pd(t_z + i)));
#endif
    for (; i < t_n; ++i) t_out[i] = std::erfc(t_z[i]);
}

// "%.6E" of t_p into t_buf (which must hold at least 32 bytes). Returns the
// length; the buffer is NOT null terminated.
inline int mtVecFormatE6(double t_p, char* t_buf)
{
#ifdef SAIGE_VEC_TO_CHARS
    const std::to_chars_result r =
        std::to_chars(t_buf, t_buf + 32, t_p, std::chars_format::scientific, 6);
    for (char* q = t_buf; q != r.ptr; ++q)
        if (*q == 'e') { *q = 'E'; break; }
    return static_cast<int>(r.ptr - t_buf);
#else
    return std::snprintf(t_buf, 32, "%.6E", t_p);
#endif
}

}  // namespace SAIGE

#endif  // SAIGE_SCORE_VEC_HPP
