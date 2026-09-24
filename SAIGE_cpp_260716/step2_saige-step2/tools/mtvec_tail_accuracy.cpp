// mtvec_tail_accuracy.cpp -- how far score_vec.hpp's chi-square(1) upper tail
// is from the boost call format_score_result makes, and how much each costs.
//
//   make tools/mtvec_tail_accuracy
//   tools/mtvec_tail_accuracy                 synthetic sweep + timings
//   tools/mtvec_tail_accuracy stat.f64        also the same table over the
//                                             stat values of a real run
//                                             (tools/mtvec_cmp_sgs.py --dump-stat)
//
// The reference is exactly what the scalar path evaluates:
//   boost::math::cdf(complement(boost::math::chi_squared(1), stat))
// with boost's default policy, which promotes double to long double
// internally. The candidate is score_vec.hpp's erfc(sqrt(stat/2)) in double.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>

#include "../score_vec.hpp"

static inline double boostTail(double stat)
{
    boost::math::chi_squared d(1);
    return boost::math::cdf(complement(d, stat));
}

static inline double vecTail(double stat)
{
    double z = std::sqrt(stat * 0.5), p;
    SAIGE::mtVecErfc(&z, &p, 1);
    return p;
}

static double ulps(double a, double b)
{
    if (a == b) return 0.0;
    if (!(std::isfinite(a) && std::isfinite(b))) return 1e300;
    int64_t ia, ib;
    std::memcpy(&ia, &a, 8);
    std::memcpy(&ib, &b, 8);
    return std::fabs(static_cast<double>(ia - ib));
}

struct Row { double hi, lo; long n, nfall; double relmax, ulpmax, dlpmax; };

static void band(const std::vector<double>& t_stats, const char* t_title)
{
    const double cut = SAIGE::mtVecStatCutoff();
    std::vector<Row> rows;
    for (int e = 0; e >= -320; e -= 10)
        rows.push_back({std::pow(10.0, e), std::pow(10.0, e - 10), 0, 0, 0, 0, 0});
    Row below{0, 0, 0, 0, 0, 0, 0};   // p == 0 exactly (the "%.1fE%d" branch)

    for (double s : t_stats) {
        if (!(s >= 0.0) || !std::isfinite(s)) continue;
        const double ref = boostTail(s);
        const bool   fall = !(s < cut);
        Row* r = nullptr;
        if (ref == 0.0) r = &below;
        else for (auto& b : rows) if (ref <= b.hi && ref > b.lo) { r = &b; break; }
        if (!r) continue;
        r->n++;
        if (fall) { r->nfall++; continue; }   // boost recomputes it: delta is 0
        const double got = vecTail(s);
        const double rel = ref > 0 ? std::fabs(got - ref) / ref : 0.0;
        r->relmax = std::max(r->relmax, rel);
        r->ulpmax = std::max(r->ulpmax, ulps(got, ref));
        const double dlp = (ref > 0 && got > 0)
                               ? std::fabs(std::log10(got) - std::log10(ref)) : 0.0;
        r->dlpmax = std::max(r->dlpmax, dlp);
    }
    printf("\n# %s   (fallback cutoff: stat >= %.5f, p < %.0e)\n",
           t_title, cut, SAIGE::MT_VEC_EXACT_BELOW_P);
    printf("%-10s %-10s %10s %10s %12s %10s %12s\n",
           "p <=", "p >", "n", "->boost", "max rel err", "max ulp", "max |dlog10p|");
    for (auto& b : rows) {
        if (b.n == 0) continue;
        printf("%-10.0e %-10.0e %10ld %10ld %12.3e %10.0f %12.3e\n",
               b.hi, b.lo, b.n, b.nfall, b.relmax, b.ulpmax, b.dlpmax);
    }
    if (below.n)
        printf("%-10s %-10s %10ld %10ld %12s %10s %12s\n",
               "p == 0", "", below.n, below.nfall, "-", "-", "-");
}

int main(int argc, char** argv)
{
#ifdef SAIGE_USE_MVEC_ERFC
    printf("erfc: libmvec _ZGVdN4v_erfc (4 wide)\n");
#else
    printf("erfc: scalar std::erfc (libmvec probe did not fire)\n");
#endif
#ifdef SAIGE_VEC_TO_CHARS
    printf("\"%%.6E\": std::to_chars(scientific, 6)\n");
#else
    printf("\"%%.6E\": snprintf (no <charconv> float support)\n");
#endif

    // ---- synthetic sweep: stat chosen so p lands on a log grid ----
    std::vector<double> syn;
    for (double lp = 0.0; lp >= -330.0; lp -= 0.002) {
        const double p = std::pow(10.0, lp);
        double stat;
        if (p > 1e-300) {
            boost::math::chi_squared d(1);
            stat = boost::math::quantile(complement(d, p));
        } else {
            const double L = -lp * std::log(10.0);
            double z = std::sqrt(L);
            for (int it = 0; it < 60; ++it)
                z = std::sqrt(L - std::log(z) - 0.5 * std::log(M_PI));
            stat = 2.0 * z * z;
        }
        if (std::isfinite(stat) && stat > 0) syn.push_back(stat);
    }
    std::sort(syn.begin(), syn.end());
    syn.erase(std::unique(syn.begin(), syn.end()), syn.end());
    band(syn, "synthetic sweep, 500 stat values per decade of p");

    // ---- the stat values of a real run, if one was handed over ----
    if (argc > 1) {
        FILE* f = std::fopen(argv[1], "rb");
        if (!f) { std::perror(argv[1]); return 1; }
        std::fseek(f, 0, SEEK_END);
        const long nb = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        std::vector<double> v(nb / 8);
        if (std::fread(v.data(), 8, v.size(), f) != v.size()) { std::fclose(f); return 1; }
        std::fclose(f);
        band(v, (std::string("run stat values from ") + argv[1]).c_str());
    }

    // ---- what each piece costs ----
    const long N = 5000000;
    std::mt19937_64 rng(12345);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> in(N), z(N), out(N);
    for (long i = 0; i < N; ++i) { const double g = nd(rng); in[i] = g * g; }
    double sink = 0.0;
    auto t = [&](const char* name, auto fn) {
        const auto t0 = std::chrono::steady_clock::now();
        fn();
        const auto t1 = std::chrono::steady_clock::now();
        for (long i = 0; i < N; i += 9973) sink += out[i];
        printf("  %-40s %8.2f ns/elem\n", name,
               std::chrono::duration<double, std::nano>(t1 - t0).count() / N);
    };
    printf("\n# cost, %ld draws of stat ~ chi-square(1), single thread\n", N);
    t("boost cdf(complement(chi_squared(1)))",
      [&] { for (long i = 0; i < N; ++i) out[i] = boostTail(in[i]); });
    t("sqrt(stat/2) then mtVecErfc over the array", [&] {
        for (long i = 0; i < N; ++i) z[i] = std::sqrt(in[i] * 0.5);
        SAIGE::mtVecErfc(z.data(), out.data(), static_cast<int>(N));
    });
    {
        char buf[64];
        long acc = 0;
        const auto t0 = std::chrono::steady_clock::now();
        for (long i = 0; i < N; ++i) acc += std::snprintf(buf, sizeof buf, "%.6E", out[i]);
        const auto t1 = std::chrono::steady_clock::now();
        for (long i = 0; i < N; ++i) acc += SAIGE::mtVecFormatE6(out[i], buf);
        const auto t2 = std::chrono::steady_clock::now();
        printf("  %-40s %8.2f ns/elem\n", "snprintf(\"%.6E\")",
               std::chrono::duration<double, std::nano>(t1 - t0).count() / N);
        printf("  %-40s %8.2f ns/elem\n", "mtVecFormatE6",
               std::chrono::duration<double, std::nano>(t2 - t1).count() / N);
        sink += acc;
    }
    printf("  (sink %.6g)\n", sink);

    // ---- mtVecFormatE6 must be snprintf byte for byte ----
    printf("\n# mtVecFormatE6 vs snprintf(\"%%.6E\")\n");
    std::mt19937_64 r2(99);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    long bad = 0, ntry = 0;
    char a[64], b[64];
    for (long i = 0; i < 50000000; ++i) {
        const double x = std::pow(10.0, -330.0 * u(r2));
        std::snprintf(a, sizeof a, "%.6E", x);
        const int len = SAIGE::mtVecFormatE6(x, b);
        b[len] = 0;
        ntry++;
        if (std::strcmp(a, b) != 0) {
            if (bad < 5) printf("  MISMATCH %.20g  snprintf=%s  vec=%s\n", x, a, b);
            bad++;
        }
    }
    const double fixed[] = {0.0, 1.0, 0.1, 1e-5, 9.9999999e-6, 5e-324, 1e-308, 1e-300,
                            0.9999999999, 1.0000005e-3, 1.4999995e-7};
    for (double x : fixed) {
        std::snprintf(a, sizeof a, "%.6E", x);
        const int len = SAIGE::mtVecFormatE6(x, b);
        b[len] = 0;
        ntry++;
        if (std::strcmp(a, b) != 0) { printf("  MISMATCH %.20g  %s vs %s\n", x, a, b); bad++; }
    }
    printf("  %ld values, %ld mismatches\n", ntry, bad);
    return bad ? 1 : 0;
}
