// erfc_port_hosttest -- compile erfc_boost53.cuh as host C++ and compare both
// ports with boost::math::erfc on the same machine, same libm. Boost evaluates
// erfc(double) in long double (promote_double), so bit-identity is not
// expected; what is checked:
//   * normal range: max ulp difference and the share of bit-identical values;
//   * subnormal range (26.55 <= z < 27.226): values differ by how many
//     subnormal grid units, and how often;
//   * the zero cutoff: arguments where exactly one side is 0.
// On the device only exp() / ldexp() can additionally differ.
//
//   g++ -std=c++17 -O2 -ffp-contract=off -I. -isystem $CONDA_PREFIX/include erfc_port_hosttest.cpp -o build/erfc_port_hosttest
#include <boost/math/special_functions/erf.hpp>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "erfc_boost53.cuh"

using saige::spa_gpu::erfc53::erfImp53;
using saige::spa_gpu::erfc53::erfImp53s;

static double ulps(double a, double b)
{
    if (a == b) return 0;
    if (!std::isfinite(a) || !std::isfinite(b) || a == 0 || b == 0) return std::numeric_limits<double>::infinity();
    return std::fabs(a - b) / std::fabs(std::nexttoward(a, b) - a);
}
static double gridUnits(double a, double b)   // in units of the smallest subnormal
{
    return std::fabs(a - b) / std::numeric_limits<double>::denorm_min();
}

struct Stat {
    const char* name; long n = 0, exact = 0, zeroMis = 0; double maxUlp = 0, zAtMaxUlp = 0;
    long subN = 0, subExact = 0; double subMaxGrid = 0; double firstZero = NAN;
};

static void feed(Stat& s, double z, double host, double port)
{
    s.n++;
    const bool sub = std::isfinite(host) && host != 0 && std::fabs(host) < std::numeric_limits<double>::min();
    if (host == port) s.exact++;
    if ((host == 0) != (port == 0)) s.zeroMis++;
    if (sub || (port != 0 && std::fabs(port) < std::numeric_limits<double>::min())) {
        s.subN++;
        if (host == port) s.subExact++;
        s.subMaxGrid = std::max(s.subMaxGrid, gridUnits(host, port));
    } else if (host != port && host != 0 && port != 0) {
        const double u = ulps(host, port);
        if (u > s.maxUlp) { s.maxUlp = u; s.zAtMaxUlp = z; }
    }
    if (std::isnan(s.firstZero) && port == 0 && z > 0) s.firstZero = z;
}

int main()
{
    std::vector<double> z;
    for (double x = -6.0; x <= 26.5; x += 1e-4) z.push_back(x);
    for (double x = 26.5; x <= 28.2; x += 1e-6) z.push_back(x);      // the subnormal zone, finely
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> U(-6, 28);
    for (int i = 0; i < 1000000; ++i) z.push_back(U(rng));
    for (double b : {0.0, 1e-10, 0.5, 1.5, 2.5, 4.5, 5.93, (double)5.93f, 28.0, -0.5}) {
        z.push_back(b); z.push_back(std::nextafter(b, 100.0)); z.push_back(std::nextafter(b, -100.0));
    }
    Stat sS{"scaled port (mode 1)"}, sL{"literal port (mode 2)"}, eS{"erf scaled"}, eL{"erf literal"};
    double hostFirstZero = NAN;
    for (double x : z) {
        const double hc = boost::math::erfc(x), he = boost::math::erf(x);
        feed(sS, x, hc, erfImp53s(x, true));
        feed(sL, x, hc, erfImp53(x, true));
        feed(eS, x, he, erfImp53s(x, false));
        feed(eL, x, he, erfImp53(x, false));
        if (std::isnan(hostFirstZero) && hc == 0 && x > 0) hostFirstZero = x;
    }
    std::printf("erfc ports, host build vs boost::math::erfc (Boost 1.85; it evaluates in long double and narrows): %zu arguments in [-6, 28.2]\n", z.size());
    std::printf("  host erfc first 0 at z = %.7f  (true erfc there ~ %.3Le; half the smallest subnormal is %.3Le)\n",
                hostFirstZero, std::exp(-(long double)hostFirstZero * hostFirstZero) / (hostFirstZero * 1.7724538509055159L), (long double)std::numeric_limits<double>::denorm_min() / 2);
    for (const Stat* s : {&sS, &sL, &eS, &eL})
        std::printf("  %-22s bit-identical %ld/%ld; normal range max ulp %.3g (z=%.6f); subnormal zone n=%ld identical=%ld max grid units off=%.0f; zero/non-zero disagreements=%ld; first 0 at z=%.7f\n",
                    s->name, s->exact, s->n, s->maxUlp, s->zAtMaxUlp, s->subN, s->subExact, s->subMaxGrid, s->zeroMis, s->firstZero);
    const double inf = std::numeric_limits<double>::infinity();
    std::printf("  erfc(+inf) host %g scaled %g; erfc(-inf) host %g scaled %g; erfc(27.2) host %.6g scaled %.6g literal %.6g; erfc(27.226) host %.6g scaled %.6g\n",
                boost::math::erfc(inf), erfImp53s(inf, true), boost::math::erfc(-inf), erfImp53s(-inf, true),
                boost::math::erfc(27.2), erfImp53s(27.2, true), erfImp53(27.2, true), boost::math::erfc(27.226), erfImp53s(27.226, true));
    return 0;
}
