// er_gpu_check.cpp — the device ER (gpu/gpu_er.hpp) against er_binary.cpp.
//
//   tools/er_gpu_check math [nRandom]        device exp / log vs the host's, bit by bit
//   tools/er_gpu_check er   [nPairs] [seed]  device ER p vs SKATExactBin_Work with the
//                                            production arguments, bit by bit
//
// Exit status 0 only when every compared value is bit-identical.
#include "er_binary.hpp"
#include "gpu_er.hpp"

#include <armadillo>
#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

using saige::gpu2::ErPairIn;
using saige::gpu2::ErPairOut;

static uint64_t bits(double d) { uint64_t u; std::memcpy(&u, &d, 8); return u; }
static double frombits(uint64_t u) { double d; std::memcpy(&d, &u, 8); return d; }
static bool same(double a, double b) { return (std::isnan(a) && std::isnan(b)) || bits(a) == bits(b); }

static int runMath(long long nRand)
{
    std::mt19937_64 rng(12345);
    std::vector<double> x;
    // edge values
    const double ed[] = {0.0, -0.0, 1.0, -1.0, 2.0, 0.5, 1e-300, -1e-300, 5e-324, -5e-324,
                         std::numeric_limits<double>::min(), std::numeric_limits<double>::max(),
                         std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN(), 709.782712893384, 709.79, -708.39641853226,
                         -745.1332191019411, -745.14, -744.4400719213812, 1024.0, -1024.0, 512.0, -512.0,
                         0x1p-54, -0x1p-54, 0x1p-55, 1.0 - 0x1p-4, 1.0 + 0x1.09p-4, 1.0 + 0x1.08fffffffffffp-4};
    for (double d : ed) { x.push_back(d); x.push_back(std::nextafter(d, 1e308)); x.push_back(std::nextafter(d, -1e308)); }
    std::uniform_real_distribution<double> uExp(-800.0, 720.0), uExp2(-60.0, 5.0), uNear(1.0 - 0.0625, 1.0 + 0.0647),
        uLogw(0.0, 50.0), uSub(0.0, 1.0);
    std::uniform_int_distribution<uint64_t> uBits;
    for (long long i = 0; i < nRand; i++) {
        switch (i % 7) {
        case 0: x.push_back(uExp(rng)); break;
        case 1: x.push_back(uExp2(rng)); break;
        case 2: x.push_back(uNear(rng)); break;
        case 3: x.push_back(std::exp(uLogw(rng) - 25.0)); break;               // positive, wide range
        case 4: x.push_back(frombits(uBits(rng))); break;                        // any bit pattern
        case 5: x.push_back(frombits(uBits(rng) & 0x000fffffffffffffULL)); break; // subnormal
        default: x.push_back(std::ldexp(uSub(rng) + 1.0, (int)(uBits(rng) % 2046) - 1022)); break;
        }
    }
    const long long n = (long long)x.size();
    std::vector<double> de(n), dl(n);
    if (!saige::gpu2::erMathCheck(0, x.data(), n, de.data(), dl.data())) {
        std::fprintf(stderr, "erMathCheck failed: %s\n", saige::gpu2::erLastError());
        return 2;
    }
    long long be = 0, bl = 0;
    for (long long i = 0; i < n; i++) {
        const double he = std::exp(x[i]), hl = std::log(x[i]);
        if (!same(he, de[i])) { if (be < 10) std::printf("exp mismatch x=%a host=%a dev=%a\n", x[i], he, de[i]); be++; }
        if (!same(hl, dl[i])) { if (bl < 10) std::printf("log mismatch x=%a host=%a dev=%a\n", x[i], hl, dl[i]); bl++; }
    }
    std::printf("math: %lld inputs; exp differing bits %lld, log differing bits %lld\n", n, be, bl);
    return (be == 0 && bl == 0) ? 0 : 1;
}

struct Trait { int n; std::vector<double> mu; std::vector<double> y; int ncase; };

static int runEr(int nPairs, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    // traits: sizes and prevalences that cover small/large n, rare/common cases
    std::vector<Trait> T;
    const int    ns[]   = {50000, 50000, 50000, 2000, 300, 50000, 40};
    const double prev[] = {0.01, 0.10, 0.50, 0.05, 0.30, 0.90, 0.5};
    for (int t = 0; t < 7; t++) {
        Trait tr; tr.n = ns[t];
        tr.mu.resize(tr.n); tr.y.resize(tr.n);
        const double a0 = std::log(prev[t] / (1 - prev[t]));
        int nc = 0;
        for (int i = 0; i < tr.n; i++) {
            double eta = a0 + 1.2 * nd(rng);
            double m = 1.0 / (1.0 + std::exp(-eta));
            if (i % 997 == 0) m = 0.1 * (1 + (i / 997) % 9);     // exactly on a bin edge
            tr.mu[i] = m;
            tr.y[i] = (ud(rng) < m) ? 1.0 : 0.0;
            nc += (int)tr.y[i];
        }
        tr.ncase = nc;
        T.push_back(std::move(tr));
    }
    int maxN = 0;
    for (auto& tr : T) maxN = std::max(maxN, tr.n);
    std::vector<double> LT(maxN + 1);
    for (int i = 0; i <= maxN; i++) LT[i] = std::log((double)i);
    std::vector<saige::gpu2::ErTraitArgs> ta(T.size());
    for (size_t t = 0; t < T.size(); t++) { ta[t].mu = T[t].mu.data(); ta[t].n = T[t].n; ta[t].ncase = T[t].ncase; }
    saige::gpu2::ErCreateArgs ca;
    ca.device = 0; ca.nTraits = (int)T.size(); ca.traits = ta.data(); ca.logTable = LT.data(); ca.maxN = maxN;
    saige::gpu2::Er* E = saige::gpu2::erCreate(ca);
    if (!E) { std::fprintf(stderr, "erCreate failed: %s\n", saige::gpu2::erLastError()); return 2; }

    std::vector<ErPairIn> in;
    std::vector<uint32_t> cIdx; std::vector<double> cG; std::vector<unsigned char> cCase;
    std::vector<double> cpu;
    std::vector<int> kk;
    for (int q = 0; q < nPairs; q++) {
        const int t = (int)(rng() % T.size());
        const Trait& tr = T[t];
        int k;
        const int r = (int)(rng() % 100);
        if (r < 60) k = 1 + (int)(rng() % 4);          // the default MACCutoffforER range
        else if (r < 95) k = 5 + (int)(rng() % 12);
        else k = 17 + (int)(rng() % 4);                // up to 20
        k = std::min(k, tr.n - 1);
        // k distinct carriers, sorted
        std::vector<uint32_t> pos;
        while ((int)pos.size() < k) {
            uint32_t p = (uint32_t)(rng() % tr.n);
            if (q % 13 == 0) p = (uint32_t)((rng() % 3 == 0) ? 0 : (rng() % 2 ? tr.n - 1 : p));  // ends
            bool dup = false;
            for (auto v : pos) dup |= (v == p);
            if (!dup) pos.push_back(p);
        }
        std::sort(pos.begin(), pos.end());
        // genotype vector
        arma::vec G(tr.n, arma::fill::zeros);
        for (int i = 0; i < k; i++) {
            const int s = (int)(rng() % 10);
            double gv = (s < 7) ? 1.0 : (s < 9 ? 2.0 : 0.25 + 1.5 * ud(rng));   // hard calls + dosages
            G(pos[i]) = gv;
        }
        arma::uvec iIndex(k), iComp(tr.n - k);
        {
            int a = 0, b = 0;
            for (int i = 0; i < tr.n; i++) { if (G(i) != 0.0) iIndex(a++) = i; else iComp(b++) = i; }
        }
        arma::vec mu(const_cast<double*>(tr.mu.data()), tr.n, false, true);
        arma::vec yv(const_cast<double*>(tr.y.data()), tr.n, false, true);
        arma::vec res = yv - mu;
        if (q % 5 == 0) {   // carriers all cases (a strong signal, small p)
            for (int i = 0; i < k; i++) res(pos[i]) = 1.0 - mu(pos[i]);
        }
        arma::mat Z(tr.n, 1); Z.col(0) = G;
        arma::vec pi1 = mu;
        arma::mat resout;
        ER::SL_set_stream((uint64_t)q + 1);
        const double p = ER::SKATExactBin_Work(Z, res, pi1, (uint32_t)tr.ncase, iIndex, iComp, resout,
                                               2e+6, 1e+4, 1e-6, 1);
        cpu.push_back(p); kk.push_back(k);
        ErPairIn e; e.trait = t; e.k = k; e.off = (long long)cIdx.size();
        in.push_back(e);
        for (int i = 0; i < k; i++) {
            cIdx.push_back(pos[i]); cG.push_back(G(pos[i])); cCase.push_back(res(pos[i]) > 0 ? 1 : 0);
        }
    }
    std::vector<ErPairOut> out(in.size());
    if (!saige::gpu2::erRun(E, in.data(), (int)in.size(), cIdx.data(), cG.data(), cCase.data(),
                            (long long)cIdx.size(), out.data())) {
        std::fprintf(stderr, "erRun failed: %s\n", saige::gpu2::erLastError());
        return 2;
    }
    double tk = 0; long long np = 0;
    saige::gpu2::erTimings(E, &tk, &np);
    long long bad = 0, small = 0;
    double minp = 1;
    for (size_t q = 0; q < in.size(); q++) {
        if (cpu[q] < 1e-3) small++;
        minp = std::min(minp, cpu[q]);
        if (!same(cpu[q], out[q].pval)) {
            if (bad < 20) std::printf("pair %zu trait %d k %d: cpu %.17g dev %.17g\n", q, in[q].trait, kk[q], cpu[q], out[q].pval);
            bad++;
        }
    }
    std::printf("er: %zu pairs (k 1..20), %lld with p < 1e-3, min p %.3g; differing bits %lld; kernel %.4f s\n",
                in.size(), small, minp, bad, tk);
    saige::gpu2::erDestroy(E);
    return bad == 0 ? 0 : 1;
}

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: er_gpu_check math [n] | er [nPairs] [seed]\n"); return 2; }
    const std::string m = argv[1];
    if (m == "math") return runMath(argc > 2 ? std::atoll(argv[2]) : 20000000LL);
    if (m == "er") return runEr(argc > 2 ? std::atoi(argv[2]) : 3000, argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 1);
    return 2;
}
