// spa_gpos_test -- is spaGposGneg (spa_binary.cpp) bit-identical to the
// accu(g.elem(find(g > 0))) / accu(g.elem(find(g < 0))) pair it replaces?
//
// Armadillo materialises the compacted vector and sums it with two interleaved
// accumulators (op_accu_meat.hpp); spaGposGneg reproduces that without the
// index vector or the copy. This checks the claim on random vectors shaped
// like the covariate-adjusted genotype g~ the SPA actually sees (almost no
// exact zeros, positives and negatives of very different magnitude), on
// vectors with many exact zeros, and on tiny vectors, at many lengths.
//
//   g++ -O3 -march=native -std=c++17 tools/spa_gpos_test.cpp spa_binary.o ... -o spa_gpos_test
//   ./spa_gpos_test [trials]
#include <armadillo>
#include <cstdio>
#include <cstdlib>
#include <random>
#include "../spa_binary.hpp"

int main(int argc, char** argv)
{
    const long trials = argc > 1 ? std::atol(argv[1]) : 2000;
    std::mt19937_64 rng(20261002);
    std::uniform_int_distribution<int> lenPick(0, 6);
    const arma::uword lens[7] = {1, 2, 3, 17, 1000, 50000, 50001};
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> Z(0.0, 1.0);
    long bad = 0, done = 0;
    for (long t = 0; t < trials; ++t) {
        const arma::uword N = lens[lenPick(rng)];
        arma::vec g(N);
        const int shape = (int)(t % 4);
        for (arma::uword i = 0; i < N; ++i) {
            double v;
            if (shape == 0) {
                // g~ = g - X b: genotype in {0,1,2} minus a small projection
                const double u = U(rng);
                const double gi = (u < 0.81) ? 0.0 : (u < 0.99 ? 1.0 : 2.0);
                v = gi - (0.19 + 0.01 * Z(rng));
            } else if (shape == 1) {
                v = Z(rng) * std::pow(10.0, (int)(U(rng) * 12) - 6);   // wild magnitudes
            } else if (shape == 2) {
                const double u = U(rng);
                v = (u < 0.5) ? 0.0 : Z(rng);                           // half exact zeros
            } else {
                v = (U(rng) < 0.5) ? 1.0 : -1.0;                        // +-1 ties everywhere
            }
            g[i] = v;
        }
        const double gposRef = arma::accu(g.elem(arma::find(g > 0)));
        const double gnegRef = arma::accu(g.elem(arma::find(g < 0)));
        double gpos, gneg;
        spaGposGneg(g, gpos, gneg);
        // bit comparison, not ==: distinguishes -0.0 from 0.0 too
        if (std::memcmp(&gpos, &gposRef, 8) != 0 || std::memcmp(&gneg, &gnegRef, 8) != 0) {
            ++bad;
            if (bad <= 5)
                std::printf("MISMATCH N=%llu shape=%d gpos %.17g vs %.17g  gneg %.17g vs %.17g\n",
                            (unsigned long long)N, shape, gpos, gposRef, gneg, gnegRef);
        }
        ++done;
    }
    std::printf("spa_gpos_test: %ld trials, %ld mismatches -> %s\n", done, bad, bad ? "FAIL" : "PASS");
    return bad ? 1 : 0;
}
