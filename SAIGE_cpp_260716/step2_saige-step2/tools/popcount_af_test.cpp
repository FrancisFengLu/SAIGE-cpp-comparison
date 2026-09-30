// popcount_af_test -- pcSumPair against the literal gather of main.cpp's
// finalize, bit for bit, on random columns / masks / dosage tables.
//   g++ -O2 -std=c++17 -march=native -I.. tools/popcount_af_test.cpp -o tools/popcount_af_test
// Each trial: N samples (random, odd sizes included), a missing rate, a case
// rate, a table fd = post-flip {0,1,2} plus a missing dosage that is an
// integer, 2*AF (mean imputation), or a random double; the gather and
// pcSumPair must agree exactly on both sums and the four hom/het counts.
#include "../popcount_af.hpp"
#include <cstdio>
#include <random>
#include <vector>
#include <cstring>

using namespace SAIGE;

int main(int argc, char** argv)
{
    const long trials = argc > 1 ? atol(argv[1]) : 20000;
    const int  maxN   = argc > 2 ? atoi(argv[2]) : 3000;
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> U(0, 1);
    long nTrial = 0, nCounts = 0, nReplay = 0, nRefuse = 0, nBad = 0;
    long nCross = 0;   // replay trials whose running sum crossed a power of two after a fraction appeared
    for (long tr = 0; tr < trials; ++tr) {
        const int N = 1 + (int)(U(rng) * (tr % 5 == 0 ? 70 : maxN));
        const double missRate = (tr % 3 == 0) ? 0.0 : (tr % 3 == 1 ? 0.01 : 0.2);
        const double maf = std::exp(std::log(0.001) + U(rng) * (std::log(0.5) - std::log(0.001)));
        const double caseRate = 0.1 + 0.8 * U(rng);
        const bool flip = U(rng) < 0.5;
        const int nW = pcWords(N);
        std::vector<uint64_t> col(nW, 0), cm(nW, 0), om(nW, 0);
        std::vector<uint8_t> code(N);
        std::vector<uint32_t> ci, oi;
        for (int i = 0; i < N; ++i) {
            uint8_t c;
            if (U(rng) < missRate) c = PC_MISSING;
            else {
                const int g = (U(rng) < maf) + (U(rng) < maf);   // alt count
                c = g == 0 ? PC_HOM_REF : (g == 1 ? PC_HET : PC_HOM_ALT);
            }
            code[i] = c;
            col[i >> 5] |= (uint64_t)c << (2 * (i & 31));
            (U(rng) < caseRate ? ci : oi).push_back((uint32_t)i);
        }
        pcBuildMask(ci.data(), ci.size(), cm.data());
        pcBuildMask(oi.data(), oi.size(), om.data());
        // the table: alt-first codes, dmap = {2 (00), miss, 1 (10), 0 (11)}
        double fd[4];
        fd[PC_HOM_ALT] = flip ? 0.0 : 2.0;
        fd[PC_HET]     = 1.0;
        fd[PC_HOM_REF] = flip ? 2.0 : 0.0;
        switch (tr % 4) {
            case 0: fd[PC_MISSING] = std::round(2 * maf); break;       // best_guess
            case 1: fd[PC_MISSING] = 2 * maf; break;                    // mean
            case 2: fd[PC_MISSING] = 0.0; break;                        // minor
            default: fd[PC_MISSING] = U(rng) * 2; break;                // anything
        }
        // the gather as main.cpp writes it
        std::vector<double> g(N);
        for (int i = 0; i < N; ++i) g[i] = fd[code[i]];
        double sc = 0, so = 0; uint32_t chom = 0, chet = 0, ohom = 0, ohet = 0;
        for (uint32_t k : ci) { const double d = g[k]; sc += d; if (d >= 1.5 && d <= 2.0) chom++; else if (d >= 0.5 && d < 1.5) chet++; }
        for (uint32_t k : oi) { const double d = g[k]; so += d; if (d >= 1.5 && d <= 2.0) ohom++; else if (d >= 0.5 && d < 1.5) ohet++; }
        // crossing diagnostic: did the sequential sum ever cross a power of two while non-integer?
        {
            double s = 0; bool frac = false;
            for (uint32_t k : ci) { const double s2 = s + g[k]; if (frac && s > 0 && pcNextPow2Above(s) <= s2) { nCross++; break; } s = s2; if (s != std::floor(s)) frac = true; }
        }
        double pc = -1, po = -1; uint32_t pch = 0, pce = 0, poh = 0, poe = 0;
        const int rc = pcSumPair(col.data(), cm.data(), om.data(), nW, ci.size(), oi.size(), fd, true,
                                 pc, po, true, pch, pce, poh, poe);
        nTrial++;
        if (rc == 0) { nRefuse++; continue; }
        if (rc == 1) nCounts++; else nReplay++;
        if (std::memcmp(&pc, &sc, 8) != 0 || std::memcmp(&po, &so, 8) != 0 ||
            pch != chom || pce != chet || poh != ohom || poe != ohet) {
            nBad++;
            if (nBad <= 10)
                std::printf("MISMATCH trial %ld N=%d miss=%.2f fdMiss=%.17g rc=%d  case %.17g vs %.17g  ctrl %.17g vs %.17g  hom/het %u/%u vs %u/%u\n",
                            tr, N, missRate, fd[PC_MISSING], rc, pc, sc, po, so, pch, pce, chom, chet);
        }
    }
    std::printf("trials %ld  counts-path %ld  replay-path %ld  refused %ld  crossings-seen %ld  mismatches %ld\n",
                nTrial, nCounts, nReplay, nRefuse, nCross, nBad);
    return nBad == 0 ? 0 : 1;
}
