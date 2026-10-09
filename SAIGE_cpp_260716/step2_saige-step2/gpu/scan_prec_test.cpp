// scan_prec_test.cpp — direct accuracy test of the step-2 reducer
// (gpu_step2.hpp) in each scan precision against a long-double reference.
// Simulated data only. Build: make USE_CUDA=1 gpu/scan_prec_test (not part of
// `all`). Run: gpu/scan_prec_test [N K1 K2 nSlots maxSlots missRate slices]
//
// Covers what the gate data sets do not: N not a multiple of 4 or 16, K2 = 0,
// markers whose missing call is mean-imputed (int8 indicator columns, and
// passes split by the indicator cap), cleaned markers (all-zero table),
// flipped tables, and B columns of very different scales. For every mode it
// prints, per result block, max |C - ref| / sqrt(sum_i (g_i b_ik)^2) -- the
// error in units of the sum's own random-walk scale -- and the same over
// sum_i |g_i| max_i |b_ik|, the scale of int8's per-column bound. Exits
// non-zero if fp64 exceeds 1e-12 or fp32 1e-4 on the first, or int8 its
// bound 2^(1 - 7 slices) (x4, + 1e-14) on the second.
#include "gpu_step2.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

using namespace saige::gpu2;

int main(int argc, char** argv)
{
    const int N        = argc > 1 ? std::atoi(argv[1]) : 10007;
    const int K1       = argc > 2 ? std::atoi(argv[2]) : 37;
    const int K2       = argc > 3 ? std::atoi(argv[3]) : 5;
    const int nSlots   = argc > 4 ? std::atoi(argv[4]) : 700;
    const int maxSlots = argc > 5 ? std::atoi(argv[5]) : 1024;
    const double miss  = argc > 6 ? std::atof(argv[6]) : 0.03;
    const int slices   = argc > 7 ? std::atoi(argv[7]) : kInt8SlicesDefault;
    const int fracAll  = argc > 8 ? std::atoi(argv[8]) : 0;   // every table fractional (int8 pass splitting)
    std::string why;
    if (!available(0, &why)) { std::printf("no device: %s\n", why.c_str()); return 0; }
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> Z(0.0, 1.0);

    // B1: column kinds by k % 5 -- residual-like (mean 0), positive around a
    // mean (mu2-like), tiny scale, large scale, skewed sparse.
    std::vector<double> B1((std::size_t)N * K1), B2((std::size_t)N * (K2 > 0 ? K2 : 1));
    for (int k = 0; k < K1; ++k)
        for (int i = 0; i < N; ++i) {
            double v;
            switch (k % 5) {
                case 0:  v = (U(rng) < 0.1 ? 0.9 : -0.1) + 0.01 * Z(rng); break;
                case 1:  v = 0.01 + 0.002 * Z(rng); break;
                case 2:  v = 2e-5 + 2e-5 * Z(rng); break;
                case 3:  v = 1e3 * Z(rng); break;
                default: v = U(rng) < 0.01 ? Z(rng) : 0.0; break;
            }
            B1[(std::size_t)k * N + i] = v;
        }
    for (int k = 0; k < K2; ++k)
        for (int i = 0; i < N; ++i) B2[(std::size_t)k * N + i] = 0.02 * U(rng) + 1e-4;

    // markers: AF spread over (0, 1), a few cleaned, half flipped
    const std::size_t nb = (std::size_t)(N + 3) / 4;
    std::vector<unsigned char> codes((std::size_t)nSlots * N);
    std::vector<double> lutv((std::size_t)nSlots * 4);
    for (int s = 0; s < nSlots; ++s) {
        const double af = (s % 7 == 0) ? 1.0 / N : std::pow(U(rng), 2.0);
        const double ms = (s % 3 == 0) ? 0.0 : miss;
        double sum = 0; int n = 0;
        for (int i = 0; i < N; ++i) {
            unsigned c;
            if (U(rng) < ms) c = 1;
            else { const int g = (U(rng) < af) + (U(rng) < af); c = g == 0 ? 3u : (g == 1 ? 2u : 0u); sum += g; ++n; }
            codes[(std::size_t)s * N + i] = (unsigned char)c;
        }
        // PLINK codes: 0 HOM_ALT(2), 1 MISSING, 2 HET(1), 3 HOM_REF(0)
        const double imp = n > 0 ? sum / n : 0.0;
        double* L = &lutv[(std::size_t)s * 4];
        L[0] = 2; L[1] = imp; L[2] = 1; L[3] = 0;
        if (s % 2 == 1) for (int c = 0; c < 4; ++c) L[c] = 2.0 - L[c];
        if (s % 11 == 5) for (int c = 0; c < 4; ++c) L[c] = 0.0;
        if (fracAll || s % 13 == 4) { L[0] = 1.75; L[2] = 0.625; L[3] = 0.125 + 0.001 * (s % 5); }
    }

    // long-double reference
    std::vector<long double> R1((std::size_t)nSlots * K1), R2((std::size_t)nSlots * (K2 > 0 ? K2 : 1));
    std::vector<long double> S1((std::size_t)nSlots * K1), S2((std::size_t)nSlots * (K2 > 0 ? K2 : 1));
    // int8's bound is per column and absolute: |err| <= sum_i |g_i| 2^(E_k - 7 slices),
    // max_i |b_ik| >= 2^(E_k - 1); A = sum_i |g_i| max_i |b_ik|
    std::vector<long double> A1((std::size_t)nSlots * K1), A2((std::size_t)nSlots * (K2 > 0 ? K2 : 1));
    std::vector<double> mx1(K1, 0.0), mx2(K2 > 0 ? K2 : 1, 0.0);
    for (int k = 0; k < K1; ++k) for (int i = 0; i < N; ++i) mx1[k] = std::max(mx1[k], std::fabs(B1[(std::size_t)k * N + i]));
    for (int k = 0; k < K2; ++k) for (int i = 0; i < N; ++i) mx2[k] = std::max(mx2[k], std::fabs(B2[(std::size_t)k * N + i]));
    for (int s = 0; s < nSlots; ++s) {
        const double* L = &lutv[(std::size_t)s * 4];
        long double ga = 0, g2a = 0;
        for (int i = 0; i < N; ++i) {
            const double g = L[codes[(std::size_t)s * N + i]];
            ga += std::fabs(g); g2a += g * g;
        }
        for (int k = 0; k < K1; ++k) A1[(std::size_t)s * K1 + k] = ga * mx1[k];
        for (int k = 0; k < K2; ++k) A2[(std::size_t)s * K2 + k] = g2a * mx2[k];
        for (int k = 0; k < K1; ++k) {
            long double a = 0, q = 0;
            for (int i = 0; i < N; ++i) {
                const long double t = (long double)L[codes[(std::size_t)s * N + i]] * B1[(std::size_t)k * N + i];
                a += t; q += t * t;
            }
            R1[(std::size_t)s * K1 + k] = a; S1[(std::size_t)s * K1 + k] = std::sqrt(q);
        }
        for (int k = 0; k < K2; ++k) {
            long double a = 0, q = 0;
            for (int i = 0; i < N; ++i) {
                const long double g = L[codes[(std::size_t)s * N + i]];
                const long double t = (long double)(double)(g * g) * B2[(std::size_t)k * N + i];
                a += t; q += t * t;
            }
            R2[(std::size_t)s * K2 + k] = a; S2[(std::size_t)s * K2 + k] = std::sqrt(q);
        }
    }

    int bad = 0;
    for (Prec p : {Prec::FP64, Prec::FP32, Prec::INT8}) {
        CreateArgs a;
        a.N = N; a.K1 = K1; a.B1 = B1.data(); a.K2 = K2; a.B2 = K2 > 0 ? B2.data() : nullptr;
        a.maxSlots = maxSlots; a.precision = p; a.int8Slices = slices;
        Reducer* r = create(a);
        if (!r) { std::printf("%s: create failed\n", precName(p)); ++bad; continue; }
        double e1 = 0, e2 = 0, f1 = 0, f2 = 0;
        for (int s0 = 0; s0 < nSlots; s0 += maxSlots) {
            const int ns = std::min(maxSlots, nSlots - s0);
            unsigned char* pk = packed(r);
            double* lu = lut(r);
            const std::size_t bpv = bytesPerSlot(r);
            for (int j = 0; j < ns; ++j) {
                unsigned char* row = pk + (std::size_t)j * bpv;
                std::memset(row, 0, nb);
                for (int i = 0; i < N; ++i)
                    row[i >> 2] |= (unsigned char)(codes[(std::size_t)(s0 + j) * N + i] << ((i & 3) * 2));
                std::memcpy(lu + (std::size_t)j * 4, &lutv[(std::size_t)(s0 + j) * 4], 4 * sizeof(double));
            }
            if (!reduce(r, ns)) { std::printf("%s: reduce failed\n", precName(p)); ++bad; break; }
            const double* C1 = outCd(r);
            const double* C2 = outC2d(r);
            const std::size_t ld = ldC(r);
            for (int j = 0; j < ns; ++j) {
                for (int k = 0; k < K1; ++k) {
                    const long double sc = S1[(std::size_t)(s0 + j) * K1 + k];
                    const long double d = std::fabs((long double)C1[(std::size_t)k * ld + j] - R1[(std::size_t)(s0 + j) * K1 + k]);
                    if (sc > 0) e1 = std::max(e1, (double)(d / sc)); else if (d > 0) e1 = INFINITY;
                    const long double as = A1[(std::size_t)(s0 + j) * K1 + k];
                    if (as > 0) f1 = std::max(f1, (double)(d / as)); else if (d > 0) f1 = INFINITY;
                }
                for (int k = 0; k < K2; ++k) {
                    const long double sc = S2[(std::size_t)(s0 + j) * K2 + k];
                    const long double d = std::fabs((long double)C2[(std::size_t)k * ld + j] - R2[(std::size_t)(s0 + j) * K2 + k]);
                    if (sc > 0) e2 = std::max(e2, (double)(d / sc)); else if (d > 0) e2 = INFINITY;
                    const long double as = A2[(std::size_t)(s0 + j) * K2 + k];
                    if (as > 0) f2 = std::max(f2, (double)(d / as)); else if (d > 0) f2 = INFINITY;
                }
            }
        }
        // fp64: 1e-12 of the scaled error; fp32: 1e-4; int8: the split's bound
        // 2^(1 - 7 slices) of A (plus fp64 rounding of the recombination).
        bool ok;
        if (p == Prec::INT8) {
            const double tol = 4.0 * std::ldexp(1.0, 1 - 7 * slices) + 1e-14;
            ok = f1 <= tol && f2 <= tol;
        } else {
            const double tol = (p == Prec::FP32) ? 1e-4 : 1e-12;
            ok = e1 <= tol && e2 <= tol;
        }
        std::printf("%-5s N=%d K1=%d K2=%d slots=%d/%d miss=%g: C1 err %.3g (abs-scaled %.3g), "
                    "C2 err %.3g (abs-scaled %.3g) %s\n",
                    precName(p), N, K1, K2, nSlots, maxSlots, miss, e1, f1, e2, f2, ok ? "ok" : "FAIL");
        if (!ok) ++bad;
        destroy(r);
    }
    return bad ? 1 : 0;
}
