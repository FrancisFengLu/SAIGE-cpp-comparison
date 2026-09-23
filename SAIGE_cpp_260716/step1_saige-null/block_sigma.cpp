#include "block_sigma.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace blocksigma {

namespace {

bool  g_enabled = false;
bool  g_verify  = false;

// Verification accumulators (fit.block_sparse_sigma_verify).
long long g_vcalls = 0;
double    g_vmaxabs = 0.0;
double    g_vmaxrel = 0.0;
long long g_vexact = 0;        // calls where every fp32 element matched bitwise

struct UnionFind {
    std::vector<int> p;
    explicit UnionFind(int n) : p((size_t)n) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) {
        while (p[(size_t)x] != x) { p[(size_t)x] = p[(size_t)p[(size_t)x]]; x = p[(size_t)x]; }
        return x;
    }
    void join(int a, int b) { a = find(a); b = find(b); if (a != b) p[(size_t)a] = b; }
};

}  // namespace

BlockSigma& instance() { static BlockSigma s; return s; }

void enable(bool on)  { g_enabled = on; }
bool enabled()        { return g_enabled; }
void enableVerify(bool on) { g_verify = on; }
bool verifyEnabled()  { return g_verify; }

void resetVerify() { g_vcalls = 0; g_vmaxabs = 0.0; g_vmaxrel = 0.0; g_vexact = 0; }

void recordVerify(const arma::fvec& ref, const arma::fvec& got) {
    if (ref.n_elem != got.n_elem) return;
    g_vcalls++;
    bool bitwise = true;
    double maxabs = 0.0, maxrel = 0.0;
    for (arma::uword i = 0; i < ref.n_elem; i++) {
        const float a = ref(i), b = got(i);
        if (a != b) bitwise = false;
        const double d = std::fabs((double)a - (double)b);
        if (d > maxabs) maxabs = d;
        const double den = std::fabs((double)a);
        if (den > 0.0) { const double r = d / den; if (r > maxrel) maxrel = r; }
    }
    if (bitwise) g_vexact++;
    if (maxabs > g_vmaxabs) g_vmaxabs = maxabs;
    if (maxrel > g_vmaxrel) g_vmaxrel = maxrel;
}

void reportVerify(const char* what) {
    if (!g_verify || g_vcalls == 0) return;
    printf("[blocksigma] verify %s: %lld calls, %lld bitwise-identical (%.1f%%), "
           "max |abs diff| %.3e, max |rel diff| %.3e\n",
           what, g_vcalls, g_vexact, 100.0 * (double)g_vexact / (double)g_vcalls,
           g_vmaxabs, g_vmaxrel);
    fflush(stdout);
}

bool BlockSigma::build(const arma::umat& loc, const arma::vec& val, int n) {
    part_ = Partition();
    inv_.clear(); scratch_.clear(); refreshed_ = false; nFloored_ = 0;
    if (n <= 0 || loc.n_cols == 0 || loc.n_rows < 2 || val.n_elem != loc.n_cols)
        return false;

    // --- components over the off-diagonal entries -----------------------
    UnionFind uf(n);
    for (arma::uword k = 0; k < loc.n_cols; k++) {
        const int a = (int)loc(0, k), b = (int)loc(1, k);
        if (a < 0 || b < 0 || a >= n || b >= n) return false;
        if (a != b) uf.join(a, b);
    }
    std::vector<int> root((size_t)n), id((size_t)n, -1);
    int nb = 0;
    for (int i = 0; i < n; i++) {
        root[(size_t)i] = uf.find(i);
        if (id[(size_t)root[(size_t)i]] < 0) id[(size_t)root[(size_t)i]] = nb++;
    }
    part_.n = n; part_.nblocks = nb;
    part_.blockOf.assign((size_t)n, -1);
    part_.localOf.assign((size_t)n, -1);
    part_.size.assign((size_t)nb, 0);
    for (int i = 0; i < n; i++) {
        const int b = id[(size_t)root[(size_t)i]];
        part_.blockOf[(size_t)i] = b;
        part_.localOf[(size_t)i] = part_.size[(size_t)b]++;
    }
    part_.start.assign((size_t)nb, 0);
    part_.invStart.assign((size_t)nb, 0);
    long long off = 0, ioff = 0;
    for (int b = 0; b < nb; b++) {
        part_.start[(size_t)b]    = (int)off;   off  += part_.size[(size_t)b];
        part_.invStart[(size_t)b] = (int)ioff;  ioff += (long long)part_.size[(size_t)b]
                                                        * part_.size[(size_t)b];
        part_.maxBlock = std::max(part_.maxBlock, part_.size[(size_t)b]);
    }
    part_.member.assign((size_t)n, -1);
    for (int i = 0; i < n; i++)
        part_.member[(size_t)(part_.start[(size_t)part_.blockOf[(size_t)i]]
                              + part_.localOf[(size_t)i])] = i;

    // --- Psi in per-block local coordinates ------------------------------
    std::vector<int> cnt((size_t)nb, 0);
    for (arma::uword k = 0; k < loc.n_cols; k++) cnt[(size_t)part_.blockOf[(size_t)loc(0, k)]]++;
    part_.psiStart.assign((size_t)nb + 1, 0);
    for (int b = 0; b < nb; b++) part_.psiStart[(size_t)b + 1] = part_.psiStart[(size_t)b] + cnt[(size_t)b];
    const int nnz = part_.psiStart[(size_t)nb];
    part_.psiR.assign((size_t)nnz, 0);
    part_.psiC.assign((size_t)nnz, 0);
    part_.psiV.assign((size_t)nnz, 0.0);
    std::vector<int> fill(part_.psiStart.begin(), part_.psiStart.end() - 1);
    for (arma::uword k = 0; k < loc.n_cols; k++) {
        const int r = (int)loc(0, k), c = (int)loc(1, k);
        const int b = part_.blockOf[(size_t)r];
        if (part_.blockOf[(size_t)c] != b) return false;   // pattern is not block diagonal
        const int at = fill[(size_t)b]++;
        part_.psiR[(size_t)at] = part_.localOf[(size_t)r];
        part_.psiC[(size_t)at] = part_.localOf[(size_t)c];
        part_.psiV[(size_t)at] = val(k);
    }

    inv_.assign((size_t)ioff, 0.0);
    part_.built = true;
    return true;
}

bool BlockSigma::upToDate(const arma::fvec& w, const arma::fvec& tau) const {
    if (!refreshed_ || lastW_.n_elem != w.n_elem || lastTau_.n_elem != tau.n_elem)
        return false;
    for (arma::uword i = 0; i < tau.n_elem; i++) if (lastTau_(i) != tau(i)) return false;
    for (arma::uword i = 0; i < w.n_elem;   i++) if (lastW_(i)   != w(i))   return false;
    return true;
}

void BlockSigma::refresh(const arma::fvec& w, const arma::fvec& tau) {
    if (!part_.built) return;
    if (upToDate(w, tau)) { nReuse_++; return; }
    nRefresh_++;
    const int nb = part_.nblocks;
    const float tau0 = tau(0), tau1f = tau(1);
    const double tau1 = (double)tau1f;

    // gen_sp_Sigma computes dtVec = (1/wVec) * tau0 in fp32 and adds it to an
    // fp64 accumulator, so the fp32 rounding of that product is part of the
    // matrix being solved. Reproduce it rather than computing in fp64.
    long long floored = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:floored)
#endif
    for (int b = 0; b < nb; b++) {
        const int s   = part_.size[(size_t)b];
        const int mo  = part_.start[(size_t)b];
        double* A = &inv_[(size_t)part_.invStart[(size_t)b]];
        for (int i = 0; i < s * s; i++) A[i] = 0.0;

        for (int k = part_.psiStart[(size_t)b]; k < part_.psiStart[(size_t)b + 1]; k++) {
            const int r = part_.psiR[(size_t)k], c = part_.psiC[(size_t)k];
            double v = part_.psiV[(size_t)k] * tau1;
            if (r == c) {
                const int sample = part_.member[(size_t)(mo + r)];
                const float dt = (1.0f / w((arma::uword)sample)) * tau0;
                v += (double)dt;
                if (v < 1e-4) { v = 1e-4; floored++; }
            }
            A[(size_t)c * s + r] = v;     // column-major
        }

        // Invert in place. Sizes here are tiny (1-10 on the benchmark GRM);
        // 1 and 2 are special-cased because they are 97% of the blocks.
        if (s == 1) {
            A[0] = 1.0 / A[0];
        } else if (s == 2) {
            const double a = A[0], c = A[1], bb = A[2], d = A[3];
            const double det = a * d - bb * c;
            A[0] =  d / det; A[1] = -c / det; A[2] = -bb / det; A[3] =  a / det;
        } else {
            arma::mat M((double*)A, s, s, false, true);   // alias, no copy
            arma::mat Mi;
            if (!arma::inv_sympd(Mi, M)) {
                if (!arma::inv(Mi, M)) Mi = arma::pinv(M);
            }
            std::copy(Mi.memptr(), Mi.memptr() + (size_t)s * s, A);
        }
    }
    nFloored_ = floored;
    lastW_ = w; lastTau_ = tau;
    refreshed_ = true;
}

arma::fvec BlockSigma::solve(const arma::fvec& b) const {
    arma::fvec out(b.n_elem, arma::fill::zeros);
    if (!ready() || (int)b.n_elem != part_.n) return out;
    const int nb = part_.nblocks;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int blk = 0; blk < nb; blk++) {
        const int s  = part_.size[(size_t)blk];
        const int mo = part_.start[(size_t)blk];
        const double* A = &inv_[(size_t)part_.invStart[(size_t)blk]];
        if (s == 1) {
            const int i = part_.member[(size_t)mo];
            out((arma::uword)i) = (float)(A[0] * (double)b((arma::uword)i));
            continue;
        }
        // Accumulate in fp64, round once on store, matching gen_spsolve_v4,
        // which solves in fp64 and converts the result to fp32 at the end.
        for (int r = 0; r < s; r++) {
            double acc = 0.0;
            for (int c = 0; c < s; c++)
                acc += A[(size_t)c * s + r] * (double)b((arma::uword)part_.member[(size_t)(mo + c)]);
            out((arma::uword)part_.member[(size_t)(mo + r)]) = (float)acc;
        }
    }
    return out;
}

}  // namespace blocksigma
