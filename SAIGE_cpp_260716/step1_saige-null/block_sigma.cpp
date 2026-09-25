#include "block_sigma.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>

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

// ---------------------------------------------------------------------------
// Small dense symmetric inverse, in place, column-major, zero allocation.
//
// arma::inv_sympd on a 3x3 costs a heap allocation for the result matrix, a
// LAPACK dispatch and a copy back; measured on this machine it is 5.7x the
// arithmetic at s=3 and still 2.3x at s=10. The crossover where LAPACK's
// blocking wins is s ~= 36 (bench_small2.cpp: ratio 1.07 at s=32, 0.97 at
// s=40), so refresh() sends [3, kCholMax] here and everything larger to
// arma::inv_sympd.
//
// Cholesky A = L L^T, then L^-1, then A^-1 = L^-T L^-1 -- the same three steps
// as dpotrf + dtrtri + dlauum, written out so nothing allocates. Returns false
// if A is not positive definite, in which case the caller falls back exactly as
// before (inv_sympd -> inv -> pinv).
static bool chol_inv(double* A, int n) {
    for (int j = 0; j < n; j++) {
        double d = A[(size_t)j * n + j];
        for (int k = 0; k < j; k++) { const double v = A[(size_t)k * n + j]; d -= v * v; }
        if (!(d > 0.0)) return false;
        d = std::sqrt(d);
        A[(size_t)j * n + j] = d;
        const double id = 1.0 / d;
        for (int i = j + 1; i < n; i++) {
            double t = A[(size_t)j * n + i];
            for (int k = 0; k < j; k++) t -= A[(size_t)k * n + i] * A[(size_t)k * n + j];
            A[(size_t)j * n + i] = t * id;
        }
    }
    // L := L^-1, in place. Column j only ever reads columns <= j.
    for (int j = 0; j < n; j++) {
        A[(size_t)j * n + j] = 1.0 / A[(size_t)j * n + j];
        for (int i = j + 1; i < n; i++) {
            double t = 0.0;
            for (int k = j; k < i; k++) t += A[(size_t)k * n + i] * A[(size_t)j * n + k];
            A[(size_t)j * n + i] = -t / A[(size_t)i * n + i];
        }
    }
    // A := L^-T L^-1. Entry (i,j) reads column j only at rows >= i, and column
    // i only at rows >= i, so writing (i,j) with i ascending inside a column is
    // safe in place.
    for (int j = 0; j < n; j++)
        for (int i = j; i < n; i++) {
            double t = 0.0;
            for (int k = i; k < n; k++) t += A[(size_t)i * n + k] * A[(size_t)j * n + k];
            A[(size_t)i * n + j] = t;
            A[(size_t)j * n + i] = t;
        }
    return true;
}

// Closed-form symmetric 3x3 and 4x4 by cofactors. Measured on the benchmark
// GRM's histogram they take the s>=3 refresh from 3.20x to 3.61x over
// inv_sympd. det == 0 falls through to chol_inv, which falls through to
// armadillo, so the failure ladder is unchanged.
static inline bool sym3_inv(double* A) {
    const double a = A[0], b = A[1], c = A[2], d = A[4], e = A[5], f = A[8];
    const double C0 = d * f - e * e, C1 = c * e - b * f, C2 = b * e - c * d;
    const double det = a * C0 + b * C1 + c * C2;
    if (det == 0.0 || !std::isfinite(det)) return false;
    const double id = 1.0 / det;
    const double m00 = C0 * id, m01 = C1 * id, m02 = C2 * id;
    const double m11 = (a * f - c * c) * id, m12 = (b * c - a * e) * id,
                 m22 = (a * d - b * b) * id;
    A[0] = m00; A[1] = m01; A[2] = m02;
    A[3] = m01; A[4] = m11; A[5] = m12;
    A[6] = m02; A[7] = m12; A[8] = m22;
    return true;
}

static inline bool sym4_inv(double* A) {
    const int L = 4;
    auto m = [&](int i, int j) { return A[(size_t)j * L + i]; };
    auto det3 = [&](int r0, int r1, int r2, int c0, int c1, int c2) {
        return m(r0, c0) * (m(r1, c1) * m(r2, c2) - m(r1, c2) * m(r2, c1))
             - m(r0, c1) * (m(r1, c0) * m(r2, c2) - m(r1, c2) * m(r2, c0))
             + m(r0, c2) * (m(r1, c0) * m(r2, c1) - m(r1, c1) * m(r2, c0));
    };
    static const int oth[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
    double C[4][4];
    for (int i = 0; i < 4; i++)
        for (int j = i; j < 4; j++) {
            const double mm = det3(oth[i][0], oth[i][1], oth[i][2],
                                   oth[j][0], oth[j][1], oth[j][2]);
            C[i][j] = C[j][i] = ((i + j) & 1) ? -mm : mm;
        }
    double det = 0.0;
    for (int j = 0; j < 4; j++) det += m(0, j) * C[0][j];
    if (det == 0.0 || !std::isfinite(det)) return false;
    const double id = 1.0 / det;
    for (int j = 0; j < 4; j++)
        for (int i = 0; i < 4; i++) A[(size_t)j * L + i] = C[i][j] * id;
    return true;
}

// Largest block still sent to chol_inv. Above this arma::inv_sympd (LAPACK's
// blocked dpotrf/dpotri) wins; see the comment on chol_inv.
static const int kCholMax = 32;

// One block's inverse, in place, with the dispatch refresh() uses. Split out so
// the cost calibration in build() times exactly what refresh() will run.
static void invert_block(double* A, int s) {
    if (s == 1) { A[0] = 1.0 / A[0]; return; }
    if (s == 2) {
        const double a = A[0], c = A[1], bb = A[2], d = A[3];
        const double det = a * d - bb * c;
        A[0] = d / det; A[1] = -c / det; A[2] = -bb / det; A[3] = a / det;
        return;
    }
    if (s == 3 && sym3_inv(A)) return;
    if (s == 4 && sym4_inv(A)) return;
    if (s <= kCholMax && chol_inv(A, s)) return;
    arma::mat M((double*)A, s, s, false, true);   // alias, no copy
    arma::mat Mi;
    if (!arma::inv_sympd(Mi, M)) {
        if (!arma::inv(Mi, M)) Mi = arma::pinv(M);
    }
    std::copy(Mi.memptr(), Mi.memptr() + (size_t)s * s, A);
}

static double now_s() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

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
    if (buildFailed_) return false;          // decided once, not per call
    part_ = Partition();
    inv_.clear(); scratch_.clear(); refreshed_ = false; nFloored_ = 0;
    if (n <= 0 || loc.n_cols == 0 || loc.n_rows < 2 || val.n_elem != loc.n_cols) {
        buildFailed_ = true;
        return false;
    }

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
        if (part_.blockOf[(size_t)c] != b) { buildFailed_ = true; return false; }  // not block diagonal
        const int at = fill[(size_t)b]++;
        part_.psiR[(size_t)at] = part_.localOf[(size_t)r];
        part_.psiC[(size_t)at] = part_.localOf[(size_t)c];
        part_.psiV[(size_t)at] = val(k);
    }

    // --- size-1 fast lane -----------------------------------------------
    // A size-1 block whose Psi contribution is exactly one diagonal entry is
    // one reciprocal; nothing about size/start/psiStart needs re-reading.
    // Anything else (no entry at all, or a duplicate) stays on the general
    // path, so the arithmetic is bit-for-bit what it was.
    part_.one_sample.clear(); part_.one_psi.clear(); part_.one_slot.clear();
    part_.genBlock.clear();
    for (int b = 0; b < nb; b++) {
        const int k0 = part_.psiStart[(size_t)b], k1 = part_.psiStart[(size_t)b + 1];
        if (part_.size[(size_t)b] == 1 && k1 - k0 == 1 &&
            part_.psiR[(size_t)k0] == 0 && part_.psiC[(size_t)k0] == 0) {
            part_.one_sample.push_back(part_.member[(size_t)part_.start[(size_t)b]]);
            part_.one_psi.push_back(part_.psiV[(size_t)k0]);
            part_.one_slot.push_back(part_.invStart[(size_t)b]);
        } else {
            part_.genBlock.push_back(b);
        }
    }

    // --- cost gate --------------------------------------------------------
    flops_ = 0.0; bytes_ = 0.0;
    for (int b = 0; b < nb; b++) {
        const double sz = (double)part_.size[(size_t)b];
        flops_ += sz * sz * sz;
        bytes_ += sz * sz * 8.0;
    }
    // Estimate the seconds one refresh costs by timing a real inverse at the
    // largest block size, with the dispatch refresh() will use, and scaling by
    // sum(b^3). Intentionally a SERIAL estimate -- refresh() is OpenMP over
    // blocks, so the measured cost is at most this and usually several times
    // less. The measured value is printed after the first refresh so a bad
    // estimate shows up instead of hiding.
    estSecs_ = estimateRefreshSeconds();
    const bool overTime  = (refreshBudget_ > 0.0 && estSecs_ > refreshBudget_);
    const bool overFlops = (flopBudget_    > 0.0 && flops_   > flopBudget_);
    const bool overBytes = (byteBudget_    > 0.0 && bytes_   > byteBudget_);
    if (overTime || overFlops || overBytes) {
        printf("[blocksigma] refusing: %d blocks, max %d, sum(b^3) = %.3g, "
               "block inverses %.1f MB, estimated %.3f s per refresh -- over the "
               "budget (%.3f s%s, %.0f MB). Falling back to spsolve.\n",
               nb, part_.maxBlock, flops_, bytes_ / 1048576.0, estSecs_,
               refreshBudget_, overFlops ? ", raw flop gate hit" : "",
               byteBudget_ / 1048576.0);
        fflush(stdout);
        part_ = Partition();
        buildFailed_ = true;
        return false;
    }

    inv_.assign((size_t)ioff, 0.0);
    part_.built = true;
    printf("[blocksigma] partition: %d blocks, max %d (%zu on the size-1 fast "
           "lane), sum(b^3) = %.3g, block inverses %.2f MB, estimated %.4f s "
           "per refresh (budget %.3f s)\n",
           nb, part_.maxBlock, part_.one_slot.size(), flops_, bytes_ / 1048576.0,
           estSecs_, refreshBudget_);
    fflush(stdout);
    return true;
}

// Seconds one refresh would cost if it ran on a single thread. Times
// invert_block at the largest block size present -- the same code path
// refresh() takes -- and extrapolates linearly in sum(b^3). Costs one dense
// inverse of the largest block, i.e. less than one refresh of it.
double BlockSigma::estimateRefreshSeconds() const {
    const int smax = part_.maxBlock;
    if (smax <= 2 || flops_ <= 0.0) return 0.0;
    const int sc = std::min(smax, 2048);
    std::vector<double> ref((size_t)sc * sc), work((size_t)sc * sc);
    std::mt19937_64 rng(20260925u);
    std::uniform_real_distribution<double> ud(-1.0, 1.0);
    for (int j = 0; j < sc; j++)
        for (int i = j; i < sc; i++) {
            const double v = (i == j) ? (1.0 + 0.5 * std::fabs(ud(rng)))
                                      : (0.2 * ud(rng) / sc);
            ref[(size_t)j * sc + i] = v;
            ref[(size_t)i * sc + j] = v;
        }
    double best = 1e30;
    int reps = 0;
    const double t_start = now_s();
    while (reps < 200 && (now_s() - t_start) < 0.05) {
        work = ref;
        const double t0 = now_s();
        invert_block(work.data(), sc);
        best = std::min(best, now_s() - t0);
        reps++;
    }
    if (!(best < 1e30) || best <= 0.0) return 0.0;
    const double per_cube = best / ((double)sc * sc * sc);
    return per_cube * flops_;
}

void BlockSigma::setBudget(double flopBudget, double byteBudget) {
    // <= 0 disables that gate. The raw flop gate is off by default now; the
    // seconds budget is the operational one.
    if (flopBudget >= 0) flopBudget_ = flopBudget;
    if (byteBudget >  0) byteBudget_ = byteBudget;
}

void BlockSigma::setRefreshBudget(double seconds) {
    if (seconds >= 0) refreshBudget_ = seconds;
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
    const double t_begin = now_s();
    const float tau0 = tau(0), tau1f = tau(1);
    const double tau1 = (double)tau1f;

    // gen_sp_Sigma computes dtVec = (1/wVec) * tau0 in fp32 and adds it to an
    // fp64 accumulator, so the fp32 rounding of that product is part of the
    // matrix being solved. Reproduce it rather than computing in fp64.
    long long floored = 0;

    // Fast lane: one reciprocal per single-sample block, over flat arrays.
    const int n1 = (int)part_.one_slot.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:floored)
#endif
    for (int k = 0; k < n1; k++) {
        double v = part_.one_psi[(size_t)k] * tau1;
        const float dt = (1.0f / w((arma::uword)part_.one_sample[(size_t)k])) * tau0;
        v += (double)dt;
        if (v < 1e-4) { v = 1e-4; floored++; }
        inv_[(size_t)part_.one_slot[(size_t)k]] = 1.0 / v;
    }

    const int ng = (int)part_.genBlock.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) reduction(+:floored)
#endif
    for (int bi = 0; bi < ng; bi++) {
        const int b   = part_.genBlock[(size_t)bi];
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

        // Invert in place. 1 and 2 closed form, 3 and 4 closed form, [5,32] a
        // hand-rolled allocation-free Cholesky inverse, above that LAPACK via
        // arma::inv_sympd -- see invert_block and the crossover measured there.
        invert_block(A, s);
    }
    nFloored_ = floored;
    lastW_ = w; lastTau_ = tau;
    refreshed_ = true;
    refreshSecs_ += now_s() - t_begin;
}

void BlockSigma::traces(double* trSigmaInvPsi, double* trSigmaInv) const {
    double tp = 0.0, ti = 0.0;
    if (!ready()) { if (trSigmaInvPsi) *trSigmaInvPsi = 0.0;
                    if (trSigmaInv)    *trSigmaInv    = 0.0; return; }
    const int nb = part_.nblocks;

    // tr(AB) = sum_ij A_ij B_ji. Psi is symmetric and both triangles are
    // stored, so every non-zero (i,j) appears once and B_ji = B_ij = v.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:tp,ti)
#endif
    for (int b = 0; b < nb; b++) {
        const int s = part_.size[(size_t)b];
        const double* A = &inv_[(size_t)part_.invStart[(size_t)b]];
        for (int k = part_.psiStart[(size_t)b]; k < part_.psiStart[(size_t)b + 1]; k++)
            tp += A[(size_t)part_.psiC[(size_t)k] * s + part_.psiR[(size_t)k]]
                  * part_.psiV[(size_t)k];
        for (int r = 0; r < s; r++) ti += A[(size_t)r * s + r];
    }
    if (trSigmaInvPsi) *trSigmaInvPsi = tp;
    if (trSigmaInv)    *trSigmaInv    = ti;
}

arma::fmat BlockSigma::psiMultiply(const arma::fmat& X) const {
    arma::fmat out(X.n_rows, X.n_cols, arma::fill::zeros);
    if (!part_.built || (int)X.n_rows != part_.n) return out;
    const int nb = part_.nblocks;
    const int p  = (int)X.n_cols;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int b = 0; b < nb; b++) {
        const int mo = part_.start[(size_t)b];
        for (int k = part_.psiStart[(size_t)b]; k < part_.psiStart[(size_t)b + 1]; k++) {
            const int gi = part_.member[(size_t)(mo + part_.psiR[(size_t)k])];
            const int gj = part_.member[(size_t)(mo + part_.psiC[(size_t)k])];
            const double v = part_.psiV[(size_t)k];
            for (int c = 0; c < p; c++)
                out((arma::uword)gi, (arma::uword)c) +=
                    (float)(v * (double)X((arma::uword)gj, (arma::uword)c));
        }
    }
    return out;
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
