#include "s2_block_solve.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <mutex>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace s2blk {

namespace {

Config g_cfg;

double now_s() {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}

int thread_id() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// One cache line per thread so the marker loop's threads never share a
// counter. Same idiom as phase_timing.hpp.
struct alignas(128) TimeSlot {
    double    tBlock = 0.0, tPcg = 0.0, tVerifyPcg = 0.0;
    long long nBlock = 0,   nPcg = 0,   nVerifyPcg = 0;
    char      pad[128 - 6 * 8];
};
TimeSlot g_time[512];

std::mutex                                g_reg_mutex;
std::vector<std::shared_ptr<const Solver>> g_solvers;

}  // namespace

void          configure(const Config& c) { g_cfg = c; }
const Config& config()                   { return g_cfg; }

void noteSolve(bool block, double secs) {
    TimeSlot& s = g_time[thread_id() & 511];
    if (block) { s.tBlock += secs; s.nBlock++; }
    else       { s.tPcg   += secs; s.nPcg++;   }
}

void noteVerifyPcg(double secs) {
    TimeSlot& s = g_time[thread_id() & 511];
    s.tVerifyPcg += secs; s.nVerifyPcg++;
}

void registerSolver(std::shared_ptr<const Solver> s) {
    if (!s) return;
    std::lock_guard<std::mutex> lk(g_reg_mutex);
    g_solvers.push_back(std::move(s));
}

bool Solver::init(const arma::umat& loc, const arma::vec& sigmaVal, int n,
                  const std::string& label, const PcgFn& pcg) {
    label_ = label;
    ready_ = false;
    const double t0 = now_s();

    // Every sample must carry its diagonal term, or its block is singular
    // (a size-1 block with no entry inverts to 1/0). The loader warns when
    // relatednessCutoff removes diagonal entries; refuse here rather than
    // hand back inf.
    if (loc.n_rows != 2 || loc.n_cols != sigmaVal.n_elem || n <= 0) {
        printf("[s2blocksolve] %s: sparse Sigma is malformed (%llu x %llu locations, "
               "%llu values); staying on PCG.\n", label_.c_str(),
               (unsigned long long)loc.n_rows, (unsigned long long)loc.n_cols,
               (unsigned long long)sigmaVal.n_elem);
        fflush(stdout);
        return false;
    }
    {
        std::vector<char> hasDiag((size_t)n, 0);
        long long nonfinite = 0;
        for (arma::uword k = 0; k < loc.n_cols; k++) {
            const arma::uword r = loc(0, k), c = loc(1, k);
            if (r >= (arma::uword)n || c >= (arma::uword)n) {
                printf("[s2blocksolve] %s: sparse Sigma has an index >= n; staying on PCG.\n",
                       label_.c_str());
                fflush(stdout);
                return false;
            }
            if (r == c) hasDiag[(size_t)r] = 1;
            if (!std::isfinite(sigmaVal(k))) nonfinite++;
        }
        long long missing = 0;
        for (int i = 0; i < n; i++) if (!hasDiag[(size_t)i]) missing++;
        if (missing > 0 || nonfinite > 0) {
            printf("[s2blocksolve] %s: %lld samples have no diagonal entry in the sparse "
                   "Sigma and %lld entries are not finite; a block inverse would be "
                   "singular there. Staying on PCG.\n", label_.c_str(), missing, nonfinite);
            fflush(stdout);
            return false;
        }
    }

    bs_.setRefreshBudget(g_cfg.refreshBudgetS);
    if (!bs_.build(loc, sigmaVal, n)) {
        // build() has printed the [blocksigma] refusal (block count, largest
        // block, estimated seconds against the budget). Its last words are
        // step 1's "Falling back to spsolve"; here that means PCG.
        printf("[s2blocksolve] %s: block inverse refused by the cost gate "
               "(blockSparseSigmaRefreshBudget_s = %.3f); every sparse-Sigma "
               "solve for this trait stays on PCG.\n",
               label_.c_str(), g_cfg.refreshBudgetS);
        fflush(stdout);
        return false;
    }

    // Sigma is already formed. tau = (0, 1) and w = 1 make refresh() compute
    // v = psiV * 1.0 + (double)((1.0f / 1.0f) * 0.0f) = psiV exactly.
    arma::fvec w((arma::uword)n, arma::fill::ones);
    arma::fvec tau(2);
    tau(0) = 0.0f; tau(1) = 1.0f;
    bs_.refresh(w, tau);
    if (bs_.flooredDiagonals() > 0) {
        printf("[s2blocksolve] %s: the 1e-4 diagonal floor fired on %lld entries -- "
               "the block inverse would be of a different matrix than PCG solves. "
               "Staying on PCG.\n", label_.c_str(), bs_.flooredDiagonals());
        fflush(stdout);
        bs_.reset();
        return false;
    }
    if (!bs_.ready()) {
        printf("[s2blocksolve] %s: block inverse not ready after refresh; staying on PCG.\n",
               label_.c_str());
        fflush(stdout);
        return false;
    }
    const blocksigma::Partition& P = bs_.partition();

    // --- the solve's own layout: flat lane for size-1, transposed copies ---
    oneSample_.clear(); oneInv_.clear(); genBlk_.clear(); genOff_.clear(); invT_.clear();
    {
        long long total = 0;
        for (int blk = 0; blk < P.nblocks; blk++) {
            const int s = P.size[(size_t)blk];
            if (s == 1) {
                oneSample_.push_back(P.member[(size_t)P.start[(size_t)blk]]);
                oneInv_.push_back(bs_.blockInverse(blk)[0]);
            } else {
                genBlk_.push_back(blk);
                genOff_.push_back(total);
                total += (long long)s * s;
            }
        }
        invT_.resize((size_t)total);
        for (size_t j = 0; j < genBlk_.size(); j++) {
            const int s = P.size[(size_t)genBlk_[j]];
            const double* A = bs_.blockInverse(genBlk_[j]);
            double* T = &invT_[(size_t)genOff_[j]];
            for (int r = 0; r < s; r++)
                for (int c = 0; c < s; c++) T[(size_t)r * s + c] = A[(size_t)c * s + r];
        }
    }
    buildSecs_ = now_s() - t0;
    ready_ = true;

    // --- solve gate: one block solve against one production PCG ----------
    // Probe shaped like the right-hand side this path sees: gtilde for a
    // MAC-10 marker, i.e. ten unit entries on a small negative background
    // (centering). Deterministic, no RNG touched. PCG's stop is an absolute
    // ||r||^2 < 0.02, so the probe's scale matters and this one matches the
    // production scale; the block solve does not care. Both timed alone,
    // cache-hot, single thread, so the numbers printed at the end of the run
    // (inside the 8-thread marker loop) can be read against them.
    {
        arma::vec probe((arma::uword)n);
        probe.fill(-10.0 / (double)n);
        for (int k = 0; k < 10; k++) probe((arma::uword)((long long)k * n / 10)) = 1.0;
        arma::vec y = solve(probe);              // warm
        const int reps = 20;
        double ts = now_s();
        for (int r = 0; r < reps; r++) { y = solve(probe); probe(0) += y(0) * 1e-300; }
        soloMs_ = 1e3 * (now_s() - ts) / reps;
        if (pcg) {
            arma::vec z = pcg(probe);            // warm
            const int preps = (soloMs_ > 5.0) ? 3 : 5;
            ts = now_s();
            for (int r = 0; r < preps; r++) { z = pcg(probe); probe(0) += z(0) * 1e-300; }
            pcgMs_ = 1e3 * (now_s() - ts) / preps;
        }
    }
    if (pcg && g_cfg.solveMargin > 0.0 && soloMs_ * g_cfg.solveMargin >= pcgMs_) {
        printf("[s2blocksolve] %s: block inverse built (%d blocks, largest %d) but "
               "refused by the solve gate: one block solve %.3f ms vs one PCG %.3f ms "
               "on the same probe, margin %.1f required (blockSparseSigmaSolveMargin). "
               "Every sparse-Sigma solve for this trait stays on PCG.\n",
               label_.c_str(), P.nblocks, P.maxBlock, soloMs_, pcgMs_, g_cfg.solveMargin);
        fflush(stdout);
        ready_ = false;
        bs_.reset();
        invT_.clear(); oneInv_.clear();
        return false;
    }
    printf("[s2blocksolve] %s: sparse-Sigma solves use the block inverse: %d blocks, "
           "largest %d, %zu single-sample, built in %.1f ms (union-find + gate "
           "estimate + inverse; the inverse alone %.1f ms); one solve alone, "
           "cache-hot, single thread: block %.3f ms vs PCG %.3f ms (margin %.1f).\n",
           label_.c_str(), P.nblocks, P.maxBlock, oneSample_.size(),
           1e3 * buildSecs_, 1e3 * bs_.refreshSeconds(), soloMs_, pcgMs_, g_cfg.solveMargin);
    fflush(stdout);
    return true;
}

arma::vec Solver::solve(const arma::vec& b) const {
    const blocksigma::Partition& P = bs_.partition();
    if (!ready_ || (int)b.n_elem != P.n)
        return arma::vec(b.n_elem, arma::fill::zeros);
    arma::vec out(b.n_elem);      // every element is written: the blocks
                                  // partition all n samples
    const double* __restrict bp = b.memptr();
    double*       __restrict op = out.memptr();

    // Single-sample blocks: one reciprocal each, flat arrays.
    const int    n1 = (int)oneSample_.size();
    const int*    __restrict os = oneSample_.data();
    const double* __restrict oi = oneInv_.data();
    for (int k = 0; k < n1; k++) { const int i = os[k]; op[i] = oi[k] * bp[i]; }

    // Everything else: dense (Sigma^-1)_block times the block's entries of b,
    // fp64 throughout, inner sum over c ascending, rows of the transposed
    // copy so the inner loop walks memory contiguously.
    const int* mem = P.member.data();
    const int  ng  = (int)genBlk_.size();
    for (int j = 0; j < ng; j++) {
        const int blk = genBlk_[(size_t)j];
        const int s   = P.size[(size_t)blk];
        const int* m  = mem + P.start[(size_t)blk];
        const double* __restrict T = &invT_[(size_t)genOff_[(size_t)j]];
        for (int r = 0; r < s; r++) {
            const double* __restrict Tr = T + (size_t)r * s;
            double acc = 0.0;
            for (int c = 0; c < s; c++) acc += Tr[c] * bp[m[c]];
            op[m[r]] = acc;
        }
    }
    return out;
}

void Solver::recordVerify(const arma::sp_mat& Sigma, const arma::vec& b,
                          const arma::vec& pcg, const arma::vec& blk,
                          const arma::vec& ref, bool pcgConverged) {
    const double nb = arma::norm(b);
    if (!(nb > 0.0)) return;
    const double resPcg = arma::norm(Sigma * pcg - b) / nb;
    const double resBlk = arma::norm(Sigma * blk - b) / nb;
    const double resRef = arma::norm(Sigma * ref - b) / nb;
    const double nref   = arma::norm(ref);
    const double errPcg = (nref > 0.0) ? arma::norm(pcg - ref) / nref : 0.0;
    const double errBlk = (nref > 0.0) ? arma::norm(blk - ref) / nref : 0.0;
    const double maxAbs = arma::abs(pcg - blk).max();
    const double nblk   = arma::norm(blk);
    const double relNrm = (nblk > 0.0) ? arma::norm(pcg - blk) / nblk : 0.0;
    const double vP = arma::dot(pcg, b), vB = arma::dot(blk, b), vR = arma::dot(ref, b);
    const double dPB = (vB != 0.0) ? std::fabs(vP - vB) / std::fabs(vB) : 0.0;
    const double dPR = (vR != 0.0) ? std::fabs(vP - vR) / std::fabs(vR) : 0.0;
    const double dBR = (vR != 0.0) ? std::fabs(vB - vR) / std::fabs(vR) : 0.0;
#ifdef _OPENMP
#pragma omp critical(s2blk_verify)
#endif
    {
        vCalls_++;
        if (!pcgConverged) vPcgNonConv_++;
        if (errBlk < errPcg) vBlockCloser_++;
        if (maxAbs > vMaxAbs_)     vMaxAbs_     = maxAbs;
        if (relNrm > vMaxRelNorm_) vMaxRelNorm_ = relNrm;
        if (resPcg > vResPcgMax_)  vResPcgMax_  = resPcg;  vResPcgSum_ += resPcg;
        if (resBlk > vResBlkMax_)  vResBlkMax_  = resBlk;  vResBlkSum_ += resBlk;
        if (resRef > vResRefMax_)  vResRefMax_  = resRef;  vResRefSum_ += resRef;
        if (errPcg > vErrPcgMax_)  vErrPcgMax_  = errPcg;  vErrPcgSum_ += errPcg;
        if (errBlk > vErrBlkMax_)  vErrBlkMax_  = errBlk;  vErrBlkSum_ += errBlk;
        if (dPB > vVarPcgVsBlkMax_) vVarPcgVsBlkMax_ = dPB;
        if (dPR > vVarPcgVsRefMax_) vVarPcgVsRefMax_ = dPR;
        if (dBR > vVarBlkVsRefMax_) vVarBlkVsRefMax_ = dBR;
    }
}

void Solver::report() const {
    if (vCalls_ == 0) return;
    const double c = (double)vCalls_;
    printf("[s2blocksolve] verify %s: %lld solves, PCG not converged on %lld\n",
           label_.c_str(), vCalls_, vPcgNonConv_);
    printf("  PCG vs block:      max |x_pcg - x_blk| %.3e, max ||diff||/||x_blk|| %.3e, "
           "max |var2 diff|/var2 %.3e\n", vMaxAbs_, vMaxRelNorm_, vVarPcgVsBlkMax_);
    printf("  residual ||Sx-b||/||b||:  PCG max %.3e mean %.3e | block max %.3e mean %.3e "
           "| reference max %.3e mean %.3e\n",
           vResPcgMax_, vResPcgSum_ / c, vResBlkMax_, vResBlkSum_ / c,
           vResRefMax_, vResRefSum_ / c);
    printf("  vs reference ||x-x_ref||/||x_ref||:  PCG max %.3e mean %.3e | block max %.3e "
           "mean %.3e | block closer on %lld/%lld\n",
           vErrPcgMax_, vErrPcgSum_ / c, vErrBlkMax_, vErrBlkSum_ / c,
           vBlockCloser_, vCalls_);
    printf("  var2 = b'S^-1 b vs reference:  PCG max rel %.3e | block max rel %.3e\n",
           vVarPcgVsRefMax_, vVarBlkVsRefMax_);
    fflush(stdout);
}

void reportAll() {
    double tB = 0, tP = 0, tV = 0;
    long long nB = 0, nP = 0, nV = 0;
    for (int i = 0; i < 512; i++) {
        tB += g_time[i].tBlock; nB += g_time[i].nBlock;
        tP += g_time[i].tPcg;   nP += g_time[i].nPcg;
        tV += g_time[i].tVerifyPcg; nV += g_time[i].nVerifyPcg;
    }
    if (nB == 0 && nP == 0) return;
    printf("[s2blocksolve] sparse-Sigma solves this run (thread-time, summed over threads):\n");
    if (nP > 0)
        printf("  PCG:           %lld solves, %.3f s total, %.3f ms each\n", nP, tP, 1e3 * tP / (double)nP);
    if (nB > 0)
        printf("  block inverse: %lld solves, %.3f s total, %.4f ms each\n", nB, tB, 1e3 * tB / (double)nB);
    if (nV > 0)
        printf("  PCG (verify, same right-hand sides): %lld solves, %.3f s total, %.3f ms each\n",
               nV, tV, 1e3 * tV / (double)nV);
    std::lock_guard<std::mutex> lk(g_reg_mutex);
    for (const auto& s : g_solvers) {
        printf("  %s: %d blocks, largest %d, build %.1f ms\n",
               s->label().c_str(), s->nblocks(), s->maxBlock(), 1e3 * s->buildSecs());
        s->report();
    }
    fflush(stdout);
}

}  // namespace s2blk
