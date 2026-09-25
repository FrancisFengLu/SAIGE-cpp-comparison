// block_sigma.hpp -- explicit block-diagonal Sigma^-1 for the sparse-GRM path.
//
// Why this exists
// ---------------
// On the sparse-GRM path Sigma = tau0*diag(1/w) + tau1*Psi is solved by
// gen_spsolve_v4, which on every call rebuilds the whole n x n sp_mat from all
// nnz and then runs a fresh SuperLU symbolic + numeric factorisation. Measured
// (fit.profile_spsolve, N=50,000, 114,878 nnz): 30-34 ms per call, of which
// ~8 ms is the rebuild and ~25 ms the factor+solve, and a single-trait fit
// makes 199 (quantitative) or 261 (binary) such calls. That is 93% of the time
// GLMM fitting plus variance-ratio estimation take.
//
// None of that work depends on the sparsity pattern, which never changes.
// Sigma has exactly Psi's pattern -- tau0*diag(1/w) is diagonal -- so after
// permuting by the connected components of Psi it is block diagonal. On the
// benchmark GRM that is 34,687 blocks with a largest size of 10: inverting
// every block densely is 745,610 flops and 1.1 MB of storage, against 30 ms
// per SuperLU solve.
//
// What this class gives that a factorisation does not
// ---------------------------------------------------
// Sigma^-1 itself, entry by entry, at exactly the positions where Psi is
// non-zero. That is what turns the AI-REML trace from 30 Hutchinson probes
// (each paying a full solve) into an exact elementwise contraction against
// Psi. This header only exposes the solver; the trace is a separate change
// with a separate switch, because it changes results and the solver must not.
//
// Binary traits are not a special case: W^-1 is diagonal, so the partition is
// identical. The only difference is that W changes at every IRLS iteration, so
// refresh() runs more often -- which is the case block inversion handles well
// and a global factorisation handles badly.
#ifndef BLOCK_SIGMA_HPP
#define BLOCK_SIGMA_HPP

#include <armadillo>
#include <vector>

namespace blocksigma {

// Connected-component partition of the sparse GRM. Depends only on the
// sparsity pattern, so it is built once per run and survives every change of
// w and tau.
struct Partition {
    int n = 0;                       // matrix dimension
    int nblocks = 0;
    int maxBlock = 0;
    std::vector<int> blockOf;        // sample -> block id
    std::vector<int> localOf;        // sample -> index within its block
    std::vector<int> start;          // block -> offset into member
    std::vector<int> size;           // block -> number of members
    std::vector<int> member;         // members, grouped by block
    // Psi in per-block local coordinates. Both triangles are kept, exactly as
    // they appear in locationMat/valueVec, so refresh() can walk them without
    // deciding which triangle it has.
    std::vector<int>    psiStart;    // block -> offset into psiR/psiC/psiV
    std::vector<int>    psiR, psiC;  // local row/col inside the block
    std::vector<double> psiV;
    std::vector<int>    invStart;    // block -> offset into the dense inverse
    // Size-1 fast lane. 86% of the blocks on the benchmark GRM hold a single
    // sample, and for those the whole "form the block and invert it" reduces to
    // one reciprocal. Precomputing (sample, Psi_ii, slot) lets refresh() walk a
    // flat array instead of re-reading size/start/psiStart per block. Only
    // size-1 blocks whose Psi contribution is exactly one diagonal entry are
    // listed here; anything odd (no entry, or a duplicate) stays on the general
    // path so the arithmetic cannot drift.
    std::vector<int>    one_sample;  // sample index
    std::vector<double> one_psi;     // Psi_ii
    std::vector<int>    one_slot;    // offset into inv_
    std::vector<int>    genBlock;    // block ids NOT on the fast lane
    bool built = false;
};

class BlockSigma {
public:
    // Build the partition from SAIGE's sparse-GRM globals. Returns false (and
    // leaves the object unusable) if the input is empty or inconsistent, in
    // which case the caller must fall back to gen_spsolve_v4.
    // Cost gate. The old gate compared sum(b^3) against a bare flop budget whose
    // default (5e8) refused a single 800-sample block -- a number with no
    // operational meaning, and wrong besides: FastSparseGRM's blocks are dense
    // cliques (BLOCK_INVERSE_ALGOS.md S1), so a large block is a genuinely dense
    // inverse, not a sparse matrix being mishandled. The gate is now a WALL-CLOCK
    // budget for one refresh: build() times a real inverse at the largest block
    // size with the same dispatch refresh() will use, extrapolates by sum(b^3),
    // divides by the threads refresh() can actually use, and refuses when the
    // estimate exceeds fit.block_sparse_sigma_refresh_budget_s. The storage
    // budget (sum(b^2) bytes) still applies.
    bool build(const arma::umat& loc, const arma::vec& val, int n);
    // True once build() has refused; the caller must not retry. Without this
    // every subsequent solve re-ran the union-find over all nnz and reprinted
    // the refusal -- 200+ times in a single fit.
    bool buildRefused() const { return buildFailed_; }

    // Drop the partition and the cached inverse. Must be called whenever the
    // sample set changes: the partition indexes samples by position in THAT
    // sample set, so reusing it across a sample-set group is a crash waiting to
    // happen (refresh() has no dimension guard, unlike solve()).
    void reset() { *this = BlockSigma(); }
    void setBudget(double flopBudget, double byteBudget);
    // Seconds allowed for ONE refresh (fit.block_sparse_sigma_refresh_budget_s).
    void setRefreshBudget(double seconds);
    double lastFlops() const { return flops_; }
    double lastBytes() const { return bytes_; }
    // Estimated and measured cost of one refresh, in seconds. The estimate is
    // what the gate used; the measurement is what it actually cost. Both are
    // printed so a bad estimate is visible rather than silent.
    double estRefreshSeconds() const { return estSecs_; }
    double refreshSeconds() const { return refreshSecs_; }

    // Form Sigma for this (w, tau) and invert it block by block. Mirrors
    // gen_sp_Sigma's arithmetic exactly, including the 1e-4 floor that is
    // applied to a diagonal entry after tau0/w_i has been added to it, and the
    // fp32 computation of tau0/w_i.
    void refresh(const arma::fvec& w, const arma::fvec& tau);

    // Sigma^-1 b. Requires a refresh() with the same (w, tau) first.
    arma::fvec solve(const arma::fvec& b) const;

    // True when the last refresh() already used this (w, tau). The trace
    // probes are 30 solves against one unchanged Sigma, so skipping the
    // re-inversion there is most of what makes this cheap.
    bool upToDate(const arma::fvec& w, const arma::fvec& tau) const;

    // Exact traces from the block inverse. tr(Sigma^-1 Psi) needs Sigma^-1
    // only where Psi is non-zero, which for a block-diagonal Sigma is exactly
    // the within-block entries we already hold -- no selected inverse needed.
    // Accumulated in fp64 over ~10^5 terms.
    void traces(double* trSigmaInvPsi, double* trSigmaInv) const;

    // Psi * X, from the same per-block storage, so the trace correction term
    // does not depend on getCrossprodMatAndKin's scaling conventions.
    arma::fmat psiMultiply(const arma::fmat& X) const;

    bool ready() const { return part_.built && refreshed_; }
    const Partition& partition() const { return part_; }

    // Raw access to the dense inverse of one block, column-major, size x size.
    // Present because the exact-trace work needs Sigma^-1 at Psi's non-zero
    // positions, which are exactly the within-block positions.
    const double* blockInverse(int b) const { return &inv_[part_.invStart[b]]; }

    // True when the last refresh() hit the 1e-4 diagonal floor anywhere. The
    // floor is a guard against a non-positive diagonal; if it fires, the block
    // inverse and SuperLU are solving genuinely different matrices only if the
    // floor changes the value, so this is reported rather than hidden.
    long long flooredDiagonals() const { return nFloored_; }

private:
    Partition part_;
    std::vector<double> inv_;        // dense inverse of each block, concatenated
    // The size-1 blocks' reciprocals again, contiguous and in one_sample order,
    // so solve() reads them sequentially instead of gathering them out of inv_
    // (where they sit interleaved with the real blocks). 8 bytes per size-1
    // block -- 240 KB on the benchmark GRM.
    std::vector<double> oneInv_;
    std::vector<double> scratch_;    // per-block dense Sigma during refresh
    bool refreshed_ = false;
    long long nFloored_ = 0;
    arma::fvec lastW_, lastTau_;
    long long nRefresh_ = 0, nReuse_ = 0;
    // Measured scaling of the whole sparse-solve stage, same data, same call
    // count, block inverse off vs on (N=50,000, P=1, quantitative):
    //   max block  10  sum(b^3) 7.5e5   56x
    //   max block  49  sum(b^3) 3.1e7   9.6x   (1.757 -> 0.184 s)
    //   max block 199  sum(b^3) 4.6e8   1.6x   (1.786 -> 1.126 s)
    //   max block 999  sum(b^3) 1.2e10  a loss
    // The gain decays smoothly and turns negative between the last two, which is
    // what the seconds budget is calibrated against.
    // flopBudget_ <= 0 means "no raw flop gate"; the seconds gate below is the
    // real one. Kept as an escape hatch for reproducing the old behaviour.
    double flopBudget_ = 0.0;
    double byteBudget_ = 2e9;        // 2 GB of block inverses
    // Seconds per refresh. Calibrated against the measured cost of the four test
    // GRMs; see SMALL_BLOCK_INVERSE.md S4.
    double refreshBudget_ = 0.25;
    double flops_ = 0.0, bytes_ = 0.0;
    // sum(b^2) over the blocks that are NOT on the size-1 lane: the work the
    // general pass of solve() does, and what decides whether it is worth a
    // parallel region.
    double genFlops_ = 0.0;
    // Times one real inverse at the largest block size and scales by sum(b^3).
    double estimateRefreshSeconds() const;
    double estSecs_ = 0.0;           // gate's estimate for one refresh
    double refreshSecs_ = 0.0;       // measured, summed over all refreshes
    bool buildFailed_ = false;

public:
    long long refreshCount() const { return nRefresh_; }
    long long reuseCount()   const { return nReuse_; }
};

// The process-wide instance used by gen_spsolve_v4 when
// fit.block_sparse_sigma is on. Lives here rather than in the caller so the
// partition is built once even though gen_spsolve_v4 is called from three
// different places.
BlockSigma& instance();

// fit.block_sparse_sigma: use the block inverse instead of SuperLU.
void enable(bool on);
bool enabled();

// fit.block_sparse_sigma_verify: run BOTH paths on every call and accumulate
// the difference. Slow by construction -- it is the correctness gate, not a
// mode anyone should benchmark.
void enableVerify(bool on);
bool verifyEnabled();
void recordVerify(const arma::fvec& ref, const arma::fvec& got);
void reportVerify(const char* what);
void resetVerify();

}  // namespace blocksigma

#endif  // BLOCK_SIGMA_HPP
