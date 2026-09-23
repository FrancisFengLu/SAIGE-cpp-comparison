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
    bool built = false;
};

class BlockSigma {
public:
    // Build the partition from SAIGE's sparse-GRM globals. Returns false (and
    // leaves the object unusable) if the input is empty or inconsistent, in
    // which case the caller must fall back to gen_spsolve_v4.
    // A connected component is not a clique: on a stress GRM with components of
    // 1,000 members only 0.90% of the within-block entries are non-zero, so a
    // dense inverse of such a block spends 2e10 flops on a matrix that is 99%
    // zeros. The right gate is therefore the total sum(b^3), not the largest
    // block. build() refuses when that budget (or the sum(b^2) storage budget)
    // is exceeded, and the caller falls back to gen_spsolve_v4.
    bool build(const arma::umat& loc, const arma::vec& val, int n);
    void setBudget(double flopBudget, double byteBudget);
    double lastFlops() const { return flops_; }
    double lastBytes() const { return bytes_; }

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
    std::vector<double> scratch_;    // per-block dense Sigma during refresh
    bool refreshed_ = false;
    long long nFloored_ = 0;
    arma::fvec lastW_, lastTau_;
    long long nRefresh_ = 0, nReuse_ = 0;
    double flopBudget_ = 1e9;        // ~0.3 s per refresh on this machine
    double byteBudget_ = 2e9;        // 2 GB of block inverses
    double flops_ = 0.0, bytes_ = 0.0;

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
