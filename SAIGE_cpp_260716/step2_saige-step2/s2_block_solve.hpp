// s2_block_solve.hpp -- step 2: explicit block-diagonal Sigma^-1 for the
// sparse-GRM variance path, in place of the per-marker PCG.
//
// Where this sits
// ---------------
// With a sparse GRM loaded, step 2 computes an exact var(T) for the markers
// whose first pass gives p < pval_cutoff_for_fastTest and whose MAC is at or
// below cateVarRatioMinMACVecExclude.back() (20.5 by default): scoreTest calls
// SAIGEClass::getPCG1ofSigmaAndGtilde, a Jacobi-preconditioned PCG on
// Sigma = tau1*K + diag(1/mu2 | tau0) to an absolute tolerance ||r||^2 < 0.02,
// once per such marker. Measured on a 20,000-marker MAC 5-20 set at N=50,000
// (GPU_WHOLE_PIPELINE.md S5): 21.9 ms thread-time per solve (quantitative),
// 38.5 ms (binary, includes the SPA rerun), on 4-6% of the markers.
//
// Sigma has exactly K's sparsity pattern, so after permuting by the connected
// components of K it is block diagonal, and on a FastSparseGRM-style GRM the
// blocks are small (34,687 blocks, largest 10, on the benchmark GRM). Step 1
// already solves this very matrix with a per-block dense inverse
// (block_sigma.hpp, 0.26-0.28 ms per solve); step 2 only lacked the
// partition. This file recovers it with the same class -- block_sigma.{hpp,cpp}
// here are a byte-for-byte copy of step 1's, see the note in block_sigma.cpp's
// commit -- and puts the solve behind the PCG's own entry point.
//
// What is different from step 1's use of the class
// -------------------------------------------------
// * Step 2 already holds Sigma, not Psi: null_model_loader scales K by tau1
//   and adds the diagonal term. So the class is fed Sigma as its "Psi" with
//   tau = (0, 1) and w = 1, which makes refresh() add exactly 0.0 to every
//   diagonal and invert Sigma itself. Its 1e-4 diagonal floor cannot fire on a
//   real Sigma (the diagonal is >= tau0 or >= 1/mu2 >= 4); if it does, init()
//   refuses, because the inverse would then be of a different matrix than the
//   one PCG solves.
// * The solve is fp64 end to end. getPCG1ofSigmaAndGtilde is arma::vec, and
//   BlockSigma::solve is fp32 in and out (step 1's Sigma^-1 v lives in
//   arma::fvec), so solve() below walks the class's public per-block inverse
//   (partition() / blockInverse()) in double instead. Nothing in the copied
//   file is changed.
// * One refresh per run per trait. There is no (w, tau) to cache across; the
//   cost gate (a wall-clock budget for that one refresh) is what refuses a GRM
//   with a component too large to invert densely, and then every solve stays
//   on PCG with the reason printed once.
// * Thread-safe by construction: solve() is const and serial. Step 2 calls it
//   from inside the marker loop's OpenMP region, one marker per thread, so a
//   parallel region inside the solve would only nest.
//
// Correctness expectation
// -----------------------
// PCG stops at ||r||^2 < 0.02; the block inverse is direct. The two are NOT
// expected to agree bitwise, and the block inverse is the more accurate. The
// verify switch (blockSparseSigmaVerify) runs both on every call, plus a PCG
// driven to ||r||^2 < 1e-18 as the reference, and reports each path's residual
// and distance to the reference, so that claim is measured rather than assumed.
#ifndef S2_BLOCK_SOLVE_HPP
#define S2_BLOCK_SOLVE_HPP

#include <armadillo>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "block_sigma.hpp"

namespace s2blk {

// YAML keys (top level, like every other step-2 switch):
//   blockSparseSigma:                 true/false   (default false)
//   blockSparseSigmaVerify:           true/false   (default false)
//   blockSparseSigmaRefreshBudget_s:  seconds      (default 0.25, step 1's)
//   blockSparseSigmaSolveMargin:      factor       (default 2, 0 disables)
//
// Two gates, because step 2 pays for a solve per marker, not per fit:
//   1. the class's own refresh gate -- the one-time inversion must fit the
//      seconds budget (a GRM with a 7,378-sample component is refused here);
//   2. a solve gate -- one block solve and one production PCG (100 iterations,
//      tol 0.02) are timed on the same gtilde-like probe at init, and the block
//      inverse is used only if it is at least `solveMargin` times faster. The
//      per-solve cost grows with sum(b^2), which gate 1 does not look at: a
//      1,616-sample component passes gate 1 (0.17 s estimated refresh) but its
//      solve is 4.4-4.8 ms against PCG's 3.8-4.1 ms on the probe, so it is
//      refused here. Measured (S2_BLOCK_SOLVE.md S4): on the benchmark GRM
//      (largest block 10) the solo numbers are 0.27 ms vs 5.1 ms and inside the
//      8-thread marker loop on this 4-core box they become 2.3 ms vs 19 ms --
//      the block solve inflates 9x under that contention, PCG 4x -- which is
//      why a solo tie is not trusted. The margin is conservative in the
//      middle: with the gate disabled (margin 0) the 1,616-sample case still
//      ran 1.2x faster end to end (7.3 ms vs 12-19 ms per solve in the loop).
struct Config {
    bool   enabled        = false;
    bool   verify         = false;
    double refreshBudgetS = 0.25;
    double solveMargin    = 2.0;
};
void          configure(const Config& c);
const Config& config();

class Solver {
public:
    // The production PCG on this Sigma, for the solve gate's timing probe.
    using PcgFn = std::function<arma::vec(const arma::vec&)>;

    // (loc, val) is Sigma exactly as SAIGEClass builds m_spSigmaMat from it:
    // tau1*K with the diagonal term already added. Prints what it built or why
    // it refused. Returns true when solve() may be used.
    bool init(const arma::umat& loc, const arma::vec& sigmaVal, int n,
              const std::string& label, const PcgFn& pcg);
    bool ready() const { return ready_; }

    // Sigma^-1 b in double, one block at a time. Serial, const, reentrant.
    arma::vec solve(const arma::vec& b) const;

    // Verify-mode accumulation for one call. `ref` is the tight-tolerance PCG
    // answer; residuals are computed here against Sigma.
    void recordVerify(const arma::sp_mat& Sigma, const arma::vec& b,
                      const arma::vec& pcg, const arma::vec& blk,
                      const arma::vec& ref, bool pcgConverged);
    void report() const;

    const std::string& label() const { return label_; }
    int    nblocks()   const { return bs_.partition().nblocks; }
    int    maxBlock()  const { return bs_.partition().maxBlock; }
    double buildSecs() const { return buildSecs_; }

private:
    blocksigma::BlockSigma bs_;
    std::string label_;
    bool   ready_     = false;
    double buildSecs_ = 0.0;   // union-find + gate estimate + one refresh
    double soloMs_ = 0.0, pcgMs_ = 0.0;   // the solve gate's two timings

    // The inverse again, laid out for the solve. BlockSigma stores each block
    // column-major and its solve reads A[c*s + r] with c innermost -- a stride
    // of s doubles, harmless at s <= 10, one cache miss per multiply at
    // s = 1,616 (28 ms per solve, all of it misses). Each general block is
    // copied here transposed so the inner loop is contiguous, and the
    // single-sample blocks become one flat (sample, reciprocal) lane. Same
    // values, same summation order, so solve() is bitwise what the direct
    // walk of blockInverse() gave; the copy costs sum(b^2) doubles once.
    std::vector<int>       oneSample_;
    std::vector<double>    oneInv_;
    std::vector<int>       genBlk_;
    std::vector<long long> genOff_;
    std::vector<double>    invT_;

    // Verify accumulators. Written under a critical section: verify is the
    // correctness gate, not a mode anyone benchmarks.
    long long vCalls_ = 0, vPcgNonConv_ = 0, vBlockCloser_ = 0;
    double vMaxAbs_ = 0.0, vMaxRelNorm_ = 0.0;
    double vResPcgMax_ = 0.0, vResPcgSum_ = 0.0;
    double vResBlkMax_ = 0.0, vResBlkSum_ = 0.0;
    double vResRefMax_ = 0.0, vResRefSum_ = 0.0;
    double vErrPcgMax_ = 0.0, vErrPcgSum_ = 0.0;
    double vErrBlkMax_ = 0.0, vErrBlkSum_ = 0.0;
    // var2 = b' Sigma^-1 b is the number that reaches the test statistic.
    double vVarPcgVsBlkMax_ = 0.0, vVarPcgVsRefMax_ = 0.0, vVarBlkVsRefMax_ = 0.0;
};

// Solve-site timing, per OpenMP thread, always on: one steady_clock pair per
// solve against a solve that costs 0.2-40 ms.
void noteSolve(bool block, double secs);
void noteVerifyPcg(double secs);   // the PCG run only for comparison

// The SAIGEClass that built a solver registers it once ready; reportAll()
// prints every registered solver and the timing table at the end of the run
// (silent when nothing took the path). The registry holds a shared_ptr
// because main() deletes the SAIGEClass instances -- and with them their
// m_blockSolver -- before it reaches the final report.
void registerSolver(std::shared_ptr<const Solver> s);
void reportAll();

}  // namespace s2blk

#endif  // S2_BLOCK_SOLVE_HPP
