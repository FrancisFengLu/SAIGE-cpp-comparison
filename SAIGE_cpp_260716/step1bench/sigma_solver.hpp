// sigma_solver.hpp -- one interface, three implementations of "invert Sigma".
//
// This is a benchmark instrument. The whole point is that the AI-REML driver
// above it makes exactly the same calls in exactly the same order whichever
// solver is plugged in, so the only thing that varies between runs is HOW
// Sigma is inverted. Everything that could differ and must not -- how Sigma is
// formed, which right-hand sides are solved, which probe vectors are drawn,
// how the clock is read -- lives here or in the driver, never in a solver.
#ifndef STEP1BENCH_SIGMA_SOLVER_HPP
#define STEP1BENCH_SIGMA_SOLVER_HPP

#include <armadillo>
#include <memory>
#include <string>

namespace step1bench {

// ---------------------------------------------------------------------------
// Sigma, formed exactly as SAIGE's gen_sp_Sigma forms it
// (step1_saige-null/SAIGE_step1_fast.cpp:6256).
//
//   dtVec        = (1/wVec) * tau0          <- computed in fp32
//   valueVecNew  = valueVec * tau1          <- fp64, tau1 promoted from fp32
//   diagonal     += dtVec(i)                <- fp32 value added into an fp64 slot
//   if (diagonal < 1e-4) diagonal = 1e-4    <- floor applied AFTER the add
//
// The fp32 rounding of tau0/w_i and the post-add floor are both part of the
// matrix being solved, not incidental. block_sigma.cpp::refresh() mirrors them
// (block_sigma.cpp:413-478); this function is what the other two solvers use,
// and the startup self-check compares all three.
//
// `loc` must carry BOTH triangles, as SAIGE's locationMat does.
arma::sp_mat gen_sp_Sigma(const arma::umat& loc, const arma::vec& val, int n,
                          const arma::fvec& w, const arma::fvec& tau);

// ---------------------------------------------------------------------------
// Stage clocks. Every solver is wrapped in TimedSolver, so the timing code is
// literally the same object for all three -- a solver cannot time itself
// better or worse than its rivals.
enum Stage { ST_BUILD = 0, ST_REFRESH, ST_SOLVE, ST_TRACE, ST_N };

struct StageAcc {
    long long calls[ST_N] = {0, 0, 0, 0};
    double    secs [ST_N] = {0.0, 0.0, 0.0, 0.0};
};

double now_s();

// Which stage solve() calls are charged to. The driver sets it; a solver never
// touches it.
void  set_stage(Stage s);
Stage cur_stage();

struct StageScope {
    Stage prev;
    explicit StageScope(Stage s) : prev(cur_stage()) { set_stage(s); }
    ~StageScope() { set_stage(prev); }
};

// ---------------------------------------------------------------------------
struct SigmaSolver {
    virtual ~SigmaSolver() = default;

    // Once per run. `loc` carries both triangles; `val` the matching values.
    virtual bool build(const arma::umat& loc, const arma::vec& val, int n) = 0;

    // Per (w, tau). Called before EVERY solve, exactly as SAIGE's
    // gen_spsolve_v4 calls blocksigma::refresh on every one of its 199-261
    // calls. Solvers that cache return immediately when (w, tau) repeat and
    // report that as a hit; solvers that do not, do not.
    virtual void refresh(const arma::fvec& w, const arma::fvec& tau) = 0;

    virtual arma::fvec solve(const arma::fvec& b) = 0;

    // True when the solver holds Sigma^-1 entry by entry and can therefore
    // produce tr(Sigma^-1 Psi) and tr(Sigma^-1) without probing.
    virtual bool hasExplicitInverse() const = 0;

    // Only meaningful when hasExplicitInverse(). The SuperLU implementation
    // aborts with a message rather than returning a number it does not have.
    virtual void traces(double* trSigmaInvPsi, double* trSigmaInv) = 0;

    virtual const char* name() const = 0;

    // --- instrument hooks, not part of the mathematics --------------------

    // The Sigma this solver would form at (w, tau), as a sparse matrix. Used
    // once at startup by the fairness self-check. Not on any timed path.
    virtual arma::sp_mat debugSigma(const arma::fvec& w,
                                    const arma::fvec& tau) const = 0;

    // (w, tau) cache accounting. Solvers with no cache report 0/0.
    virtual long long refreshCount() const { return 0; }
    virtual long long reuseCount()   const { return 0; }

    // Free-form extra line for the report (e.g. SuperLU's per-solve rebuild
    // cost, which lives inside solve() and would otherwise be invisible).
    virtual std::string extraReport() const { return std::string(); }
};

// Wraps any solver and charges its time to the stage clocks. Construction and
// destruction are the only places timing logic exists in this program.
class TimedSolver : public SigmaSolver {
public:
    explicit TimedSolver(std::unique_ptr<SigmaSolver> inner)
        : in_(std::move(inner)) {}

    bool build(const arma::umat& loc, const arma::vec& val, int n) override {
        const double t = now_s();
        const bool ok = in_->build(loc, val, n);
        acc_.calls[ST_BUILD]++; acc_.secs[ST_BUILD] += now_s() - t;
        return ok;
    }
    void refresh(const arma::fvec& w, const arma::fvec& tau) override {
        const double t = now_s();
        in_->refresh(w, tau);
        acc_.calls[ST_REFRESH]++; acc_.secs[ST_REFRESH] += now_s() - t;
    }
    arma::fvec solve(const arma::fvec& b) override {
        const Stage s = cur_stage() == ST_TRACE ? ST_TRACE : ST_SOLVE;
        const double t = now_s();
        arma::fvec x = in_->solve(b);
        acc_.calls[s]++; acc_.secs[s] += now_s() - t;
        return x;
    }
    bool hasExplicitInverse() const override { return in_->hasExplicitInverse(); }
    void traces(double* a, double* b) override {
        const double t = now_s();
        in_->traces(a, b);
        acc_.calls[ST_TRACE]++; acc_.secs[ST_TRACE] += now_s() - t;
    }
    const char* name() const override { return in_->name(); }
    arma::sp_mat debugSigma(const arma::fvec& w, const arma::fvec& tau) const override {
        return in_->debugSigma(w, tau);
    }
    long long refreshCount() const override { return in_->refreshCount(); }
    long long reuseCount()   const override { return in_->reuseCount(); }
    std::string extraReport() const override { return in_->extraReport(); }

    const StageAcc& acc() const { return acc_; }
    SigmaSolver& inner() { return *in_; }

private:
    std::unique_ptr<SigmaSolver> in_;
    StageAcc acc_;
};

// Connected-component summary of the sparsity pattern. Reported for every
// solver (it is a property of the GRM, not of the solver), because it is the
// thing that decides whether the block path is applicable at all.
struct BlockSummary {
    int nblocks = 0;
    int maxBlock = 0;
    std::vector<std::pair<int, long long>> hist;   // (size, count), ascending
};
BlockSummary summarize_blocks(const arma::umat& loc, int n);

// name -> solver. Returns nullptr for an unknown name.
std::unique_ptr<SigmaSolver> make_solver(const std::string& name);

}  // namespace step1bench

#endif  // STEP1BENCH_SIGMA_SOLVER_HPP
