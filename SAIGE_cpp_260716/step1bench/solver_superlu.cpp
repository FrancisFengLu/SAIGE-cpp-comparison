// solver_superlu.cpp -- SAIGE's shape, and the baseline every other number in
// this harness is measured against.
//
// Mirrors gen_spsolve_v4 (step1_saige-null/SAIGE_step1_fast.cpp:6581) with its
// caches switched off, which is what upstream SAIGE-GPU 1.3.3 and stock SAIGE
// do: every call rebuilds the whole n x n sp_mat from all nnz via gen_sp_Sigma
// and runs a fresh arma::spsolve -- a complete SuperLU symbolic AND numeric
// factorisation -- for one right-hand side. Nothing is kept between calls.
//
// refresh() therefore does nothing but record (w, tau): there is no state to
// refresh. The cost shows up in solve(), 199-261 times a fit, and the report
// splits out how much of it was the rebuild.
//
// hasExplicitInverse() is false. traces() aborts rather than returning a
// number this solver cannot produce.
#include "sigma_solver.hpp"

#include <cstdio>
#include <cstdlib>

namespace step1bench {

namespace {

class SuperLUSolver : public SigmaSolver {
public:
    bool build(const arma::umat& loc, const arma::vec& val, int n) override {
        loc_ = loc; val_ = val; n_ = n;
        return true;
    }

    // No state, no cache -- exactly gen_spsolve_v4 with fit.cache_sparse_solve
    // and fit.block_sparse_sigma off.
    void refresh(const arma::fvec& w, const arma::fvec& tau) override {
        w_ = w; tau_ = tau;
    }

    arma::fvec solve(const arma::fvec& b) override {
        const double t0 = now_s();
        arma::sp_mat A = gen_sp_Sigma(loc_, val_, n_, w_, tau_);
        tBuild_ += now_s() - t0;

        arma::vec bd = arma::conv_to<arma::vec>::from(b);
        const double t1 = now_s();
        arma::vec x = arma::spsolve(A, bd);
        tFactSolve_ += now_s() - t1;
        nSolve_++;
        return arma::conv_to<arma::fvec>::from(x);
    }

    bool hasExplicitInverse() const override { return false; }

    void traces(double*, double*) override {
        fprintf(stderr,
                "step1bench: --trace exact is not available with --solver superlu.\n"
                "  arma::spsolve returns Sigma^-1 b, never Sigma^-1 itself, so there is\n"
                "  no entry of the inverse to contract against Psi. Use --trace hutchinson\n"
                "  (that combination IS the SAIGE baseline), or a solver with an explicit\n"
                "  inverse (block, cholmod, cholmod_reuse).\n");
        std::abort();
    }

    const char* name() const override { return "superlu"; }

    arma::sp_mat debugSigma(const arma::fvec& w, const arma::fvec& tau) const override {
        return gen_sp_Sigma(loc_, val_, n_, w, tau);
    }

    std::string extraReport() const override {
        char buf[256];
        snprintf(buf, sizeof buf,
                 "inside solve(): rebuild %.4fs, SuperLU factor+solve %.4fs over "
                 "%lld solves (no (w,tau) cache, by design)",
                 tBuild_, tFactSolve_, nSolve_);
        return std::string(buf);
    }

private:
    arma::umat loc_;
    arma::vec  val_;
    int n_ = 0;
    arma::fvec w_, tau_;
    double tBuild_ = 0.0, tFactSolve_ = 0.0;
    long long nSolve_ = 0;
};

}  // namespace

std::unique_ptr<SigmaSolver> make_superlu_solver() {
    return std::unique_ptr<SigmaSolver>(new SuperLUSolver());
}

}  // namespace step1bench
