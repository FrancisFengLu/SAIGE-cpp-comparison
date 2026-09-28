// solver_block.cpp -- the block-diagonal explicit inverse.
//
// Nothing here implements an algorithm; it is a thin adapter onto
// blocksigma::BlockSigma, copied verbatim from ../step1_saige-null/ so that
// what this benchmark measures is the production class and not a
// reimplementation of it.
#include "sigma_solver.hpp"
#include "block_sigma.hpp"

#include <cstdio>
#include <cstdlib>

namespace step1bench {

namespace {

class BlockSolver : public SigmaSolver {
public:
    bool build(const arma::umat& loc, const arma::vec& val, int n) override {
        // The wall-clock gate exists upstream to refuse GRMs whose blocks are
        // too large to invert densely. In the harness it is left at the
        // production default: if it refuses, that is a real answer about this
        // GRM and must not be papered over.
        return bs_.build(loc, val, n);
    }

    void refresh(const arma::fvec& w, const arma::fvec& tau) override {
        bs_.refresh(w, tau);            // no-op when (w, tau) repeat
    }

    arma::fvec solve(const arma::fvec& b) override { return bs_.solve(b); }

    bool hasExplicitInverse() const override { return true; }

    void traces(double* trSigmaInvPsi, double* trSigmaInv) override {
        bs_.traces(trSigmaInvPsi, trSigmaInv);
    }

    const char* name() const override { return "block"; }

    long long refreshCount() const override { return bs_.refreshCount(); }
    long long reuseCount()   const override { return bs_.reuseCount(); }

    std::string extraReport() const override {
        char buf[256];
        snprintf(buf, sizeof buf,
                 "partition %d blocks, max %d; floored diagonals last refresh %lld",
                 bs_.partition().nblocks, bs_.partition().maxBlock,
                 bs_.flooredDiagonals());
        return std::string(buf);
    }

    // Reconstruct the Sigma that refresh() forms, from the partition's own
    // stored copy of Psi. This deliberately repeats the arithmetic of
    // block_sigma.cpp:428-462 rather than calling it, because the thing being
    // checked IS that arithmetic plus the partition data it reads: if build()
    // dropped, duplicated or misplaced an entry, this picks it up. The startup
    // self-check also does a numeric residual check on the real inverse, which
    // covers the other half (that refresh() actually produced what this says).
    arma::sp_mat debugSigma(const arma::fvec& w, const arma::fvec& tau) const override {
        const blocksigma::Partition& P = bs_.partition();
        if (!P.built) return arma::sp_mat((arma::uword)0, (arma::uword)0);
        const float tau0 = tau(0), tau1f = tau(1);
        const double tau1 = (double)tau1f;

        const int nnz = P.psiStart.empty() ? 0 : P.psiStart.back();
        arma::umat loc(2, (arma::uword)nnz);
        arma::vec  v((arma::uword)nnz);
        arma::uword at = 0;
        for (int b = 0; b < P.nblocks; b++) {
            const int mo = P.start[(size_t)b];
            for (int k = P.psiStart[(size_t)b]; k < P.psiStart[(size_t)b + 1]; k++) {
                const int r = P.psiR[(size_t)k], c = P.psiC[(size_t)k];
                const int gi = P.member[(size_t)(mo + r)];
                const int gj = P.member[(size_t)(mo + c)];
                double val = P.psiV[(size_t)k] * tau1;
                if (r == c) {
                    const float dt = (1.0f / w((arma::uword)gi)) * tau0;
                    val += (double)dt;
                    if (val < 1e-4) val = 1e-4;
                }
                loc(0, at) = (arma::uword)gi;
                loc(1, at) = (arma::uword)gj;
                v(at) = val;
                at++;
            }
        }
        return arma::sp_mat(loc, v, (arma::uword)P.n, (arma::uword)P.n);
    }

    // Used by the self-check's residual test and nowhere else.
    blocksigma::BlockSigma& raw() { return bs_; }

private:
    blocksigma::BlockSigma bs_;
};

}  // namespace

std::unique_ptr<SigmaSolver> make_block_solver() {
    return std::unique_ptr<SigmaSolver>(new BlockSolver());
}

}  // namespace step1bench
