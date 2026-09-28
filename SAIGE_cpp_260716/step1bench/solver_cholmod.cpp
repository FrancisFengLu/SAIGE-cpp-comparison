// solver_cholmod.cpp -- TorchGWAS2's shape: a simplicial LL' with natural
// ordering, inverted explicitly by sparse-solving against an identity.
//
// Ported from, and deliberately kept in the same order as:
//   /opt/saige/logs/priorart_mt/TorchGWAS2/src/SparseInverse.cpp:388  inv_spamat
//   /opt/saige/logs/priorart_mt/TorchGWAS2/src/GMMAT.cpp:1112-1135    per-iteration call
//   /opt/saige/logs/priorart_mt/TorchGWAS2/src/GMMAT.cpp:1161-1164    traces
//
// Their settings, reproduced exactly:
//   cm.supernodal          = CHOLMOD_SIMPLICIAL
//   cm.method[0].ordering  = CHOLMOD_NATURAL
//   cm.postorder           = 0
//   cm.nmethods            = 1
//   cm.final_ll            = 1
//   cm.final_pack          = 1
//   cm.final_asis          = 0, cm.final_monotonic = 1   (set after analyze)
//
// Their inverse is not cholmod_spsolve against a dense identity; it is
// cholmod_factor_to_sparse -> CXSparse -> one cs_spsolve per column against a
// sparse identity, giving the columns of L^-1, then Sigma^-1 = (L^-1)' (L^-1).
// That is the cost being measured, so it is what is implemented.
//
// Two variants share this file:
//   "cholmod"       -- cholmod_analyze inside every refresh, as they do.
//   "cholmod_reuse" -- cholmod_analyze hoisted to build(). Same numbers, and
//                      the difference between the two IS the answer to "how
//                      much of their cost is the redundant symbolic step".
#include "sigma_solver.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

extern "C" {
#include <suitesparse/cholmod.h>
#include <suitesparse/cs.h>
}

namespace step1bench {

namespace {

// cholmod_sparse from an arma::sp_mat, stype = 1 (upper triangle used, lower
// ignored) -- the same call TorchGWAS2 makes on their Eigen CSC matrix.
cholmod_sparse* to_cholmod(const arma::sp_mat& A, cholmod_common* cm) {
    const arma::uword n   = A.n_rows;
    const arma::uword nnz = A.n_nonzero;
    cholmod_sparse* ssm = cholmod_allocate_sparse(
        (size_t)n, (size_t)A.n_cols, (size_t)nnz,
        /*sorted*/ 1, /*packed*/ 1, /*stype*/ 1, CHOLMOD_REAL, cm);
    if (!ssm) throw std::runtime_error("cholmod_allocate_sparse failed");

    int*    p = (int*)ssm->p;
    int*    i = (int*)ssm->i;
    double* x = (double*)ssm->x;
    const arma::uword* cp = A.col_ptrs;
    const arma::uword* ri = A.row_indices;
    for (arma::uword c = 0; c <= A.n_cols; c++) p[c] = (int)cp[c];
    for (arma::uword k = 0; k < nnz; k++) { i[k] = (int)ri[k]; x[k] = A.values[k]; }
    return ssm;
}

// cholmodSparseToCsparse, verbatim in effect (SparseInverse.cpp:325).
cs_di* to_csparse(cholmod_sparse* L) {
    const int n = (int)L->nrow;
    const int nzmax = (int)L->nzmax;
    cs_di* M = cs_di_spalloc(n, n, nzmax, 1, 0);
    if (!M) throw std::runtime_error("cs_di_spalloc failed");
    std::memcpy(M->p, L->p, (size_t)(n + 1) * sizeof(int));
    std::memcpy(M->i, L->i, (size_t)nzmax * sizeof(int));
    std::memcpy(M->x, L->x, (size_t)nzmax * sizeof(double));
    M->nz = -1;                      // compressed column form
    return M;
}

cs_di* sparse_identity(int m) {
    cs_di* I = cs_di_spalloc(m, m, m, 1, 0);
    if (!I) throw std::runtime_error("cs_di_spalloc(identity) failed");
    for (int i = 0; i < m; ++i) { I->p[i] = i; I->i[i] = i; I->x[i] = 1.0; }
    I->p[m] = m;
    I->nz = -1;
    return I;
}

class CholmodSolver : public SigmaSolver {
public:
    explicit CholmodSolver(bool reuse) : reuse_(reuse) {}

    ~CholmodSolver() override {
        if (started_) {
            if (factor_) cholmod_free_factor(&factor_, &cm_);
            cholmod_finish(&cm_);
        }
    }

    bool build(const arma::umat& loc, const arma::vec& val, int n) override {
        loc_ = loc; val_ = val; n_ = n;
        psi_ = arma::sp_mat(loc, val, (arma::uword)n, (arma::uword)n);
        if (reuse_) {
            start_common();
            // Symbolic factorisation on the pattern Sigma will have. The
            // pattern is Psi's, for every (w, tau) with tau1 != 0; refresh()
            // re-analyses (and says so) if it ever differs.
            arma::fvec w1((arma::uword)n, arma::fill::ones);
            arma::fvec t0(2); t0(0) = 1.0f; t0(1) = 1.0f;
            arma::sp_mat S = gen_sp_Sigma(loc_, val_, n_, w1, t0);
            analyze_on(S);
        }
        return true;
    }

    void refresh(const arma::fvec& w, const arma::fvec& tau) override {
        if (upToDate(w, tau)) { nReuse_++; return; }
        nRefresh_++;

        arma::sp_mat S = gen_sp_Sigma(loc_, val_, n_, w, tau);

        if (!reuse_) start_common();          // cholmod_start, as they do
        cholmod_sparse* spMatrix = to_cholmod(S, &cm_);

        if (reuse_) {
            if (!patternMatches(S)) {         // never expected; counted, not hidden
                nReanalyze_++;
                if (factor_) cholmod_free_factor(&factor_, &cm_);
                analyze_on(S);
            }
        } else {
            const double ta = now_s();
            factor_ = cholmod_analyze(spMatrix, &cm_);
            tAnalyze_ += now_s() - ta;
            if (!factor_) throw std::runtime_error("cholmod_analyze failed");
        }

        cm_.final_asis = 0;
        cm_.final_monotonic = 1;
        const double tf = now_s();
        const int status = cholmod_factorize(spMatrix, factor_, &cm_);
        tFactorize_ += now_s() - tf;
        if (!status) throw std::runtime_error("cholmod_factorize failed");

        // --- explicit inverse, their way -------------------------------
        const double ti = now_s();
        cholmod_sparse* L = cholmod_factor_to_sparse(factor_, &cm_);
        cs_di* csFactor = to_csparse(L);
        cholmod_free_sparse(&L, &cm_);

        const int nc = n_;
        std::vector<int>    xi((size_t)(2 * nc));
        std::vector<double> xw((size_t)nc);
        cs_di* csb = sparse_identity(nc);

        std::vector<arma::uword> rr, cc;
        std::vector<double> vv;
        rr.reserve((size_t)nc * 4); cc.reserve((size_t)nc * 4); vv.reserve((size_t)nc * 4);
        for (int col = 0; col < nc; ++col) {
            const int top = cs_di_spsolve(csFactor, csb, col, xi.data(), xw.data(),
                                          nullptr, 1);
            if (top == -1) throw std::runtime_error("cs_di_spsolve failed");
            for (int i = top; i < nc; ++i) {
                rr.push_back((arma::uword)xi[(size_t)i]);
                cc.push_back((arma::uword)col);
                vv.push_back(xw[(size_t)xi[(size_t)i]]);
            }
        }
        cs_di_spfree(csb);
        cs_di_spfree(csFactor);

        arma::umat locs(2, rr.size());
        for (size_t k = 0; k < rr.size(); k++) { locs(0, k) = rr[k]; locs(1, k) = cc[k]; }
        arma::vec vals(vv.data(), vv.size());
        arma::sp_mat Linv(locs, vals, (arma::uword)nc, (arma::uword)nc);
        sinv_ = Linv.t() * Linv;                 // Sigma^-1 = L^-T L^-1
        tInverse_ += now_s() - ti;

        cholmod_free_sparse(&spMatrix, &cm_);
        if (!reuse_) {
            cholmod_free_factor(&factor_, &cm_);
            cholmod_finish(&cm_);                // ... and cholmod_finish
            started_ = false;
        }
        lastW_ = w; lastTau_ = tau; valid_ = true;
    }

    arma::fvec solve(const arma::fvec& b) override {
        arma::vec bd = arma::conv_to<arma::vec>::from(b);
        arma::vec x  = sinv_ * bd;
        return arma::conv_to<arma::fvec>::from(x);
    }

    bool hasExplicitInverse() const override { return true; }

    // GMMAT.cpp:1161-1164: sigma_i.cwiseProduct(kin).sum().
    void traces(double* trSigmaInvPsi, double* trSigmaInv) override {
        if (trSigmaInvPsi) *trSigmaInvPsi = arma::accu(sinv_ % psi_);
        if (trSigmaInv)    *trSigmaInv    = arma::trace(sinv_);
    }

    const char* name() const override { return reuse_ ? "cholmod_reuse" : "cholmod"; }

    arma::sp_mat debugSigma(const arma::fvec& w, const arma::fvec& tau) const override {
        return gen_sp_Sigma(loc_, val_, n_, w, tau);
    }

    long long refreshCount() const override { return nRefresh_; }
    long long reuseCount()   const override { return nReuse_; }

    std::string extraReport() const override {
        char buf[320];
        snprintf(buf, sizeof buf,
                 "cholmod_analyze %.4fs (%s), cholmod_factorize %.4fs, "
                 "explicit inverse %.4fs, Sigma^-1 nnz %llu%s",
                 tAnalyze_, reuse_ ? "hoisted to build()" : "inside every refresh",
                 tFactorize_, tInverse_,
                 (unsigned long long)sinv_.n_nonzero,
                 nReanalyze_ ? ", pattern changed mid-run (re-analysed)" : "");
        return std::string(buf);
    }

private:
    void start_common() {
        cholmod_start(&cm_);
        started_ = true;
        cm_.supernodal         = CHOLMOD_SIMPLICIAL;
        cm_.method[0].ordering = CHOLMOD_NATURAL;
        cm_.postorder          = 0;
        cm_.nmethods           = 1;
        cm_.final_ll           = 1;
        cm_.final_pack         = 1;
    }

    void analyze_on(const arma::sp_mat& S) {
        cholmod_sparse* sp = to_cholmod(S, &cm_);
        const double ta = now_s();
        factor_ = cholmod_analyze(sp, &cm_);
        tAnalyze_ += now_s() - ta;
        patN_ = S.n_nonzero;
        patCol_.assign(S.col_ptrs, S.col_ptrs + S.n_cols + 1);
        patRow_.assign(S.row_indices, S.row_indices + S.n_nonzero);
        cholmod_free_sparse(&sp, &cm_);
        if (!factor_) throw std::runtime_error("cholmod_analyze failed");
    }

    bool patternMatches(const arma::sp_mat& S) const {
        if (S.n_nonzero != patN_) return false;
        if (std::memcmp(S.col_ptrs, patCol_.data(),
                        (size_t)(S.n_cols + 1) * sizeof(arma::uword)) != 0) return false;
        return std::memcmp(S.row_indices, patRow_.data(),
                           (size_t)S.n_nonzero * sizeof(arma::uword)) == 0;
    }

    bool upToDate(const arma::fvec& w, const arma::fvec& tau) const {
        if (!valid_ || lastW_.n_elem != w.n_elem || lastTau_.n_elem != tau.n_elem)
            return false;
        if (std::memcmp(lastTau_.memptr(), tau.memptr(), tau.n_elem * sizeof(float)) != 0)
            return false;
        return std::memcmp(lastW_.memptr(), w.memptr(), w.n_elem * sizeof(float)) == 0;
    }

    bool reuse_;
    bool started_ = false;
    cholmod_common cm_{};
    cholmod_factor* factor_ = nullptr;

    arma::umat loc_;
    arma::vec  val_;
    int n_ = 0;
    arma::sp_mat psi_;
    arma::sp_mat sinv_;

    arma::uword patN_ = 0;
    std::vector<arma::uword> patCol_, patRow_;

    arma::fvec lastW_, lastTau_;
    bool valid_ = false;
    long long nRefresh_ = 0, nReuse_ = 0, nReanalyze_ = 0;
    double tAnalyze_ = 0.0, tFactorize_ = 0.0, tInverse_ = 0.0;
};

}  // namespace

std::unique_ptr<SigmaSolver> make_cholmod_solver(bool reuse_symbolic) {
    return std::unique_ptr<SigmaSolver>(new CholmodSolver(reuse_symbolic));
}

}  // namespace step1bench
