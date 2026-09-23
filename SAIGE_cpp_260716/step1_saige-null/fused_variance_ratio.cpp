// fused_variance_ratio.cpp  -- see fused_variance_ratio.hpp
#include "fused_variance_ratio.hpp"

#include <chrono>
#include <numeric>
#include <unordered_map>
#include <vector>

// Sparse GRM currently installed in SAIGE_step1_fast.cpp (subset to this
// group's samples by main.cpp before setupSparseGRM()).  Stored as a full
// symmetric COO: both triangles present, diagonal explicit.
arma::umat export_sparse_grm_locations();
arma::vec  export_sparse_grm_values();
int        export_sparse_grm_dim();
int        export_sparse_grm_nnz();

namespace saige {

namespace {

struct DSU {
    std::vector<int> p;
    explicit DSU(int n) : p(n) { std::iota(p.begin(), p.end(), 0); }
    int find(int a) { while (p[a] != a) { p[a] = p[p[a]]; a = p[a]; } return a; }
    void join(int a, int b) { a = find(a); b = find(b); if (a != b) p[a] = b; }
};

} // namespace

FusedVrAnchor compute_fused_vr_anchor(const arma::vec& W,
                                      const arma::vec& tau,
                                      const arma::mat& X,
                                      int max_block)
{
    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    FusedVrAnchor R;

    const int n = export_sparse_grm_dim();
    if (n <= 0) { R.why = "no sparse GRM is installed (dimNum == 0)"; return R; }
    if ((int)X.n_rows != n) {
        R.why = "sparse GRM dimension " + std::to_string(n) +
                " does not match the design (" + std::to_string(X.n_rows) + " rows)";
        return R;
    }
    if (W.n_elem != (arma::uword)n || tau.n_elem < 2) {
        R.why = "W / tau have the wrong length";
        return R;
    }
    R.n = n;

    const arma::umat loc = export_sparse_grm_locations();
    const arma::vec  val = export_sparse_grm_values();
    const arma::uword nnz = val.n_elem;

    // ---- connected components of Psi -----------------------------------------
    DSU dsu(n);
    for (arma::uword k = 0; k < nnz; ++k) {
        const int i = (int)loc(0, k), j = (int)loc(1, k);
        if (i != j) dsu.join(i, j);
    }
    std::vector<int> root(n), comp(n, -1);
    int nb = 0;
    for (int i = 0; i < n; ++i) { root[i] = dsu.find(i); }
    for (int i = 0; i < n; ++i) if (root[i] == i) comp[i] = nb++;
    std::vector<int> cid(n);
    for (int i = 0; i < n; ++i) cid[i] = comp[root[i]];
    std::vector<int> csize(nb, 0);
    for (int i = 0; i < n; ++i) csize[cid[i]]++;
    R.nblocks = nb;
    R.maxblock = csize.empty() ? 0 : *std::max_element(csize.begin(), csize.end());
    if (R.maxblock > max_block) {
        R.why = "the sparse GRM has a connected component of " +
                std::to_string(R.maxblock) + " samples (cap " +
                std::to_string(max_block) + "); the exact block inverse would be "
                "O(block^3) and the closed form is not worth it";
        return R;
    }

    // members of each block, and each sample's position inside its block
    std::vector<int> bstart(nb + 1, 0);
    for (int b = 0; b < nb; ++b) bstart[b + 1] = bstart[b] + csize[b];
    std::vector<int> members(n), fill(bstart.begin(), bstart.end() - 1), slot(n);
    for (int i = 0; i < n; ++i) { slot[i] = fill[cid[i]] - bstart[cid[i]]; members[fill[cid[i]]++] = i; }

    // Psi entries grouped by block (indices are block-local)
    std::vector<std::vector<std::tuple<int,int,double>>> bent(nb);
    for (arma::uword k = 0; k < nnz; ++k) {
        const int i = (int)loc(0, k), j = (int)loc(1, k);
        bent[cid[i]].emplace_back(slot[i], slot[j], val(k));
    }

    // ---- per-block dense inverse: tr(Sigma^-1), tr(Sigma^-1 Psi), Sigma^-1 X ---
    const double t0tau = tau(0), t1tau = tau(1);
    const int p = (int)X.n_cols;
    arma::mat SiX(n, p, arma::fill::zeros);
    double trSi = 0.0, trSiPsi = 0.0;
    int floor_hits = 0;

    arma::mat S, Si, Xb;
    for (int b = 0; b < nb; ++b) {
        const int k = csize[b];
        S.zeros(k, k);
        for (const auto& e : bent[b]) S(std::get<0>(e), std::get<1>(e)) += t1tau * std::get<2>(e);
        for (int a = 0; a < k; ++a) {
            const int gi = members[bstart[b] + a];
            S(a, a) += t0tau / W(gi);
            // gen_sp_Sigma() floors the assembled diagonal at 1e-4; mirror it so
            // the anchor describes the same Sigma the PCG/spsolve path uses.
            if (S(a, a) < 1e-4) { S(a, a) = 1e-4; ++floor_hits; }
        }
        if (k == 1) { Si.set_size(1, 1); Si(0, 0) = 1.0 / S(0, 0); }
        else if (!arma::inv_sympd(Si, arma::symmatu(S))) {
            if (!arma::inv(Si, S)) {
                R.why = "a Sigma block of size " + std::to_string(k) + " is singular";
                return R;
            }
        }
        trSi += arma::trace(Si);
        for (const auto& e : bent[b]) trSiPsi += Si(std::get<0>(e), std::get<1>(e)) * std::get<2>(e);

        Xb.set_size(k, p);
        for (int a = 0; a < k; ++a) Xb.row(a) = X.row(members[bstart[b] + a]);
        const arma::mat SiXb = Si * Xb;
        for (int a = 0; a < k; ++a) SiX.row(members[bstart[b] + a]) = SiXb.row(a);
    }
    R.floor_hits = floor_hits;

    // ---- covariate corrections ------------------------------------------------
    const arma::mat XtSiX = X.t() * SiX;
    arma::mat XtSiX_inv;
    if (!arma::inv_sympd(XtSiX_inv, arma::symmatu(XtSiX))) {
        R.why = "X' Sigma^-1 X is singular";
        return R;
    }
    // Psi * SiX
    arma::mat PsiSiX(n, p, arma::fill::zeros);
    for (arma::uword k2 = 0; k2 < nnz; ++k2) {
        const arma::uword i = loc(0, k2), j = loc(1, k2);
        PsiSiX.row(i) += val(k2) * SiX.row(j);
    }

    const double trPsi = [&]{ double s = 0; for (arma::uword k2 = 0; k2 < nnz; ++k2)
                                 if (loc(0, k2) == loc(1, k2)) s += val(k2); return s; }();
    // tr(H Psi) = tr((X'X)^-1 X' Psi X)
    arma::mat PsiX(n, p, arma::fill::zeros);
    for (arma::uword k2 = 0; k2 < nnz; ++k2) PsiX.row(loc(0, k2)) += val(k2) * X.row(loc(1, k2));
    arma::mat XtX_inv;
    if (!arma::inv_sympd(XtX_inv, arma::symmatu(X.t() * X))) {
        R.why = "X'X is singular";
        return R;
    }
    const double trHPsi = arma::accu(XtX_inv % (X.t() * PsiX).t());

    // tr((I - 11'/n) Psi) = tr(Psi) - (1' Psi 1)/n, the denominator the
    // "null_noXadj" row uses (that row mean-centres the genotype instead of
    // projecting the covariates out).
    double sumPsi = 0.0;
    for (arma::uword k2 = 0; k2 < nnz; ++k2) sumPsi += val(k2);

    // ---- binary denominator (see the header for why only one term survives)
    // W == 1 reduces every line below to its quantitative counterpart, which is
    // the self-check: trPsi_Wproj must then equal trPsi_proj.
    arma::mat WX(n, p);
    for (int i = 0; i < n; ++i) WX.row(i) = W(i) * X.row(i);
    arma::mat PsiWX(n, p, arma::fill::zeros);
    for (arma::uword k2 = 0; k2 < nnz; ++k2)
        PsiWX.row(loc(0, k2)) += val(k2) * WX.row(loc(1, k2));
    double trWPsi = 0.0;
    for (arma::uword k2 = 0; k2 < nnz; ++k2)
        if (loc(0, k2) == loc(1, k2)) trWPsi += W(loc(0, k2)) * val(k2);
    arma::mat XtWX_inv;
    double trPsi_Wproj = std::numeric_limits<double>::quiet_NaN();
    if (arma::inv_sympd(XtWX_inv, arma::symmatu(X.t() * WX)))
        trPsi_Wproj = trWPsi - arma::accu(XtWX_inv % (WX.t() * PsiWX).t());

    R.trPsi         = trPsi;
    R.trPsi_mean    = trPsi - sumPsi / n;
    R.trPsi_proj    = trPsi - trHPsi;
    R.trSigmaInvPsi = trSiPsi;
    R.trPPsi        = trSiPsi - arma::accu(XtSiX_inv % (SiX.t() * PsiSiX).t());
    R.trSigmaInv    = trSi;
    R.trP           = trSi - arma::accu(XtSiX_inv % (SiX.t() * SiX).t());
    R.trWPsi        = trWPsi;
    R.trPsi_Wproj   = trPsi_Wproj;
    R.anchor_binary = R.trPPsi / trPsi_Wproj;
    R.anchor        = R.trPPsi / R.trPsi_proj;
    R.anchor_noXadj = R.trPPsi / R.trPsi_mean;
    R.anchor_raw    = R.trPPsi / R.trPsi;
    R.trP_over_N    = R.trP / n;
    R.ok            = std::isfinite(R.anchor) && R.anchor > 0.0;
    if (!R.ok) R.why = "the anchor came out non-finite or non-positive";
    R.seconds = std::chrono::duration<double>(clk::now() - t0).count();
    return R;
}

} // namespace saige
