// s2_gpu_sparse.cpp -- see s2_gpu_sparse.hpp.
#include "s2_gpu_sparse.hpp"

#include <chrono>
#include <cstdint>
#include <unordered_map>

#include "block_sigma.hpp"

namespace s2gs {

std::string build(const SAIGE::MTContext& ctx,
                  const std::vector<SAIGE::SAIGEClass*>& objs,
                  double refreshBudget_s, long long maxPairs, Plan& out)
{
    const auto t0 = std::chrono::steady_clock::now();
    const int P = ctx.P, N = ctx.N, nBin = ctx.nBin;
    out = Plan();
    out.N = N; out.nBin = nBin;
    out.on.assign(std::max(nBin, 1), 0);
    out.YBY.assign(P, arma::mat());
    out.Yres.assign(P, arma::vec());
    if (nBin <= 0) return "no binary traits";
    out.XVt.zeros(N, ctx.sumPbin);
    out.BY.zeros(N, ctx.sumPbin);
    out.Bdiag.zeros(N, nBin);

    // Per-trait pair weights, keyed by (i << 32 | j), merged into one list.
    std::unordered_map<uint64_t, long long> pairIdx;
    std::vector<std::vector<std::pair<long long, double>>> tw(nBin);

    for (int t = 0; t < P; t++) {
        const SAIGE::TraitMeta& M = ctx.meta[t];
        if (M.kind != SAIGE::TraitKind::Binary) continue;
        SAIGE::SAIGEClass* obj = objs[t];
        if (!obj->m_flagSparseGRM) continue;
        const int b = M.binIdx;
        const arma::sp_mat& S = obj->m_spSigmaMat;
        if ((int)S.n_rows != N || (int)S.n_cols != N)
            return "trait '" + M.name + "': sparse Sigma is not N x N";
        if (M.p > 64) return "trait '" + M.name + "' has more than 64 covariates";
        if (obj->m_isVarPsadj)
            return "trait '" + M.name + "': the model carries Sigma_iXXSigma_iX (variance P-adjustment), "
                   "which the device sparse variance does not implement";
        if (!ctx.samp.empty() && !ctx.samp[t].sameAsUnion)
            return "trait '" + M.name + "' has its own sample list";

        // Sigma exactly as SAIGEClass holds it.
        arma::umat loc(2, S.n_nonzero);
        arma::vec  val(S.n_nonzero);
        {
            arma::uword k = 0;
            for (auto it = S.begin(); it != S.end(); ++it, ++k) {
                loc(0, k) = it.row(); loc(1, k) = it.col(); val(k) = *it;
            }
        }
        std::vector<char> hasDiag((size_t)N, 0);
        for (arma::uword k = 0; k < loc.n_cols; k++)
            if (loc(0, k) == loc(1, k)) hasDiag[(size_t)loc(0, k)] = 1;
        for (int i = 0; i < N; i++)
            if (!hasDiag[(size_t)i])
                return "trait '" + M.name + "': a sample has no diagonal entry in the sparse Sigma";

        blocksigma::BlockSigma bs;
        bs.setRefreshBudget(refreshBudget_s);
        if (!bs.build(loc, val, N))
            return "trait '" + M.name + "': block inverse refused by the cost gate "
                   "(blockSparseSigmaRefreshBudget_s)";
        // tau = (0, 1), w = 1: refresh() inverts the given matrix itself (s2_block_solve.cpp).
        arma::fvec wv((arma::uword)N, arma::fill::ones);
        arma::fvec tau(2); tau(0) = 0.0f; tau(1) = 1.0f;
        bs.refresh(wv, tau);
        if (bs.flooredDiagonals() > 0 || !bs.ready())
            return "trait '" + M.name + "': block inverse not usable (diagonal floor fired)";
        const blocksigma::Partition& Pt = bs.partition();
        out.nBlocksMax = std::max(out.nBlocksMax, Pt.nblocks);
        out.maxBlock = std::max(out.maxBlock, Pt.maxBlock);

        double* bd = out.Bdiag.colptr((arma::uword)b);
        for (int blk = 0; blk < Pt.nblocks; blk++) {
            const int s = Pt.size[(size_t)blk];
            const int* m = &Pt.member[(size_t)Pt.start[(size_t)blk]];
            const double* A = bs.blockInverse(blk);       // s x s column-major
            for (int r = 0; r < s; r++) bd[m[r]] = A[(size_t)r * s + r];
            for (int r = 0; r < s; r++)
                for (int c = r + 1; c < s; c++) {
                    int i = m[r], j = m[c];
                    if (i > j) std::swap(i, j);
                    const uint64_t key = ((uint64_t)(uint32_t)i << 32) | (uint32_t)j;
                    auto ins = pairIdx.emplace(key, (long long)pairIdx.size());
                    if (ins.second) { out.pi.push_back(i); out.pj.push_back(j); }
                    tw[b].push_back({ins.first->second, A[(size_t)c * s + r] + A[(size_t)r * s + c]});
                }
        }
        if ((long long)pairIdx.size() > maxPairs)
            return "within-block pairs exceed " + std::to_string(maxPairs) + " (largest block " +
                   std::to_string(Pt.maxBlock) + ")";

        // The covariate-side columns, as scoreTest's g~ = g - XXVX_inv (XV g)
        // uses them: XV' as is, and B XXVX_inv block by block in double.
        const arma::mat& Y  = obj->m_XXVX_inv;     // N x p
        const arma::mat& XV = obj->m_XV;           // p x N
        const int p = M.p;
        if ((int)Y.n_rows != N || (int)Y.n_cols != p || (int)XV.n_rows != p || (int)XV.n_cols != N)
            return "trait '" + M.name + "': XV / XXVX_inv have unexpected shapes";
        const arma::uword w0 = (arma::uword)M.binOff;
        for (int c = 0; c < p; c++) {
            out.XVt.col(w0 + c) = XV.row((arma::uword)c).t();
            const double* x = Y.colptr((arma::uword)c);
            double* y = out.BY.colptr(w0 + c);
            for (int blk = 0; blk < Pt.nblocks; blk++) {
                const int s = Pt.size[(size_t)blk];
                const int* m = &Pt.member[(size_t)Pt.start[(size_t)blk]];
                const double* A = bs.blockInverse(blk);
                for (int r = 0; r < s; r++) {
                    double acc = 0.0;
                    for (int q = 0; q < s; q++) acc += A[(size_t)q * s + r] * x[m[q]];
                    y[m[r]] = acc;
                }
            }
        }
        out.YBY[t]  = Y.t() * out.BY.cols(w0, w0 + p - 1);
        out.Yres[t] = Y.t() * ctx.RES.col((arma::uword)t);
        out.on[b] = 1;
        out.nOn++;
    }
    if (out.nOn == 0) return "no binary trait carries a sparse GRM";
    out.nPairs = (long long)pairIdx.size();
    out.w.assign((size_t)nBin * (size_t)out.nPairs, 0.0);
    for (int b = 0; b < nBin; b++)
        for (const auto& e : tw[b]) out.w[(size_t)b * out.nPairs + e.first] += e.second;
    out.secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return "";
}

}  // namespace s2gs
