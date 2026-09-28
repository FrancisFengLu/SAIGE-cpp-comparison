#include "sigma_solver.hpp"

#include <chrono>
#include <functional>
#include <map>

namespace step1bench {

// ---------------------------------------------------------------------------
// Verbatim arithmetic of SAIGE's gen_sp_Sigma. Every deviation here would void
// the benchmark, so it is deliberately written to look like the original
// instead of like better C++.
arma::sp_mat gen_sp_Sigma(const arma::umat& loc, const arma::vec& val, int n,
                          const arma::fvec& w, const arma::fvec& tau) {
    arma::fvec dtVec = (1 / w) * (tau(0));          // fp32, as upstream
    arma::vec  valueVecNew = val * tau(1);          // fp64 * (double)tau1

    const arma::uword nnonzero = val.n_elem;
    for (arma::uword i = 0; i < nnonzero; i++) {
        if (loc(0, i) == loc(1, i)) {
            valueVecNew(i) = valueVecNew(i) + dtVec(loc(0, i));
            if (valueVecNew(i) < 1e-4) valueVecNew(i) = 1e-4;
        }
    }
    return arma::sp_mat(loc, valueVecNew, (arma::uword)n, (arma::uword)n);
}

// ---------------------------------------------------------------------------
double now_s() {
    using clk = std::chrono::steady_clock;
    return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}

static Stage g_stage = ST_SOLVE;
void  set_stage(Stage s) { g_stage = s; }
Stage cur_stage()        { return g_stage; }

// ---------------------------------------------------------------------------
// Union-find over the off-diagonal entries. Same rule as BlockSigma::build, so
// the histogram reported next to a superlu or cholmod run describes the same
// partition the block path would have used.
BlockSummary summarize_blocks(const arma::umat& loc, int n) {
    std::vector<int> p((size_t)n);
    for (int i = 0; i < n; i++) p[(size_t)i] = i;
    std::function<int(int)> find = [&](int x) {
        while (p[(size_t)x] != x) { p[(size_t)x] = p[(size_t)p[(size_t)x]]; x = p[(size_t)x]; }
        return x;
    };
    for (arma::uword k = 0; k < loc.n_cols; k++) {
        const int a = (int)loc(0, k), b = (int)loc(1, k);
        if (a == b) continue;
        const int ra = find(a), rb = find(b);
        if (ra != rb) p[(size_t)ra] = rb;
    }
    std::map<int, int> sz;
    for (int i = 0; i < n; i++) sz[find(i)]++;

    BlockSummary out;
    std::map<int, long long> h;
    for (auto& kv : sz) { h[kv.second]++; out.maxBlock = std::max(out.maxBlock, kv.second); }
    out.nblocks = (int)sz.size();
    for (auto& kv : h) out.hist.emplace_back(kv.first, kv.second);
    return out;
}

// Defined in the three solver translation units.
std::unique_ptr<SigmaSolver> make_block_solver();
std::unique_ptr<SigmaSolver> make_cholmod_solver(bool reuse_symbolic);
std::unique_ptr<SigmaSolver> make_superlu_solver();

std::unique_ptr<SigmaSolver> make_solver(const std::string& name) {
    if (name == "block")         return make_block_solver();
    if (name == "cholmod")       return make_cholmod_solver(false);
    if (name == "cholmod_reuse") return make_cholmod_solver(true);
    if (name == "superlu")       return make_superlu_solver();
    return nullptr;
}

}  // namespace step1bench
