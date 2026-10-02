#ifndef SPA_BINARY_HPP
#define SPA_BINARY_HPP

#include <armadillo>

// Struct replacements for Rcpp::List returns

struct RootResult {
    double root;
    int niter;
    bool Isconverge;
};

struct SaddleResult {
    double pval;
    bool isSaddle;
};

struct SPAResult {
    double pvalue;
    bool Isconverge;
};

// Config key spaScratch (default false; parsed in main.cpp). With it set the
// SPA block does exactly the same arithmetic in the same order but allocates
// nothing per call:
//   * gpos / gneg -- the sums of the positive and of the negative entries of
//     g~ that bound the root -- are computed ONCE per SPA call by spaGposGneg
//     (below) instead of once per root-find by accu(g.elem(find(g > 0))),
//     which built a ~25k-entry index vector and a compacted copy, twice per
//     root, four times per pair;
//   * the four carrier / non-carrier subset vectors of getMarkerPval and the
//     N-vector temporary of getadjGFast live in thread_local grow-only
//     buffers (saige_test.cpp).
// Default off keeps every byte of every output as it was; on, the outputs are
// byte-identical too (tests/run_phase0_gate.sh), the switch exists so that can
// be checked rather than assumed. S2_REMAINDER.md section 6.2 has the cost.
extern bool g_spaScratch;

// gpos = accu(g.elem(find(g > 0))), gneg = accu(g.elem(find(g < 0))), bit for
// bit. Armadillo's accu over the compacted vector runs two interleaved
// accumulators (even positions into one, odd into the other) and adds them at
// the end (op_accu_meat.hpp, apply_proxy_linear), so this walks g once in
// index order and does the same. tools/spa_gpos_test.cpp checks it.
void spaGposGneg(const arma::vec& g, double& gpos, double& gneg);

// --- Standard (non-fast) functions ---

double Korg_Binom(double t1, arma::vec & mu, arma::vec & g);
double K1_adj_Binom(double t1, arma::vec & mu, arma::vec & g, double q);
double K2_Binom(double t1, arma::vec & mu, arma::vec & g);

RootResult getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter = 1000);
// The same root-find with gpos / gneg supplied by the caller (spaScratch).
RootResult getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter, double gpos, double gneg);
SaddleResult Get_Saddle_Prob_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, bool logp = false);
SPAResult SPA_binary(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp = false);

// --- Fast functions (normal approximation for zero-genotype subset) ---

double Korg_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Binom(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K2_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);

RootResult getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol, int maxiter = 1000);
RootResult getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol, int maxiter, double gpos, double gneg);
SaddleResult Get_Saddle_Prob_fast_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, bool logp = false);
SPAResult SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol);

#endif
