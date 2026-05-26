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

// --- Standard (non-fast) functions ---

double Korg_Binom(double t1, arma::vec & mu, arma::vec & g);
double K1_adj_Binom(double t1, arma::vec & mu, arma::vec & g, double q);
double K2_Binom(double t1, arma::vec & mu, arma::vec & g);

RootResult getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter = 1000);
SaddleResult Get_Saddle_Prob_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, bool logp = false);
SPAResult SPA_binary(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp = false);

// --- Fast functions (normal approximation for zero-genotype subset) ---

double Korg_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K1_adj_fast_Binom(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);
double K2_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma);

RootResult getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol, int maxiter = 1000);
SaddleResult Get_Saddle_Prob_fast_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, bool logp = false);
SPAResult SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol);

#endif
