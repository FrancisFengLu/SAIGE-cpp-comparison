// Standalone port of SKAT/Burden/SKAT-O tests for SAIGE Step 2
// Replaces: SKAT:::Met_SKAT_Get_Pvalue (R) with pure C++ implementation
//
// Key algorithms:
//   1. Davies method (qfc) for P(sum(lambda_i * chi2_1) > q)
//   2. Liu moment-matching method (fallback when Davies fails)
//   3. SKAT-O optimal.adj method (grid over rho, integration correction)
//
// Dependencies:
//   - armadillo (eigendecomposition)
//   - boost::math (chi_squared, normal distributions)
//   - No R, no Rcpp

#ifndef SKAT_HPP
#define SKAT_HPP

#include <armadillo>
#include <string>

// ============================================================
// Result struct matching get_SKAT_pvalue() R function output
// ============================================================
struct SKATResult {
    double pvalue_SKATO;    // SKAT-O p-value (optimal rho)
    double pvalue_Burden;   // Burden test p-value (rho=1)
    double pvalue_SKAT;     // SKAT test p-value (rho=0)
    double beta_Burden;     // Burden effect size estimate
    double se_Burden;       // Burden standard error
};

// ============================================================
// Main entry point — matches get_SKAT_pvalue from R
// ============================================================
// Score: weighted score vector (m x 1), already multiplied by weights
// Phi: weighted variance-covariance matrix (m x m), already w*Phi*w
// r_corr: correlation parameter grid (usually {0, 0.01, 0.04, ..., 0.81, 1})
// regionTestType: "SKAT", "BURDEN", or "SKATO"
SKATResult get_SKAT_pvalue(const arma::vec& Score,
                            const arma::mat& Phi,
                            const arma::vec& r_corr,
                            const std::string& regionTestType);

// ============================================================
// Davies method: P(Q > q) where Q = sum(lambda_i * chi2_1)
// ============================================================
// q: observed test statistic
// lambda: eigenvalues (mixture weights)
// Returns: tail probability P(Q > q), or NaN on failure
// Uses the qfc algorithm (Davies 1980, AS 155)
double davies_pvalue(double q, const arma::vec& lambda);

// ============================================================
// Liu moment-matching method (fallback when Davies fails)
// ============================================================
// Approximates the distribution of sum(lambda_i * chi2_1)
// by a scaled noncentral chi-squared distribution
// matching the first 4 moments
double liu_pvalue(double q, const arma::vec& lambda);

// ============================================================
// Internal: SKAT-O optimal.adj computation
// ============================================================
// Implements SKAT:::Met_SKAT_Get_Pvalue with method="optimal.adj"
// Returns the SKAT-O p-value given per-rho p-values and test statistics
struct SKATOParam {
    arma::vec rho;          // rho grid values
    arma::vec p_val_each;   // p-value at each rho
    arma::vec q_val_each;   // Q statistic at each rho
};

// Compute SKAT-O p-value from per-rho p-values
// Uses integration approach from Lee et al. (2012) SKAT-O paper
double SKATO_optimal_pvalue(const arma::vec& Score,
                             const arma::mat& Phi,
                             const arma::vec& rho_vec,
                             const arma::vec& pval_each,
                             const arma::vec& q_each);

#endif
