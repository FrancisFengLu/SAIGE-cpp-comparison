#pragma once
#include "saige_null.hpp"
#include <string>
#include <vector>

namespace saige {

struct CovariateOffsetResult {
  std::vector<double> beta;       // length p
  std::vector<double> offset;     // length n (X[:,1:] * beta[1:])
  std::vector<double> eta;        // length n (X * beta)
  std::vector<double> mu;         // length n
  std::vector<double> dev_trace;  // deviance per IRLS step (binary only)
  int  iterations{0};
  bool converged{true};
};

// Fit initial GLM on the full design (X, y) — the covariate_offset path.
//   trait = "binary"       : logistic IRLS with deviance check + step halving
//   trait = "quantitative" : ridge-regularized normal equations (closed form)
// If `prior_offset` is non-null, the linear predictor is η = X·β + prior_offset
//   (matches R's glm.fit(family=binomial, offset=...)). For binary IRLS the
//   working response becomes z = η + (y-μ)/W − prior_offset.
// If dump_prefix is non-empty, writes
//   <dump_prefix>_cov_input.csv  (iid, y, X0..Xp-1)
//   <dump_prefix>_cov_beta.csv   (col, beta)
//   <dump_prefix>_cov_fit.csv    (sample, eta, mu, offset)
//   <dump_prefix>_cov_trace.csv  (iter, deviance)  [binary only]
CovariateOffsetResult fit_covariate_offset(const Design& design,
                                           const std::string& trait,
                                           int max_iter,
                                           double tol,
                                           const std::string& dump_prefix,
                                           const std::vector<double>* prior_offset = nullptr);

} // namespace saige
