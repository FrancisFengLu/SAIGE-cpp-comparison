#include "covariate_offset.hpp"
#include <RcppArmadillo.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace saige {

static bool ieq_local(const std::string& a, const char* b) {
  const size_t bl = std::string(b).size();
  if (a.size() != bl) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i])) return false;
  return true;
}

static double logistic_deviance(const arma::vec& y, const arma::vec& mu_clipped) {
  return -2.0 * arma::sum(y % arma::log(mu_clipped) + (1.0 - y) % arma::log(1.0 - mu_clipped));
}

CovariateOffsetResult fit_covariate_offset(const Design& design,
                                           const std::string& trait,
                                           int max_iter,
                                           double tol,
                                           const std::string& dump_prefix,
                                           const std::vector<double>* prior_offset) {
  CovariateOffsetResult res;
  if (design.p <= 0 || design.n <= 0) return res;

  const int n = design.n;
  const int p = design.p;

  // Build X (n x p) from row-major design.X
  arma::mat X(n, p);
  arma::vec yv(n);
  arma::vec off0(n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) yv(i) = design.y[i];
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      X(i, j) = design.X[(size_t)i * (size_t)p + (size_t)j];
  if (prior_offset && (int)prior_offset->size() == n)
    for (int i = 0; i < n; ++i) off0(i) = (*prior_offset)[i];
  const bool have_off = (prior_offset && (int)prior_offset->size() == n);

  arma::vec beta(p, arma::fill::zeros);
  const bool is_binary = ieq_local(trait, "binary");

  if (!is_binary) {
    // Quantitative with optional offset:  β = (X'X + λI)⁻¹ X'(y − prior_offset)
    const double lam = 1e-8;
    beta = arma::solve(X.t() * X + lam * arma::eye(p, p), X.t() * (yv - off0));
    res.iterations = 1;
    res.converged = true;
  } else {
    // Binary: logistic IRLS with deviance check + step halving, with offset
    arma::vec b(p, arma::fill::zeros);
    arma::vec eta, mu, w, z, b_trial;

    for (int it = 0; it < max_iter; ++it) {
      eta = X * b + off0;
      mu = 1.0 / (1.0 + arma::exp(-eta));
      arma::vec mu_c = arma::clamp(mu, 1e-8, 1.0 - 1e-8);
      w = mu_c % (1.0 - mu_c);
      // Working response in terms of Xβ only (subtract offset contribution)
      z = (eta - off0) + (yv - mu_c) / w;
      arma::mat XtW = X.t() * arma::diagmat(w);

      try {
        b_trial = arma::solve(XtW * X, XtW * z);
      } catch (const std::exception& e) {
        std::cerr << "[cov_offset] solve failed at iter " << it << ": " << e.what() << "\n";
        res.converged = false;
        break;
      }

      double dev_old = logistic_deviance(yv, mu_c);
      res.dev_trace.push_back(dev_old);

      // Step halving: require new deviance to not increase (within small tol)
      int halve = 0;
      while (halve < 20) {
        arma::vec eta_t = X * b_trial + off0;
        arma::vec mu_t = arma::clamp(1.0 / (1.0 + arma::exp(-eta_t)), 1e-8, 1.0 - 1e-8);
        double dev_new = logistic_deviance(yv, mu_t);
        if (std::isfinite(dev_new) && dev_new <= dev_old + 1e-8) break;
        b_trial = 0.5 * (b_trial + b);
        halve++;
      }

      double delta = arma::norm(b_trial - b, "inf");
      b = b_trial;
      res.iterations = it + 1;

      if (delta < tol) {
        res.converged = true;
        break;
      }
    }
    beta = b;
    if (res.iterations == max_iter && res.dev_trace.size() >= 2) {
      double last  = res.dev_trace.back();
      double prev  = res.dev_trace[res.dev_trace.size() - 2];
      if (std::abs(last - prev) > tol) res.converged = false;
    }
  }

  // Final eta, mu, offset
  arma::vec eta_f = X * beta + off0;
  arma::vec mu_f;
  if (is_binary) mu_f = 1.0 / (1.0 + arma::exp(-eta_f));
  else           mu_f = eta_f;

  arma::vec off(n, arma::fill::zeros);
  if (p > 1) off = X.cols(1, p - 1) * beta.subvec(1, p - 1);

  res.beta.assign(beta.begin(), beta.end());
  res.eta.assign(eta_f.begin(), eta_f.end());
  res.mu.assign(mu_f.begin(), mu_f.end());
  res.offset.assign(off.begin(), off.end());

  // CSV dumps
  if (!dump_prefix.empty()) {
    {
      std::ofstream f(dump_prefix + "_cov_input.csv");
      f << "iid,y";
      for (int j = 0; j < p; ++j) f << ",X" << j;
      f << "\n";
      for (int i = 0; i < n; ++i) {
        f << (i < (int)design.iid.size() ? design.iid[i] : std::to_string(i))
          << "," << yv(i);
        for (int j = 0; j < p; ++j) f << "," << X(i, j);
        f << "\n";
      }
    }
    {
      std::ofstream f(dump_prefix + "_cov_beta.csv");
      f << "col,beta\n";
      for (int j = 0; j < p; ++j) f << j << "," << beta(j) << "\n";
    }
    {
      std::ofstream f(dump_prefix + "_cov_fit.csv");
      f << "sample,eta,mu,offset\n";
      for (int i = 0; i < n; ++i)
        f << i << "," << eta_f(i) << "," << mu_f(i) << "," << off(i) << "\n";
    }
    if (is_binary) {
      std::ofstream f(dump_prefix + "_cov_trace.csv");
      f << "iter,deviance\n";
      for (size_t it = 0; it < res.dev_trace.size(); ++it)
        f << it << "," << res.dev_trace[it] << "\n";
    }
    std::cout << "[cov_offset] dumped CSVs to: " << dump_prefix
              << "_cov_{input,beta,fit" << (is_binary ? ",trace" : "") << "}.csv\n";
  }

  std::cout << "[cov_offset] trait=" << trait << " n=" << n << " p=" << p
            << " iters=" << res.iterations << " converged=" << (res.converged ? "true" : "false")
            << " beta[0]=" << beta(0);
  if (p > 1) std::cout << " beta[1]=" << beta(1);
  std::cout << " |offset|=" << arma::norm(off);
  if (is_binary && !res.dev_trace.empty())
    std::cout << " final_dev=" << res.dev_trace.back();
  std::cout << "\n";

  return res;
}

} // namespace saige
