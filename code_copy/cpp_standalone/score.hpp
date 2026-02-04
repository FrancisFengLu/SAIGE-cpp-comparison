#pragma once
#include <armadillo>
#include <vector>

namespace saige {

// ---------- Lazy projector: applies P1·v without forming P1 ----------
struct P1Projector {
  // Σ^{-1}Y and Σ^{-1}X come from getCoefficients_cpp(...)
  arma::fvec Sigma_iY;        // unused here but kept for completeness
  arma::fmat Sigma_iX;        // n x p
  arma::fmat cov;             // p x p = (X' Σ^{-1} X)^{-1}
  // PCG apply for Σ^{-1}·v
  arma::fvec (*SigmaInvApply)(const arma::fvec& w, const arma::fvec& tau,
                              const arma::fvec& v, int iters, float tol);
  // params captured for Σ^{-1}
  arma::fvec w;               // working weights (n)
  arma::fvec tau;             // variance comps (size 2)
  int        maxiterPCG{500};
  float      tolPCG{1e-5f};

  // P1·v = Σ^{-1}v - Σ^{-1}X cov X' Σ^{-1}v
  arma::fvec apply(const arma::fvec& v) const;
};

// ---------- Binary/quant: score-null pack (matches your R function needs) ----------
struct ScoreNull {
  arma::fmat XV;          // p x n  (== t(X * V)), V = mu2 (binary) or 1/tau0 (quant)
  arma::fmat XVX;         // p x p  (== X' V X)
  arma::fmat XVX_inv;     // p x p
  arma::fmat XXVX_inv;    // n x p  (== X * XVX_inv)
  arma::fmat XVX_inv_XV;  // n x n  implicit carrier as (XXVX_inv * V), we store as n x p with V applied: (X * XVX_inv) ⊙ Vrow
  arma::fvec V;           // n
  arma::fvec res;         // n
  arma::fvec S_a;         // p, == colSums(X ⊙ res)
};

// ---------- Survival variant (mirrors your R survival score-null) ----------
struct ScoreNullSurv {
  arma::fvec y, mu, res, V;      // V = mu for survival
  arma::fmat X1;                 // original X (possibly transformed)
  arma::fmat XV, XVX, XXVX_inv, XVX_inv, XVX_inv_XV;
  // “fg” = with intercept column appended, as your R code did:
  arma::fmat X1_fg, XV_fg, XVX_fg, XXVX_inv_fg, XVX_inv_fg, XVX_inv_XV_fg;
  arma::fvec S_a;                // on X1_fg
};

// ---------- Builders (no Rcpp) ----------
ScoreNull      build_score_null_binary(const arma::fmat& X,
                                       const arma::fvec& y,
                                       const arma::fvec& mu);

ScoreNull      build_score_null_quant (const arma::fmat& X,
                                       const arma::fvec& y,
                                       const arma::fvec& mu,
                                       float tau0_inv);  // 1/tau[0]

ScoreNullSurv  build_score_null_survival(const arma::fmat& X1,
                                         const arma::fvec& y,
                                         const arma::fvec& mu);

} // namespace saige
