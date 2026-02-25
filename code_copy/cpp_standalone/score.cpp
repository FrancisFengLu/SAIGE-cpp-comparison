#include "score.hpp"
#include "saige_null.hpp"
#include "SAIGE_step1_fast.hpp"
namespace saige {

// ===== P1 lazy apply =====
arma::fvec P1Projector::apply(const arma::fvec& v) const {
  // s = Σ^{-1} v
  arma::fvec s = SigmaInvApply(w, tau, v, maxiterPCG, tolPCG);
  // t = X' s
  arma::fvec t = Sigma_iX.t() * v;             // note: v is already Σ^{-1}·something in typical use
  // proj = Σ^{-1}X cov (X' Σ^{-1} v)
  arma::fvec proj = Sigma_iX * (cov * (Sigma_iX.t() * v));
  // P1 v = s - proj
  return s - proj;
}

// ===== Score-null (binary/quant) =====
static ScoreNull build_score_null_generic(const arma::fmat& X,
                                          const arma::fvec& y,
                                          const arma::fvec& mu,
                                          const arma::fvec& V)
{
  const arma::uword n = X.n_rows, p = X.n_cols;

  ScoreNull out;
  out.V   = V;
  out.res = y - mu;

  // XV = (X ⊙ V)'  i.e. each column of X multiplied by V elementwise, then transposed
  out.XV.set_size(p, n);
  for (arma::uword j = 0; j < p; ++j) out.XV.row(j) = (X.col(j) % V).t();

  out.XVX      = X.t() * (X.each_col() % V);
  out.XVX_inv  = arma::inv_sympd(arma::symmatu(out.XVX));
  out.XXVX_inv = X * out.XVX_inv;

  // XVX_inv_XV = (X * XVX_inv) ⊙ Vrow; we keep as n x p matrix that equals
  //              [X * XVX_inv] row-wise scaled by V (so multiplying by a vector reproduces original)
  out.XVX_inv_XV = out.XXVX_inv.each_col() % V;

  // S_a = colSums(X ⊙ res)
  out.S_a.set_size(p);
  for (arma::uword j = 0; j < p; ++j) out.S_a(j) = arma::sum(X.col(j) % out.res);

  return out;
}

ScoreNull build_score_null_binary(const arma::fmat& X,
                                  const arma::fvec& y,
                                  const arma::fvec& mu)
{
  arma::fvec V = mu % (1.0f - mu);   // mu2
  return build_score_null_generic(X, y, mu, V);
}

ScoreNull build_score_null_quant(const arma::fmat& X,
                                 const arma::fvec& y,
                                 const arma::fvec& mu,
                                 float tau0_inv)
{
  arma::fvec V(X.n_rows, arma::fill::value(tau0_inv)); // constant
  return build_score_null_generic(X, y, mu, V);
}

// ===== Score-null (survival) =====
ScoreNullSurv build_score_null_survival(const arma::fmat& X1,
                                        const arma::fvec& y,
                                        const arma::fvec& mu)
{
  // This mirrors your R ‘ScoreTest_NULL_Model_survival’ shape:
  // V = mu, res = y - mu, and an “fg” design with appended 1-column.
  ScoreNullSurv out;
  out.y = y; out.mu = mu; out.res = y - mu;
  out.V = mu;
  out.X1 = X1;

  // XV, XVX on X1
  out.XV  = (X1.each_col() % out.V).t();
  out.XVX = X1.t() * (X1.each_col() % out.V);
  out.XVX_inv = arma::inv_sympd(arma::symmatu(out.XVX));
  out.XXVX_inv = X1 * out.XVX_inv;
  out.XVX_inv_XV = out.XXVX_inv.each_col() % out.V;

  // fg: append intercept
  arma::fmat X1_fg(X1.n_rows, X1.n_cols + 1);
  X1_fg.cols(0, X1.n_cols - 1) = X1;
  X1_fg.col(X1.n_cols).ones();
  out.X1_fg = X1_fg;

  out.XV_fg  = (X1_fg.each_col() % out.V).t();
  out.XVX_fg = X1_fg.t() * (X1_fg.each_col() % out.V);
  out.XVX_inv_fg = arma::inv_sympd(arma::symmatu(out.XVX_fg));
  out.XXVX_inv_fg = X1_fg * out.XVX_inv_fg;
  out.XVX_inv_XV_fg = out.XXVX_inv_fg.each_col() % out.V;

  // S_a on X1_fg
  out.S_a.set_size(out.X1_fg.n_cols);
  for (arma::uword j = 0; j < out.X1_fg.n_cols; ++j) out.S_a(j) = arma::sum(out.X1_fg.col(j) % out.res);
  return out;
}

// ===== Convert ScoreNull → ScoreNullPack for serialization =====
ScoreNullPack to_pack(const ScoreNull& sn,
                       const arma::fmat& X_orig,
                       const arma::fvec& y_vec,
                       const arma::fvec& mu_vec,
                       const std::string& trait)
{
  ScoreNullPack pk;
  pk.n = sn.V.n_elem;
  pk.p = sn.S_a.n_elem;
  pk.trait_type = trait;

  // N-length vectors (float → double)
  pk.V   = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(sn.V));
  pk.mu  = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(mu_vec));
  pk.res = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(sn.res));
  pk.y   = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(y_vec));

  // p-length vector
  pk.S_a = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(sn.S_a));

  // Flatten matrices row-major (arma is col-major, so transpose first)
  auto flatten = [](const arma::fmat& m) -> std::vector<double> {
    arma::mat md = arma::conv_to<arma::mat>::from(m);
    arma::mat mt = md.t();  // row-major = transpose of col-major
    return arma::conv_to<std::vector<double>>::from(arma::vectorise(mt));
  };

  pk.XV         = flatten(sn.XV);
  pk.XVX        = flatten(sn.XVX);
  pk.XVX_inv    = flatten(sn.XVX_inv);
  pk.XXVX_inv   = flatten(sn.XXVX_inv);
  pk.XVX_inv_XV = flatten(sn.XVX_inv_XV);
  pk.X          = flatten(X_orig);

  return pk;
}

} // namespace saige
