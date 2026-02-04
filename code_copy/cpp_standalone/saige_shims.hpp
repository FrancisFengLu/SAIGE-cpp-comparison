#pragma once
#include <armadillo>

// Declare the existing global definitions from SAIGE_step1_fast.cpp:
extern arma::fvec getCrossprodMatAndKin(arma::fvec& bVec);
extern arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w, const arma::fvec& tau, const arma::fvec& v,
                                          int maxiterPCG, float tolPCG);
extern arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& w, const arma::fvec& tau, const arma::fvec& v,
                                               int maxiterPCG, float tolPCG);
extern float      GetTrace(const arma::fmat Sigma_iX, const arma::fmat& X, const arma::fvec& w, const arma::fvec& tau,
                           const arma::fmat& cov, int nrun, int maxiterPCG, float tolPCG,
                           float traceCVcutoff);
extern arma::fvec GetTrace_q(arma::fmat Sigma_iX, arma::fmat& X, arma::fvec& w, arma::fvec& tau,
                             arma::fmat& cov, int nrun, int maxiterPCG, float tolPCG,
                             float traceCVcutoff);

// Provide the namespaced, const-friendly wrappers that saige_ai.hpp promised:
namespace saige {

arma::fvec getCrossprodMatAndKin(const arma::fvec& v) {
  arma::fvec tmp = v;
  return ::getCrossprodMatAndKin(tmp);
}
extern arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w, const arma::fvec& tau, const arma::fvec& v,
                                          int maxiterPCG, float tolPCG);
extern arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& w, const arma::fvec& tau, const arma::fvec& v,
                                               int maxiterPCG, float tolPCG);

// float GetTrace(const arma::fmat& Sigma_iX, const arma::fmat& X, const arma::fvec& w,
//                const arma::fvec& tau, const arma::fmat& cov, int nrun,
//                int maxiterPCG, float tolPCG, float traceCVcutoff) {
//   arma::fmat Sigma_iXc = Sigma_iX, Xc = X, covc = cov;
//   arma::fvec wc = w, tc = tau;
//   return ::GetTrace(Sigma_iXc, Xc, wc, tc, covc, nrun, maxiterPCG, tolPCG, traceCVcutoff);
// }

arma::fvec GetTrace_q(const arma::fmat& Sigma_iX, arma::fmat& X, arma::fvec& w, arma::fvec& tau,
                      arma::fmat& cov, int nrun, int maxiterPCG, float tolPCG,
                      float traceCVcutoff) {
  arma::fmat Sigma_iXc = Sigma_iX;
  return ::GetTrace_q(Sigma_iXc, X, w, tau, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
}

} // namespace saige
