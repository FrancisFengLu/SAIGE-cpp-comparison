#pragma once
#include <armadillo>

namespace saige {

// ---------- Plain C++ return types ----------
struct CoefficientsOut {
  arma::fvec Sigma_iY;   // n
  arma::fmat Sigma_iX;   // n x p
  arma::fmat cov;        // p x p  = (X' Σ^{-1} X)^{-1}
  arma::fvec alpha;      // p
  arma::fvec eta;        // n      = Y - τ0*(Σ^{-1}Y - Σ^{-1}X α) ./ w
};

struct AIScoreOut {
  float      YPAPY;      // scalar score term
  float      Trace;      // scalar trace
  arma::fvec PY;         // n      = Σ^{-1}Y - Σ^{-1}X cov X'Σ^{-1}Y
  arma::fvec APY;        // n      = K * PY (GRM times PY)
  float      AI;         // scalar average-information entry
};

struct AIScoreQOut {
  float       YPA0PY;    // quantitative-only
  float       YPAPY;
  arma::fvec  Trace;     // length 2: [Trace(A0P), Trace(AP)]
  arma::fmat  AI;        // 2x2
};

struct FitAIOut {
  arma::fvec tau;   // variance components (size 2 in your code)
  arma::fmat cov;   // p x p
  arma::fvec alpha; // p
  arma::fvec eta;   // n
};

// ---------- Forward decls: implemented in your existing TU ----------
arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w,
                                   const arma::fvec& tau,
                                   const arma::fvec& v,
                                   int maxiterPCG, float tolPCG);

arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& w,
                                        const arma::fvec& tau,
                                        const arma::fvec& v,
                                        int maxiterPCG, float tolPCG);

arma::fvec getCrossprodMatAndKin(const arma::fvec& v);
float      GetTrace(const arma::fmat& Sigma_iX, const arma::fmat& X,
                    const arma::fvec& w, const arma::fvec& tau,
                    const arma::fmat& cov, int nrun,
                    int maxiterPCG, float tolPCG, float traceCVcutoff);

arma::fvec GetTrace_q(const arma::fmat& Sigma_iX, arma::fmat& X,
                      arma::fvec& w, arma::fvec& tau,
                      arma::fmat& cov, int nrun,
                      int maxiterPCG, float tolPCG, float traceCVcutoff);

// ---------- Non-LOCO ----------
CoefficientsOut getCoefficients_cpp(const arma::fvec& Y,
                                    const arma::fmat& X,
                                    const arma::fvec& w,
                                    const arma::fvec& tau,
                                    int maxiterPCG, float tolPCG);

AIScoreOut      getAIScore_cpp(const arma::fvec& Y,
                               const arma::fmat& X,
                               const arma::fvec& w,
                               const arma::fvec& tau,
                               const arma::fvec& Sigma_iY,
                               const arma::fmat& Sigma_iX,
                               const arma::fmat& cov,
                               int nrun, int maxiterPCG,
                               float tolPCG, float traceCVcutoff);

FitAIOut        fitglmmaiRPCG_cpp(const arma::fvec& Y,
                                  const arma::fmat& X,
                                  const arma::fvec& w,
                                  arma::fvec tau,
                                  const arma::fvec& Sigma_iY,
                                  const arma::fmat& Sigma_iX,
                                  const arma::fmat& cov,
                                  int nrun, int maxiterPCG,
                                  float tolPCG, float tol,
                                  float traceCVcutoff);

// ---------- Quantitative ----------
AIScoreQOut     getAIScore_q_cpp(const arma::fvec& Y,
                                 arma::fmat& X,
                                 arma::fvec& w,
                                 arma::fvec& tau,
                                 const arma::fvec& Sigma_iY,
                                 const arma::fmat& Sigma_iX,
                                 arma::fmat& cov,
                                 int nrun, int maxiterPCG,
                                 float tolPCG, float traceCVcutoff);

FitAIOut        fitglmmaiRPCG_q_cpp(const arma::fvec& Y,
                                    arma::fmat& X,
                                    arma::fvec& w,
                                    arma::fvec  tau,
                                    const arma::fvec& Sigma_iY,
                                    const arma::fmat& Sigma_iX,
                                    arma::fmat& cov,
                                    int nrun, int maxiterPCG,
                                    float tolPCG, float tol,
                                    float traceCVcutoff);

// ---------- LOCO variants ----------
CoefficientsOut getCoefficients_LOCO_cpp(const arma::fvec& Y,
                                         const arma::fmat& X,
                                         const arma::fvec& w,
                                         const arma::fvec& tau,
                                         int maxiterPCG, float tolPCG);

// If you need LOCO quantitative AI-score/fits, add symmetrical *_q_LOCO wrappers.

} // namespace saige
