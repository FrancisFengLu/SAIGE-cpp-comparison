#include "saige_ai.hpp"
#include <stdexcept>
#include <iostream>
#include "SAIGE_step1_fast.hpp"

namespace saige {

// ----- helpers: safe inverse (mirror your original try/catch) -----
static inline arma::fmat inv_psd_or_pinv(const arma::fmat& A) {
  try {
    return arma::inv_sympd(arma::symmatu(A));
  } catch (const std::exception&) {
    std::cout << "inv_sympd failed, inverted with pinv\n";
    return arma::pinv(arma::symmatu(A));
  }
}

// ================= Non-LOCO =================

static int g_getCoefficients_cpp_call_count = 0;

CoefficientsOut getCoefficients_cpp(const arma::fvec& Y,
                                    const arma::fmat& X,
                                    const arma::fvec& w,
                                    const arma::fvec& tau,
                                    int maxiterPCG, float tolPCG)
{
  g_getCoefficients_cpp_call_count++;
  std::cout << "ENTRY_getCoefficients_cpp_CALLED" << std::endl;
  std::cout << "getCoefficients_cpp call #" << g_getCoefficients_cpp_call_count
            << " |Y|=" << arma::norm(Y) << " |W|=" << arma::norm(w)
            << " tau=[" << tau(0) << "," << tau(1) << "]" << std::endl;

  std::cout << "\n===== C++ getCoefficients_cpp INPUT (call #" << g_getCoefficients_cpp_call_count << ") =====" << std::endl;
  std::cout << "Yvec size: " << Y.n_elem << ", |Yvec|: " << arma::norm(Y) << std::endl;
  std::cout << "Xmat size: " << X.n_rows << "x" << X.n_cols << std::endl;
  std::cout << "wVec size: " << w.n_elem << ", |wVec|: " << arma::norm(w) << std::endl;
  std::cout << "tauVec: [" << tau(0) << ", " << tau(1) << "]" << std::endl;
  std::cout << "Yvec[0:5]: " << Y(0) << " " << Y(1) << " " << Y(2) << " " << Y(3) << " " << Y(4) << std::endl;
  std::cout << "wVec[0:5]: " << w(0) << " " << w(1) << " " << w(2) << " " << w(3) << " " << w(4) << std::endl;
  std::cout << "============================================" << std::endl;
  std::cout << std::flush;
  // Sigma^{-1}Y and Sigma^{-1}X via PCG
  std::cout << "[DEBUG] About to call getPCG1ofSigmaAndVector for Y..." << std::endl << std::flush;
  arma::fvec Sigma_iY = getPCG1ofSigmaAndVector(w, tau, Y, maxiterPCG, tolPCG);
  std::cout << "[DEBUG] getPCG1ofSigmaAndVector for Y done" << std::endl << std::flush;                // :contentReference[oaicite:0]{index=0}
  arma::fmat Sigma_iX(Y.n_rows, X.n_cols);
  for (int j = 0; j < static_cast<int>(X.n_cols); ++j) {
    std::cout << "[DEBUG] getPCG for X col " << j << "..." << std::endl << std::flush;
    Sigma_iX.col(j) = getPCG1ofSigmaAndVector(w, tau, X.col(j), maxiterPCG, tolPCG);
  }
  std::cout << "[DEBUG] All Sigma_iX done" << std::endl << std::flush;

  // cov = (X' Σ^{-1} X)^{-1} with PSD fallback
  arma::fmat cov = inv_psd_or_pinv(X.t() * Sigma_iX);
  std::cout << "[DEBUG] cov done" << std::endl << std::flush;

  // alpha = cov * X' Σ^{-1} Y
  arma::fvec alpha = cov * (Sigma_iX.t() * Y);
  std::cout << "[DEBUG] alpha done" << std::endl << std::flush;

  // eta = Y - τ0 * (Σ^{-1}Y - Σ^{-1}X α) ./ w
  arma::fvec eta = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / w;
  std::cout << "[DEBUG] eta done, returning from getCoefficients_cpp" << std::endl << std::flush;

  return {Sigma_iY, Sigma_iX, cov, alpha, eta};
}

AIScoreOut getAIScore_cpp(const arma::fvec& Y,
                          const arma::fmat& X,
                          const arma::fvec& w,
                          const arma::fvec& tau,
                          const arma::fvec& Sigma_iY,
                          const arma::fmat& Sigma_iX,
                          const arma::fmat& cov,
                          int nrun, int maxiterPCG,
                          float tolPCG, float traceCVcutoff)
{
  // ===== DEBUG: Input values before Trace calculation =====
  std::cout << "\n===== C++ getAIScore_cpp DEBUG =====" << std::endl;
  std::cout << "tau: [" << tau(0) << ", " << tau(1) << "]" << std::endl;
  std::cout << "|Y| (norm): " << arma::norm(Y) << std::endl;
  std::cout << "|w| (norm): " << arma::norm(w) << std::endl;
  std::cout << "|Sigma_iY| (norm): " << arma::norm(Sigma_iY) << std::endl;
  std::cout << "|Sigma_iX| (frobenius): " << arma::norm(Sigma_iX, "fro") << std::endl;
  std::cout << "|cov| (frobenius): " << arma::norm(cov, "fro") << std::endl;
  std::cout << "Y[0:5]: " << Y(0) << " " << Y(1) << " " << Y(2) << " " << Y(3) << " " << Y(4) << std::endl;
  std::cout << "w[0:5]: " << w(0) << " " << w(1) << " " << w(2) << " " << w(3) << " " << w(4) << std::endl;
  std::cout << "Sigma_iY[0:5]: " << Sigma_iY(0) << " " << Sigma_iY(1) << " " << Sigma_iY(2) << " " << Sigma_iY(3) << " " << Sigma_iY(4) << std::endl;

  arma::fmat Sigma_iXt = Sigma_iX.t();
  arma::fvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));                                // :contentReference[oaicite:5]{index=5}
  arma::fvec APY = getCrossprodMatAndKin(PY);                                                    // :contentReference[oaicite:6]{index=6}
  float YPAPY = arma::dot(PY, APY);

  std::cout << "|PY| (norm): " << arma::norm(PY) << std::endl;
  std::cout << "|APY| (norm): " << arma::norm(APY) << std::endl;
  std::cout << "PY[0:5]: " << PY(0) << " " << PY(1) << " " << PY(2) << " " << PY(3) << " " << PY(4) << std::endl;
  std::cout << "APY[0:5]: " << APY(0) << " " << APY(1) << " " << APY(2) << " " << APY(3) << " " << APY(4) << std::endl;
  std::cout << "YPAPY: " << YPAPY << std::endl;
  std::cout << "=====================================" << std::endl;

  float Trace  = GetTrace(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);    // :contentReference[oaicite:8]{index=8}
  arma::fvec PAPY_1 = getPCG1ofSigmaAndVector(w, tau, APY, maxiterPCG, tolPCG);                  // :contentReference[oaicite:9]{index=9}
  arma::fvec PAPY   = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
  float AI = arma::dot(APY, PAPY);                                                               // :contentReference[oaicite:10]{index=10}
  return {YPAPY, Trace, PY, APY, AI};
}

FitAIOut fitglmmaiRPCG_cpp(const arma::fvec& Y,
                           const arma::fmat& X,
                           const arma::fvec& w,
                           arma::fvec tau,
                           const arma::fvec& Sigma_iY,
                           const arma::fmat& Sigma_iX,
                           const arma::fmat& cov,
                           int nrun, int maxiterPCG,
                           float tolPCG, float tol,
                           float traceCVcutoff)
{
  // Single AI step (caller typically loops + checks convergence in higher level)
  AIScoreOut re = getAIScore_cpp(Y, X, w, tau, Sigma_iY, Sigma_iX, cov,
                                 nrun, maxiterPCG, tolPCG, traceCVcutoff);                        // :contentReference[oaicite:11]{index=11}
  // In your R path you solve AI * Δτ = score; here AI is scalar (binary/surv, 1 VC shown)
  // Caller usually computes updated tau outside; we keep tau unchanged here and only
  // return fixed-effect updates like your R fitglmmaiRPCG did after AI step.

  arma::fmat cov_upd = inv_psd_or_pinv(X.t() * Sigma_iX);
  arma::fvec alpha   = cov_upd * (Sigma_iX.t() * Y);
  arma::fvec eta     = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / w;

  FitAIOut out{tau, cov_upd, alpha, eta};
  return out;
}

// ================= Quantitative =================

AIScoreQOut getAIScore_q_cpp(const arma::fvec& Y,
                             arma::fmat& X,
                             arma::fvec& w,
                             arma::fvec& tau,
                             const arma::fvec& Sigma_iY,
                             const arma::fmat& Sigma_iX,
                             arma::fmat& cov,
                             int nrun, int maxiterPCG,
                             float tolPCG, float traceCVcutoff)
{
  arma::fmat Sigma_iXt = Sigma_iX.t();
  arma::fmat Xmatt     = X.t();
  arma::fmat cov1      = inv_psd_or_pinv(Xmatt * Sigma_iX);                                      // :contentReference[oaicite:12]{index=12}

  arma::fvec PY    = Sigma_iY - Sigma_iX * (cov1 * (Sigma_iXt * Y));
  std::cout << "[DEBUG-Q] PY computed, |PY|=" << arma::norm(PY) << std::endl << std::flush;
  arma::fvec APY   = getCrossprodMatAndKin(PY);
  std::cout << "[DEBUG-Q] APY computed, |APY|=" << arma::norm(APY) << std::endl << std::flush;
  float YPAPY      = arma::dot(PY, APY);
  arma::fvec A0PY  = PY;
  float YPA0PY     = arma::dot(PY, A0PY);

  std::cout << "[DEBUG-Q] About to call GetTrace_q..." << std::endl << std::flush;
  arma::fvec Trace = GetTrace_q(Sigma_iX, X, w, tau, cov1, nrun, maxiterPCG, tolPCG, traceCVcutoff);
  std::cout << "[DEBUG-Q] GetTrace_q done" << std::endl << std::flush;

  arma::fvec PA0PY_1 = getPCG1ofSigmaAndVector(w, tau, A0PY, maxiterPCG, tolPCG);
  arma::fvec PA0PY   = PA0PY_1 - Sigma_iX * (cov1 * (Sigma_iXt * PA0PY_1));
  arma::fvec PAPY_1  = getPCG1ofSigmaAndVector(w, tau, APY,   maxiterPCG, tolPCG);
  arma::fvec PAPY    = PAPY_1  - Sigma_iX * (cov1 * (Sigma_iXt * PAPY_1));
  arma::fmat AI(2,2, arma::fill::zeros);
  AI(0,0) = arma::dot(A0PY, PA0PY);
  AI(1,1) = arma::dot(APY,  PAPY);
  AI(0,1) = arma::dot(A0PY, PAPY);
  AI(1,0) = AI(0,1);                                                                              // :contentReference[oaicite:17]{index=17}

  return {YPA0PY, YPAPY, Trace, AI};
}

FitAIOut fitglmmaiRPCG_q_cpp(const arma::fvec& Y,
                              arma::fmat& X,
                              arma::fvec& w,
                              arma::fvec  tau,
                              const arma::fvec& Sigma_iY,
                              const arma::fmat& Sigma_iX,
                              arma::fmat& cov,
                              int nrun, int maxiterPCG,
                              float tolPCG, float /*tol*/,
                              float traceCVcutoff)
{
  AIScoreQOut re = getAIScore_q_cpp(Y, X, w, tau, Sigma_iY, Sigma_iX, cov,
                                    nrun, maxiterPCG, tolPCG, traceCVcutoff);                     // :contentReference[oaicite:18]{index=18}
  // As in your R path, do a fixed-effect update here (tau update strategy is up to caller).
  arma::fmat cov_upd = inv_psd_or_pinv(X.t() * Sigma_iX);
  arma::fvec alpha   = cov_upd * (Sigma_iX.t() * Y);
  arma::fvec eta     = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / w;

  FitAIOut out{tau, cov_upd, alpha, eta};
  return out;
}

// ================= LOCO =================

CoefficientsOut getCoefficients_LOCO_cpp(const arma::fvec& Y,
                                         const arma::fmat& X,
                                         const arma::fvec& w,
                                         const arma::fvec& tau,
                                         int maxiterPCG, float tolPCG)
{
  arma::fvec Sigma_iY = getPCG1ofSigmaAndVector_LOCO(w, tau, Y, maxiterPCG, tolPCG);             // :contentReference[oaicite:19]{index=19}
  arma::fmat Sigma_iX(Y.n_rows, X.n_cols);
  for (int j = 0; j < static_cast<int>(X.n_cols); ++j) {
    Sigma_iX.col(j) = getPCG1ofSigmaAndVector_LOCO(w, tau, X.col(j), maxiterPCG, tolPCG);        // :contentReference[oaicite:20]{index=20}
  }
  arma::fmat cov = inv_psd_or_pinv(X.t() * Sigma_iX);                                            // :contentReference[oaicite:21]{index=21}
  arma::fvec alpha = cov * (Sigma_iX.t() * Y);
  arma::fvec eta   = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / w;                              // :contentReference[oaicite:22]{index=22}
  return {Sigma_iY, Sigma_iX, cov, alpha, eta};
}

} // namespace saige
