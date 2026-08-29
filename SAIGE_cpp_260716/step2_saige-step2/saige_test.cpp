// Standalone port of SAIGE/src/SAIGE_test.cpp
// SAIGEClass: score test, SPA, variance ratio, Firth correction

#define ARMA_USE_SUPERLU 1

#include <armadillo>

#include "saige_test.hpp"
#include "spa.hpp"
#include "er_binary.hpp"
#include "UTIL.hpp"
#include "getMem.hpp"

// A3 globals (declared extern in saige_test.hpp)
thread_local bool g_firthDefer = false;
std::atomic<std::uint64_t> g_firthFitCalls{0};
#include <thread>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

namespace SAIGE {

// P1 fix (2026-05-09): direct log P(X > stat) for chi-square(df=1) without
// underflow. The previous code computed log(boost::math::cdf(...)) which
// underflows to -inf when the cdf itself underflows below ~1e-300 (i.e. at
// genuinely significant markers). This matches R's pchisq(q, 1, lower.tail=FALSE,
// log.p=TRUE).
//
// chi-square(df=1) upper tail: P(X > stat) = erfc(sqrt(stat/2))
// Boost provides erfc; for very large argument it underflows too, so use
// the asymptotic expansion log(erfc(z)) ≈ -z² - log(z) - 0.5*log(π) for z > 6.
static inline double log_chisq1_uppertail(double stat) {
    if (!std::isfinite(stat) || stat <= 0.0) return 0.0;  // log(p=1) = 0
    const double z = std::sqrt(stat / 2.0);
    if (z > 6.0) {
        // first-order asymptotic for erfc(z), accurate to ~1e-3 already at z=4
        // and machine precision at z=10+.
        return -z * z - std::log(z) - 0.5 * std::log(boost::math::constants::pi<double>());
    }
    const double e = boost::math::erfc(z);
    if (e <= 0.0 || !std::isfinite(e)) {
        // fallback: same asymptotic if boost's erfc itself returned 0
        return -z * z - std::log(z) - 0.5 * std::log(boost::math::constants::pi<double>());
    }
    return std::log(e);
}

// Fused-kernel mode + A/B validation accumulators (Pillar 1). See saige_test.hpp.
int    g_fusedMode        = 0;
double g_fusedMaxRelTstat = 0.0;
double g_fusedMaxRelVar1  = 0.0;
long   g_fusedNCompared   = 0;
long   g_fusedNExceed     = 0;


SAIGEClass::SAIGEClass(
	arma::mat & t_XVX,
	arma::mat  t_XXVX_inv,
	arma::mat & t_XV,
	arma::mat & t_XVX_inv_XV,
	arma::mat & t_Sigma_iXXSigma_iX,
	arma::mat & t_X,
	arma::vec &  t_S_a,
	arma::vec & t_res,
	arma::vec & t_mu2,
	arma::vec & t_mu,
	arma::vec & t_varRatio_sparse,
	arma::vec & t_varRatio_null,
	arma::vec & t_varRatio_null_noXadj,
	arma::vec & t_cateVarRatioMinMACVecExclude,
        arma::vec & t_cateVarRatioMaxMACVecInclude,
	double t_SPA_Cutoff,
	arma::vec & t_tauvec,
	std::string t_traitType,
	arma::vec & t_y,
	std::string t_impute_method,
	bool t_flagSparseGRM,
	bool t_isFastTest,
	bool t_isnoadjCov,
	double t_pval_cutoff_for_fastTest,
	arma::umat & t_locationMat,
	arma::vec & t_valueVec,
        int t_dimNum,
	bool t_isCondition,
        std::vector<uint32_t> & t_condition_genoIndex,
	bool t_is_Firth_beta,
        double t_pCutoffforFirth,
        arma::vec & t_offset,
	arma::vec & t_resout){


    m_XVX = t_XVX;
    m_XV = t_XV;
    m_XXVX_inv = t_XXVX_inv;
    m_XVX_inv_XV = t_XVX_inv_XV;
    m_Sigma_iXXSigma_iX = t_Sigma_iXXSigma_iX;
    m_isVarPsadj = false;
    if(m_Sigma_iXXSigma_iX.n_cols == 1 && m_Sigma_iXXSigma_iX.n_rows == 1){
	m_isVarPsadj = false;
    }else{
	m_isVarPsadj = true;
    }
    m_X = t_X;
    // Pillar 1: transposed copies (p x N) so each carrier's covariate row is a
    // contiguous column in scoreTestFast_fused. Read-only after construction.
    m_Xt = m_X.t();
    m_At = m_XVX_inv_XV.t();
    m_S_a = t_S_a;
    m_res = t_res;
    m_resout = t_resout;
    m_mu2 = t_mu2;
    // Pillar 1 (noadjCov variant): precompute the marker-independent full-N sums
    // once (scoreTestFast_noadjCov recomputes arma::accu each marker). Same value.
    m_sum_mu2 = arma::accu(m_mu2);
    m_sum_res = arma::accu(m_res);
    m_mu = t_mu;
    m_varRatio_sparse = t_varRatio_sparse;
    m_varRatio_null = t_varRatio_null;
    m_varRatio_null_noXadj = t_varRatio_null_noXadj;
    m_cateVarRatioMinMACVecExclude = t_cateVarRatioMinMACVecExclude;
    m_cateVarRatioMaxMACVecInclude = t_cateVarRatioMaxMACVecInclude;
    m_tauvec = t_tauvec;
    m_traitType = t_traitType;
    m_y = t_y;

    m_case_indices = arma::find(m_y == 1);
    m_ctrl_indices = arma::find(m_y == 0);

    m_n = t_y.size();
    m_p = t_XV.n_rows;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_impute_method =  t_impute_method;
    m_isCondition = t_isCondition;
    m_condition_genoIndex = t_condition_genoIndex;
    if(m_isCondition){
	        m_numMarker_cond = t_condition_genoIndex.size();
    }else{
		m_numMarker_cond = 0;
    }

    if(m_traitType == "binary"){
        m_case_indices = arma::find(m_y == 1);
        m_ctrl_indices = arma::find(m_y == 0);
	m_n_case = m_case_indices.n_elem;
	m_n_ctrl = m_ctrl_indices.n_elem;
	m_is_Firth_beta = t_is_Firth_beta;
	m_pCutoffforFirth = t_pCutoffforFirth;
	m_offset = t_offset;
    }
    m_flagSparseGRM = t_flagSparseGRM;
    m_isFastTest = t_isFastTest;
    m_isnoadjCov = t_isnoadjCov;
    m_pval_cutoff_for_fastTest = t_pval_cutoff_for_fastTest;
    if(t_dimNum != 0){
        m_spSigmaMat = arma::sp_mat(t_locationMat, t_valueVec, t_dimNum, t_dimNum);
	m_diagSigma = arma::vec(m_spSigmaMat.diag());
    }
}


// Replaced R's set.seed with std::mt19937
void SAIGEClass::set_seed(unsigned int seed){
  m_rng_engine.seed(seed);
}

void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
		     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     arma::vec & t_P2Vec,
		     double& t_gy,
		     bool t_is_region,
		     arma::uvec & t_indexForNonZero){
    PerMarkerCtx ctx{m_flagSparseGRM_cur, m_isnoadjCov_cur, m_varRatioVal};
    scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_str, t_pval, t_islogp,
              t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec,
              t_gy, t_is_region, t_indexForNonZero, ctx);
}

// Ctx overload (Phase A).
void SAIGEClass::scoreTest(arma::vec & t_GVec,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
		     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     arma::vec & t_gtilde,
		     arma::vec & t_P2Vec,
		     double& t_gy,
		     bool t_is_region,
		     arma::uvec & t_indexForNonZero,
                     const PerMarkerCtx& ctx){
    arma::vec Sm, var2m;
    double S, var2;
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);

    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
     }
    S = dot(t_gtilde, m_res);
    S = S/m_tauvec[0];

    if(!ctx.flagSparseGRM_cur){
      t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
      var2m = dot(t_P2Vec , t_gtilde);
    }else{
      t_P2Vec = getPCG1ofSigmaAndGtilde(t_gtilde, 100, 0.02);
      var2m = dot(t_P2Vec , t_gtilde);
      if(m_isVarPsadj){
	var2m = var2m - t_gtilde.t() * m_Sigma_iXXSigma_iX * m_X.t() * t_P2Vec;
      }
    }
    var2 = var2m(0,0);
    double var1 = var2 * ctx.varRatioVal;

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
	  stat = 0.0;
      }
    }
    char pValueBuf[100];

    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
	t_islogp = false;
    }else{
        // R::pchisq(stat,1,false,true) = log(P(X > stat)) for chi-squared(1)
        // P1 fix: direct log(upper-tail) avoiding cdf underflow at p<1e-300
        double logp = log_chisq1_uppertail(stat);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
        }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
	t_pval = logp;
	t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){
    PerMarkerCtx ctx{m_flagSparseGRM_cur, m_isnoadjCov_cur, m_varRatioVal};
    scoreTestFast(t_GVec, t_indexForNonZero, t_Beta, t_seBeta, t_pval_str,
                  t_pval, t_islogp, t_altFreq, t_Tstat, t_var1, t_var2, ctx);
}

// Ctx overload (Phase A).
void SAIGEClass::scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
		     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     const PerMarkerCtx& ctx){

    // --- thread_local persistent scratch (2026-07-15) ---------------------
    // The carrier-subset gathers below (.elem/.rows) previously allocated a
    // fresh arma vec/mat every marker; at whole-cohort N those exceed the
    // mallopt M_MMAP_THRESHOLD and mmap/munmap on each call -> page faults +
    // per-process mmap_lock contention (the "score-gather" + kernel buckets in
    // profiling). Grow-only thread_local buffers, reused via exact-size aliases,
    // remove the per-marker allocation. Arithmetic is unchanged -> bit-identical.
    const arma::uword nnz = t_indexForNonZero.n_elem;
    const arma::uword p   = m_X.n_cols;
    thread_local arma::vec tl_g1, tl_res1, tl_mu21, tl_B, tl_g1t;
    thread_local arma::mat tl_X1, tl_A1;
    if(tl_g1.n_elem < nnz){ tl_g1.set_size(nnz); tl_res1.set_size(nnz); tl_mu21.set_size(nnz); tl_B.set_size(nnz); tl_g1t.set_size(nnz); }
    if(tl_X1.n_rows < nnz){ tl_X1.set_size(nnz, p); tl_A1.set_size(nnz, p); }
    // Exact-size (nnz) aliases over the persistent buffers; no allocation. The
    // physical buffers may be larger (grow-only) but each alias is written and
    // read with the same nnz-stride, so the layout is self-consistent.
    arma::vec g1(tl_g1.memptr(), nnz, false, false);
    arma::vec res1(tl_res1.memptr(), nnz, false, false);
    arma::vec B(tl_B.memptr(), nnz, false, false);
    arma::vec g1_tilde(tl_g1t.memptr(), nnz, false, false);
    arma::mat X1(tl_X1.memptr(), nnz, p, false, false);
    arma::mat A1(tl_A1.memptr(), nnz, p, false, false);
    // Manual gather (no arma temporaries): restrict G, res, and the X / XVX_inv_XV
    // rows to the carrier indices. Same values arma's .elem/.rows would produce.
    {
      const arma::uword* idxp = t_indexForNonZero.memptr();
      const double* Gp   = t_GVec.memptr();
      const double* resp = m_res.memptr();
      for(arma::uword k=0;k<nnz;++k){ const arma::uword s=idxp[k]; g1[k]=Gp[s]; res1[k]=resp[s]; }
      for(arma::uword j=0;j<p;++j){
        const double* Xc = m_X.colptr(j);
        const double* Ac = m_XVX_inv_XV.colptr(j);
        double* X1c = X1.colptr(j);
        double* A1c = A1.colptr(j);
        for(arma::uword k=0;k<nnz;++k){ const arma::uword s=idxp[k]; X1c[k]=Xc[s]; A1c[k]=Ac[s]; }
      }
    }
    arma::vec Z = A1.t() * g1;
    B = X1 * Z;
    g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary" || m_traitType == "survival"){
      arma::vec mu21(tl_mu21.memptr(), nnz, false, false);
      { const arma::uword* idxp=t_indexForNonZero.memptr(); const double* mu2p=m_mu2.memptr();
        for(arma::uword k=0;k<nnz;++k) mu21[k]=mu2p[idxp[k]]; }
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }

    var1 = var2 * ctx.varRatioVal;
    S1 = dot(res1, g1_tilde);
    arma::mat res1X1_temp = (res1.t()) * X1;
    arma::vec res1X1 = res1X1_temp.t();
    S_a2 = m_S_a - res1X1;
    S2 = - arma::dot(S_a2,  Z);
    S = S1 + S2;
    S = S/m_tauvec[0];

    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));

      }else{
          t_pval = 1;
	  stat = 0.0;
      }
    }
    char pValueBuf[100];
    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
	t_islogp = false;
    }else {
        // P1 fix: direct log(upper-tail) avoiding cdf underflow at p<1e-300
	double logp = log_chisq1_uppertail(stat);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
         }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
	t_pval = logp;
	t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}


void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
		     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2){
    PerMarkerCtx ctx{m_flagSparseGRM_cur, m_isnoadjCov_cur, m_varRatioVal};
    scoreTestFast_noadjCov(t_GVec, t_indexForNonZero, t_Beta, t_seBeta,
                           t_pval_str, t_pval, t_islogp, t_altFreq,
                           t_Tstat, t_var1, t_var2, ctx);
}

// Ctx overload (Phase A).
void SAIGEClass::scoreTestFast_noadjCov(arma::vec & t_GVec,
		     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2,
                     const PerMarkerCtx& ctx){

      arma::vec g1 = t_GVec.elem(t_indexForNonZero);
      arma::vec m_mu21 = m_mu2.elem(t_indexForNonZero);
      arma::vec m_res1 = m_res.elem(t_indexForNonZero);
    double S, var2;

    double var2_a = dot(m_mu21,pow(g1,2));
    double var2_b = dot(m_mu21, 2*2*t_altFreq*g1);
    double var2_c = arma::accu(m_mu2)*pow(2*t_altFreq, 2);
    var2 = var2_a - var2_b + var2_c;
    var2 = var2 *(m_tauvec[0]);

    S = dot(g1, m_res1)  - arma::accu(m_res)*(2*t_altFreq);
    S = S/m_tauvec[0];

    double var1 = var2 * ctx.varRatioVal;
    double stat = S*S/var1;
    if (var1 <= std::numeric_limits<double>::min()){
          t_pval = 1;
    }else{
      if(!std::isnan(stat) && std::isfinite(stat)){
          boost::math::chi_squared chisq_dist(1);
          t_pval = boost::math::cdf(complement(chisq_dist, stat));
      }else{
          t_pval = 1;
          stat = 0.0;
      }
    }
    char pValueBuf[100];

    if (t_pval != 0){
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    }else{
        // P1 fix: direct log(upper-tail) avoiding cdf underflow at p<1e-300
        double logp = log_chisq1_uppertail(stat);
        double log10p = logp/(log(10));
        int exponent = floor(log10p);
        double fraction = pow(10.0, log10p - exponent);
        if (fraction >= 9.95) {
          fraction = 1;
           exponent++;
        }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp;
        t_islogp = true;
    }
    std::string buffAsStdStr = pValueBuf;
    t_pval_str = buffAsStdStr;
    t_Beta = S/var1;
    t_seBeta = fabs(t_Beta) / sqrt(fabs(stat));
    t_Tstat = S;
    t_var1 = var1;
    t_var2 = var2;
}





// ============================================================
// Phase C: matrix-level first pass via single BLAS-3 GEMM block.
// For a block of B markers (G is N x B dense), compute Tstat / var1 / var2
// using the same math as scoreTestFast (covar-adjusted, non-sparse GRM)
// but vectorized over columns.
//
// Math (per column j):
//   tildeG_j = G_j - X * (X^TWX)^{-1} X^TW * G_j
//   S_j = tildeG_j^T (y - mu) / tau_0
//   var2_j = tildeG_j^T W tildeG_j   (binary/survival)
//          = tildeG_j^T tildeG_j     (quantitative, with tau_0 scale on X-part)
//   var1_j = var2_j * varRatio_j
// Vectorized:
//   Z = (X^TW) * G          [p x B]
//   B_X = X * (X^TWX)^{-1} * Z  -- but we use the cached identity
//                                  m_XXVX_inv * (X^TW) = X (X^TWX)^{-1} X^TW
//                                  expressed via Z and m_XVX_inv
//   In scoreTestFast the math uses:
//     Z = m_XVX_inv_XV.rows(idx).t() * g1  (= (X^TWX)^{-1} X^TW g restricted)
//     B = X * Z = m_XXVX_inv.rows(idx) * Z (= X (X^TWX)^{-1} X^TW g restricted)
//   Here for dense block over all N samples (no idx restriction since G is full
//   length, dense, with zeros explicit), we compute:
//     Z = m_XV * G            [p x B]  // X^TW * G
//     tildeG = G - m_XXVX_inv * (m_XVX_inv * Z)  -- BUT m_XXVX_inv already
//   absorbs (X^TWX)^{-1}: m_XXVX_inv == X * (X^TWX)^{-1}. And X^TW is m_XV.
//   So tildeG = G - m_XXVX_inv * (m_XV * G) -- one GEMM for Z, one for subtract.
//   This is exactly getadjG done for all columns at once.
//
// For var2 in binary/survival, we use the column-by-column formula
//   var2_j = Z_j^T * m_XVX * Z_j - sum_i mu2_i * B_{ij}^2 + sum_i mu2_i * tildeG_{ij}^2
// where B = m_XXVX_inv * Z (without the X subtraction, just the projected part).
// This matches scoreTestFast (which uses g1, X1, A1 restricted to nonzero idx,
// but algebraically is equivalent when G is full-length dense and we account
// for zero rows correctly).
//
// For quantitative:
//   var2_j = Z_j^T * m_XVX * Z_j * tau_0 + G_j^T G_j - 2 * G_j^T B_j
// (matches the scalar formula with g1 = G_j, dot products over full N since
// zero rows contribute 0).
//
// S (Tstat) for all columns:
//   S_j = tildeG_j^T m_res / tau_0
//        = ( G_j^T m_res - Z_j^T (X^T m_res) ) / tau_0
// We use the equivalent:
//   S_j = G_j^T m_res / tau_0 - Z_j^T * m_S_a
// Note scoreTestFast uses: S = S1 + S2; S1 = dot(res1, g1_tilde); the
// algebraic equivalent on full N is dot(m_res, tildeG_j) = dot(m_res, G_j) -
// dot(m_res, B_j) where B_j = m_XXVX_inv * Z_j. Then divide by tauvec[0].
void SAIGEClass::scoreTestFast_block(const arma::mat& G,
                                      const arma::vec& varRatioVec,
                                      const std::vector<bool>& validMask,
                                      arma::vec& Beta,
                                      arma::vec& seBeta,
                                      arma::vec& Tstat,
                                      arma::vec& var1,
                                      arma::vec& var2,
                                      arma::vec& StdStat,
                                      arma::vec& pvalNoadj,
                                      std::vector<bool>& pvalIsLog,
                                      std::vector<std::string>& pvalStr) const {
    const unsigned int B = G.n_cols;
    const double tau0 = m_tauvec[0];

    // Matrix layout notes (from null_model_engine score builder):
    //   m_X           : N x p   = covariate design
    //   m_XV          : p x N   = X^TW
    //   m_XXVX_inv    : N x p   = X (X^TWX)^{-1}
    //   m_XVX_inv_XV  : N x p   = X (X^TWX)^{-1} per row scaled by V
    //   m_XVX         : p x p   = X^TWX
    //
    // Per scalar scoreTestFast (over idx of nonzero g rows):
    //   Z   = m_XVX_inv_XV.rows(idx).t() * g1   = (X^TWX)^{-1} X^TW g   (p-vec)
    //   B   = X1 * Z                            = X * Z   restricted to idx
    //         (Z already absorbed the (X^TWX)^{-1} factor, so the projection
    //         from latent space back to data space is plain X * Z, NOT
    //         m_XXVX_inv * Z which would re-apply (X^TWX)^{-1}.)
    // Block analogue:
    //   Z   = m_XVX_inv_XV.t() * G                                      [p x B]
    //   BX  = m_X * Z                                                   [N x B]
    // Since G has explicit zeros for the "zero" rows (post-impute), the full-N
    // dot products coincide with the |idx|-restricted ones.
    arma::mat Z = m_XVX_inv_XV.t() * G;          // p x B
    arma::mat BX = m_X * Z;                       // N x B
    // tildeG only needed for var2 computation (column-wise squared sums).
    // tildeG = G - BX
    // We avoid materializing tildeG separately by computing the column squared
    // sums directly: but for clarity and simplicity allocate once.
    arma::mat tildeG = G - BX;                    // N x B

    // ---- S = (G^T m_res - Z^T X^T m_res) / tau0 -----------------------------
    // Use S = tildeG^T m_res / tau0   (matches scalar exactly)
    arma::vec S_vec = (tildeG.t() * m_res) / tau0;        // length B

    // ---- var2 ---------------------------------------------------------------
    arma::vec var2_vec(B, arma::fill::zeros);
    if (m_traitType == "binary" || m_traitType == "survival") {
        // var2_j = Z_j^T * m_XVX * Z_j  - sum_i mu2_i * B_{ij}^2  + sum_i mu2_i * tildeG_{ij}^2
        // ZtXVXZ diag: diag(Z^T m_XVX Z) -- compute via (m_XVX * Z) elem mul Z.
        arma::mat XVXZ = m_XVX * Z;                          // p x B
        arma::vec ZtXVXZ_diag = arma::sum(Z % XVXZ, 0).t();  // length B
        // Bmu2_j = sum_i mu2_i * B_{ij}^2
        arma::vec Bmu2 = (arma::square(BX).t()) * m_mu2;     // length B
        arma::vec g1tildemu2 = (arma::square(tildeG).t()) * m_mu2;
        var2_vec = ZtXVXZ_diag - Bmu2 + g1tildemu2;
    } else if (m_traitType == "quantitative") {
        // var2_j = ZtXVXZ_diag * tau0 + G_j^T G_j - 2 * G_j^T B_j
        arma::mat XVXZ = m_XVX * Z;
        arma::vec ZtXVXZ_diag = arma::sum(Z % XVXZ, 0).t();
        arma::vec GtG = arma::sum(arma::square(G), 0).t();
        arma::vec GtB = arma::sum(G % BX, 0).t();
        var2_vec = ZtXVXZ_diag * tau0 + GtG - 2.0 * GtB;
    }

    // ---- per-column finalize (cheap scalar work) ----------------------------
    Beta.set_size(B);
    seBeta.set_size(B);
    Tstat.set_size(B);
    var1.set_size(B);
    var2.set_size(B);
    StdStat.set_size(B);
    pvalNoadj.set_size(B);
    if (pvalIsLog.size() != B) pvalIsLog.assign(B, false);
    if (pvalStr.size() != B) pvalStr.assign(B, std::string());

    for (unsigned int j = 0; j < B; j++) {
        if (!validMask[j]) {
            Beta[j] = arma::datum::nan;
            seBeta[j] = arma::datum::nan;
            Tstat[j] = arma::datum::nan;
            var1[j] = arma::datum::nan;
            var2[j] = arma::datum::nan;
            StdStat[j] = arma::datum::nan;
            pvalNoadj[j] = 1.0;
            pvalIsLog[j] = false;
            pvalStr[j].clear();
            continue;
        }
        double S = S_vec[j];
        double v2 = var2_vec[j];
        double v1 = v2 * varRatioVec[j];
        double stat = S * S / v1;
        double pval;
        bool islogp = false;
        if (v1 <= std::numeric_limits<double>::min()) {
            pval = 1.0;
        } else if (!std::isnan(stat) && std::isfinite(stat)) {
            boost::math::chi_squared chisq_dist(1);
            pval = boost::math::cdf(complement(chisq_dist, stat));
        } else {
            pval = 1.0;
            stat = 0.0;
        }
        char pValueBuf[100];
        if (pval != 0) {
            sprintf(pValueBuf, "%.6E", pval);
            islogp = false;
        } else {
            double logp = log_chisq1_uppertail(stat);
            double log10p = logp / std::log(10);
            int exponent = (int)std::floor(log10p);
            double fraction = std::pow(10.0, log10p - exponent);
            if (fraction >= 9.95) { fraction = 1; exponent++; }
            sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
            pval = logp;
            islogp = true;
        }
        Beta[j] = S / v1;
        seBeta[j] = std::fabs(Beta[j]) / std::sqrt(std::fabs(stat));
        Tstat[j] = S;
        var1[j] = v1;
        var2[j] = v2;
        StdStat[j] = std::fabs(S) / std::sqrt(v1);
        pvalNoadj[j] = pval;
        pvalIsLog[j] = islogp;
        pvalStr[j] = std::string(pValueBuf);
    }
}

// ============================================================
// scoreTestFast_fused (Pillar 1) — single fused pass over carrier samples.
// Algebraic reduction of scoreTestFast (covar-adjusted, non-sparse GRM). No
// dense N-vector, no .elem()/.rows() temporaries. With A = m_XVX_inv_XV,
// X = m_X, carriers k, i = idx[k], g = G[i]:
//   Z[j]  = Sum_k A(i,j) g                     (= A1.t() * g1)
//   rg    = Sum_k res(i) g
//   binary/survival: gsq = Sum_k g^2 mu2(i) ;  W[j] = Sum_k g mu2(i) X(i,j)
//   quantitative   : gsq = Sum_k g^2         ;  W[j] = Sum_k g       X(i,j)
// After the rX and Zt(Sum mu2 XiXi')Z cancellations (see MEMORY_LATENCY doc):
//   S    = ( rg - m_S_a . Z ) / tau0
//   var2 = Z' m_XVX Z            + gsq - 2 Z.W     (binary/survival)
//        = Z' m_XVX Z * tau0     + gsq - 2 Z.W     (quantitative)
//   var1 = var2 * varRatioVal
// NOT bit-identical to scoreTestFast (hand loops + rearranged sums); matches
// to ~1e-10 rel. p-value formatting identical to scoreTestFast.
void SAIGEClass::scoreTestFast_fused(const arma::vec& t_GVec,
                                     const arma::uvec& t_indexForNonZero,
                                     double varRatioVal,
                                     double& t_Beta,
                                     double& t_seBeta,
                                     std::string& t_pval_str,
                                     double& t_pval,
                                     bool& t_islogp,
                                     double& t_Tstat,
                                     double& t_var1,
                                     double& t_var2,
                                     double& t_StdStat) const {
    const arma::uword nz = t_indexForNonZero.n_elem;
    const arma::uword p  = m_At.n_rows;                 // # covariates (small)
    const double tau0    = m_tauvec[0];
    const bool   isBin   = (m_traitType == "binary" || m_traitType == "survival");

    std::vector<double> Zacc(p, 0.0), Wacc(p, 0.0);
    double rg = 0.0, gsq = 0.0;

    const double*      gp   = t_GVec.memptr();
    const double*      resp = m_res.memptr();
    const double*      mu2p = m_mu2.memptr();
    const arma::uword* ip   = t_indexForNonZero.memptr();

    // ---- the only O(nz*p) work: fused carrier pass over contiguous columns ----
    for (arma::uword k = 0; k < nz; ++k) {
        const arma::uword i  = ip[k];
        const double      g  = gp[i];
        const double*     Ai = m_At.colptr(i);          // contiguous p doubles
        const double*     Xi = m_Xt.colptr(i);          // contiguous p doubles
        rg += resp[i] * g;
        if (isBin) {
            const double gmu = g * mu2p[i];
            gsq += g * gmu;                             // Sum g^2 mu2
            for (arma::uword j = 0; j < p; ++j) { Zacc[j] += Ai[j] * g; Wacc[j] += Xi[j] * gmu; }
        } else {
            gsq += g * g;                              // Sum g^2
            for (arma::uword j = 0; j < p; ++j) { Zacc[j] += Ai[j] * g; Wacc[j] += Xi[j] * g; }
        }
    }

    // ---- O(p)/O(p^2) finalize (tiny) via small arma ops ----
    arma::vec Zv(Zacc.data(), p, false, false);
    arma::vec Wv(Wacc.data(), p, false, false);
    const double ZtXVXZ = arma::as_scalar(Zv.t() * m_XVX * Zv);
    const double ZdotW  = arma::dot(Zv, Wv);
    const double SaZ    = arma::dot(m_S_a, Zv);

    const double var2 = isBin ? (ZtXVXZ        + gsq - 2.0 * ZdotW)
                              : (ZtXVXZ * tau0  + gsq - 2.0 * ZdotW);
    const double var1 = var2 * varRatioVal;
    const double S    = (rg - SaZ) / tau0;

    // ---- p-value: identical formatting to scoreTestFast ----
    double stat = S * S / var1;
    if (var1 <= std::numeric_limits<double>::min()) {
        t_pval = 1.0; stat = 0.0;
    } else if (!std::isnan(stat) && std::isfinite(stat)) {
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    } else {
        t_pval = 1.0; stat = 0.0;
    }
    char pValueBuf[100];
    if (t_pval != 0) {
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    } else {
        double logp = log_chisq1_uppertail(stat);
        double log10p = logp / std::log(10);
        int exponent = (int)std::floor(log10p);
        double fraction = std::pow(10.0, log10p - exponent);
        if (fraction >= 9.95) { fraction = 1; exponent++; }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp;
        t_islogp = true;
    }
    t_pval_str = pValueBuf;
    t_Beta    = S / var1;
    t_seBeta  = std::fabs(t_Beta) / std::sqrt(std::fabs(stat));
    t_Tstat   = S;
    t_var1    = var1;
    t_var2    = var2;
    t_StdStat = std::fabs(S) / std::sqrt(var1);
}

// ============================================================
// scoreTestFast_noadjCov_fused (Pillar 1, no-covariate-adjustment variant).
// Fused carrier pass for the isnoadjCov score test. Over carriers k, i=idx[k],
// g=G[i]:  g2mu = Sum mu2(i) g^2 ; gmu = Sum mu2(i) g ; rg = Sum res(i) g.
//   var2 = ( g2mu - 4*altFreq*gmu + m_sum_mu2*(2*altFreq)^2 ) * tau0
//   S    = ( rg - m_sum_res*2*altFreq ) / tau0 ;  var1 = var2 * varRatioVal
// No X/A gathers. m_sum_mu2/m_sum_res are precomputed (== scalar's per-call
// arma::accu). Matches scoreTestFast_noadjCov to ~1e-13.
void SAIGEClass::scoreTestFast_noadjCov_fused(const arma::vec& t_GVec,
                                              const arma::uvec& t_indexForNonZero,
                                              double t_altFreq,
                                              double varRatioVal,
                                              double& t_Beta,
                                              double& t_seBeta,
                                              std::string& t_pval_str,
                                              double& t_pval,
                                              bool& t_islogp,
                                              double& t_Tstat,
                                              double& t_var1,
                                              double& t_var2,
                                              double& t_StdStat) const {
    const arma::uword nz = t_indexForNonZero.n_elem;
    const double tau0    = m_tauvec[0];

    double g2mu = 0.0, gmu = 0.0, rg = 0.0;
    const double*      gp   = t_GVec.memptr();
    const double*      resp = m_res.memptr();
    const double*      mu2p = m_mu2.memptr();
    const arma::uword* ip   = t_indexForNonZero.memptr();

    for (arma::uword k = 0; k < nz; ++k) {
        const arma::uword i  = ip[k];
        const double      g  = gp[i];
        const double      m2 = mu2p[i];
        g2mu += m2 * g * g;
        gmu  += m2 * g;
        rg   += resp[i] * g;
    }

    const double af2  = 2.0 * t_altFreq;
    const double var2 = (g2mu - 2.0 * af2 * gmu + m_sum_mu2 * af2 * af2) * tau0;
    const double var1 = var2 * varRatioVal;
    const double S    = (rg - m_sum_res * af2) / tau0;

    double stat = S * S / var1;
    if (var1 <= std::numeric_limits<double>::min()) {
        t_pval = 1.0; stat = 0.0;
    } else if (!std::isnan(stat) && std::isfinite(stat)) {
        boost::math::chi_squared chisq_dist(1);
        t_pval = boost::math::cdf(complement(chisq_dist, stat));
    } else {
        t_pval = 1.0; stat = 0.0;
    }
    char pValueBuf[100];
    if (t_pval != 0) {
        sprintf(pValueBuf, "%.6E", t_pval);
        t_islogp = false;
    } else {
        double logp = log_chisq1_uppertail(stat);
        double log10p = logp / std::log(10);
        int exponent = (int)std::floor(log10p);
        double fraction = std::pow(10.0, log10p - exponent);
        if (fraction >= 9.95) { fraction = 1; exponent++; }
        sprintf(pValueBuf, "%.1fE%d", fraction, exponent);
        t_pval = logp; t_islogp = true;
    }
    t_pval_str = pValueBuf;
    t_Beta    = S / var1;
    t_seBeta  = std::fabs(t_Beta) / std::sqrt(std::fabs(stat));
    t_Tstat   = S;
    t_var1    = var1;
    t_var2    = var2;
    t_StdStat = std::fabs(S) / std::sqrt(var1);
}


void SAIGEClass::getadjG(arma::vec & t_GVec, arma::vec & g){
   g = m_XV * t_GVec;
    g = t_GVec - m_XXVX_inv * g;
}


void SAIGEClass::getadjGFast(arma::vec & t_GVec, arma::vec & g, arma::uvec & iIndex)
{
  // To increase computational efficiency when lots of GVec elements are 0
 arma::vec m_XVG(m_p, arma::fill::zeros);
  for(unsigned int i = 0; i < iIndex.n_elem; i++){
      m_XVG += m_XV.col(iIndex(i)) * t_GVec(iIndex(i));
  }
  g = t_GVec - m_XXVX_inv * m_XVG;
}


void SAIGEClass::get_mu(arma::vec & t_mu){
    t_mu = m_mu;
}

void SAIGEClass::getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices){
     t_case_indices = m_case_indices;
     t_ctrl_indices = m_ctrl_indices;
  }


// revised
// need to add sparse Sigma version
// This function only uses variance ratio and does not use sparse GRM
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta,
                               double& t_seBeta,
			       std::string& t_pval,
			       std::string& t_pval_noSPA,
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1,
			       bool & t_isSPAConverge,
			       arma::vec & t_gtilde,
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
			       std::string& t_pval_c,
                               std::string& t_pval_noSPA_c,
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde,
				bool & t_isFirth,
				bool & t_isFirthConverge,
				bool t_isER,
				bool t_isnoadjCov,
				bool t_isSparseGRM)
{
    PerMarkerCtx ctx{m_flagSparseGRM_cur, m_isnoadjCov_cur, m_varRatioVal};
    getMarkerPval(t_GVec, iIndex, iIndexComVec, t_Beta, t_seBeta, t_pval,
                  t_pval_noSPA, t_altFreq, t_Tstat, t_gy, t_var1,
                  t_isSPAConverge, t_gtilde, is_gtilde, is_region, t_P2Vec,
                  t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c,
                  t_pval_noSPA_c, t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde,
                  t_isFirth, t_isFirthConverge, t_isER, t_isnoadjCov,
                  t_isSparseGRM, ctx);
}

// Ctx overload (Phase A): same body as the member-state version, but reads
// ctx.flagSparseGRM_cur / ctx.varRatioVal instead of class members. Nested
// scoreTest*/scoreTestFast* calls are routed to their ctx overloads so they
// too use ctx rather than mutable class state.
void SAIGEClass::getMarkerPval(arma::vec & t_GVec,
			       arma::uvec & iIndex,
			       arma::uvec & iIndexComVec,
                               double& t_Beta,
                               double& t_seBeta,
			       std::string& t_pval,
			       std::string& t_pval_noSPA,
                               double t_altFreq,
                               double& t_Tstat,
			       double& t_gy,
			       double& t_var1,
			       bool & t_isSPAConverge,
			       arma::vec & t_gtilde,
			       bool & is_gtilde,
			       bool  is_region,
                               arma::vec & t_P2Vec,
			       bool t_isCondition,
			       double& t_Beta_c,
                           	double& t_seBeta_c,
			       std::string& t_pval_c,
                               std::string& t_pval_noSPA_c,
                           	double& t_Tstat_c,
                           	double& t_varT_c,
			   	arma::rowvec & t_G1tilde_P_G2tilde,
				bool & t_isFirth,
				bool & t_isFirthConverge,
				bool t_isER,
				bool t_isnoadjCov,
				bool t_isSparseGRM,
				const PerMarkerCtx& ctx)
{



  t_isFirth = false;
  std::string t_pval_str;
  double t_var2, t_SPApval;
  bool isScoreFast = true;
  if(ctx.flagSparseGRM_cur){
    isScoreFast = false;
  }


  double pval_noadj, pval, t_qval_Firth;
  bool ispvallog;

if(!ctx.flagSparseGRM_cur && t_isnoadjCov){
	is_gtilde = false;
        isScoreFast = true;
        // Pillar 1 fused-kernel wiring (no-covariate-adjustment variant).
        if (g_fusedMode == 2) {
            double sd_f;
            scoreTestFast_noadjCov_fused(t_GVec, iIndex, t_altFreq, ctx.varRatioVal,
                                         t_Beta, t_seBeta, t_pval_noSPA, pval_noadj,
                                         ispvallog, t_Tstat, t_var1, t_var2, sd_f);
        } else if (g_fusedMode == 1) {
            scoreTestFast_noadjCov(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq,t_Tstat, t_var1, t_var2, ctx);
            double fB, fsB, fpv, fT, fv1, fv2, fsd; bool flog; std::string fstr;
            scoreTestFast_noadjCov_fused(t_GVec, iIndex, t_altFreq, ctx.varRatioVal,
                                         fB, fsB, fstr, fpv, flog, fT, fv1, fv2, fsd);
            double relT = std::fabs(fT - t_Tstat) / (std::fabs(t_Tstat) + 1e-300);
            double relV = std::fabs(fv1 - t_var1) / (std::fabs(t_var1) + 1e-300);
            #pragma omp critical(fusedAB)
            {
                g_fusedNCompared++;
                if (relT > g_fusedMaxRelTstat) g_fusedMaxRelTstat = relT;
                if (relV > g_fusedMaxRelVar1)  g_fusedMaxRelVar1  = relV;
                if (relT > 1e-6 || relV > 1e-6) g_fusedNExceed++;
            }
        } else {
            scoreTestFast_noadjCov(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq,t_Tstat, t_var1, t_var2, ctx);
        }

}else if(ctx.flagSparseGRM_cur){
	is_gtilde = true;
        isScoreFast = false;
        scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex, ctx);
}else{
	is_gtilde = false;
        isScoreFast = true;
        // Pillar 1 fused-kernel wiring. This is exactly the covar-adjusted,
        // non-sparse-GRM scoreTestFast path the fused kernel replaces; SPA/ER/
        // Firth below rebuild t_gtilde on demand (is_gtilde stays false), so the
        // swap is transparent to everything downstream.
        if (g_fusedMode == 2) {
            double sd_f;
            scoreTestFast_fused(t_GVec, iIndex, ctx.varRatioVal, t_Beta, t_seBeta,
                                t_pval_noSPA, pval_noadj, ispvallog,
                                t_Tstat, t_var1, t_var2, sd_f);
        } else if (g_fusedMode == 1) {
            // A/B: scalar is authoritative (drives output); fused is compared.
            scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, ctx);
            double fB, fsB, fpv, fT, fv1, fv2, fsd; bool flog; std::string fstr;
            scoreTestFast_fused(t_GVec, iIndex, ctx.varRatioVal, fB, fsB, fstr,
                                fpv, flog, fT, fv1, fv2, fsd);
            double relT = std::fabs(fT - t_Tstat) / (std::fabs(t_Tstat) + 1e-300);
            double relV = std::fabs(fv1 - t_var1) / (std::fabs(t_var1) + 1e-300);
            #pragma omp critical(fusedAB)
            {
                g_fusedNCompared++;
                if (relT > g_fusedMaxRelTstat) g_fusedMaxRelTstat = relT;
                if (relV > g_fusedMaxRelVar1)  g_fusedMaxRelVar1  = relV;
                if (relT > 1e-6 || relV > 1e-6) g_fusedNExceed++;
            }
        } else {
            scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, ctx);
        }
}


  double StdStat = std::abs(t_Tstat) / sqrt(t_var1);

  t_isSPAConverge = false;

  double q, qinv, m1, NAmu, NAsigma, tol1, p_iIndexComVecSize;


  unsigned int iIndexComVecSize = iIndexComVec.n_elem;
  unsigned int iIndexSize = iIndex.n_elem;
  arma::vec gNB(iIndexSize, arma::fill::none);
  arma::vec gNA(iIndexComVecSize, arma::fill::none);
  arma::vec muNB(iIndexSize, arma::fill::none);
  arma::vec muNA(iIndexComVecSize, arma::fill::none);


  double gmuNB;


if((StdStat > m_SPA_Cutoff || std::isnan(StdStat)) && m_traitType != "quantitative" && t_isER){
	t_isER = true;
}else{
	t_isER = false;
}


if(!t_isER){


  if(!std::isnan(StdStat) && (StdStat > m_SPA_Cutoff) && m_traitType != "quantitative"){

       if(!is_gtilde){
          t_gtilde.resize(m_n);
          getadjGFast(t_GVec, t_gtilde, iIndex);
	  is_gtilde = true;
       }
        p_iIndexComVecSize = double(iIndexComVecSize)/m_n;
   	m1 = dot(m_mu, t_gtilde);

	if(p_iIndexComVecSize >= 0.5){

	gNB = t_gtilde(iIndex);
	gNA = t_gtilde(iIndexComVec);
   	muNB = m_mu(iIndex);
   	muNA = m_mu(iIndexComVec);

  	gmuNB = dot(gNB,muNB);
   	NAmu= m1-gmuNB;

   }

   	if(m_traitType == "binary"){
                q = t_Tstat/sqrt(t_var1/t_var2) + m1;

                if((q-m1) > 0){
                        qinv = -1 * std::abs(q-m1) + m1;
                }else if ((q-m1) == 0){
                        qinv =  m1;
                }else{
                        qinv = std::abs(q-m1) + m1;
                }
		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % (1-muNB) % arma::pow(gNB,2));
		}
        }else if(m_traitType == "survival"){
                q = t_Tstat/sqrt(t_var1/t_var2);
                qinv = -q;
  		if(p_iIndexComVecSize >= 0.5){
           		NAsigma = t_var2 - arma::sum(muNB % arma::pow(gNB,2));
		}
  }
	double tol0 = std::numeric_limits<double>::epsilon();
	tol1 = std::pow(tol0, 0.25);
	if(p_iIndexComVecSize >= 0.5 && !ctx.flagSparseGRM_cur){
        	SPA_fast(m_mu, t_gtilde, q, qinv, pval_noadj, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, t_SPApval, t_isSPAConverge);

	}else{
		SPA(m_mu, t_gtilde, q, qinv, pval_noadj, tol1, ispvallog, m_traitType, t_SPApval, t_isSPAConverge);
	}

    boost::math::normal ns;
    double t_qval;



    if(t_isSPAConverge){
        try {
          if(ispvallog){
              // R::qnorm(t_SPApval/2, 0, 1, false, true) -- log-scale quantile
              // t_SPApval is already log(p), so t_SPApval/2 = log(p)/2
              // This is log(sqrt(p)), and we need quantile of complement at that scale
              // TODO: Verify this log-scale qnorm conversion matches R exactly
              // For now, use exp to convert back from log scale, then take quantile
              double half_pval_log = t_SPApval / 2.0;
              double half_pval = std::exp(half_pval_log);
              if(half_pval > 0 && half_pval < 1){
                  t_qval = boost::math::quantile(complement(ns, half_pval));
              } else {
                  t_qval = std::numeric_limits<double>::infinity();
                  t_isSPAConverge = false;
              }
          }else{
              t_qval = boost::math::quantile(complement(ns, t_SPApval/2));
          }
          t_qval = fabs(t_qval);
          t_seBeta = fabs(t_Beta)/t_qval;
        }catch (const std::overflow_error&) {
          t_qval = std::numeric_limits<double>::infinity();
	  t_isSPAConverge = false;
        }
    }


    if(!ispvallog && t_SPApval == 0){
          t_isSPAConverge = false;
    }


  }
   char pValueBuf_SPA[100];
   if(m_traitType!="quantitative"){
        if(t_isSPAConverge){
		if(!ispvallog){
		  sprintf(pValueBuf_SPA, "%.6E", t_SPApval);
		}else{
        	  double t_SPApval_log10 = t_SPApval/(log(10));
        	  int exponent = floor(t_SPApval_log10);
        	  double fraction = pow(10.0, t_SPApval_log10 - exponent);
        	  if (fraction >= 9.95) {
          	    fraction = 1;
           	    exponent++;
         	  }
        	  sprintf(pValueBuf_SPA, "%.1fE%d", fraction, exponent);
		}
		std::string buffAsStdStr_SPA = pValueBuf_SPA;
	        t_pval = buffAsStdStr_SPA;
		pval = t_SPApval;
        }else{
                t_pval = t_pval_noSPA;
		pval = pval_noadj;
        }
     if(m_traitType=="binary"){
	if(!ispvallog){
		if(m_is_Firth_beta && pval <= m_pCutoffforFirth){
			t_isFirth = true;
                        boost::math::normal ns_firth;
			t_qval_Firth = boost::math::quantile(complement(ns_firth, pval/2));
		}

	}else{
		if(m_is_Firth_beta && pval <= std::log(m_pCutoffforFirth)){
			t_isFirth = true;
                        // R::qnorm(pval/2, 0, 1, false, true) -- log-scale
                        double half_pval = std::exp(pval / 2.0);
                        boost::math::normal ns_firth;
                        if(half_pval > 0 && half_pval < 1){
                            t_qval_Firth = boost::math::quantile(complement(ns_firth, half_pval));
                        } else {
                            t_qval_Firth = std::numeric_limits<double>::infinity();
                        }
		}
	}
     }

   }else{
        t_pval = t_pval_noSPA;
	pval = pval_noadj;
   }

}else{ //if(!t_isER){

    arma::mat Z_er(t_GVec.n_elem, 1);
    Z_er.col(0) = t_GVec;
    arma::vec res_er = m_res;
    arma::vec pi1_er = m_mu;
    arma::mat resout_er;
    if(m_resout.n_elem > 0){
        resout_er = arma::mat(m_resout.n_elem, 1);
        resout_er.col(0) = m_resout;
    }
    // W1-4: seed the ER resampling RNG from the marker's own stream id so the
    // result does not depend on which thread the marker landed on.
    ER::SL_set_stream(ctx.erSeedStream);
    double pval_ER = ER::SKATExactBin_Work(Z_er, res_er, pi1_er, m_n_case, iIndex, iIndexComVec, resout_er, 2e+6, 1e+4, 1e-6, 1);
    // P2 fix (2026-05-09): ER can return NaN/inf on degenerate ultra-rare configs
    // (e.g. all carriers in same case/control group → variance 0 in resampling).
    // Without this guard, the NaN propagates into boost::math::quantile which
    // throws an uncaught domain_error (the existing catch only handles
    // overflow_error).
    if (!std::isfinite(pval_ER) || pval_ER < 0.0 || pval_ER > 1.0) {
        pval_ER = 1.0;  // treat as non-significant; let caller re-route to score-test
    }
    char pValueBuf_ER[100];
    sprintf(pValueBuf_ER, "%.6E", pval_ER);
    std::string buffAsStdStr_ER = pValueBuf_ER;
    t_pval = pValueBuf_ER;

    pval = pval_ER;
    boost::math::normal ns;
    double t_qval_ER;
    try{
      t_qval_ER = boost::math::quantile(ns, pval_ER/2);
      t_qval_ER = fabs(t_qval_ER);
      if (t_qval_ER == 0.0 || !std::isfinite(t_qval_ER)) {
          t_seBeta = 0;
      } else {
          t_seBeta = fabs(t_Beta)/t_qval_ER;
      }
      t_isSPAConverge = true;
    }catch (const std::exception&) {  // widen catch beyond overflow_error
      t_qval_ER = std::numeric_limits<double>::infinity();
      t_seBeta = 0;
    }

    if(m_is_Firth_beta && pval <= m_pCutoffforFirth){
	t_isFirth = true;
	t_qval_Firth = t_qval_ER;
    }
}
   if(t_isFirth && !g_firthDefer){
	if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	// A3: persistent N x 2 design matrix (col0 = intercept ones, col1 = gtilde).
	// Constant N x 2 => set_size reuses storage after first call (no fault).
	thread_local arma::mat x;
	x.set_size(t_GVec.n_elem, 2);
	x.col(0).ones();
	x.col(1) = t_gtilde;
	arma::vec init(2, arma::fill::zeros);
	g_firthFitCalls.fetch_add(1, std::memory_order_relaxed);
	fast_logistf_fit_simple(x, m_y, m_offset, true, init, 50, 15, 15, 1e-5, 1e-5, 1e-5, t_Beta ,t_seBeta, t_isFirthConverge);
	//back calculates se based on beta from firth adjustion and the p-value that accounts for case-control imbalance
	t_seBeta = fabs(t_Beta)/fabs(t_qval_Firth);
   }


   //condition
   if(t_isCondition){
	if(!is_gtilde){
        	getadjGFast(t_GVec, t_gtilde, iIndex);
        	is_gtilde = true;
        }
        t_G1tilde_P_G2tilde = sqrt(ctx.varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
        arma::vec t_Tstat_ctemp =  t_G1tilde_P_G2tilde * m_VarInvMat_cond * m_Tstat_cond;
	arma::mat tempgP2 = t_gtilde.t() * m_P2Mat_cond;

    	t_Tstat_c = t_Tstat - t_Tstat_ctemp(0);
    	arma::vec t_varT_ctemp = t_G1tilde_P_G2tilde * m_VarInvMat_cond * (t_G1tilde_P_G2tilde.t());
    	t_varT_c = t_var1 - t_varT_ctemp(0);

    double S_c = t_Tstat_c;

    double stat_c = S_c*S_c/t_varT_c;

    double pval_noSPA_c;

     if (t_varT_c <= std::numeric_limits<double>::min()){
        pval_noSPA_c = 1;
	stat_c = 0;
     }else{
       if(!std::isnan(stat_c) && std::isfinite(stat_c)){
        boost::math::chi_squared chisq_dist(1);
        pval_noSPA_c = boost::math::cdf(complement(chisq_dist, stat_c));
       }else{
        pval_noSPA_c = 1;
	stat_c = 0;
       }
     }

    char pValueBuf_c[100];

    if (pval_noSPA_c != 0){
        sprintf(pValueBuf_c, "%.6E", pval_noSPA_c);
        ispvallog = false;
    }else {
        // P1 fix: direct log(upper-tail) avoiding cdf underflow at p<1e-300
	double logp_c = log_chisq1_uppertail(stat_c);
	double log10p_c = logp_c/(log(10));
	int exponent_c = floor(log10p_c);
        double fraction_c = pow(10.0, log10p_c - exponent_c);
        if (fraction_c >= 9.95) {
          fraction_c = 1;
           exponent_c++;
         }
        sprintf(pValueBuf_c, "%.1fE%d", fraction_c, exponent_c);
	ispvallog = true;
	pval_noSPA_c = logp_c;
    }
    std::string buffAsStdStr_c = pValueBuf_c;
    t_pval_noSPA_c = buffAsStdStr_c;
    t_Beta_c = S_c/t_varT_c;
    t_seBeta_c = fabs(t_Beta_c) / sqrt(stat_c);
    t_Tstat_c = S_c;


    bool t_isSPAConverge_c;
    if(m_traitType != "quantitative" && stat_c > std::pow(m_SPA_Cutoff,2)){
	double q_c, qinv_c, pval_noadj_c, SPApval_c;
	if(m_traitType == "binary"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2) + m1;
                if((q_c-m1) > 0){
                        qinv_c = -1 * std::abs(q_c-m1) + m1;
                }else if ((q_c-m1) == 0){
                        qinv_c =  m1;
                }else{
                        qinv_c = std::abs(q_c-m1) + m1;
                }
        }else if(m_traitType == "survival"){
                q_c = t_Tstat_c/sqrt(t_varT_c/t_var2);
                qinv = -q_c;
        }


        if(p_iIndexComVecSize >= 0.5 && !ctx.flagSparseGRM_cur){
                SPA_fast(m_mu, t_gtilde, q_c, qinv_c, pval_noSPA_c, ispvallog, gNA, gNB, muNA, muNB, NAmu, NAsigma, tol1, m_traitType, SPApval_c, t_isSPAConverge_c);
        }else{
                SPA(m_mu, t_gtilde, q_c, qinv_c, pval_noSPA_c, tol1, ispvallog, m_traitType, SPApval_c, t_isSPAConverge_c);
        }


        boost::math::normal ns;
        double t_qval_c;
     if(t_isSPAConverge_c){
        try {
          if(!ispvallog){
           t_qval_c = boost::math::quantile(ns, SPApval_c/2);
          }else{
              // R::qnorm(SPApval_c/2, 0, 1, false, true) -- log-scale
              double half_pval = std::exp(SPApval_c / 2.0);
              if(half_pval > 0 && half_pval < 1){
                  t_qval_c = boost::math::quantile(complement(ns, half_pval));
              } else {
                  t_qval_c = std::numeric_limits<double>::infinity();
                  t_isSPAConverge_c = false;
              }
          }
           t_qval_c = fabs(t_qval_c);
           t_seBeta_c = fabs(t_Beta_c)/t_qval_c;
        }catch (const std::overflow_error&) {
          t_qval_c = std::numeric_limits<double>::infinity();
	  t_isSPAConverge_c = false;
        }
      }

	if(!ispvallog && SPApval_c == 0){
		t_isSPAConverge_c = false;
	}


      char pValueBuf_SPA_c[100];
        if(t_isSPAConverge_c){
                if(!ispvallog){
                  sprintf(pValueBuf_SPA_c, "%.6E", SPApval_c/2);
                }else{
                  double SPApval_c_log10 = SPApval_c/(log(10));
                  int exponent = floor(SPApval_c_log10);
                  double fraction = pow(10.0, SPApval_c_log10 - exponent);
                  if (fraction >= 9.95) {
                    fraction = 1;
                    exponent++;
                  }
                  sprintf(pValueBuf_SPA_c, "%.1fE%d", fraction, exponent);
                }
                std::string buffAsStdStr_SPA_c = pValueBuf_SPA_c;
                t_pval_c = buffAsStdStr_SPA_c;
        }else{
                t_pval_c = t_pval_noSPA_c;
        }

    }else{
	t_isSPAConverge_c = false;
    	t_pval_c = t_pval_noSPA_c;
    }
 }


    gNA.clear();
    gNB.clear();
    muNA.clear();
    gNB.clear();


    if(is_region && !is_gtilde){
	getadjGFast(t_GVec, t_gtilde, iIndex);
	is_gtilde = true;
    }

    if(is_region && isScoreFast){

      t_gy = dot(t_gtilde, m_y);
      if(!ctx.flagSparseGRM_cur){
        t_P2Vec = t_gtilde % m_mu2 *m_tauvec[0];
      }else{
	t_P2Vec = getPCG1ofSigmaAndGtilde(t_gtilde, 100, 0.02);
      }
    }
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
	m_varRatio = m_varRatio_sparse;
    }else{
	m_varRatio = m_varRatio_null;
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
		m_varRatioVal = m_varRatio(i);
		hasVarRatio = true;
	}
    }

    if(!hasVarRatio){
	if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
		m_varRatioVal = m_varRatio(0);
		hasVarRatio = true;
	}
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
		m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

    return(hasVarRatio);
}


bool SAIGEClass::assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj){
    bool hasVarRatio = false;
    arma::vec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        if(!isnoXadj){
            m_varRatio = m_varRatio_null;
        }else{
            m_varRatio = m_varRatio_null_noXadj;
        }
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++)
    {
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
                m_varRatioVal = m_varRatio(i);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
                m_varRatioVal = m_varRatio(0);
                hasVarRatio = true;
        }
    }

    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
                m_varRatioVal = m_varRatio.back();
                hasVarRatio = true;
        }
    }

    return(hasVarRatio);
}





void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR){
    arma::vec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        m_varRatio = m_varRatio_null;
    }
    m_varRatioVal = m_varRatio(0);
}


void SAIGEClass::assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj){
    arma::rowvec m_varRatio;
    if(issparseforVR){
        m_varRatio = m_varRatio_sparse;
    }else{
        if(isnoXadj){
            m_varRatio = m_varRatio_null_noXadj;
        }else{
            m_varRatio = m_varRatio_null;
        }
    }
    m_varRatioVal = m_varRatio(0);
}



void SAIGEClass::assignSingleVarianceRatio_withinput(double t_varRatioVal){
        m_varRatioVal = t_varRatioVal;
}


// ============================================================
// Pure computeVarianceRatio* helpers (Phase A of step2 parallelism plan).
// These return the value that the corresponding assignVarianceRatio* would
// have written into m_varRatioVal, without mutating any class member. They
// let per-marker call-sites build a PerMarkerCtx without touching shared
// state, enabling thread-safe parallel marker processing in Phase B.
// ============================================================

double SAIGEClass::computeVarianceRatio(double MAC, bool issparseforVR, bool& hasVarRatio) const {
    hasVarRatio = false;
    double out = 1.0;
    arma::vec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        varRatio = m_varRatio_null;
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++){
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
            out = varRatio(i);
            hasVarRatio = true;
        }
    }
    if(!hasVarRatio){
        if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
            out = varRatio(0);
            hasVarRatio = true;
        }
    }
    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
            out = varRatio.back();
            hasVarRatio = true;
        }
    }
    return out;
}

double SAIGEClass::computeVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj, bool& hasVarRatio) const {
    hasVarRatio = false;
    double out = 1.0;
    arma::vec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        if(!isnoXadj){
            varRatio = m_varRatio_null;
        }else{
            varRatio = m_varRatio_null_noXadj;
        }
    }
    for(unsigned int i = 0; i < m_cateVarRatioMaxMACVecInclude.n_elem; i++){
        if(MAC <= m_cateVarRatioMaxMACVecInclude(i) && MAC > m_cateVarRatioMinMACVecExclude(i)){
            out = varRatio(i);
            hasVarRatio = true;
        }
    }
    if(!hasVarRatio){
        if(MAC <= m_cateVarRatioMinMACVecExclude(0)){
            out = varRatio(0);
            hasVarRatio = true;
        }
    }
    if(!hasVarRatio){
        if(MAC > m_cateVarRatioMaxMACVecInclude.back()){
            out = varRatio.back();
            hasVarRatio = true;
        }
    }
    return out;
}

double SAIGEClass::computeSingleVarianceRatio(bool issparseforVR) const {
    arma::vec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        varRatio = m_varRatio_null;
    }
    return varRatio(0);
}

double SAIGEClass::computeSingleVarianceRatio(bool issparseforVR, bool isnoXadj) const {
    arma::rowvec varRatio;
    if(issparseforVR){
        varRatio = m_varRatio_sparse;
    }else{
        if(isnoXadj){
            varRatio = m_varRatio_null_noXadj;
        }else{
            varRatio = m_varRatio_null;
        }
    }
    return varRatio(0);
}





void SAIGEClass::assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      std::vector<std::string> & t_p_cond
      ){
	m_P2Mat_cond = t_P2Mat_cond;
	m_VarInvMat_cond = t_VarInvMat_cond;
	m_VarMat_cond = t_VarMat_cond;
	m_Tstat_cond = t_Tstat_cond;
	m_MAF_cond = t_MAF_cond;
	m_qsum_cond = t_qsum_cond;
	m_gsum_cond = t_gsum_cond;
	m_G2_Weight_cond = t_G2_Weight_cond;
	m_p_cond = t_p_cond;
}

void SAIGEClass::assignConditionFactors_scalefactor(
	arma::vec & t_scalefactor_G2_cond
		){
	m_scalefactor_G2_cond = t_scalefactor_G2_cond;
	arma::mat scalefactor_G2_cond_Mat = arma::diagmat(arma::sqrt(m_scalefactor_G2_cond));
	arma::mat weightMat_G2_G2 = m_G2_Weight_cond * m_G2_Weight_cond.t();
	arma::mat VarMat_cond_scaled = scalefactor_G2_cond_Mat * m_VarMat_cond * scalefactor_G2_cond_Mat;
	arma::mat VarMat_cond_scaled_weighted = VarMat_cond_scaled % weightMat_G2_G2;
	m_VarInvMat_cond_scaled_weighted = arma::pinv(VarMat_cond_scaled_weighted);
}

void SAIGEClass::extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv){
	t_XV = m_XV;
	t_XXVX_inv = m_XXVX_inv;
}



void SAIGEClass::fast_logistf_fit_simple(arma::mat & x,
                arma::vec & y,
                arma::vec & offset,
                bool firth,
        arma::vec init,
        int maxit,
        int maxstep,
        int maxhs,
        double lconv,
        double gconv,
        double xconv,
        double & beta_G,
        double & sebeta_G,
	bool & isfirthconverge){
  isfirthconverge = false;
  int n = x.n_rows;
  int k = x.n_cols;
  arma::vec beta = init;
  int iter = 0;
  // A2 (2026-07-15): thread_local persistent workspace for the N-sized and N x k
  // objects allocated every Newton iteration. Firth runs on the FULL sample
  // (n == N constant across all markers/iterations), so set_size() reuses
  // storage after the first call -> no per-iteration mmap/fault (the ~54% minor-
  // fault source in the significant region). Math and evaluation order are
  // unchanged -> bit-identical. XW2/XX_XW2 are fully overwritten column-by-column
  // so set_size (no zero-fill) is equivalent to the old fill::zeros and also
  // drops the dscal zeroing faults.
  thread_local arma::vec ws_pi_0, ws_pi, ws_wpi, ws_W2, ws_h, ws_ypih, ws_xcol;
  thread_local arma::mat ws_XW2, ws_Q, ws_XX_XW2;
  ws_pi_0 = -x * beta - offset;
  ws_pi_0 = arma::exp(ws_pi_0) + 1;
  ws_pi = 1/ws_pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
        beta_old = beta;
        ws_wpi = ws_pi % (1 - ws_pi);
        ws_W2 = arma::sqrt(ws_wpi);
        ws_XW2.set_size(n, k);
        for(int j = 0; j < k; j++){
                ws_XW2.col(j) = x.col(j) % ws_W2;
        }

        arma::mat R;
        arma::qr_econ(ws_Q, R, ws_XW2);
        ws_h = ws_Q % ws_Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        if(firth){
                ws_ypih = (y - ws_pi) + (ws_h % (0.5 - ws_pi));
        }else{
                ws_ypih = (y - ws_pi);
        }
        U_star = x.t() * ws_ypih;

        ws_XX_XW2.set_size(n, k);
        for(int j = 0; j < k; j++){
                ws_xcol = x.col(j);
                ws_XX_XW2.col(j) = ws_xcol % ws_W2;
        }
        arma::mat XX_Fisher = ws_XX_XW2.t() * (ws_XX_XW2);
        bool isinv = arma::inv_sympd (XX_covs, XX_Fisher);

	if(!isinv){
                break;
        }
        arma::vec delta = XX_covs * U_star;

        double mx = arma::max(arma::abs(delta))/maxstep;
        if(mx > 1){
                delta = delta/mx;
        }
        evals = evals + 1;
        iter = iter + 1;
        beta = beta + delta;
        ws_pi_0 = -x * beta - offset;
        ws_pi_0 = arma::exp(ws_pi_0) + 1;
        ws_pi = 1/ws_pi_0;
        if((iter == maxit) || ( (arma::max(arma::abs(delta)) <= xconv) & (abs(U_star).is_zero(gconv)))){
		isfirthconverge = true;
                break;
        }
  }
        arma::mat var;
        if(XX_covs.has_nan()){
                var = XX_covs;
                beta_G = arma::datum::nan;
                sebeta_G = arma::datum::nan;
        }else{
                beta_G = beta(1);
                sebeta_G = sqrt(XX_covs(1,1));
        }
}

void SAIGEClass::set_flagSparseGRM_cur(bool t_flagSparseGRM_cur){
	m_flagSparseGRM_cur = t_flagSparseGRM_cur;
}


void SAIGEClass::set_isnoadjCov_cur(bool t_isnoadjCov_cur){
        m_isnoadjCov_cur = t_isnoadjCov_cur;
}

arma::vec SAIGEClass::getPCG1ofSigmaAndGtilde(arma::vec& bVec, int maxiterPCG, double tolPCG) {
    int Nnomissing = m_spSigmaMat.n_rows;
    arma::vec xVec(Nnomissing, arma::fill::zeros); // Initialize xVec to zeros
    arma::vec rVec = bVec; // Residual vector
    arma::vec zVec(Nnomissing);
    arma::vec minvVec = 1.0 / m_diagSigma; // Convert diagonal view to dense vector
    zVec = minvVec % rVec; // Apply preconditioner
    double sumr2 = arma::dot(rVec, rVec); // Initial residual norm
    arma::vec pVec = zVec; // Search direction

    int iter = 0;
    while (sumr2 > tolPCG && iter < maxiterPCG) {
        iter++;
        arma::vec ApVec = m_spSigmaMat * pVec; // Sparse matrix-vector multiplication
        double alpha = arma::dot(rVec, zVec) / arma::dot(pVec, ApVec); // Step size
        xVec = xVec + alpha * pVec;

	arma::vec r1Vec = rVec - alpha * ApVec; // Update residual

        arma::vec z1Vec = minvVec % r1Vec; // Apply preconditioner to new residual
        double beta = arma::dot(z1Vec, r1Vec) / arma::dot(zVec, rVec); // Update beta

        pVec = z1Vec + beta * pVec; // Update search direction
        zVec = z1Vec; // Update preconditioned residual
        rVec = r1Vec; // Update residual
        sumr2 = arma::dot(rVec, rVec); // Update residual norm
    }

    if (iter >= maxiterPCG) {
        std::cout << "PCG in getPCG1ofSigmaAndGtilde did not converge. You may increase maxiter number." << std::endl;
    }

    return xVec;
}



}
