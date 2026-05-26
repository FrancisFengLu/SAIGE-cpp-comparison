// Standalone port of SAIGE/src/SAIGE_test.cpp
// SAIGEClass: score test, SPA, variance ratio, Firth correction

#define ARMA_USE_SUPERLU 1

#include <armadillo>

#include "saige_test.hpp"
#include "spa.hpp"
#include "er_binary.hpp"
#include "UTIL.hpp"
#include "getMem.hpp"
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
    m_S_a = t_S_a;
    m_res = t_res;
    m_resout = t_resout;
    m_mu2 = t_mu2;
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
    arma::vec Sm, var2m;
    double S, var2;
    getadjGFast(t_GVec, t_gtilde, t_indexForNonZero);

    if(t_is_region && m_traitType == "binary"){
      t_gy = dot(t_gtilde, m_y);
     }
    S = dot(t_gtilde, m_res);
    S = S/m_tauvec[0];

    if(!m_flagSparseGRM_cur){
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
    double var1 = var2 * m_varRatioVal;

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

    arma::vec g1 = t_GVec.elem(t_indexForNonZero);
    arma::mat X1 = m_X.rows(t_indexForNonZero);
    arma::mat A1 = m_XVX_inv_XV.rows(t_indexForNonZero);
    arma::vec mu21;
    arma::vec res1 = m_res.elem(t_indexForNonZero);
    arma::vec Z = A1.t() * g1;
    arma::vec B = X1 * Z;
    arma::vec g1_tilde = g1 - B;
    double var1, var2, S, S1, S2, g1tildemu2;
    arma::vec S_a2;
    double Bmu2;
    arma::mat  ZtXVXZ = Z.t() * m_XVX * Z;
    if(m_traitType == "binary" || m_traitType == "survival"){
      mu21  = m_mu2.elem(t_indexForNonZero);
      g1tildemu2 = dot(square(g1_tilde), mu21);
      Bmu2 = arma::dot(square(B),  mu21);
      var2 = ZtXVXZ(0,0) - Bmu2 + g1tildemu2;
    }else if(m_traitType == "quantitative"){
      Bmu2 = dot(g1, B);
      var2 = ZtXVXZ(0,0)*m_tauvec[0] +  dot(g1,g1) - 2*Bmu2;
    }

    var1 = var2 * m_varRatioVal;
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

    double var1 = var2 * m_varRatioVal;
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



  t_isFirth = false;
  std::string t_pval_str;
  double t_var2, t_SPApval;
  bool isScoreFast = true;
  if(m_flagSparseGRM_cur){
    isScoreFast = false;
  }


  double pval_noadj, pval, t_qval_Firth;
  bool ispvallog;

if(!m_flagSparseGRM_cur && t_isnoadjCov){
	is_gtilde = false;
        isScoreFast = true;
	scoreTestFast_noadjCov(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq,t_Tstat, t_var1, t_var2);

}else if(m_flagSparseGRM_cur){
	is_gtilde = true;
        isScoreFast = false;
        scoreTest(t_GVec, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2, t_gtilde, t_P2Vec, t_gy, is_region, iIndex);
}else{
	is_gtilde = false;
        isScoreFast = true;
        scoreTestFast(t_GVec, iIndex, t_Beta, t_seBeta, t_pval_noSPA, pval_noadj, ispvallog, t_altFreq, t_Tstat, t_var1, t_var2);
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
	if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
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
   if(t_isFirth){
	if(!is_gtilde){
                getadjGFast(t_GVec, t_gtilde, iIndex);
                is_gtilde = true;
        }
	arma::mat x(t_GVec.n_elem, 2, arma::fill::ones);
	x.col(1) = t_gtilde;
	arma::vec init(2, arma::fill::zeros);
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
        t_G1tilde_P_G2tilde = sqrt(m_varRatioVal) * t_gtilde.t() * m_P2Mat_cond;
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


        if(p_iIndexComVecSize >= 0.5 && !m_flagSparseGRM_cur){
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
      if(!m_flagSparseGRM_cur){
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
  arma::vec pi_0 = -x * beta - offset;
  pi_0 = arma::exp(pi_0) + 1;
  arma::vec pi = 1/pi_0;
  int evals = 1;
  arma::vec beta_old;
  arma::mat oneVec(k, 1 , arma::fill::ones);
  arma::mat XX_covs(k, k, arma::fill::zeros);
  while(iter <= maxit){
        beta_old = beta;
        arma::vec wpi = pi % (1 - pi);
        arma::vec W2 = arma::sqrt(wpi);
        arma::mat XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                XW2.col(j) = x.col(j) % W2;
        }

        arma::mat Q;
        arma::mat R;
        arma::qr_econ(Q, R, XW2);
        arma::vec h = Q % Q * oneVec;
        arma::vec U_star(2, arma::fill::zeros);
        arma::vec ypih;
        if(firth){
                ypih = (y - pi) + (h % (0.5 - pi));
        }else{
                ypih = (y - pi);
        }
        arma::vec xcol(n, arma::fill::zeros);
        U_star = x.t() * ypih;

        arma::mat XX_XW2(n, k, arma::fill::zeros);
        for(int j = 0; j < k; j++){
                xcol = x.col(j);
                XX_XW2.col(j) = xcol % W2;
        }
        arma::mat XX_Fisher = XX_XW2.t() * (XX_XW2);
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
        pi_0 = -x * beta - offset;
        pi_0 = arma::exp(pi_0) + 1;
        pi = 1/pi_0;
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
