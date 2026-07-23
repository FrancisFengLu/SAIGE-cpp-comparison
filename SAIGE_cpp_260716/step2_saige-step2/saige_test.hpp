// Standalone port of SAIGE/src/SAIGE_test.hpp
// SAIGEClass: core association test class

#ifndef SAIGE_HPP
#define SAIGE_HPP

#include <armadillo>
#include <random>
#include <atomic>
#include <cstdint>

// A3 (2026-07-15): two-stage-lite Firth. g_firthDefer (thread_local) makes the
// first getMarkerPval pass record the Firth candidate (sets t_isFirth) but skip
// execution; Firth runs once, in the fast-test recompute. Eliminates duplicate
// Firth (firthCut 0.01 < fastCut 0.05 => every Firth marker recomputes).
// g_firthFitCalls counts real fast_logistf_fit_simple invocations.
extern thread_local bool g_firthDefer;
extern std::atomic<std::uint64_t> g_firthFitCalls;


namespace SAIGE{

// Per-marker context (Phase A of step2 parallelism plan).
// Holds the small set of per-marker scalars that were previously mutated on
// SAIGEClass between each compute call. By passing this struct as a parameter
// to the new scoreTest*/getMarkerPval overloads, multiple markers can be
// processed concurrently in Phase B without racing on shared scalars.
struct PerMarkerCtx {
    bool flagSparseGRM_cur = false;
    bool isnoadjCov_cur = false;
    double varRatioVal = 1.0;
};

// Fused-kernel mode + A/B validation accumulators (Pillar 1).
//   g_fusedMode: 0 = off (scalar scoreTestFast, default)
//                1 = A/B check (run BOTH, compare, output scalar — validation)
//                2 = fused production (use scoreTestFast_fused for output)
// The accumulators are updated under critical(fusedAB) in mode 1; main() prints
// a summary and the pass/fail (<1e-6) verdict after the marker loop.
extern int    g_fusedMode;
extern double g_fusedMaxRelTstat;
extern double g_fusedMaxRelVar1;
extern long   g_fusedNCompared;
extern long   g_fusedNExceed;

class SAIGEClass
{
    private:
      arma::mat m_XVX;
      arma::mat m_XVX_inv_XV;
      arma::mat m_X;
      // Fused-kernel (Pillar 1) transposed copies: p x N, so each sample's
      // covariate row X(i,:) / A(i,:) is a CONTIGUOUS column (cache-friendly,
      // prefetchable) for scoreTestFast_fused's single-pass carrier loop.
      arma::mat m_Xt;   // == m_X.t()
      arma::mat m_At;   // == m_XVX_inv_XV.t()
      // Precomputed full-N sums (marker-independent) for scoreTestFast_noadjCov_fused.
      double m_sum_mu2 = 0.0;   // == arma::accu(m_mu2)
      double m_sum_res = 0.0;   // == arma::accu(m_res)
      arma::mat m_Sigma_iXXSigma_iX;
      arma::vec m_res;
      arma::vec m_resout;
      arma::vec m_mu;
      arma::vec m_mu2;
      arma::vec m_tauvec;
      arma::vec  m_S_a;
      std::string m_impute_method;
      // RNG engine replacing R's set.seed
      std::mt19937 m_rng_engine;

    public:
      std::vector<uint32_t> m_condition_genoIndex;
           std::string m_traitType;
      arma::mat m_XXVX_inv;
      arma::mat m_XV;
      int m_n, m_p; //MAIN Dimensions: sample size, number of covariates
      double m_varRatioVal;
      arma::vec m_varRatio_sparse;
      arma::vec m_varRatio_null;
      arma::vec m_varRatio_null_noXadj;
      arma::vec m_y;

      bool m_isOutputAFinCaseCtrl;
      bool m_isOutputNinCaseCtrl;
      bool m_isOutputHetHomCountsinCaseCtrl;
      arma::uvec m_case_indices;
      arma::uvec m_ctrl_indices;
      arma::uvec m_case_hom_indices;
      arma::uvec m_case_het_indices;
      arma::uvec m_ctrl_hom_indices;
      arma::uvec m_ctrl_het_indices;
      int m_n_case;
      int m_n_ctrl;
      arma::sp_mat m_SigmaMat_sp;
      bool m_flagSparseGRM;
      bool m_flagSparseGRM_cur;
      bool m_isFastTest;
      bool m_isnoadjCov;
      bool m_isnoadjCov_cur;
      double m_pval_cutoff_for_fastTest;
      double m_SPA_Cutoff;
      arma::vec m_cateVarRatioMinMACVecExclude;
      arma::vec m_cateVarRatioMaxMACVecInclude;
      arma::mat m_P2Mat_cond;
      int m_numMarker_cond;
      arma::mat m_VarInvMat_cond;
      arma::mat m_VarMat_cond;
      arma::vec m_Tstat_cond;
      arma::vec m_G2_Weight_cond;
      arma::vec m_MAF_cond;
      double  m_qsum_cond;
      arma::vec m_gsum_cond;
      std::vector<std::string> m_p_cond;
      arma::vec m_scalefactor_G2_cond;
      arma::mat m_VarInvMat_cond_scaled_weighted;
      bool m_isCondition;
      bool m_is_Firth_beta;
      double m_pCutoffforFirth;
     arma::vec  m_offset;
      bool m_isVarPsadj;
      bool m_islog10p;
	arma::sp_mat m_spSigmaMat; // Declare the sparse matrix attribute
	arma::vec m_diagSigma;    // Precompute diagonal


  ////////////////////// -------------------- functions ----------------
  //------------------ //////////////////////

  SAIGEClass(
        arma::mat & t_XVX,
        arma::mat  t_XXVX_inv,
        arma::mat & t_XV,
        arma::mat & t_XVX_inv_XV,
	arma::mat & t_sigmainvX_XsignmainXtXinv,
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
	arma::vec & t_resout);

   // Replaced R's set.seed with std::mt19937
   void set_seed(unsigned int seed);

   void scoreTest(arma::vec & t_GVec,
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
		     arma::uvec & t_indexForNonZero);

   // Ctx-based overload (Phase A): uses ctx.flagSparseGRM_cur / ctx.varRatioVal
   // instead of the class-member singletons.
   void scoreTest(arma::vec & t_GVec,
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
                     const PerMarkerCtx& ctx);


void scoreTestFast(arma::vec & t_GVec,
                     arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);

// Ctx-based overload (Phase A): reads ctx.varRatioVal instead of m_varRatioVal.
void scoreTestFast(arma::vec & t_GVec,
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
                     const PerMarkerCtx& ctx);

// Pillar 1: single fused pass over the carrier (nonzero) samples — no dense
// N-vector, no .elem()/.rows() temporaries. Algebraic reduction of
// scoreTestFast (covar-adjusted, non-sparse GRM). NOT bit-identical to
// scoreTestFast (hand loops + term cancellations); algebraically equal,
// matches to ~1e-10 rel. Requires m_Xt / m_At.
void scoreTestFast_fused(const arma::vec& t_GVec,
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
                     double& t_StdStat) const;

// Pillar 1 (no-covariate-adjustment variant): fused carrier pass for the
// isnoadjCov score test. Needs only res/mu2 at carriers + precomputed full-N
// sums m_sum_mu2 / m_sum_res. No X/A, no m_Xt/m_At. Matches
// scoreTestFast_noadjCov to ~1e-13.
void scoreTestFast_noadjCov_fused(const arma::vec& t_GVec,
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
                     double& t_StdStat) const;

void scoreTestFast_noadjCov(arma::vec & t_GVec,
                arma::uvec & t_indexForNonZero,
                     double& t_Beta,
                     double& t_seBeta,
                     std::string& t_pval_str,
                     double& t_pval,
                     bool& t_islogp,
                     double t_altFreq,
                     double &t_Tstat,
                     double &t_var1,
                     double &t_var2);

// Ctx-based overload (Phase A): reads ctx.varRatioVal instead of m_varRatioVal.
void scoreTestFast_noadjCov(arma::vec & t_GVec,
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
                     const PerMarkerCtx& ctx);

// Phase C: matrix-level first pass. Computes Tstat / var1 / var2 for a block
// of B markers via a single BLAS-3 GEMM, replacing B BLAS-2 dot products.
// This is the dense-block analogue of scoreTestFast (with covariate adjustment,
// non-sparse GRM). Caller is responsible for routing columns that need a
// different path (sparse-GRM, noadjCov, SPA/ER/Firth recompute, conditional)
// to the scalar code; this block method does only the cheap first-pass score
// test math.
//
// G: N x B dense block of genotypes (already imputed / QC-passed).
// varRatioVec: length B; per-column variance ratio (NaN entries skipped).
// Outputs are length B; entries for skipped columns are left untouched.
// validMask: length B, true for columns that were computed.
void scoreTestFast_block(const arma::mat& G,
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
                          std::vector<std::string>& pvalStr) const;



     void set_flagSparseGRM_cur(bool t_flagSparseGRM_cur);

     void get_mu(arma::vec & t_mu);

     void getadjG(arma::vec & t_GVec, arma::vec & g);
     void getadjGFast(arma::vec & t_GVec, arma::vec & g,  arma::uvec & iIndex);

     void getMarkerPval(arma::vec & t_GVec,
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
                                bool t_isSparseGRM);

   // Ctx-based overload (Phase A): uses ctx.* instead of class-member singletons.
   void getMarkerPval(arma::vec & t_GVec,
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
                                const PerMarkerCtx& ctx);

    void getindices(arma::uvec & t_case_indices,
      arma::uvec & t_ctrl_indices);

    bool assignVarianceRatio(double MAC, bool issparseforVR);

    // Pure (Phase A): returns the variance ratio that assignVarianceRatio would
    // have written into m_varRatioVal, without mutating any member.
    // hasVarRatio is set to the same flag the mutating version returns.
    double computeVarianceRatio(double MAC, bool issparseforVR, bool& hasVarRatio) const;

    void assignSingleVarianceRatio(bool issparseforVR);

    // Pure (Phase A): same value assignSingleVarianceRatio would write.
    double computeSingleVarianceRatio(bool issparseforVR) const;

    void assignSingleVarianceRatio_withinput(double t_varRatioVal);

    void assignConditionFactors(
      arma::mat & t_P2Mat_cond,
      arma::mat & t_VarInvMat_cond,
      arma::mat & t_VarMat_cond,
      arma::vec & t_Tstat_cond,
      arma::vec & t_G2_Weight_cond,
      arma::vec & t_MAF_cond,
      double t_qsum_cond,
      arma::vec & t_gsum_cond,
      std::vector<std::string> & t_p_cond);

     void assignConditionFactors_scalefactor(
        arma::vec & t_scalefactor_G2_cond);


    void extract_XV_XXVX_inv(arma::mat & t_XV, arma::mat & t_XXVX_inv);

    void fast_logistf_fit_simple(arma::mat & x,
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
	bool & isfirthconverge);

     void set_isnoadjCov_cur(bool t_isnoadjCov_cur);
     bool  assignVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj);

     // Pure (Phase A): see computeVarianceRatio above.
     double computeVarianceRatio(double MAC, bool issparseforVR, bool isnoXadj, bool& hasVarRatio) const;

     void assignSingleVarianceRatio(bool issparseforVR, bool isnoXadj);

     // Pure (Phase A): same value assignSingleVarianceRatio (3-arg) would write.
     double computeSingleVarianceRatio(bool issparseforVR, bool isnoXadj) const;

     arma::vec getPCG1ofSigmaAndGtilde(arma::vec& bVec, int maxiterPCG, double tolPCG);

};
}
#endif
