// Standalone C++ port of SAIGE Step 2 main entry point
// Ported from: SAIGE/src/Main.cpp
//
// Phase 3: CLI entry point (YAML config), setAssocTest_GlobalVarsInCPP,
//           Unified_getMarkerPval, openOutfile_single, writeOutfile_single,
//           mainMarkerInCPP
// Phase 6: Region/gene-based association testing (mainRegionInCPP),
//           openOutfile, writeOutfile_BURDEN, openOutfile_singleinGroup,
//           writeOutfile_singleInGroup, setRegion_GlobalVarsInCPP,
//           SPA Phi adjustment, SKAT/BURDEN/SKAT-O dispatch, CCT combination
//
// Conversions from original SAIGE Main.cpp:
//   1. #include <RcppArmadillo.h> --> #include <armadillo>
//   2. Rcpp::stop(...)            --> throw std::runtime_error(...)
//   3. Rcpp::List returns         --> C++ structs (not needed here)
//   4. // [[Rcpp::export]]        --> removed
//   5. R statistical functions    --> boost::math
//   6. Rcpp::Environment set.seed --> std::mt19937
//   7. Global variables preserved exactly as in SAIGE Main.cpp
//   8. R orchestration (SAIGE_SPATest_Region.R) merged into C++

#include <armadillo>

#include <vector>
#include <thread>
#include <memory>
#include <chrono>
#include <omp.h>

// OpenBLAS thread-count knob (no public header in some installs); declared here
// so we can clamp BLAS to 1 thread per OMP worker and avoid nested oversubscription.
extern "C" void openblas_set_num_threads(int);

#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <numeric>

#include <yaml-cpp/yaml.h>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

#include "null_model_loader.hpp"
#include "genotype_reader.hpp"
#include "saige_test.hpp"
#include "UTIL.hpp"
#include "cct.hpp"
#include "spa.hpp"
#include "getMem.hpp"
#include "group_file.hpp"
#include "skat.hpp"
#include "ldmat.hpp"

// ============================================================
// Global variables (exact match from SAIGE/src/Main.cpp)
// ============================================================

// Global objects for analysis methods
static SAIGE::SAIGEClass* ptr_gSAIGEobj = NULL;

// Global variables for analysis
std::string g_impute_method;        // "mean", "minor", or "drop"
double g_missingRate_cutoff;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;      // for RVs whose MAC < this value
double g_marker_minINFO_cutoff;
arma::vec g_region_maxMAF_cutoff;
double g_min_gourpmac_for_burdenonly;
double g_maxMAFLimit;
unsigned int g_region_maxMarkers_cutoff;
bool g_isOutputMoreDetails;
int g_marker_chunksize;
int g_nThreads = 1;
// Phase D (Wave 1.3): number of BGEN decoder threads in BgenStreamer.
// Only used by mainMarkerInCPP when t_genoType=="bgen". Set via YAML key
// `bgenDecoders` (default 4). Region/group path (mainRegionInCPP) is not
// wired to use the streamer.
int g_bgenDecoders = 4;

std::string g_method_to_CollapseUltraRare;
double g_DosageCutoff_for_UltraRarePresence;

double g_dosage_zerod_MAC_cutoff;
double g_dosage_zerod_cutoff;
bool g_markerTestEnd = false;
arma::vec g_weights_beta(2);

bool g_is_Firth_beta;
double g_pCutoffforFirth;
double g_MACCutoffforER;
bool g_is_rewrite_XnonPAR_forMales = false;

arma::uvec g_indexInModel_male;
arma::umat g_X_PARregion_mat;

// Output file streams
std::ofstream OutFile;
std::ofstream OutFile_singleInGroup;
std::ofstream OutFile_single;
std::ofstream OutFile_singleInGroup_temp;

// Output file prefix strings
std::string g_outputFilePrefixGroup;
std::string g_outputFilePrefixSingleInGroup;
std::string g_outputFilePrefixSingleInGroup_temp;
std::string g_outputFilePrefixSingle;


// ============================================================
// setAssocTest_GlobalVarsInCPP
// Direct port from SAIGE/src/Main.cpp lines 142-167
// ============================================================
void setAssocTest_GlobalVarsInCPP(std::string t_impute_method,
                                   double t_missing_cutoff,
                                   double t_min_maf_marker,
                                   double t_min_mac_marker,
                                   double t_min_info_marker,
                                   double t_dosage_zerod_cutoff,
                                   double t_dosage_zerod_MAC_cutoff,
                                   arma::vec & t_weights_beta,
                                   std::string t_outputFilePrefix,
                                   double t_MACCutoffforER)
{
    g_impute_method = t_impute_method;
    g_missingRate_cutoff = t_missing_cutoff;
    g_marker_minMAF_cutoff = t_min_maf_marker;
    g_marker_minMAC_cutoff = t_min_mac_marker;
    g_marker_minINFO_cutoff = t_min_info_marker;
    g_dosage_zerod_cutoff = t_dosage_zerod_cutoff;
    g_dosage_zerod_MAC_cutoff = t_dosage_zerod_MAC_cutoff;
    g_weights_beta = t_weights_beta;
    g_outputFilePrefixGroup = t_outputFilePrefix;
    g_outputFilePrefixSingleInGroup = t_outputFilePrefix + ".singleAssoc.txt";
    g_outputFilePrefixSingleInGroup_temp = t_outputFilePrefix + ".singleAssoc.txt_temp";
    g_outputFilePrefixSingle = t_outputFilePrefix;
    g_MACCutoffforER = t_MACCutoffforER;
}


// ============================================================
// setMarker_GlobalVarsInCPP
// Direct port from SAIGE/src/Main.cpp lines 183-191
// ============================================================
void setMarker_GlobalVarsInCPP(bool t_isOutputMoreDetails,
                                int t_marker_chunksize)
{
    g_isOutputMoreDetails = t_isOutputMoreDetails;
    g_marker_chunksize = t_marker_chunksize;
}


// ============================================================
// Unified_getMarkerPval
// Direct port from SAIGE/src/Main.cpp lines 807-847
// A unified function to get marker-level p-value
// ============================================================
void Unified_getMarkerPval(
    arma::vec & t_GVec,
    bool t_isOnlyOutputNonZero,
    arma::uvec & t_indexForNonZero_vec,
    arma::uvec & t_indexForZero_vec,
    double& t_Beta,
    double& t_seBeta,
    std::string& t_pval,
    std::string& t_pval_noSPA,
    double& t_Tstat,
    double& t_gy,
    double& t_varT,
    double t_altFreq,
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
    arma::rowvec & t_G1tilde_P_G2tilde_Vec,
    bool & t_isFirth,
    bool & t_isFirthConverge,
    bool t_isER,
    bool t_isnoadjCov,
    bool t_isSparseGRM)
{
    if (t_isOnlyOutputNonZero == true)
        throw std::runtime_error(
            "When using SAIGE method to calculate marker-level p-values, "
            "'t_isOnlyOutputNonZero' should be false.");

    ptr_gSAIGEobj->getMarkerPval(
        t_GVec, t_indexForNonZero_vec, t_indexForZero_vec,
        t_Beta, t_seBeta, t_pval, t_pval_noSPA,
        t_altFreq, t_Tstat, t_gy, t_varT,
        t_isSPAConverge, t_gtilde, is_gtilde,
        is_region, t_P2Vec,
        t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c,
        t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec,
        t_isFirth, t_isFirthConverge, t_isER,
        t_isnoadjCov, t_isSparseGRM);
}


// Ctx-based overload (Phase A of step2 parallelism plan): same as above but
// routes to the SAIGEClass::getMarkerPval ctx overload so per-marker scalars
// (flagSparseGRM_cur / isnoadjCov_cur / varRatioVal) come from `ctx` rather
// than mutable class state.
void Unified_getMarkerPval(
    arma::vec & t_GVec,
    bool t_isOnlyOutputNonZero,
    arma::uvec & t_indexForNonZero_vec,
    arma::uvec & t_indexForZero_vec,
    double& t_Beta,
    double& t_seBeta,
    std::string& t_pval,
    std::string& t_pval_noSPA,
    double& t_Tstat,
    double& t_gy,
    double& t_varT,
    double t_altFreq,
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
    arma::rowvec & t_G1tilde_P_G2tilde_Vec,
    bool & t_isFirth,
    bool & t_isFirthConverge,
    bool t_isER,
    bool t_isnoadjCov,
    bool t_isSparseGRM,
    const SAIGE::PerMarkerCtx& ctx)
{
    if (t_isOnlyOutputNonZero == true)
        throw std::runtime_error(
            "When using SAIGE method to calculate marker-level p-values, "
            "'t_isOnlyOutputNonZero' should be false.");

    ptr_gSAIGEobj->getMarkerPval(
        t_GVec, t_indexForNonZero_vec, t_indexForZero_vec,
        t_Beta, t_seBeta, t_pval, t_pval_noSPA,
        t_altFreq, t_Tstat, t_gy, t_varT,
        t_isSPAConverge, t_gtilde, is_gtilde,
        is_region, t_P2Vec,
        t_isCondition, t_Beta_c, t_seBeta_c, t_pval_c, t_pval_noSPA_c,
        t_Tstat_c, t_varT_c, t_G1tilde_P_G2tilde_Vec,
        t_isFirth, t_isFirthConverge, t_isER,
        t_isnoadjCov, t_isSparseGRM, ctx);
}


// ============================================================
// setSAIGEobjInCPP
// Direct port from SAIGE/src/Main.cpp lines 918-989
// ============================================================
void setSAIGEobjInCPP(arma::mat & t_XVX,
                       arma::mat & t_XXVX_inv,
                       arma::mat & t_XV,
                       arma::mat & t_XVX_inv_XV,
                       arma::mat & t_Sigma_iXXSigma_iX,
                       arma::mat & t_X,
                       arma::vec & t_S_a,
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
                       arma::vec & t_resout)
{
    ptr_gSAIGEobj = new SAIGE::SAIGEClass(
        t_XVX,
        t_XXVX_inv,
        t_XV,
        t_XVX_inv_XV,
        t_Sigma_iXXSigma_iX,
        t_X,
        t_S_a,
        t_res,
        t_mu2,
        t_mu,
        t_varRatio_sparse,
        t_varRatio_null,
        t_varRatio_null_noXadj,
        t_cateVarRatioMinMACVecExclude,
        t_cateVarRatioMaxMACVecInclude,
        t_SPA_Cutoff,
        t_tauvec,
        t_traitType,
        t_y,
        t_impute_method,
        t_flagSparseGRM,
        t_isFastTest,
        t_isnoadjCov,
        t_pval_cutoff_for_fastTest,
        t_locationMat,
        t_valueVec,
        t_dimNum,
        t_isCondition,
        t_condition_genoIndex,
        t_is_Firth_beta,
        t_pCutoffforFirth,
        t_offset,
        t_resout);
}


// ============================================================
// assign_conditionMarkers_factors
// Ported from SAIGE/src/Main.cpp lines 2189-2366
// Reads conditioning marker genotypes, computes their test stats,
// and stores the results in SAIGEClass for conditional analysis.
// ============================================================
void assign_conditionMarkers_factors(
    std::string t_genoType,
    std::vector<uint32_t> & t_genoIndex,
    unsigned int t_n,
    arma::vec & t_weight_cond)
{
    ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
    bool isImpute = false;
    unsigned int q = t_genoIndex.size();
    arma::mat P1Mat(q, t_n);
    arma::mat P2Mat(t_n, q);
    arma::mat VarInvMat(q, q);
    arma::vec TstatVec(q);
    std::vector<std::string> pVec(q, "NA");
    arma::vec MAFVec(q);
    arma::vec gyVec(q);
    arma::vec w0G2_cond_Vec(q);
    arma::vec gsumVec(t_n, arma::fill::zeros);
    boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);
    arma::vec GVec(t_n);
    std::string pval, pval_noSPA;
    double Beta, seBeta, Tstat, varT, gy, w0G2_cond;
    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    arma::vec P2Vec(t_n);

    std::string pval_c, pval_noSPA_c;
    double Beta_c, seBeta_c, Tstat_c, varT_c;
    arma::rowvec G1tilde_P_G2tilde_Vec;
    bool isCondition = false;  // Don't condition on self when reading conditioning markers

    for (unsigned int i = 0; i < q; i++) {
        double altFreq, altCounts, missingRate, imputeInfo;
        std::vector<uint> indexForMissing;
        std::vector<uint> indexZeroVec;
        std::vector<uint> indexNonZeroVec;
        std::string chr, ref, alt, marker;
        uint32_t pd;
        bool flip = false;

        bool isOutputIndexForMissing = true;
        bool isOnlyOutputNonZero = false;

        uint64_t gIndex = (uint64_t)t_genoIndex[i];
        uint64_t gIndex_prev = 0;
        if (i > 0) {
            gIndex_prev = (uint64_t)t_genoIndex[i-1];
        }

        bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr,
                                                  altFreq, altCounts, missingRate, imputeInfo,
                                                  isOutputIndexForMissing,
                                                  indexForMissing,
                                                  isOnlyOutputNonZero,
                                                  indexNonZeroVec, GVec, isImpute);
        if (!isReadMarker) {
            throw std::runtime_error("Failed to read conditioning marker at index " +
                                     std::to_string(t_genoIndex[i]));
        }

        std::string info = chr + ":" + std::to_string(pd) + "_" + ref + "/" + alt;

        double MAF = std::min(altFreq, 1.0 - altFreq);
        double MAC = MAF * 2 * t_n * (1.0 - missingRate);

        bool hasVarRatio;
        if ((ptr_gSAIGEobj->m_varRatio_null).n_elem == 1) {
            ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM, false);
        } else {
            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM, false);
        }

        flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing,
                                  g_impute_method, g_dosage_zerod_cutoff,
                                  g_dosage_zerod_MAC_cutoff, MAC,
                                  indexZeroVec, indexNonZeroVec);

        arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
        indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
        indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

        MAF = std::min(altFreq, 1.0 - altFreq);
        MAC = std::min(altCounts, 2.0 * t_n - altCounts);

        arma::vec gtildeVec;

        if (MAC <= g_MACCutoffforER && ptr_gSAIGEobj->m_traitType == "binary") {
            Unified_getMarkerPval(
                GVec, false, indexNonZeroVec_arma, indexZeroVec_arma,
                Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq,
                isSPAConverge, gtildeVec, is_gtilde, true, P2Vec,
                isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                is_Firth, is_FirthConverge, true, false,
                ptr_gSAIGEobj->m_flagSparseGRM);
        } else {
            Unified_getMarkerPval(
                GVec, false, indexNonZeroVec_arma, indexZeroVec_arma,
                Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq,
                isSPAConverge, gtildeVec, is_gtilde, true, P2Vec,
                isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                is_Firth, is_FirthConverge, false, false,
                ptr_gSAIGEobj->m_flagSparseGRM);
        }

        P1Mat.row(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal) * gtildeVec.t();
        P2Mat.col(i) = sqrt(ptr_gSAIGEobj->m_varRatioVal) * P2Vec;
        MAFVec(i) = MAF;

        if (MAF == 0.0) {
            std::cerr << "ERROR: Conditioning marker " << info << " is monomorphic" << std::endl;
        }

        if (!t_weight_cond.is_zero()) {
            w0G2_cond = t_weight_cond(i);
        } else {
            w0G2_cond = boost::math::pdf(beta_dist, MAF);
        }
        w0G2_cond_Vec(i) = w0G2_cond;
        gyVec(i) = gy * w0G2_cond;
        gsumVec = gsumVec + GVec * w0G2_cond;
        TstatVec(i) = Tstat;
        pVec.at(i) = pval;

        std::cout << "    Condition marker " << i << ": " << info
                  << " MAF=" << MAF << " MAC=" << MAC
                  << " Beta=" << Beta << " Tstat=" << Tstat
                  << " pval=" << pval << std::endl;
    }

    arma::mat VarMat = P1Mat * P2Mat;
    VarInvMat = arma::pinv(VarMat);
    double qsum = arma::accu(gyVec);
    arma::vec gsumtildeVec;

    ptr_gSAIGEobj->getadjG(gsumVec, gsumtildeVec);
    ptr_gSAIGEobj->assignConditionFactors(
        P2Mat,
        VarInvMat,
        VarMat,
        TstatVec,
        w0G2_cond_Vec,
        MAFVec,
        qsum,
        gsumtildeVec,
        pVec);

    std::cout << "  Conditioning factors assigned successfully." << std::endl;
}


// ============================================================
// openOutfile_single
// Direct port from SAIGE/src/Main.cpp lines 2587-2643
// ============================================================
bool openOutfile_single(std::string t_traitType,
                         bool t_isImputation,
                         bool isappend,
                         bool t_isMoreOutput)
{
    bool isopen;
    if (!isappend) {
        OutFile_single.open(g_outputFilePrefixSingle.c_str());
        isopen = OutFile_single.is_open();
        if (isopen) {
            OutFile_single << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
            if (t_isImputation) {
                OutFile_single << "imputationInfo\t";
            } else {
                OutFile_single << "MissingRate\t";
            }
            OutFile_single << "BETA\tSE\tTstat\tvar\tp.value\t";
            if (t_traitType == "binary" || t_traitType == "survival") {
                OutFile_single << "p.value.NA\tIs.SPA\t";
            }

            if (ptr_gSAIGEobj->m_isCondition) {
                OutFile_single << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile_single << "p.value.NA_c\t";
                }
            }

            if (t_traitType == "binary") {
                OutFile_single << "AF_case\tAF_ctrl\tN_case\tN_ctrl";
                if (t_isMoreOutput) {
                    OutFile_single << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
                }
                OutFile_single << "\n";
            } else if (t_traitType == "quantitative") {
                OutFile_single << "N\n";
            } else if (t_traitType == "survival") {
                OutFile_single << "AF_event\tAF_censor\tN_event\tN_censor";
                if (t_isMoreOutput) {
                    OutFile_single << "\tN_event_hom\tN_event_het\tN_censor_hom\tN_censor_het";
                }
                OutFile_single << "\n";
            }
        }
    } else {
        OutFile_single.open(g_outputFilePrefixSingle.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_single.is_open();
    }

    return isopen;
}


// ============================================================
// writeOutfile_single
// Direct port from SAIGE/src/Main.cpp lines 2646-2779
// ============================================================
void writeOutfile_single(bool t_isMoreOutput,
                          bool t_isImputation,
                          bool t_isCondition,
                          bool t_isFirth,
                          int mFirth,
                          int mFirthConverge,
                          std::string t_traitType,
                          std::vector<std::string> & chrVec,
                          std::vector<std::string> & posVec,
                          std::vector<std::string> & markerVec,
                          std::vector<std::string> & refVec,
                          std::vector<std::string> & altVec,
                          std::vector<double> & altCountsVec,
                          std::vector<double> & altFreqVec,
                          std::vector<double> & imputationInfoVec,
                          std::vector<double> & missingRateVec,
                          std::vector<double> & BetaVec,
                          std::vector<double> & seBetaVec,
                          std::vector<double> & TstatVec,
                          std::vector<double> & varTVec,
                          std::vector<std::string> & pvalVec,
                          std::vector<std::string> & pvalNAVec,
                          std::vector<bool> & isSPAConvergeVec,
                          std::vector<double> & Beta_cVec,
                          std::vector<double> & seBeta_cVec,
                          std::vector<double> & Tstat_cVec,
                          std::vector<double> & varT_cVec,
                          std::vector<std::string> & pval_cVec,
                          std::vector<std::string> & pvalNA_cVec,
                          std::vector<double> & AF_caseVec,
                          std::vector<double> & AF_ctrlVec,
                          std::vector<uint32_t> & N_caseVec,
                          std::vector<uint32_t> & N_ctrlVec,
                          std::vector<double> & N_case_homVec,
                          std::vector<double> & N_ctrl_hetVec,
                          std::vector<double> & N_case_hetVec,
                          std::vector<double> & N_ctrl_homVec,
                          std::vector<uint32_t> & N_Vec)
{
    int numtest = 0;
    for (unsigned int k = 0; k < pvalVec.size(); k++) {
        if (pvalVec.at(k) != "NA") {
            numtest = numtest + 1;
            OutFile_single << chrVec.at(k);
            OutFile_single << "\t";
            OutFile_single << posVec.at(k);
            OutFile_single << "\t";
            OutFile_single << markerVec.at(k);
            OutFile_single << "\t";
            OutFile_single << refVec.at(k);
            OutFile_single << "\t";
            OutFile_single << altVec.at(k);
            OutFile_single << "\t";
            OutFile_single << altCountsVec.at(k);
            OutFile_single << "\t";
            OutFile_single << altFreqVec.at(k);
            OutFile_single << "\t";

            if (t_isImputation) {
                OutFile_single << imputationInfoVec.at(k);
                OutFile_single << "\t";
            } else {
                OutFile_single << missingRateVec.at(k);
                OutFile_single << "\t";
            }
            OutFile_single << BetaVec.at(k);
            OutFile_single << "\t";
            OutFile_single << seBetaVec.at(k);
            OutFile_single << "\t";
            OutFile_single << TstatVec.at(k);
            OutFile_single << "\t";
            OutFile_single << varTVec.at(k);
            OutFile_single << "\t";
            OutFile_single << pvalVec.at(k);
            OutFile_single << "\t";

            if (t_traitType == "binary" || t_traitType == "survival") {
                OutFile_single << pvalNAVec.at(k);
                OutFile_single << "\t";
                OutFile_single << std::boolalpha << isSPAConvergeVec.at(k);
                OutFile_single << "\t";
            }
            if (t_isCondition) {
                OutFile_single << Beta_cVec.at(k);
                OutFile_single << "\t";
                OutFile_single << seBeta_cVec.at(k);
                OutFile_single << "\t";
                OutFile_single << Tstat_cVec.at(k);
                OutFile_single << "\t";
                OutFile_single << varT_cVec.at(k);
                OutFile_single << "\t";
                OutFile_single << pval_cVec.at(k);
                OutFile_single << "\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile_single << pvalNA_cVec.at(k);
                    OutFile_single << "\t";
                }
            }
            if (t_traitType == "binary" || t_traitType == "survival") {
                OutFile_single << AF_caseVec.at(k);
                OutFile_single << "\t";
                OutFile_single << AF_ctrlVec.at(k);
                OutFile_single << "\t";
                OutFile_single << N_caseVec.at(k);
                OutFile_single << "\t";
                OutFile_single << N_ctrlVec.at(k);

                if (t_isMoreOutput) {
                    OutFile_single << "\t";
                    OutFile_single << N_case_homVec.at(k);
                    OutFile_single << "\t";
                    OutFile_single << N_case_hetVec.at(k);
                    OutFile_single << "\t";
                    OutFile_single << N_ctrl_homVec.at(k);
                    OutFile_single << "\t";
                    OutFile_single << N_ctrl_hetVec.at(k);
                }
                OutFile_single << "\n";
            } else if (t_traitType == "quantitative") {
                OutFile_single << N_Vec.at(k);
                OutFile_single << "\n";
            }
        }
    }
    std::cout << numtest << " markers were tested." << std::endl;
    if (t_traitType == "binary") {
        if (t_isFirth) {
            std::cout << "Firth approx was applied to " << mFirth
                      << " markers. " << mFirthConverge
                      << " successfully converged." << std::endl;
        }
    }
}


// ============================================================
// mainMarkerInCPP
// Direct port from SAIGE/src/Main.cpp lines 219-679
// Single-variant marker testing loop
// ============================================================
void mainMarkerInCPP(
    std::string & t_genoType,     // "plink", "bgen", etc.
    std::string & t_traitType,
    std::vector<std::string> & t_genoIndex_prev,
    std::vector<std::string> & t_genoIndex,
    bool & t_isMoreOutput,
    bool & t_isImputation,
    bool & t_isFirth)
{
    int q = t_genoIndex.size();  // number of markers

    // set up output vectors
    std::vector<std::string> markerVec(q);
    std::vector<std::string> chrVec(q);
    std::vector<std::string> posVec(q);
    std::vector<std::string> refVec(q);
    std::vector<std::string> altVec(q);

    std::vector<std::string> infoVec(q);
    std::vector<double> altFreqVec(q);
    std::vector<double> altCountsVec(q);
    std::vector<double> imputationInfoVec(q);
    std::vector<double> missingRateVec(q);
    std::vector<double> BetaVec(q, arma::datum::nan);
    std::vector<double> seBetaVec(q, arma::datum::nan);
    std::vector<std::string> pvalVec(q, "NA");
    std::vector<double> TstatVec(q, arma::datum::nan);
    std::vector<double> varTVec(q, arma::datum::nan);
    std::vector<std::string> pvalNAVec(q, "NA");

    bool isCondition = ptr_gSAIGEobj->m_isCondition;
    std::vector<double> Beta_cVec(q, arma::datum::nan);
    std::vector<double> seBeta_cVec(q, arma::datum::nan);
    std::vector<std::string> pval_cVec(q, "NA");
    std::vector<double> Tstat_cVec(q, arma::datum::nan);
    std::vector<double> varT_cVec(q, arma::datum::nan);
    std::vector<std::string> pvalNA_cVec(q, "NA");
    // BUG 1 fix: G1tilde_P_G2tilde_Vec was previously declared here (shared across
    // threads). Inside the OMP parallel for below, each thread calls .clear() and
    // then has it written by getMarkerPval -> SAIGEClass::getMarkerPval (saige_test.cpp
    // line 820), which assigns `sqrt(varRatio)*gtilde.t()*m_P2Mat_cond` to it.
    // Thread A's .clear() (size becomes 1x0) raced with thread B's read at the
    // matmul `t_G1tilde_P_G2tilde * m_VarInvMat_cond` (3x3), producing
    // "incompatible matrix dimensions: 1x0 and 3x3". Declaring it inside the
    // loop body (see below) gives each thread its own per-iteration scratch.
    // Cache the cond marker count once (used only for sizing inside the loop).
    const int g_numMarker_cond_local = ptr_gSAIGEobj->m_numMarker_cond;

    std::vector<bool> isSPAConvergeVec(q);
    std::vector<double> AF_caseVec(q);
    std::vector<double> AF_ctrlVec(q);
    std::vector<uint32_t> N_caseVec(q);
    std::vector<uint32_t> N_ctrlVec(q);
    std::vector<double> N_case_homVec(q);
    std::vector<double> N_ctrl_hetVec(q);
    std::vector<double> N_case_hetVec(q);
    std::vector<double> N_ctrl_homVec(q);
    std::vector<uint32_t> N_Vec(q);
    int n = ptr_gSAIGEobj->m_n;

    bool isSingleVarianceRatio = true;
    if ((ptr_gSAIGEobj->m_varRatio_null).n_elem == 1) {
        // One-time setup before parallel region; OK to mutate here because no
        // threads are running yet.
        ptr_gSAIGEobj->assignSingleVarianceRatio(
            ptr_gSAIGEobj->m_flagSparseGRM,
            ptr_gSAIGEobj->m_isnoadjCov);
    } else {
        isSingleVarianceRatio = false;
    }

    int mFirth = 0;
    int mFirthConverge = 0;
    // Phase B: end-of-stream flag (some readers signal EOF mid-iteration). With
    // OpenMP we cannot 'break' out of a parallel for, so we record the lowest
    // failing index and skip work for i >= that index.
    int firstEndIdx = q;  // q == "no end yet"

    // Phase D (Wave 1.3): BGEN block-read + decode pipeline. For t_genoType
    // == "bgen" we spin up a streamer (1 reader + N decoders + bounded queue)
    // and consume markers in original order via getNext(). Other genoTypes go
    // through the existing per-iteration Unified_getOneMarker path.
    std::unique_ptr<BGEN::BgenStreamer> bgenStreamer;
    if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "mainMarkerInCPP: BGEN object not initialized but t_genoType=='bgen'.");
        }
        int nDec = (g_bgenDecoders > 0 ? g_bgenDecoders : 4);
        bgenStreamer.reset(new BGEN::BgenStreamer(
            ptr_gBGENobj, t_genoIndex, t_isImputation, nDec, /*queueCap*/ 64));
        std::cout << "BGEN streamer: " << nDec << " decoders, queueCap=64" << std::endl;
    }

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < q; i++) {
        // Phase B: per-iteration scratch — declared inside the loop body so
        // each thread has its own copies (no false sharing, no races).
        arma::vec t_GVec(n);
        arma::vec gtildeVec(n);
        arma::vec t_P2Vec;
        std::vector<uint> indexZeroVec;
        std::vector<uint> indexNonZeroVec;
        std::vector<uint> indexForMissing;
        // BUG 1 fix: per-thread conditional-analysis scratch (was previously a
        // single shared instance above this loop; see comment near declaration
        // of g_numMarker_cond_local). getMarkerPval writes into it via
        // t_G1tilde_P_G2tilde = sqrt(varRatio)*gtilde.t()*m_P2Mat_cond.
        arma::rowvec G1tilde_P_G2tilde_Vec(g_numMarker_cond_local);

        // Skip if another thread already hit end-of-stream at an earlier index.
        if (i >= firstEndIdx) continue;

        if ((i + 1) % g_marker_chunksize == 0) {
            #pragma omp critical(progress)
            {
                std::cout << "Completed " << (i + 1) << "/" << q
                          << " markers in the chunk." << std::endl;
            }
        }

        // information of marker
        double altFreq, altCounts, missingRate, imputeInfo;
        double AF_case, AF_ctrl, N_case_hom, N_ctrl_het, N_case_het, N_ctrl_hom;
        std::string chr, ref, alt, marker;
        uint32_t pd, N_case, N_ctrl, N;

        bool flip = false;

        bool isOutputIndexForMissing = true;
        bool isOnlyOutputNonZero = false;

        // clear vectors
        indexZeroVec.clear();
        indexNonZeroVec.clear();
        indexForMissing.clear();

        bool isReadMarker;
        if (t_genoType == "bgen") {
            // Phase D: BGEN streamer is internally thread-safe (single reader
            // thread + N decoder threads + bounded queue). No critical needed.
            BGEN::BgenDecodedMarker dm;
            isReadMarker = bgenStreamer->getNext(dm);
            if (isReadMarker) {
                ref         = dm.alleles.size() > 0 ? dm.alleles[0] : "";
                alt         = dm.alleles.size() > 1 ? dm.alleles[1] : "";
                marker      = dm.rsID;
                pd          = dm.physpos;
                chr         = dm.chr;
                altFreq     = dm.altFreq;
                altCounts   = dm.altCounts;
                missingRate = dm.missingRate;
                imputeInfo  = dm.info;
                indexForMissing = std::move(dm.indexForMissing);
                if (isOnlyOutputNonZero) {
                    indexNonZeroVec = std::move(dm.indexForNonZero);
                }
                t_GVec = std::move(dm.dosages);
            }
        } else {
            std::string t_genoIndex_str = t_genoIndex.at(i);
            char* end;
            uint64_t gIndex = std::strtoull(t_genoIndex_str.c_str(), &end, 10);

            // For absolute-seek readers (PLINK, PGEN, VCF), gIndex_prev=0
            // forces SEEK_SET inside the reader — required for correctness when
            // iterations execute out of order under OpenMP.
            uint64_t gIndex_prev = 0;

            // Phase B: PLINK/PGEN/VCF readers are not thread-safe. Serialize
            // disk I/O via critical(genoread); parallelize compute below.
            #pragma omp critical(genoread)
            {
                isReadMarker = Unified_getOneMarker(
                    t_genoType, gIndex_prev, gIndex,
                    ref, alt, marker, pd, chr,
                    altFreq, altCounts, missingRate, imputeInfo,
                    isOutputIndexForMissing,
                    indexForMissing,
                    isOnlyOutputNonZero,
                    indexNonZeroVec, t_GVec, t_isImputation);
            }
        }

        if (!isReadMarker) {
            #pragma omp critical(endflag)
            {
                if (i < firstEndIdx) firstEndIdx = i;
                g_markerTestEnd = true;
            }
            // pvalVec[i] already initialized to "NA" sentinel; just stop work.
            continue;
        }

        std::string pds = std::to_string(pd);
        std::string info = chr + ":" + pds + ":" + ref + ":" + alt;

        chrVec.at(i) = chr;
        posVec.at(i) = pds;
        refVec.at(i) = ref;
        altVec.at(i) = alt;
        markerVec.at(i) = marker;
        infoVec.at(i) = info;
        altFreqVec.at(i) = altFreq;
        missingRateVec.at(i) = missingRate;
        imputationInfoVec.at(i) = imputeInfo;

        // MAF and MAC are for Quality Control (QC)
        double MAF = std::min(altFreq, 1 - altFreq);
        double MAC = MAF * n * (1 - missingRate) * 2;

        // Quality Control (QC) based on missing rate, MAF, and MAC
        if ((missingRate > g_missingRate_cutoff) ||
            (MAF < g_marker_minMAF_cutoff) ||
            (MAC < g_marker_minMAC_cutoff) ||
            (imputeInfo < g_marker_minINFO_cutoff)) {
            continue;
        } else {
            indexZeroVec.clear();
            indexNonZeroVec.clear();

            flip = imputeGenoAndFlip(
                t_GVec, altFreq, altCounts,
                indexForMissing, g_impute_method,
                g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff,
                MAC, indexZeroVec, indexNonZeroVec);

            // === CHECKPOINT OUTPUT BEGIN ===
            if (g_writeCheckpoints && i == 0 && !g_checkpointDir.empty()) {
                std::string cpFile = g_checkpointDir + "/ckpt_09_marker0_imputed.txt";
                std::ofstream ofs(cpFile);
                if (ofs.is_open()) {
                    ofs << std::setprecision(15);
                    ofs << "field\tvalue" << std::endl;
                    ofs << "altFreq\t" << altFreq << std::endl;
                    ofs << "altCounts\t" << altCounts << std::endl;
                    ofs << "MAC\t" << MAC << std::endl;
                    ofs << "MAF\t" << MAF << std::endl;
                    ofs << "flip\t" << flip << std::endl;
                    ofs.close();
                    std::cout << "[CHECKPOINT] Wrote ckpt_09_marker0_imputed.txt to " << g_checkpointDir << std::endl;
                }
            }
            // === CHECKPOINT OUTPUT END ===

            MAC = std::min(altCounts, 2.0 * n - altCounts);
            MAF = std::min(altFreq, 1 - altFreq);

            if ((MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff)) {
                continue;
            } else {
                altFreqVec.at(i) = altFreq;
                altCountsVec.at(i) = altCounts;

                // analysis results for single-marker
                double Beta, seBeta, Tstat, varT, gy;
                double Beta_c, seBeta_c, Tstat_c, varT_c;
                std::string pval, pval_noSPA, pval_c, pval_noSPA_c;
                bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;

                arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
                indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);
                indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

                indexZeroVec.clear();
                indexNonZeroVec.clear();
                t_P2Vec.clear();
                G1tilde_P_G2tilde_Vec.clear();

                // Phase B: build per-marker ctx (no mutation of SAIGEClass).
                // Previously this block called set_flagSparseGRM_cur(),
                // assignSingleVarianceRatio(), and assignVarianceRatio() which
                // wrote to ptr_gSAIGEobj fields — those are data races under
                // OpenMP. ctx now carries the per-marker scalars directly.
                SAIGE::PerMarkerCtx ctx_first;
                ctx_first.flagSparseGRM_cur = ptr_gSAIGEobj->m_isFastTest
                    ? false : ptr_gSAIGEobj->m_flagSparseGRM;
                ctx_first.isnoadjCov_cur = ptr_gSAIGEobj->m_isnoadjCov;
                {
                    bool dummyHas;
                    if (isSingleVarianceRatio) {
                        ctx_first.varRatioVal = ptr_gSAIGEobj->computeSingleVarianceRatio(
                            ctx_first.flagSparseGRM_cur, ctx_first.isnoadjCov_cur);
                    } else {
                        ctx_first.varRatioVal = ptr_gSAIGEobj->computeVarianceRatio(
                            MAC, ctx_first.flagSparseGRM_cur, ctx_first.isnoadjCov_cur,
                            dummyHas);
                    }
                }

                // === CHECKPOINT OUTPUT BEGIN ===
                if (g_writeCheckpoints && i == 0 && !g_checkpointDir.empty()) {
                    std::string cpFile = g_checkpointDir + "/ckpt_10_marker0_vr.txt";
                    std::ofstream ofs(cpFile);
                    if (ofs.is_open()) {
                        ofs << std::setprecision(15);
                        ofs << "field\tvalue" << std::endl;
                        ofs << "varRatioVal\t" << ctx_first.varRatioVal << std::endl;
                        ofs << "isSingleVarianceRatio\t" << isSingleVarianceRatio << std::endl;
                        ofs << "isFastTest\t" << ptr_gSAIGEobj->m_isFastTest << std::endl;
                        ofs.close();
                        std::cout << "[CHECKPOINT] Wrote ckpt_10_marker0_vr.txt to " << g_checkpointDir << std::endl;
                    }
                }
                // === CHECKPOINT OUTPUT END ===

                bool is_region = false;

                if (MAC <= g_MACCutoffforER && t_traitType == "binary") {
                    Unified_getMarkerPval(
                        t_GVec,
                        false,  // bool t_isOnlyOutputNonZero
                        indexNonZeroVec_arma, indexZeroVec_arma,
                        Beta, seBeta, pval, pval_noSPA,
                        Tstat, gy, varT,
                        altFreq, isSPAConverge,
                        gtildeVec, is_gtilde,
                        is_region, t_P2Vec,
                        isCondition,
                        Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                        Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                        is_Firth, is_FirthConverge,
                        true,  // t_isER
                        ctx_first.isnoadjCov_cur,
                        ctx_first.flagSparseGRM_cur,
                        ctx_first);
                } else {
                    Unified_getMarkerPval(
                        t_GVec,
                        false,  // bool t_isOnlyOutputNonZero
                        indexNonZeroVec_arma, indexZeroVec_arma,
                        Beta, seBeta, pval, pval_noSPA,
                        Tstat, gy, varT,
                        altFreq, isSPAConverge,
                        gtildeVec, is_gtilde,
                        is_region, t_P2Vec,
                        isCondition,
                        Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                        Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                        is_Firth, is_FirthConverge,
                        false,  // t_isER
                        ctx_first.isnoadjCov_cur,
                        ctx_first.flagSparseGRM_cur,
                        ctx_first);
                }

                double pval_num;
                try {
                    pval_num = std::stod(pval);
                } catch (const std::invalid_argument&) {
                    std::cerr << "Argument is invalid\n";
                    pval_num = 0;
                } catch (const std::out_of_range&) {
                    std::cerr << "Argument is out of range for a double\n";
                    pval_num = 0;
                }

                // Fast test re-evaluation (exact match of SAIGE logic)
                if ((t_traitType == "binary" && MAC > g_MACCutoffforER) ||
                    t_traitType != "binary") {

                    if (ptr_gSAIGEobj->m_isFastTest &&
                        pval_num < (ptr_gSAIGEobj->m_pval_cutoff_for_fastTest)) {
                        // Phase B: build ctx for the fast-test re-eval. The
                        // previous block called set_flagSparseGRM_cur(),
                        // set_isnoadjCov_cur(false), and assignVarianceRatio*
                        // on the shared SAIGEClass — those are races under OMP.
                        SAIGE::PerMarkerCtx ctx_fast;
                        ctx_fast.flagSparseGRM_cur =
                            (MAC > ptr_gSAIGEobj->m_cateVarRatioMinMACVecExclude.back())
                                ? false : ptr_gSAIGEobj->m_flagSparseGRM;
                        ctx_fast.isnoadjCov_cur = false;
                        {
                            bool dummyHas;
                            if (!isSingleVarianceRatio) {
                                ctx_fast.varRatioVal = ptr_gSAIGEobj->computeVarianceRatio(
                                    MAC, ctx_fast.flagSparseGRM_cur,
                                    ctx_fast.isnoadjCov_cur, dummyHas);
                            } else {
                                ctx_fast.varRatioVal = ptr_gSAIGEobj->computeSingleVarianceRatio(
                                    ctx_fast.flagSparseGRM_cur, ctx_fast.isnoadjCov_cur);
                            }
                        }

                        Unified_getMarkerPval(
                            t_GVec,
                            false,
                            indexNonZeroVec_arma, indexZeroVec_arma,
                            Beta, seBeta, pval, pval_noSPA,
                            Tstat, gy, varT,
                            altFreq, isSPAConverge,
                            gtildeVec, is_gtilde,
                            is_region, t_P2Vec,
                            isCondition,
                            Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                            Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                            is_Firth, is_FirthConverge,
                            false,
                            ctx_fast.isnoadjCov_cur,
                            ctx_fast.flagSparseGRM_cur,
                            ctx_fast);
                    }
                }  // if((t_traitType == "binary" && MAC > g_MACCutoffforER) || t_traitType != "binary")

                if (t_traitType == "binary") {
                    if (is_Firth) {
                        #pragma omp atomic
                        mFirth = mFirth + 1;
                        if (is_FirthConverge) {
                            #pragma omp atomic
                            mFirthConverge = mFirthConverge + 1;
                        }
                    }
                }

                indexNonZeroVec_arma.clear();
                indexZeroVec_arma.clear();

                BetaVec.at(i) = Beta * (1 - 2 * flip);
                seBetaVec.at(i) = seBeta;
                pvalVec.at(i) = pval;
                pvalNAVec.at(i) = pval_noSPA;
                TstatVec.at(i) = Tstat * (1 - 2 * flip);
                varTVec.at(i) = varT;

                if (isCondition) {
                    Beta_cVec.at(i) = Beta_c * (1 - 2 * flip);
                    seBeta_cVec.at(i) = seBeta_c;
                    pval_cVec.at(i) = pval_c;
                    pvalNA_cVec.at(i) = pval_noSPA_c;
                    Tstat_cVec.at(i) = Tstat_c * (1 - 2 * flip);
                    varT_cVec.at(i) = varT_c;
                }

                if (t_traitType == "binary" || t_traitType == "survival") {
                    arma::vec dosage_case = t_GVec.elem(ptr_gSAIGEobj->m_case_indices);
                    arma::vec dosage_ctrl = t_GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                    AF_case = arma::mean(dosage_case) / 2;
                    AF_ctrl = arma::mean(dosage_ctrl) / 2;
                    N_case = dosage_case.n_elem;
                    N_ctrl = dosage_ctrl.n_elem;
                    if (flip) {
                        AF_case = 1 - AF_case;
                        AF_ctrl = 1 - AF_ctrl;
                    }
                    isSPAConvergeVec.at(i) = isSPAConverge;
                    AF_caseVec.at(i) = AF_case;
                    AF_ctrlVec.at(i) = AF_ctrl;

                    N_caseVec.at(i) = N_case;
                    N_ctrlVec.at(i) = N_ctrl;

                    arma::uvec N_case_ctrl_het_hom0;
                    if (t_isMoreOutput) {
                        N_case_ctrl_het_hom0 = arma::find(dosage_case <= 2 && dosage_case >= 1.5);
                        N_case_homVec.at(i) = N_case_ctrl_het_hom0.n_elem;
                        N_case_ctrl_het_hom0 = arma::find(dosage_case < 1.5 && dosage_case >= 0.5);
                        N_case_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
                        N_case_ctrl_het_hom0 = arma::find(dosage_ctrl <= 2 && dosage_ctrl >= 1.5);
                        N_ctrl_homVec.at(i) = N_case_ctrl_het_hom0.n_elem;
                        N_case_ctrl_het_hom0 = arma::find(dosage_ctrl < 1.5 && dosage_ctrl >= 0.5);
                        N_ctrl_hetVec.at(i) = N_case_ctrl_het_hom0.n_elem;
                        if (flip) {
                            N_case_homVec.at(i) = N_case - N_case_hetVec.at(i) - N_case_homVec.at(i);
                            N_ctrl_homVec.at(i) = N_ctrl - N_ctrl_hetVec.at(i) - N_ctrl_homVec.at(i);
                        }
                    }
                } else if (t_traitType == "quantitative") {
                    N_Vec.at(i) = n;
                }

                // === CHECKPOINT OUTPUT BEGIN ===
                if (g_writeCheckpoints && i == 0 && !g_checkpointDir.empty()) {
                    // Write first marker score test results
                    std::string cpFile = g_checkpointDir + "/ckpt_11_marker0_pval.txt";
                    std::ofstream ofs(cpFile);
                    if (ofs.is_open()) {
                        ofs << std::setprecision(15);
                        ofs << "field\tvalue" << std::endl;
                        ofs << "marker\t" << markerVec.at(i) << std::endl;
                        ofs << "chr\t" << chrVec.at(i) << std::endl;
                        ofs << "pos\t" << posVec.at(i) << std::endl;
                        ofs << "ref\t" << refVec.at(i) << std::endl;
                        ofs << "alt\t" << altVec.at(i) << std::endl;
                        ofs << "altFreq\t" << altFreqVec.at(i) << std::endl;
                        ofs << "altCounts\t" << altCountsVec.at(i) << std::endl;
                        ofs << "Beta\t" << BetaVec.at(i) << std::endl;
                        ofs << "seBeta\t" << seBetaVec.at(i) << std::endl;
                        ofs << "Tstat\t" << TstatVec.at(i) << std::endl;
                        ofs << "varT\t" << varTVec.at(i) << std::endl;
                        ofs << "pval\t" << pvalVec.at(i) << std::endl;
                        ofs << "pval_noSPA\t" << pvalNAVec.at(i) << std::endl;
                        ofs << "isSPAConverge\t" << std::boolalpha << isSPAConverge << std::endl;
                        ofs << "flip\t" << flip << std::endl;
                        ofs << "MAC\t" << MAC << std::endl;
                        ofs << "MAF\t" << MAF << std::endl;
                        ofs.close();
                        std::cout << "[CHECKPOINT] Wrote ckpt_11_marker0_pval.txt to " << g_checkpointDir << std::endl;
                    }
                }
                // === CHECKPOINT OUTPUT END ===

            }  // if((MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff))
        }  // if((missingRate > g_missingRate_cutoff) || ...)
    }  // for(int i = 0; i < q; i++)

    // output
    writeOutfile_single(t_isMoreOutput,
                         t_isImputation,
                         isCondition,
                         t_isFirth,
                         mFirth,
                         mFirthConverge,
                         t_traitType,
                         chrVec,
                         posVec,
                         markerVec,
                         refVec,
                         altVec,
                         altCountsVec,
                         altFreqVec,
                         imputationInfoVec,
                         missingRateVec,
                         BetaVec,
                         seBetaVec,
                         TstatVec,
                         varTVec,
                         pvalVec,
                         pvalNAVec,
                         isSPAConvergeVec,
                         Beta_cVec,
                         seBeta_cVec,
                         Tstat_cVec,
                         varT_cVec,
                         pval_cVec,
                         pvalNA_cVec,
                         AF_caseVec,
                         AF_ctrlVec,
                         N_caseVec,
                         N_ctrlVec,
                         N_case_homVec,
                         N_ctrl_hetVec,
                         N_case_hetVec,
                         N_ctrl_homVec,
                         N_Vec);
}


// ============================================================
// setRegion_GlobalVarsInCPP
// Direct port from SAIGE/src/Main.cpp lines 200-211
// ============================================================
void setRegion_GlobalVarsInCPP(
    arma::vec t_max_maf_region,
    unsigned int t_max_markers_region,
    double t_MACCutoff_to_CollapseUltraRare,
    double t_min_gourpmac_for_burdenonly)
{
    g_region_maxMAF_cutoff = t_max_maf_region;
    g_maxMAFLimit = g_region_maxMAF_cutoff.max();
    g_region_maxMarkers_cutoff = t_max_markers_region;
    g_region_minMAC_cutoff = t_MACCutoff_to_CollapseUltraRare;
    g_min_gourpmac_for_burdenonly = t_min_gourpmac_for_burdenonly;
}


// ============================================================
// openOutfile (region/gene output)
// Direct port from SAIGE/src/Main.cpp lines 2503-2528
// ============================================================
bool openOutfile(std::string t_traitType, bool isappend) {
    bool isopen;
    if (!isappend) {
        OutFile.open(g_outputFilePrefixGroup.c_str());
        isopen = OutFile.is_open();
        if (isopen) {
            OutFile << std::setprecision(15);
            OutFile << "Region\tGroup\tmax_MAF\tPvalue_Burden\tBETA_Burden\tSE_Burden\t";
            if (ptr_gSAIGEobj->m_isCondition) {
                OutFile << "Pvalue_Burden_c\tBeta_Burden_c\tseBeta_Burden_c\t";
            }
            OutFile << "MAC\t";
            if (t_traitType == "binary") {
                OutFile << "MAC_case\tMAC_control\t";
            }
            if (t_traitType == "survival") {
                OutFile << "MAC_event\tMAC_censor\t";
            }
            OutFile << "Number_rare\tNumber_ultra_rare\n";
        }
    } else {
        OutFile.open(g_outputFilePrefixGroup.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile.is_open();
        if (isopen) {
            OutFile << std::setprecision(15);
        }
    }
    return isopen;
}


// ============================================================
// openOutfile_SKATO (region/gene output for SKAT-O/SKAT tests)
// Writes header with Pvalue, Pvalue_Burden, Pvalue_SKAT columns
// ============================================================
bool openOutfile_SKATO(std::string t_traitType, bool isappend) {
    bool isopen;
    if (!isappend) {
        OutFile.open(g_outputFilePrefixGroup.c_str());
        isopen = OutFile.is_open();
        if (isopen) {
            OutFile << std::setprecision(15);
            OutFile << "Region\tGroup\tmax_MAF\tPvalue\tPvalue_Burden\tPvalue_SKAT\tBETA_Burden\tSE_Burden\t";
            if (ptr_gSAIGEobj->m_isCondition) {
                OutFile << "Pvalue_cond\tPvalue_Burden_cond\tPvalue_SKAT_cond\tBETA_Burden_cond\tSE_Burden_cond\t";
            }
            OutFile << "MAC\t";
            if (t_traitType == "binary") {
                OutFile << "MAC_case\tMAC_control\t";
            }
            if (t_traitType == "survival") {
                OutFile << "MAC_event\tMAC_censor\t";
            }
            OutFile << "Number_rare\tNumber_ultra_rare\n";
        }
    } else {
        OutFile.open(g_outputFilePrefixGroup.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile.is_open();
        if (isopen) {
            OutFile << std::setprecision(15);
        }
    }
    return isopen;
}


// ============================================================
// openOutfile_singleinGroup
// Direct port from SAIGE/src/Main.cpp lines 2531-2583
// ============================================================
bool openOutfile_singleinGroup(std::string t_traitType, bool t_isImputation,
                                bool isappend, bool t_isMoreOutput) {
    bool isopen;
    if (!isappend) {
        OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str());
        isopen = OutFile_singleInGroup.is_open();
        if (isopen) {
            OutFile_singleInGroup << "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
            if (t_isImputation) {
                OutFile_singleInGroup << "imputationInfo\t";
            } else {
                OutFile_singleInGroup << "MissingRate\t";
            }
            OutFile_singleInGroup << "BETA\tSE\tTstat\tvar\tp.value\t";
            if (t_traitType == "binary" || t_traitType == "survival") {
                OutFile_singleInGroup << "p.value.NA\tIs.SPA\t";
            }
            if (ptr_gSAIGEobj->m_isCondition) {
                OutFile_singleInGroup << "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile_singleInGroup << "p.value.NA_c\t";
                }
            }
            if (t_traitType == "binary") {
                OutFile_singleInGroup << "AF_case\tAF_ctrl\tN_case\tN_ctrl";
                if (t_isMoreOutput) {
                    OutFile_singleInGroup << "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
                }
                OutFile_singleInGroup << "\n";
            } else if (t_traitType == "quantitative") {
                OutFile_singleInGroup << "N\n";
            } else if (t_traitType == "survival") {
                OutFile_singleInGroup << "AF_event\tAF_censor\tN_event\tN_censor";
                if (t_isMoreOutput) {
                    OutFile_singleInGroup << "\tN_event_hom\tN_event_het\tN_censor_hom\tN_censor_het";
                }
                OutFile_singleInGroup << "\n";
            }
        }
    } else {
        OutFile_singleInGroup.open(g_outputFilePrefixSingleInGroup.c_str(),
                                    std::ofstream::out | std::ofstream::app);
        isopen = OutFile_singleInGroup.is_open();
    }
    return isopen;
}


// ============================================================
// writeOutfile_BURDEN
// Direct port from SAIGE/src/Main.cpp lines 2816-2891
// ============================================================
void writeOutfile_BURDEN(std::string regionName,
                          std::vector<std::string>& BURDEN_AnnoName_Vec,
                          std::vector<std::string>& BURDEN_maxMAFName_Vec,
                          std::vector<std::string>& BURDEN_pval_Vec,
                          std::vector<double>& BURDEN_Beta_Vec,
                          std::vector<double>& BURDEN_seBeta_Vec,
                          std::vector<std::string>& BURDEN_pval_cVec,
                          std::vector<double>& BURDEN_Beta_cVec,
                          std::vector<double>& BURDEN_seBeta_cVec,
                          arma::vec& MAC_GroupVec,
                          arma::vec& MACCase_GroupVec,
                          arma::vec& MACControl_GroupVec,
                          arma::vec& NumRare_GroupVec,
                          arma::vec& NumUltraRare_GroupVec,
                          double cctpval,
                          double cctpval_cond,
                          unsigned int q_anno,
                          unsigned int q_maf,
                          bool isCondition,
                          std::string t_traitType) {
    unsigned int i;
    for (unsigned int j = 0; j < q_anno; j++) {
        for (unsigned int m = 0; m < q_maf; m++) {
            i = j * q_maf + m;
            if (BURDEN_pval_Vec.at(i) != "NA") {
                OutFile << regionName << "\t";
                OutFile << BURDEN_AnnoName_Vec.at(i) << "\t";
                OutFile << BURDEN_maxMAFName_Vec.at(i) << "\t";
                OutFile << BURDEN_pval_Vec.at(i) << "\t";
                OutFile << BURDEN_Beta_Vec.at(i) << "\t";
                OutFile << BURDEN_seBeta_Vec.at(i) << "\t";
                if (isCondition) {
                    OutFile << BURDEN_pval_cVec.at(i) << "\t";
                    OutFile << BURDEN_Beta_cVec.at(i) << "\t";
                    OutFile << BURDEN_seBeta_cVec.at(i) << "\t";
                }
                OutFile << MAC_GroupVec(i) << "\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile << MACCase_GroupVec(i) << "\t";
                    OutFile << MACControl_GroupVec(i) << "\t";
                }
                OutFile << NumRare_GroupVec(i) << "\t";
                OutFile << NumUltraRare_GroupVec(i) << "\n";
            }
        }
    }
    OutFile << regionName << "\tCauchy\tNA\t";
    OutFile << cctpval << "\tNA\tNA\t";
    if (isCondition) {
        OutFile << cctpval_cond << "\tNA\tNA\t";
    }
    OutFile << "NA\t";
    if (t_traitType == "binary" || t_traitType == "survival") {
        OutFile << "NA\t";
        OutFile << "NA\t";
    }
    OutFile << "NA\t";
    OutFile << "NA\n";
}


// ============================================================
// writeOutfile_singleInGroup
// Direct port from SAIGE/src/Main.cpp lines 2894-3026
// ============================================================
int writeOutfile_singleInGroup(bool t_isMoreOutput,
                                bool t_isImputation,
                                bool t_isCondition,
                                bool t_isFirth,
                                int mFirth,
                                int mFirthConverge,
                                std::string t_traitType,
                                std::vector<std::string>& chrVec,
                                std::vector<std::string>& posVec,
                                std::vector<std::string>& markerVec,
                                std::vector<std::string>& refVec,
                                std::vector<std::string>& altVec,
                                std::vector<double>& altCountsVec,
                                std::vector<double>& altFreqVec,
                                std::vector<double>& imputationInfoVec,
                                std::vector<double>& missingRateVec,
                                std::vector<double>& BetaVec,
                                std::vector<double>& seBetaVec,
                                std::vector<double>& TstatVec,
                                std::vector<double>& varTVec,
                                std::vector<std::string>& pvalVec,
                                std::vector<std::string>& pvalNAVec,
                                std::vector<bool>& isSPAConvergeVec,
                                std::vector<double>& Beta_cVec,
                                std::vector<double>& seBeta_cVec,
                                std::vector<double>& Tstat_cVec,
                                std::vector<double>& varT_cVec,
                                std::vector<std::string>& pval_cVec,
                                std::vector<std::string>& pvalNA_cVec,
                                std::vector<double>& AF_caseVec,
                                std::vector<double>& AF_ctrlVec,
                                std::vector<uint32_t>& N_caseVec,
                                std::vector<uint32_t>& N_ctrlVec,
                                std::vector<double>& N_case_homVec,
                                std::vector<double>& N_ctrl_hetVec,
                                std::vector<double>& N_case_hetVec,
                                std::vector<double>& N_ctrl_homVec,
                                std::vector<uint32_t>& N_Vec,
                                std::ofstream& t_OutFile_singleInGroup) {
    int numofUR = 0;
    for (unsigned int k = 0; k < pvalVec.size(); k++) {
        if (pvalVec.at(k) != "NA") {
            if (chrVec.at(k) == "UR") {
                numofUR = numofUR + 1;
            }
            t_OutFile_singleInGroup << chrVec.at(k) << "\t";
            t_OutFile_singleInGroup << posVec.at(k) << "\t";
            t_OutFile_singleInGroup << markerVec.at(k) << "\t";
            t_OutFile_singleInGroup << refVec.at(k) << "\t";
            t_OutFile_singleInGroup << altVec.at(k) << "\t";
            t_OutFile_singleInGroup << altCountsVec.at(k) << "\t";
            t_OutFile_singleInGroup << altFreqVec.at(k) << "\t";
            if (t_isImputation) {
                t_OutFile_singleInGroup << imputationInfoVec.at(k) << "\t";
            } else {
                t_OutFile_singleInGroup << missingRateVec.at(k) << "\t";
            }
            t_OutFile_singleInGroup << BetaVec.at(k) << "\t";
            t_OutFile_singleInGroup << seBetaVec.at(k) << "\t";
            t_OutFile_singleInGroup << TstatVec.at(k) << "\t";
            t_OutFile_singleInGroup << varTVec.at(k) << "\t";
            t_OutFile_singleInGroup << pvalVec.at(k) << "\t";
            if (t_traitType == "binary" || t_traitType == "survival") {
                t_OutFile_singleInGroup << pvalNAVec.at(k) << "\t";
                t_OutFile_singleInGroup << std::boolalpha << isSPAConvergeVec.at(k) << "\t";
            }
            if (t_isCondition) {
                t_OutFile_singleInGroup << Beta_cVec.at(k) << "\t";
                t_OutFile_singleInGroup << seBeta_cVec.at(k) << "\t";
                t_OutFile_singleInGroup << Tstat_cVec.at(k) << "\t";
                t_OutFile_singleInGroup << varT_cVec.at(k) << "\t";
                t_OutFile_singleInGroup << pval_cVec.at(k) << "\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    t_OutFile_singleInGroup << pvalNA_cVec.at(k) << "\t";
                }
            }
            if (t_traitType == "binary" || t_traitType == "survival") {
                t_OutFile_singleInGroup << AF_caseVec.at(k) << "\t";
                t_OutFile_singleInGroup << AF_ctrlVec.at(k) << "\t";
                t_OutFile_singleInGroup << N_caseVec.at(k) << "\t";
                t_OutFile_singleInGroup << N_ctrlVec.at(k);
                if (t_isMoreOutput) {
                    t_OutFile_singleInGroup << "\t" << N_case_homVec.at(k);
                    t_OutFile_singleInGroup << "\t" << N_case_hetVec.at(k);
                    t_OutFile_singleInGroup << "\t" << N_ctrl_homVec.at(k);
                    t_OutFile_singleInGroup << "\t" << N_ctrl_hetVec.at(k);
                }
                t_OutFile_singleInGroup << "\n";
            } else if (t_traitType == "quantitative") {
                t_OutFile_singleInGroup << N_Vec.at(k) << "\n";
            }
        }
    }
    return numofUR;
}


// ============================================================
// convert_str_to_log
// Port from SAIGE/R/Util.R: convert_str_to_log()
// Converts a p-value string in scientific notation ("1.23E-5")
// to log(p). Returns NaN if input is "NA".
// ============================================================
static double convert_str_to_log(const std::string& a) {
    if (a == "NA" || a.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    // Try direct conversion first (handles both "1.23E-5" and "0.05")
    try {
        double val = std::stod(a);
        if (val <= 0.0) {
            return -std::numeric_limits<double>::infinity();
        }
        return std::log(val);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}


// ============================================================
// SPA_ER_kernel_related_Phiadj_fast_new
// Port from SAIGE/R/SAIGE_SPATest_Region_Func.R lines 190-272
// Adjusts the variance-covariance matrix Phi using SPA results.
//
// p_new: vector of log(p-values) for each variant (from SPA)
// Score: weighted score vector
// Phi: weighted variance-covariance matrix
// p_value_burden: log(p-value) of the burden test (from SPA)
// regionTestType: "SKATO", "SKAT", or "BURDEN"
//
// Returns: adjusted Phi and scale factors
// ============================================================
struct PhiAdjResult {
    arma::mat Phi_adj;          // Adjusted variance-covariance matrix
    arma::vec scaleFactor;      // Scale factors per variant
};

static PhiAdjResult SPA_ER_kernel_related_Phiadj_fast_new(
    const arma::vec& p_new,     // log(p-values)
    const arma::vec& Score,
    const arma::mat& Phi,
    double p_value_burden,      // log(p-value) of burden
    const std::string& regionTestType)
{
    int p_m = (int)Score.n_elem;
    arma::vec VarS_org = Phi.diag();
    arma::vec VarS = Score % Score / 500.0;  // default, will be overwritten

    // Find indices where p_new is not NaN
    boost::math::chi_squared chi2_1(1.0);
    for (int i = 0; i < p_m; i++) {
        if (std::isfinite(p_new(i))) {
            // R: VarS[idx_p0] = zscore.all_0[idx_p0]^2 / qchisq(p.new[idx_p0], df=1, lower.tail=FALSE, log.p=TRUE)
            // p_new(i) is log(p), so exp(p_new(i)) is the actual p-value
            double actual_p = std::exp(p_new(i));
            if (actual_p > 0.0 && actual_p < 1.0) {
                try {
                    double qchi = boost::math::quantile(boost::math::complement(chi2_1, actual_p));
                    if (qchi > 0.0) {
                        VarS(i) = Score(i) * Score(i) / qchi;
                    }
                } catch (...) {
                    // keep default
                }
            }
        }
    }

    // Handle Inf values
    for (int i = 0; i < p_m; i++) {
        if (!std::isfinite(VarS(i))) {
            VarS(i) = 0.0;
        }
    }

    // Scale factor
    arma::vec VarStoorg(p_m);
    for (int i = 0; i < p_m; i++) {
        if (VarS_org(i) > 0.0) {
            VarStoorg(i) = VarS(i) / VarS_org(i);
        } else {
            VarStoorg(i) = 0.0;
        }
    }
    arma::vec scaleFactor = arma::sqrt(VarStoorg);

    // G2_adj_n = t(t(Phi * sqrt(VarStoorg)) * sqrt(VarStoorg))
    // = outer(sqrt(VarStoorg), sqrt(VarStoorg)) .* Phi
    arma::mat G2_adj_n(p_m, p_m);
    for (int i = 0; i < p_m; i++) {
        for (int j = 0; j < p_m; j++) {
            G2_adj_n(i, j) = Phi(i, j) * scaleFactor(i) * scaleFactor(j);
        }
    }

    // Burden adjustment
    double VarQ = arma::accu(G2_adj_n);
    double Q_b = arma::sum(Score) * arma::sum(Score);

    double VarQ_2;
    if (std::isfinite(p_value_burden)) {
        double actual_p_burden = std::exp(p_value_burden);
        if (actual_p_burden > 0.0 && actual_p_burden < 1.0) {
            try {
                double qchi_burden = boost::math::quantile(boost::math::complement(chi2_1, actual_p_burden));
                VarQ_2 = (qchi_burden > 0.0) ? Q_b / qchi_burden : 0.0;
            } catch (...) {
                VarQ_2 = 0.0;
            }
        } else {
            VarQ_2 = 0.0;
        }
    } else {
        VarQ_2 = 0.0;
    }

    double r;
    if (VarQ_2 == 0.0) {
        r = 1.0;
    } else {
        r = VarQ / VarQ_2;
    }
    r = std::min(r, 1.0);

    // Phi_ccadj = G2_adj_n / r
    arma::mat Phi_ccadj = G2_adj_n / r;
    scaleFactor = scaleFactor / std::sqrt(r);

    PhiAdjResult result;
    result.Phi_adj = Phi_ccadj;
    result.scaleFactor = scaleFactor;
    return result;
}


// ============================================================
// get_newPhi_scaleFactor_traitType
// Port from SAIGE/R/SAIGE_SPATest_Region_Func.R lines 295-315
// Computes SPA-adjusted Phi for binary/survival traits.
// ============================================================
static PhiAdjResult get_newPhi_scaleFactor_traitType(
    double q_sum,               // sum of gy * AnnoWeights for the group
    const arma::vec& mu_a,      // fitted values from null model (N-vector)
    const arma::vec& g_sum,     // weighted genotype sum column (N-vector)
    const arma::vec& p_new,     // log(p-values) for variants in group
    const arma::vec& Score,     // weighted score vector
    const arma::mat& Phi,       // weighted variance-covariance matrix
    const std::string& regionTestType,
    const std::string& traitType)
{
    // m1 = sum(mu * g.sum)
    double m1 = arma::dot(mu_a, g_sum);

    double var1;
    double qinv;
    if (traitType == "binary") {
        // var1 = sum(mu*(1-mu)*g.sum^2)
        var1 = arma::dot(mu_a % (1.0 - mu_a), g_sum % g_sum);
        qinv = -((q_sum - m1) > 0 ? 1.0 : ((q_sum - m1) < 0 ? -1.0 : 0.0))
               * std::abs(q_sum - m1) + m1;
    } else {
        // survival: var1 = sum(mu * g.sum^2)
        var1 = arma::dot(mu_a, g_sum % g_sum);
        double q_diff = q_sum - m1;
        qinv = -q_diff + m1; // = m1 - (q_sum - m1) = 2*m1 - q_sum
    }

    // pval_noadj = pchisq((q.sum - m1)^2/var1, df=1, lower.tail=FALSE, log.p=TRUE)
    double chi2_stat = (q_sum - m1) * (q_sum - m1) / var1;
    double pval_noadj;
    try {
        boost::math::chi_squared chi2_1(1.0);
        double actual_p = boost::math::cdf(boost::math::complement(chi2_1, chi2_stat));
        pval_noadj = (actual_p > 0.0) ? std::log(actual_p) : -std::numeric_limits<double>::infinity();
    } catch (...) {
        pval_noadj = -std::numeric_limits<double>::infinity();
    }

    bool isSPAConverge = true;
    double p_value_burden;

    // if abs(q.sum - m1)/sqrt(var1) < 2: use normal approx
    if (std::abs(q_sum - m1) / std::sqrt(var1) < 2.0) {
        p_value_burden = pval_noadj;
    } else {
        // SPA correction
        arma::vec mu_a_copy = mu_a;
        arma::vec g_sum_copy = g_sum;
        double eps = std::pow(std::numeric_limits<double>::epsilon(), 0.25);
        double spa_pval = SPA_pval(mu_a_copy, g_sum_copy, q_sum, qinv,
                                    pval_noadj, eps, true, traitType, isSPAConverge);
        if (isSPAConverge && std::isfinite(spa_pval)) {
            p_value_burden = spa_pval; // SPA_pval returns log(p) when logp=true
        } else {
            p_value_burden = pval_noadj;
        }
    }

    return SPA_ER_kernel_related_Phiadj_fast_new(p_new, Score, Phi, p_value_burden, regionTestType);
}


// ============================================================
// get_CCT_pvalue
// Port from SAIGE/R/SAIGE_SPATest_Region_Func.R: get_CCT_pvalue()
// Combines p-values using Cauchy combination test,
// handling NA values.
// ============================================================
static double get_CCT_pvalue(const std::vector<double>& pvals) {
    std::vector<double> validPvals;
    for (double p : pvals) {
        if (std::isfinite(p) && p >= 0.0 && p <= 1.0) {
            validPvals.push_back(p);
        }
    }
    if (validPvals.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    arma::vec pvalArma = arma::conv_to<arma::vec>::from(validPvals);
    return CCT_cpp(pvalArma);
}


// ============================================================
// mainRegionInCPP
// Unified standalone function that combines:
//   A) The C++ marker loop from SAIGE/src/Main.cpp:1032-2179
//   B) The R orchestration from SAIGE/R/SAIGE_SPATest_Region.R:376-1000
//
// This is the standalone version: marker loop, URV collapsing,
// VarMat construction, SPA Phi adjustment, SKAT/BURDEN/SKAT-O,
// CCT combination, and output writing all happen in one function.
// ============================================================
void mainRegionInCPP(
    std::string t_genoType,
    std::string t_traitType,
    RegionData& region,
    arma::vec& maxMAFVec,
    std::string t_outputFile,
    unsigned int t_n,
    arma::mat& P1Mat,
    arma::mat& P2Mat,
    std::string t_regionTestType,  // "SKATO", "SKAT", or "BURDEN"
    bool t_isImputation,
    arma::vec& t_weight,
    bool t_isSingleinGroupTest,
    bool t_isMoreOutput,
    arma::vec& r_corr,
    arma::vec& mu)
{
    std::string regionName = region.regionName;
    std::vector<std::string>& t_genoIndex = region.genoIndex;
    std::vector<std::string>& t_genoIndex_prev = region.genoIndex_prev;
    arma::imat& annoIndicatorMat_input = region.annoIndicatorMat;
    std::vector<std::string>& annoStringVec = region.annoVec;

    bool isWeightCustomized = false;
    bool isEqualWeights = false;
    unsigned int q0 = t_genoIndex.size();

    if (t_weight.n_elem == q0 && !t_weight.is_zero()) {
        isWeightCustomized = true;
        bool all_ones = arma::all(t_weight == 1.0);
        if (all_ones) {
            isEqualWeights = true;
        }
    }

    unsigned int q_anno = annoIndicatorMat_input.n_cols;
    unsigned int q_maf = maxMAFVec.n_elem;
    unsigned int q_anno_maf = q_anno * q_maf;
    arma::mat genoURMat(t_n, q_anno_maf, arma::fill::zeros);
    arma::mat genoURMat_noweights(t_n, q_anno_maf, arma::fill::zeros);
    unsigned int q = q0 + q_anno_maf;
    arma::imat annoMAFIndicatorMat(q, q_anno_maf, arma::fill::zeros);
    arma::ivec annoMAFIndicatorVec(q_anno_maf);
    annoMAFIndicatorVec.zeros();
    arma::vec maxMAFperAnno(q_anno, arma::fill::zeros);
    arma::vec MAFIndicatorVec(maxMAFVec.n_elem);
    MAFIndicatorVec.zeros();

    // Beta distribution for weights
    boost::math::beta_distribution<> beta_dist(g_weights_beta[0], g_weights_beta[1]);

    bool isCondition = ptr_gSAIGEobj->m_isCondition;

    arma::mat genoSumMat(t_n, q_anno_maf, arma::fill::zeros);
    arma::vec genoSumcount_noweight(q_anno_maf, arma::fill::zeros);

    // Group-level stats
    arma::vec MAC_GroupVec(q_anno_maf, arma::fill::zeros);
    arma::vec MACCase_GroupVec(q_anno_maf, arma::fill::zeros);
    arma::vec MACControl_GroupVec(q_anno_maf, arma::fill::zeros);
    arma::vec NumRare_GroupVec(q_anno_maf, arma::fill::zeros);
    arma::vec NumUltraRare_GroupVec(q_anno_maf, arma::fill::zeros);

    // Single-variant assoc output
    arma::uvec indicatorVec(q, arma::fill::zeros);
    std::vector<std::string> markerVec(q);
    std::vector<std::string> chrVec(q);
    std::vector<std::string> posVec(q);
    std::vector<std::string> refVec(q);
    std::vector<std::string> altVec(q);
    std::vector<std::string> infoVec(q);
    std::vector<double> altFreqVec(q, arma::datum::nan);
    std::vector<double> MACVec(q, arma::datum::nan);
    std::vector<double> MAFVec(q, arma::datum::nan);
    std::vector<double> altCountsVec(q, arma::datum::nan);
    std::vector<double> imputationInfoVec(q, arma::datum::nan);
    std::vector<double> missingRateVec(q, arma::datum::nan);
    std::vector<double> BetaVec(q, arma::datum::nan);
    std::vector<double> seBetaVec(q, arma::datum::nan);
    std::vector<std::string> pvalVec(q, "NA");
    std::vector<double> TstatVec(q, arma::datum::nan);
    std::vector<double> TstatVec_flip(q, arma::datum::nan);
    std::vector<double> gyVec(q, arma::datum::nan);
    std::vector<double> varTVec(q, arma::datum::nan);
    std::vector<std::string> pvalNAVec(q, "NA");
    std::vector<bool> isSPAConvergeVec(q);

    std::vector<double> AF_caseVec(q, arma::datum::nan);
    std::vector<double> AF_ctrlVec(q, arma::datum::nan);
    std::vector<uint32_t> N_caseVec(q);
    std::vector<uint32_t> N_ctrlVec(q);
    std::vector<double> N_case_homVec(q, arma::datum::nan);
    std::vector<double> N_ctrl_hetVec(q, arma::datum::nan);
    std::vector<double> N_case_hetVec(q, arma::datum::nan);
    std::vector<double> N_ctrl_homVec(q, arma::datum::nan);
    std::vector<uint32_t> N_Vec(q);

    // Conditional analysis vectors (structural support)
    std::vector<double> Beta_cVec(q, arma::datum::nan);
    std::vector<double> seBeta_cVec(q, arma::datum::nan);
    std::vector<std::string> pval_cVec(q, "NA");
    std::vector<double> Tstat_cVec(q, arma::datum::nan);
    std::vector<double> varT_cVec(q, arma::datum::nan);
    std::vector<std::string> pvalNA_cVec(q, "NA");
    unsigned int q_cond = (ptr_gSAIGEobj->m_VarInvMat_cond).n_rows;
    arma::rowvec G1tilde_P_G2tilde_Vec(q_cond);

    // Marker testing variables
    unsigned int m1 = g_region_maxMarkers_cutoff;
    std::vector<unsigned int> mPassCVVec;
    std::string pval, pval_noSPA, pval_c, pval_noSPA_c;
    double Beta, seBeta, Tstat, varT, gy;
    double Beta_c, seBeta_c, Tstat_c, varT_c;
    bool isSPAConverge, is_gtilde, is_Firth, is_FirthConverge;
    arma::vec P1Vec(t_n), P2Vec(t_n);
    arma::vec GVec(t_n);
    arma::vec gtildeVec;
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    double MACgroup, MACcasegroup, MACcontrolgroup, AF_case, AF_ctrl;
    uint32_t N_case, N_ctrl;

    bool hasVarRatio = true;
    bool isSingleVarianceRatio = true;
    if ((ptr_gSAIGEobj->m_varRatio_null).n_elem > 1) {
        isSingleVarianceRatio = false;
    } else {
        ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur);
    }

    unsigned int nchunks = 0;
    unsigned int ichunk = 0;
    unsigned int i1InChunk = 0;
    unsigned int i1 = 0;    // non-URV markers
    unsigned int i2 = 0;    // URV markers
    unsigned int jm;

    // ===== Marker loop =====
    for (unsigned int i = 0; i < q0; i++) {
        double altFreq, altCounts, missingRate, imputeInfo;
        std::vector<uint32_t> indexForMissing;
        std::string chr, ref, alt, marker;
        uint32_t pd;
        bool flip = false;
        bool isOutputIndexForMissing = true;
        bool isOnlyOutputNonZero = false;

        GVec.resize(t_n);
        GVec.zeros();

        std::string t_genoIndex_str = t_genoIndex.at(i);
        char* end;
        uint64_t gIndex = std::strtoull(t_genoIndex_str.c_str(), &end, 10);

        uint64_t gIndex_prev = 0;
        if (i == 0) {
            gIndex_prev = 0;
        } else {
            char* end_prev;
            std::string t_genoIndex_prev_str;
            if (t_genoType == "bgen") {
                t_genoIndex_prev_str = t_genoIndex_prev.at(i - 1);
            } else if (t_genoType == "plink" || t_genoType == "pgen" || t_genoType == "vcf") {
                t_genoIndex_prev_str = t_genoIndex.at(i - 1);
            }
            gIndex_prev = std::strtoull(t_genoIndex_prev_str.c_str(), &end_prev, 10);
        }

        bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex,
            ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
            isOutputIndexForMissing, indexForMissing,
            isOnlyOutputNonZero, indexNonZeroVec, GVec, t_isImputation);

        if (!isReadMarker) {
            std::cout << "ERROR: Reading " << i << "th marker failed." << std::endl;
            break;
        }

        std::string pds = std::to_string(pd);
        std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;

        double MAF = std::min(altFreq, 1.0 - altFreq);
        double w0;
        double MAC = MAF * 2 * t_n * (1 - missingRate);
        flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing,
            g_impute_method, g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff,
            MAC, indexZeroVec, indexNonZeroVec);

        arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
        MAF = std::min(altFreq, 1.0 - altFreq);
        MAC = std::min(altCounts, (double)t_n * 2 - altCounts);

        chrVec.at(i) = chr;
        posVec.at(i) = pds;
        refVec.at(i) = ref;
        altVec.at(i) = alt;
        markerVec.at(i) = marker;
        infoVec.at(i) = info;
        altFreqVec.at(i) = altFreq;
        missingRateVec.at(i) = missingRate;
        altCountsVec.at(i) = altCounts;
        MACVec.at(i) = MAC;
        MAFVec.at(i) = MAF;
        imputationInfoVec.at(i) = imputeInfo;

        if ((missingRate > g_missingRate_cutoff) || (MAF > g_maxMAFLimit) ||
            (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff) ||
            (imputeInfo < g_marker_minINFO_cutoff)) {
            continue;
        }

        if (isWeightCustomized) {
            w0 = t_weight(i);
        } else {
            w0 = boost::math::pdf(beta_dist, MAF);
        }

        indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
        uint nNonZero = indexNonZeroVec_arma.n_elem;

        if (MAC > g_region_minMAC_cutoff) {
            // ===== Non-URV marker =====
            indicatorVec.at(i) = 1;
            if (i1InChunk == 0) {
                std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
            }

            if (MAC > ptr_gSAIGEobj->m_cateVarRatioMinMACVecExclude.back()) {
                ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
            } else {
                ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
            }

            if (!isSingleVarianceRatio) {
                hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC, ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
            } else {
                ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
            }

            if (t_regionTestType != "BURDEN" || t_isSingleinGroupTest) {
                indexZeroVec_arma = arma::conv_to<arma::uvec>::from(indexZeroVec);

                // Phase A: build ctx mirroring the set_flagSparseGRM_cur and
                // assignVarianceRatio* mutations just above. The region path
                // hard-codes isnoXadj=false at the assignVarianceRatio call,
                // and t_isnoadjCov is passed as false at both call sites.
                SAIGE::PerMarkerCtx ctx_region;
                ctx_region.flagSparseGRM_cur =
                    (MAC > ptr_gSAIGEobj->m_cateVarRatioMinMACVecExclude.back())
                        ? false : ptr_gSAIGEobj->m_flagSparseGRM;
                ctx_region.isnoadjCov_cur = false;
                {
                    bool dummyHas;
                    if (!isSingleVarianceRatio) {
                        ctx_region.varRatioVal = ptr_gSAIGEobj->computeVarianceRatio(
                            MAC, ctx_region.flagSparseGRM_cur, false, dummyHas);
                    } else {
                        ctx_region.varRatioVal = ptr_gSAIGEobj->computeSingleVarianceRatio(
                            ctx_region.flagSparseGRM_cur, false);
                    }
                }

                if (MAC <= g_MACCutoffforER && t_traitType == "binary") {
                    Unified_getMarkerPval(GVec, false, indexNonZeroVec_arma, indexZeroVec_arma,
                        Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq,
                        isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition,
                        Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c,
                        G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge,
                        true, false, ptr_gSAIGEobj->m_flagSparseGRM_cur, ctx_region);
                } else {
                    Unified_getMarkerPval(GVec, false, indexNonZeroVec_arma, indexZeroVec_arma,
                        Beta, seBeta, pval, pval_noSPA, Tstat, gy, varT, altFreq,
                        isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition,
                        Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c,
                        G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge,
                        false, false, ptr_gSAIGEobj->m_flagSparseGRM_cur, ctx_region);
                }

                BetaVec.at(i) = Beta * (1 - 2 * flip);
                seBetaVec.at(i) = seBeta;
                pvalVec.at(i) = pval;
                pvalNAVec.at(i) = pval_noSPA;
                TstatVec.at(i) = Tstat * (1 - 2 * flip);
                TstatVec_flip.at(i) = Tstat;
                gyVec.at(i) = gy;
                varTVec.at(i) = varT;
                isSPAConvergeVec.at(i) = isSPAConverge;

                if (t_regionTestType != "BURDEN") {
                    P1Mat.row(i1InChunk) = std::sqrt(ptr_gSAIGEobj->m_varRatioVal) * gtildeVec.t();
                    P2Mat.col(i1InChunk) = std::sqrt(ptr_gSAIGEobj->m_varRatioVal) * P2Vec;
                }
            }

            i1 += 1;
            i1InChunk += 1;

            // Accumulate group statistics
            arma::vec dosage_case, dosage_ctrl;
            if (t_traitType == "binary" || t_traitType == "survival") {
                dosage_case = GVec.elem(ptr_gSAIGEobj->m_case_indices);
                dosage_ctrl = GVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                MACcasegroup = arma::accu(dosage_case);
                MACcontrolgroup = arma::accu(dosage_ctrl);
            }

            MAFIndicatorVec.zeros();
            MAFIndicatorVec.elem(arma::find(maxMAFVec >= MAF)).ones();
            annoMAFIndicatorVec.zeros();
            for (unsigned int j = 0; j < q_anno; j++) {
                if (annoIndicatorMat_input(i, j) == 1) {
                    maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
                    for (unsigned int m = 0; m < q_maf; m++) {
                        if (MAFIndicatorVec(m) == 1) {
                            jm = j * q_maf + m;
                            annoMAFIndicatorVec(jm) = 1;
                            MAC_GroupVec(jm) = MAC_GroupVec(jm) + MAC;
                            if (t_traitType == "binary" || t_traitType == "survival") {
                                MACCase_GroupVec(jm) = MACCase_GroupVec(jm) + MACcasegroup;
                                MACControl_GroupVec(jm) = MACControl_GroupVec(jm) + MACcontrolgroup;
                            }
                            for (unsigned int k = 0; k < nNonZero; k++) {
                                genoSumMat(indexNonZeroVec_arma(k), jm) += w0 * GVec(indexNonZeroVec_arma(k));
                                genoSumcount_noweight(jm) += GVec(indexNonZeroVec_arma(k));
                            }
                            NumRare_GroupVec(jm) += 1;
                        }
                    }
                }
            }
            annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();

            if (t_regionTestType != "BURDEN" || t_isSingleinGroupTest) {
                if (t_traitType == "binary" || t_traitType == "survival") {
                    AF_case = arma::mean(dosage_case) / 2;
                    AF_ctrl = arma::mean(dosage_ctrl) / 2;
                    if (flip) { AF_case = 1 - AF_case; AF_ctrl = 1 - AF_ctrl; }
                    AF_caseVec.at(i) = AF_case;
                    AF_ctrlVec.at(i) = AF_ctrl;
                    N_case = dosage_case.n_elem;
                    N_ctrl = dosage_ctrl.n_elem;
                    N_caseVec.at(i) = N_case;
                    N_ctrlVec.at(i) = N_ctrl;
                } else if (t_traitType == "quantitative") {
                    N_Vec.at(i) = t_n;
                }
            }

        } else {
            // ===== Ultra-Rare Variant (URV) =====
            indicatorVec.at(i) = 2;
            arma::vec MAFIndVec_local(maxMAFVec.n_elem, arma::fill::zeros);
            MAFIndVec_local.elem(arma::find(maxMAFVec >= MAF)).ones();
            annoMAFIndicatorVec.zeros();
            for (unsigned int j = 0; j < q_anno; j++) {
                if (annoIndicatorMat_input(i, j) == 1) {
                    maxMAFperAnno(j) = std::max(maxMAFperAnno(j), MAF);
                    for (unsigned int m = 0; m < q_maf; m++) {
                        if (MAFIndVec_local(m) == 1) {
                            jm = j * q_maf + m;
                            annoMAFIndicatorVec(jm) = 2;
                            if (!isWeightCustomized) {
                                for (unsigned int k = 0; k < nNonZero; k++) {
                                    genoURMat(indexNonZeroVec_arma(k), jm) = std::max(
                                        genoURMat(indexNonZeroVec_arma(k), jm),
                                        GVec(indexNonZeroVec_arma(k)));
                                }
                            } else {
                                for (unsigned int k = 0; k < nNonZero; k++) {
                                    genoURMat(indexNonZeroVec_arma(k), jm) = std::max(
                                        genoURMat(indexNonZeroVec_arma(k), jm),
                                        t_weight(i) * GVec(indexNonZeroVec_arma(k)));
                                    genoURMat_noweights(indexNonZeroVec_arma(k), jm) = std::max(
                                        genoURMat_noweights(indexNonZeroVec_arma(k), jm),
                                        GVec(indexNonZeroVec_arma(k)));
                                }
                            }
                            NumUltraRare_GroupVec(jm) += 1;
                        }
                    }
                }
            }
            annoMAFIndicatorMat.row(i) = annoMAFIndicatorVec.t();
            i2 += 1;
        }

        // Chunk management
        if (i1InChunk == m1) {
            std::cout << "In chunks 0-" << ichunk << ", " << i2
                      << " markers are ultra-rare and " << i1
                      << " markers are not ultra-rare." << std::endl;
            if (t_regionTestType != "BURDEN") {
                P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
                P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
            }
            mPassCVVec.push_back(m1);
            ichunk += 1;
            i1InChunk = 0;
            nchunks += 1;
        }
    } // end marker loop

    // Save last chunk of non-URV markers
    if (i1InChunk != 0) {
        std::cout << "In chunks 0-" << ichunk << ", " << i2
                  << " markers are ultra-rare and " << i1
                  << " markers are not ultra-rare." << std::endl;
        if (t_regionTestType != "BURDEN") {
            P1Mat = P1Mat.rows(0, i1InChunk - 1);
            P2Mat = P2Mat.cols(0, i1InChunk - 1);
            P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
            P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        }
        ichunk += 1;
        mPassCVVec.push_back(i1InChunk);
        nchunks += 1;
        i1InChunk = 0;
    }

    // ===== Process URV pseudo-markers =====
    if (i2 > 0) {
        int m1new = std::max(m1, q_anno_maf);
        if (t_regionTestType != "BURDEN") {
            P1Mat.resize(m1new, P1Mat.n_cols);
            P2Mat.resize(P2Mat.n_rows, m1new);
        }

        arma::mat XV, XXVX_inv;
        ptr_gSAIGEobj->extract_XV_XXVX_inv(XV, XXVX_inv);

        unsigned int i_ur;
        for (unsigned int j = 0; j < q_anno; j++) {
            for (unsigned int m = 0; m < q_maf; m++) {
                jm = j * q_maf + m;
                arma::vec genoURVec = genoURMat.col(jm);
                arma::vec genoURVec_noweights = genoURMat_noweights.col(jm);
                arma::uvec indexForNonZero = arma::find(genoURVec != 0);
                i_ur = q0 + jm;
                markerVec.at(i_ur) = "UR";

                if (indexForNonZero.n_elem > 0) {
                    double altFreq_ur = arma::mean(genoURVec) / 2;
                    double altCounts_ur = arma::accu(genoURVec);
                    double missingRate_ur = 0;
                    bool flip_ur = false;
                    double MAF_ur = std::min(altFreq_ur, 1.0 - altFreq_ur);
                    double w0_ur;
                    double MAC_ur = MAF_ur * 2 * t_n;
                    std::vector<uint32_t> indexForMissing_ur;
                    flip_ur = imputeGenoAndFlip(genoURVec, altFreq_ur, altCounts_ur,
                        indexForMissing_ur, g_impute_method,
                        g_dosage_zerod_cutoff, g_dosage_zerod_MAC_cutoff,
                        MAC_ur, indexZeroVec, indexNonZeroVec);

                    if (isWeightCustomized) {
                        for (unsigned int k = 0; k < indexForNonZero.n_elem; k++) {
                            genoSumMat(indexForNonZero(k), jm) += genoURVec(indexForNonZero(k));
                            genoSumcount_noweight(jm) += genoURVec_noweights(indexForNonZero(k));
                        }
                    } else {
                        w0_ur = boost::math::pdf(beta_dist, MAF_ur);
                        for (unsigned int k = 0; k < indexForNonZero.n_elem; k++) {
                            genoSumMat(indexForNonZero(k), jm) += genoURVec(indexForNonZero(k)) * w0_ur;
                            genoSumcount_noweight(jm) += genoURVec(indexForNonZero(k));
                        }
                    }

                    if (t_regionTestType != "BURDEN") {
                        arma::vec genoSumMatvec1 = genoSumMat.col(jm);
                        arma::vec genoSumMatvec2 = XV * genoSumMatvec1;
                        arma::vec genoSumMatvec3 = genoSumMatvec1 - XXVX_inv * genoSumMatvec2;
                        genoSumMat.col(jm) = genoSumMatvec3;
                    }

                    MAC_ur = MAF_ur * 2 * t_n;

                    if (t_regionTestType != "BURDEN" || t_isSingleinGroupTest) {
                        if (MAC_ur > ptr_gSAIGEobj->m_cateVarRatioMinMACVecExclude.back()) {
                            ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
                        } else {
                            ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
                        }
                        if (!isSingleVarianceRatio) {
                            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC_ur, ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
                        } else {
                            ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
                        }

                        annoMAFIndicatorVec.zeros();
                        annoMAFIndicatorVec(jm) = 1;
                        annoMAFIndicatorMat.row(i_ur) = annoMAFIndicatorVec.t();
                        altFreqVec.at(i_ur) = altFreq_ur;
                        altCountsVec.at(i_ur) = altCounts_ur;
                        missingRateVec.at(i_ur) = missingRate_ur;
                        MACVec.at(i_ur) = MAC_ur;
                        MAFVec.at(i_ur) = MAF_ur;

                        arma::uvec indexZeroVec_arma_ur = arma::conv_to<arma::uvec>::from(indexZeroVec);
                        arma::uvec indexNonZeroVec_arma_ur = arma::conv_to<arma::uvec>::from(indexNonZeroVec);

                        if (MAC_ur <= g_MACCutoffforER && t_traitType == "binary"
                            && (!isWeightCustomized || isEqualWeights)) {
                            ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma_ur, indexZeroVec_arma_ur,
                                Beta, seBeta, pval, pval_noSPA, altFreq_ur, Tstat, gy, varT,
                                isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition,
                                Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c,
                                G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge,
                                true, false, ptr_gSAIGEobj->m_flagSparseGRM_cur);
                        } else {
                            ptr_gSAIGEobj->getMarkerPval(genoURVec, indexNonZeroVec_arma_ur, indexZeroVec_arma_ur,
                                Beta, seBeta, pval, pval_noSPA, altFreq_ur, Tstat, gy, varT,
                                isSPAConverge, gtildeVec, is_gtilde, true, P2Vec, isCondition,
                                Beta_c, seBeta_c, pval_c, pval_noSPA_c, Tstat_c, varT_c,
                                G1tilde_P_G2tilde_Vec, is_Firth, is_FirthConverge,
                                false, false, ptr_gSAIGEobj->m_flagSparseGRM_cur);
                        }

                        BetaVec.at(i_ur) = Beta * (1 - 2 * flip_ur);
                        seBetaVec.at(i_ur) = seBeta;
                        pvalVec.at(i_ur) = pval;
                        pvalNAVec.at(i_ur) = pval_noSPA;
                        TstatVec.at(i_ur) = Tstat * (1 - 2 * flip_ur);
                        TstatVec_flip.at(i_ur) = Tstat;
                        gyVec.at(i_ur) = gy;
                        varTVec.at(i_ur) = varT;
                        isSPAConvergeVec.at(i_ur) = isSPAConverge;

                        chrVec.at(i_ur) = "UR";
                        posVec.at(i_ur) = "UR";
                        refVec.at(i_ur) = "UR";
                        altVec.at(i_ur) = "UR";

                        std::string maxMAFStr = std::to_string(maxMAFVec.at(m));
                        maxMAFStr.erase(maxMAFStr.find_last_not_of('0') + 1, std::string::npos);
                        markerVec.at(i_ur) = regionName + ":" + annoStringVec.at(j) + ":" + maxMAFStr;

                        MAC_GroupVec(jm) += MAC_ur;
                        if (t_traitType == "binary" || t_traitType == "survival") {
                            arma::vec dosage_case_ur = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                            arma::vec dosage_ctrl_ur = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                            MACCase_GroupVec(jm) += arma::accu(dosage_case_ur);
                            MACControl_GroupVec(jm) += arma::accu(dosage_ctrl_ur);
                            AF_caseVec.at(i_ur) = arma::mean(dosage_case_ur) / 2;
                            AF_ctrlVec.at(i_ur) = arma::mean(dosage_ctrl_ur) / 2;
                            N_caseVec.at(i_ur) = dosage_case_ur.n_elem;
                            N_ctrlVec.at(i_ur) = dosage_ctrl_ur.n_elem;
                        } else if (t_traitType == "quantitative") {
                            N_Vec.at(i_ur) = t_n;
                        }

                        if (t_regionTestType != "BURDEN") {
                            P1Mat.row(i1InChunk) = std::sqrt(ptr_gSAIGEobj->m_varRatioVal) * gtildeVec.t();
                            P2Mat.col(i1InChunk) = std::sqrt(ptr_gSAIGEobj->m_varRatioVal) * P2Vec;
                        }
                    } else {
                        // BURDEN-only path for URV
                        MAC_GroupVec(jm) += MAC_ur;
                        if (t_traitType == "binary" || t_traitType == "survival") {
                            arma::vec dosage_case_ur = genoURVec.elem(ptr_gSAIGEobj->m_case_indices);
                            arma::vec dosage_ctrl_ur = genoURVec.elem(ptr_gSAIGEobj->m_ctrl_indices);
                            MACCase_GroupVec(jm) += arma::accu(dosage_case_ur);
                            MACControl_GroupVec(jm) += arma::accu(dosage_ctrl_ur);
                        }
                    }

                    i1InChunk += 1;
                    i1 += 1;
                }
            }
        }

        // Save last URV chunk
        if (i1InChunk != 0) {
            nchunks += 1;
            if (t_regionTestType != "BURDEN") {
                P1Mat = P1Mat.rows(0, i1InChunk - 1);
                P2Mat = P2Mat.cols(0, i1InChunk - 1);
                P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
                P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
            }
            ichunk += 1;
            mPassCVVec.push_back(i1InChunk);
        }
    } // end if(i2 > 0)

    int mPassCVVecsize = mPassCVVec.size();
    nchunks = mPassCVVecsize;

    // ===== Build VarMat =====
    arma::mat VarMat;
    if (t_regionTestType != "BURDEN") {
        VarMat.resize(i1, i1);
        if (nchunks == 1) {
            VarMat = P1Mat * P2Mat;
        }
        if (nchunks > 1) {
            int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
            for (unsigned int index1 = 0; index1 < nchunks; index1++) {
                last_row = first_row + mPassCVVec.at(index1) - 1;
                std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
                P1Mat.load(P1MatFile);
                if (P1Mat.n_cols == 0) continue;

                for (unsigned int index2 = 0; index2 < index1; index2++) {
                    P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
                    if (P2Mat.n_cols == 0) continue;
                    arma::mat offVarMat = P1Mat * P2Mat;
                    last_col = first_col + mPassCVVec.at(index2) - 1;
                    VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
                    VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
                    first_col = last_col + 1;
                }

                last_col = first_col + mPassCVVec.at(index1) - 1;
                P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
                arma::mat diagVarMat = P1Mat * P2Mat;
                VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;
                first_row = last_row + 1;
                first_col = 0;
            }

            // Clean up chunk files
            for (unsigned int index1 = 0; index1 < nchunks; index1++) {
                std::string P1f = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
                std::string P2f = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
                std::remove(P1f.c_str());
                std::remove(P2f.c_str());
            }
        }
    }

    // ===== Check max MAF per annotation =====
    arma::uvec q_maf_for_anno(q_anno);
    for (unsigned int j = 0; j < q_anno; j++) {
        arma::uvec jtemp = arma::find(maxMAFVec >= maxMAFperAnno(j));
        q_maf_for_anno(j) = jtemp.min();
    }

    // ===== BURDEN-only path =====
    if (t_regionTestType == "BURDEN") {
        std::vector<std::string> BURDEN_pval_Vec(q_anno_maf, "NA");
        std::vector<std::string> BURDEN_pval_cVec(q_anno_maf, "NA");
        std::vector<std::string> BURDEN_AnnoName_Vec(q_anno_maf);
        std::vector<std::string> BURDEN_maxMAFName_Vec(q_anno_maf);
        std::vector<double> BURDEN_Beta_Vec(q_anno_maf);
        std::vector<double> BURDEN_seBeta_Vec(q_anno_maf);
        std::vector<double> BURDEN_Beta_cVec(q_anno_maf);
        std::vector<double> BURDEN_seBeta_cVec(q_anno_maf);

        bool isregion = ptr_gSAIGEobj->m_flagSparseGRM;
        unsigned int q_maf_m;
        bool isPolyMarker;
        std::string AnnoName;
        double maxMAFName;

        for (unsigned int j = 0; j < q_anno; j++) {
            q_maf_m = q_maf_for_anno(j);
            AnnoName = annoStringVec[j];
            isPolyMarker = true;
            for (unsigned int m = 0; m < q_maf; m++) {
                maxMAFName = maxMAFVec(m);
                jm = j * q_maf + m;
                unsigned int i_b = jm;
                if (m <= q_maf_m) {
                    arma::vec genoSumVec = genoSumMat.col(jm);
                    arma::uvec indexNZ = arma::find(genoSumVec != 0);
                    arma::uvec indexZ = arma::find(genoSumVec == 0);
                    double altCounts_b = genoSumcount_noweight(i_b);
                    double altFreq_b = altCounts_b / (2.0 * t_n);
                    double MAC_b, MAF_b;
                    if (altFreq_b > 1) { MAF_b = 1; MAC_b = t_n; }
                    else { MAF_b = std::min(altFreq_b, 1.0 - altFreq_b); MAC_b = MAF_b * 2 * t_n; }

                    if (indexNZ.n_elem > 0 && MAC_b >= g_min_gourpmac_for_burdenonly) {
                        isPolyMarker = true;
                        if (MAC_b > ptr_gSAIGEobj->m_cateVarRatioMinMACVecExclude.back()) {
                            ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
                        } else {
                            ptr_gSAIGEobj->set_flagSparseGRM_cur(ptr_gSAIGEobj->m_flagSparseGRM);
                        }
                        if (!isSingleVarianceRatio) {
                            hasVarRatio = ptr_gSAIGEobj->assignVarianceRatio(MAC_b, ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
                        } else {
                            ptr_gSAIGEobj->assignSingleVarianceRatio(ptr_gSAIGEobj->m_flagSparseGRM_cur, false);
                        }

                        if (MAC_b <= g_MACCutoffforER && t_traitType == "binary"
                            && (!isWeightCustomized || isEqualWeights)) {
                            ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNZ, indexZ,
                                Beta, seBeta, pval, pval_noSPA, altFreq_b, Tstat, gy, varT,
                                isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec,
                                isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                                Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                                is_Firth, is_FirthConverge, true, false,
                                ptr_gSAIGEobj->m_flagSparseGRM_cur);
                        } else {
                            ptr_gSAIGEobj->getMarkerPval(genoSumVec, indexNZ, indexZ,
                                Beta, seBeta, pval, pval_noSPA, altFreq_b, Tstat, gy, varT,
                                isSPAConverge, gtildeVec, is_gtilde, isregion, P2Vec,
                                isCondition, Beta_c, seBeta_c, pval_c, pval_noSPA_c,
                                Tstat_c, varT_c, G1tilde_P_G2tilde_Vec,
                                is_Firth, is_FirthConverge, false, false,
                                ptr_gSAIGEobj->m_flagSparseGRM_cur);
                        }

                        BURDEN_AnnoName_Vec.at(i_b) = AnnoName;
                        std::string str = std::to_string(maxMAFName);
                        str.erase(str.find_last_not_of('0') + 1, std::string::npos);
                        BURDEN_maxMAFName_Vec.at(i_b) = str;
                        BURDEN_pval_Vec.at(i_b) = pval;
                        BURDEN_Beta_Vec.at(i_b) = Beta;
                        BURDEN_seBeta_Vec.at(i_b) = seBeta;
                        if (isCondition) {
                            BURDEN_pval_cVec.at(i_b) = pval_c;
                            BURDEN_Beta_cVec.at(i_b) = Beta_c;
                            BURDEN_seBeta_cVec.at(i_b) = seBeta_c;
                        }
                    } else {
                        isPolyMarker = false;
                    }
                } else {
                    if (isPolyMarker) {
                        BURDEN_AnnoName_Vec.at(i_b) = AnnoName;
                        BURDEN_maxMAFName_Vec.at(i_b) = std::to_string(maxMAFName);
                        BURDEN_pval_Vec.at(i_b) = pval;
                        BURDEN_Beta_Vec.at(i_b) = Beta;
                        BURDEN_seBeta_Vec.at(i_b) = seBeta;
                        if (isCondition) {
                            BURDEN_pval_cVec.at(i_b) = pval_c;
                            BURDEN_Beta_cVec.at(i_b) = Beta_c;
                            BURDEN_seBeta_cVec.at(i_b) = seBeta_c;
                        }
                    }
                }
            }
        }

        // CCT across BURDEN results
        std::vector<double> nonMissingPvalVec_std;
        double cctpval, cctpval_cond = std::numeric_limits<double>::quiet_NaN();
        for (unsigned int ix = 0; ix < BURDEN_pval_Vec.size(); ix++) {
            if (BURDEN_pval_Vec.at(ix) != "NA") {
                try {
                    nonMissingPvalVec_std.push_back(std::stod(BURDEN_pval_Vec.at(ix)));
                } catch (...) {}
            }
        }
        arma::vec nonMissingPvalVec = arma::conv_to<arma::vec>::from(nonMissingPvalVec_std);
        cctpval = CCT_cpp(nonMissingPvalVec);

        writeOutfile_BURDEN(regionName, BURDEN_AnnoName_Vec, BURDEN_maxMAFName_Vec,
            BURDEN_pval_Vec, BURDEN_Beta_Vec, BURDEN_seBeta_Vec,
            BURDEN_pval_cVec, BURDEN_Beta_cVec, BURDEN_seBeta_cVec,
            MAC_GroupVec, MACCase_GroupVec, MACControl_GroupVec,
            NumRare_GroupVec, NumUltraRare_GroupVec,
            cctpval, cctpval_cond, q_anno, q_maf, isCondition, t_traitType);

    } else {
        // ===== SKAT-O / SKAT path =====
        // This merges the R orchestration from SAIGE_SPATest_Region.R

        q_maf_for_anno = q_maf_for_anno + 1; // R convention: 1-indexed maxMAF0

        // Build noNAIndices (indices of non-NA pvalVec entries)
        std::vector<unsigned int> noNAIndices;
        // First convert pvalVec to log(p)
        std::vector<double> pvalVec_log(q);
        for (unsigned int ix = 0; ix < q; ix++) {
            pvalVec_log[ix] = convert_str_to_log(pvalVec[ix]);
        }
        for (unsigned int ix = 0; ix < q; ix++) {
            if (std::isfinite(pvalVec_log[ix])) {
                noNAIndices.push_back(ix);
            }
        }

        if (noNAIndices.empty()) {
            std::cout << "No valid markers for region " << regionName << ", skipping." << std::endl;
            return;
        }

        // Subset to non-NA
        unsigned int nValid = noNAIndices.size();
        arma::vec StatVec(nValid);
        arma::vec adjPVec(nValid);  // log(p) values
        arma::vec MAFVecSub(nValid);

        for (unsigned int ix = 0; ix < nValid; ix++) {
            StatVec(ix) = TstatVec_flip[noNAIndices[ix]];
            adjPVec(ix) = pvalVec_log[noNAIndices[ix]];
            MAFVecSub(ix) = MAFVec[noNAIndices[ix]];
        }

        // Subset VarMat to non-NA
        arma::mat VarMatSub = VarMat;  // already of size (i1 x i1)

        // Subset annoMAFIndicatorMat
        arma::imat annoMAFIndicatorMatSub(nValid, q_anno_maf, arma::fill::zeros);
        for (unsigned int ix = 0; ix < nValid; ix++) {
            for (unsigned int jj = 0; jj < q_anno_maf; jj++) {
                annoMAFIndicatorMatSub(ix, jj) = annoMAFIndicatorMat(noNAIndices[ix], jj);
            }
        }

        // Build AnnoWeights
        arma::vec AnnoWeights(nValid);
        if (t_weight.n_elem > 0 && arma::sum(t_weight) > 0) {
            // Use custom weights, with weight=1 for URV pseudo-markers
            arma::vec fullWeights(q, arma::fill::ones);
            for (unsigned int ix = 0; ix < q0 && ix < t_weight.n_elem; ix++) {
                fullWeights(ix) = t_weight(ix);
            }
            for (unsigned int ix = 0; ix < nValid; ix++) {
                AnnoWeights(ix) = fullWeights(noNAIndices[ix]);
            }
        } else {
            for (unsigned int ix = 0; ix < nValid; ix++) {
                AnnoWeights(ix) = boost::math::pdf(beta_dist, MAFVecSub(ix));
            }
        }

        // Weighted stats
        arma::mat weightMat = AnnoWeights * AnnoWeights.t();
        arma::vec wStatVec = StatVec % AnnoWeights;
        arma::mat wadjVarSMat = VarMatSub % weightMat;

        // Subset gyVec for binary traits
        arma::vec gyVecSub;
        if (t_traitType == "binary" || t_traitType == "survival") {
            gyVecSub.resize(nValid);
            for (unsigned int ix = 0; ix < nValid; ix++) {
                gyVecSub(ix) = gyVec[noNAIndices[ix]];
            }
        }

        // ===== Loop over annotation x MAF groups =====
        std::vector<double> pval_SKATO_vec;
        std::vector<double> pval_Burden_vec;
        std::vector<double> pval_SKAT_vec;
        std::vector<double> beta_Burden_vec;
        std::vector<double> se_Burden_vec;
        std::vector<std::string> annoName_vec;
        std::vector<double> maxMAFName_vec;
        std::vector<unsigned int> annoMAFIndVec;
        std::vector<double> mac_vec;
        std::vector<double> mac_case_vec;
        std::vector<double> mac_ctrl_vec;
        std::vector<double> nrare_vec;
        std::vector<double> nultra_vec;

        for (unsigned int j = 0; j < q_anno; j++) {
            std::string AnnoName = annoStringVec[j];
            unsigned int maxMAF0 = q_maf_for_anno(j); // 1-indexed
            bool isPolyRegion = true;

            for (unsigned int m = 0; m < q_maf; m++) {
                jm = j * q_maf + m;
                double maxMAFName = maxMAFVec(m);

                if (m + 1 <= maxMAF0) {
                    // Find indices where annoMAFIndicatorMatSub[,jm] == 1
                    std::vector<unsigned int> tempPos;
                    for (unsigned int ix = 0; ix < nValid; ix++) {
                        if (annoMAFIndicatorMatSub(ix, jm) == 1) {
                            tempPos.push_back(ix);
                        }
                    }

                    if (!tempPos.empty()) {
                        isPolyRegion = true;
                        annoMAFIndVec.push_back(jm);

                        unsigned int nSub = tempPos.size();
                        arma::mat Phi(nSub, nSub);
                        arma::vec Score(nSub);
                        for (unsigned int a = 0; a < nSub; a++) {
                            Score(a) = wStatVec(tempPos[a]);
                            for (unsigned int b = 0; b < nSub; b++) {
                                Phi(a, b) = wadjVarSMat(tempPos[a], tempPos[b]);
                            }
                        }

                        // SPA Phi adjustment for binary/survival traits
                        if (t_traitType == "binary" || t_traitType == "survival") {
                            arma::vec p_new(nSub);
                            for (unsigned int a = 0; a < nSub; a++) {
                                p_new(a) = adjPVec(tempPos[a]);
                            }
                            arma::vec g_sum = genoSumMat.col(jm);
                            double q_sum = 0.0;
                            for (unsigned int a = 0; a < nSub; a++) {
                                q_sum += gyVecSub(tempPos[a]) * AnnoWeights(tempPos[a]);
                            }

                            PhiAdjResult re_phi = get_newPhi_scaleFactor_traitType(
                                q_sum, mu, g_sum, p_new, Score, Phi,
                                t_regionTestType, t_traitType);
                            Phi = re_phi.Phi_adj;
                        }

                        // Call SKAT/BURDEN/SKAT-O
                        SKATResult groupResult = get_SKAT_pvalue(Score, Phi, r_corr, t_regionTestType);

                        pval_SKATO_vec.push_back(groupResult.pvalue_SKATO);
                        pval_Burden_vec.push_back(groupResult.pvalue_Burden);
                        pval_SKAT_vec.push_back(groupResult.pvalue_SKAT);
                        beta_Burden_vec.push_back(groupResult.beta_Burden);
                        se_Burden_vec.push_back(groupResult.se_Burden);
                        annoName_vec.push_back(AnnoName);
                        maxMAFName_vec.push_back(maxMAFName);
                        mac_vec.push_back(MAC_GroupVec(jm));
                        if (t_traitType == "binary" || t_traitType == "survival") {
                            mac_case_vec.push_back(MACCase_GroupVec(jm));
                            mac_ctrl_vec.push_back(MACControl_GroupVec(jm));
                        }
                        nrare_vec.push_back(NumRare_GroupVec(jm));
                        nultra_vec.push_back(NumUltraRare_GroupVec(jm));

                    } else {
                        isPolyRegion = false;
                    }
                } else {
                    if (isPolyRegion && !pval_SKATO_vec.empty()) {
                        // Repeat last result for higher MAF threshold
                        annoMAFIndVec.push_back(jm);
                        pval_SKATO_vec.push_back(pval_SKATO_vec.back());
                        pval_Burden_vec.push_back(pval_Burden_vec.back());
                        pval_SKAT_vec.push_back(pval_SKAT_vec.back());
                        beta_Burden_vec.push_back(beta_Burden_vec.back());
                        se_Burden_vec.push_back(se_Burden_vec.back());
                        annoName_vec.push_back(AnnoName);
                        maxMAFName_vec.push_back(maxMAFName);
                        mac_vec.push_back(MAC_GroupVec(jm));
                        if (t_traitType == "binary" || t_traitType == "survival") {
                            mac_case_vec.push_back(MACCase_GroupVec(jm));
                            mac_ctrl_vec.push_back(MACControl_GroupVec(jm));
                        }
                        nrare_vec.push_back(NumRare_GroupVec(jm));
                        nultra_vec.push_back(NumUltraRare_GroupVec(jm));
                    }
                }
            }
        }

        // ===== Write region output =====
        if (!pval_SKATO_vec.empty()) {
            // Write each annotation x MAF result
            for (size_t ix = 0; ix < pval_SKATO_vec.size(); ix++) {
                OutFile << regionName << "\t";
                OutFile << annoName_vec[ix] << "\t";
                OutFile << maxMAFName_vec[ix] << "\t";
                OutFile << pval_SKATO_vec[ix] << "\t";
                OutFile << pval_Burden_vec[ix] << "\t";
                OutFile << pval_SKAT_vec[ix] << "\t";
                OutFile << beta_Burden_vec[ix] << "\t";
                OutFile << se_Burden_vec[ix] << "\t";
                OutFile << mac_vec[ix] << "\t";
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile << mac_case_vec[ix] << "\t";
                    OutFile << mac_ctrl_vec[ix] << "\t";
                }
                OutFile << nrare_vec[ix] << "\t";
                OutFile << nultra_vec[ix] << "\n";
            }

            // CCT combination across annotation x MAF groups
            if (annoStringVec.size() > 1 || q_maf > 1) {
                double cctpval_SKATO = get_CCT_pvalue(pval_SKATO_vec);
                double cctpval_Burden = get_CCT_pvalue(pval_Burden_vec);
                double cctpval_SKAT = get_CCT_pvalue(pval_SKAT_vec);

                OutFile << regionName << "\tCauchy\tNA\t";
                OutFile << cctpval_SKATO << "\t";
                OutFile << cctpval_Burden << "\t";
                OutFile << cctpval_SKAT << "\t";
                OutFile << "NA\tNA\t";  // BETA, SE
                OutFile << "NA\t";       // MAC
                if (t_traitType == "binary" || t_traitType == "survival") {
                    OutFile << "NA\tNA\t";
                }
                OutFile << "NA\tNA\n"; // Number_rare, Number_ultra_rare
            }
        }
    }

    // ===== Write singleInGroup output =====
    if (t_isSingleinGroupTest) {
        int numofUR0 = writeOutfile_singleInGroup(t_isMoreOutput, t_isImputation,
            isCondition, is_Firth, 0, 0, t_traitType,
            chrVec, posVec, markerVec, refVec, altVec,
            altCountsVec, altFreqVec, imputationInfoVec, missingRateVec,
            BetaVec, seBetaVec, TstatVec, varTVec, pvalVec, pvalNAVec,
            isSPAConvergeVec, Beta_cVec, seBeta_cVec, Tstat_cVec, varT_cVec,
            pval_cVec, pvalNA_cVec, AF_caseVec, AF_ctrlVec,
            N_caseVec, N_ctrlVec, N_case_homVec, N_ctrl_hetVec,
            N_case_hetVec, N_ctrl_homVec, N_Vec, OutFile_singleInGroup);
    }
}


// ============================================================
// main() -- CLI entry point
// Parses YAML config file, loads null model, sets up genotype
// reader, runs single-variant or region testing.
// ============================================================
int main(int argc, char* argv[])
{
    try {
        // ---- 1. Parse command-line ----
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <config.yaml>" << std::endl;
            std::cerr << std::endl;
            std::cerr << "SAIGE Step 2: Association testing (standalone C++ port)" << std::endl;
            std::cerr << "  Single-variant testing (default) or region/gene-based testing" << std::endl;
            std::cerr << std::endl;
            std::cerr << "Required YAML config keys:" << std::endl;
            std::cerr << "  modelFile:         Path to null model directory (from Step 1)" << std::endl;
            std::cerr << "  varianceRatioFile: Path to varianceRatio.txt" << std::endl;
            std::cerr << "  plinkFile:         Path to PLINK prefix (.bed/.bim/.fam) [for genoType=plink]" << std::endl;
            std::cerr << "  vcfFile:           Path to VCF/BCF/VCF.GZ file [for genoType=vcf]" << std::endl;
            std::cerr << "  bgenFile:          Path to BGEN file [for genoType=bgen]" << std::endl;
            std::cerr << "  bgenSampleFile:    Path to BGEN .sample file [for genoType=bgen]" << std::endl;
            std::cerr << "  outputFile:        Output file path" << std::endl;
            std::cerr << std::endl;
            std::cerr << "Optional YAML config keys (general):" << std::endl;
            std::cerr << "  minMAF:            Minimum MAF (default: 0)" << std::endl;
            std::cerr << "  minMAC:            Minimum MAC (default: 0.5)" << std::endl;
            std::cerr << "  maxMissRate:       Maximum missing rate (default: 0.15)" << std::endl;
            std::cerr << "  minINFO:           Minimum imputation info (default: 0)" << std::endl;
            std::cerr << "  AlleleOrder:       'alt-first' or 'ref-first' (default: 'alt-first')" << std::endl;
            std::cerr << "  dosage_zerod_cutoff:     (default: 0.2)" << std::endl;
            std::cerr << "  dosage_zerod_MAC_cutoff: (default: 10)" << std::endl;
            std::cerr << "  pgenFile:          Path to .pgen file [for genoType=pgen]" << std::endl;
            std::cerr << "  pvarFile:          Path to .pvar file [for genoType=pgen]" << std::endl;
            std::cerr << "  psamFile:          Path to .psam file [for genoType=pgen]" << std::endl;
            std::cerr << "  genoType:          'plink', 'vcf', 'bgen', or 'pgen' (default: 'plink')" << std::endl;
            std::cerr << "  vcfField:          'GT' or 'DS' (default: 'GT', for genoType=vcf)" << std::endl;
            std::cerr << "  isImputation:      true/false (default: false)" << std::endl;
            std::cerr << "  isMoreOutput:      true/false (default: false)" << std::endl;
            std::cerr << "  marker_chunksize:  markers per progress report (default: 10000)" << std::endl;
            std::cerr << "  nThreads:          OpenMP threads for marker loop (default: 1)" << std::endl;
            std::cerr << "  MACCutoffforER:    MAC cutoff for ER (default: 4)" << std::endl;
            std::cerr << "  weights_beta:      [a, b] for Beta weights (default: [1, 25])" << std::endl;
            std::cerr << "  checkpointDir:     Directory for checkpoint outputs (optional)" << std::endl;
            std::cerr << std::endl;
            std::cerr << "Region/gene-based testing (add groupFile to enable):" << std::endl;
            std::cerr << "  groupFile:         Path to group file (gene definitions)" << std::endl;
            std::cerr << "  annotationList:    [\"lof\", \"lof;missense\", ...] (required for region)" << std::endl;
            std::cerr << "  maxMAFList:        [0.0001, 0.001, 0.01] (required for region)" << std::endl;
            std::cerr << "  regionTestType:    'SKATO', 'SKAT', or 'BURDEN' (default: 'SKATO')" << std::endl;
            std::cerr << "  MACCutoff_to_CollapseUltraRare: (default: 10)" << std::endl;
            std::cerr << "  markers_per_chunk_in_groupTest:  (default: 500)" << std::endl;
            std::cerr << "  groups_per_chunk:  regions per I/O chunk (default: 100)" << std::endl;
            std::cerr << "  r_corr:            correlation parameter (0=SKAT-O, 1=BURDEN) (default: 0)" << std::endl;
            std::cerr << "  isSingleInGroupTest:   true/false (default: true)" << std::endl;
            std::cerr << "  isOutputMarkerList:    true/false (default: false)" << std::endl;
            std::cerr << "  max_markers_region:    max markers per region (default: 100000)" << std::endl;
            std::cerr << "  min_gourpmac_for_burdenonly: (default: 5)" << std::endl;
            std::cerr << std::endl;
            std::cerr << "LD matrix generation (requires groupFile):" << std::endl;
            std::cerr << "  isLDMatrix:        true/false (default: false)" << std::endl;
            std::cerr << "  ldmat_maxMAF:      max MAF for LD matrix (default: 0.5)" << std::endl;
            return 1;
        }

        std::string configFile = argv[1];
        std::cout << "==========================================" << std::endl;
        std::cout << "  SAIGE Step 2: Association Testing        " << std::endl;
        std::cout << "  Standalone C++ Port                     " << std::endl;
        std::cout << "==========================================" << std::endl;
        std::cout << std::endl;
        std::cout << "Loading config from: " << configFile << std::endl;

        // ---- 2. Read YAML config ----
        YAML::Node config = YAML::LoadFile(configFile);

        // Required keys
        if (!config["modelFile"]) {
            throw std::runtime_error("Config missing required key: modelFile");
        }
        if (!config["varianceRatioFile"]) {
            throw std::runtime_error("Config missing required key: varianceRatioFile");
        }
        if (!config["outputFile"]) {
            throw std::runtime_error("Config missing required key: outputFile");
        }

        // Determine genotype type early (needed for input file validation)
        std::string genoType_early = config["genoType"] ? config["genoType"].as<std::string>() : "plink";

        // Validate genotype input file
        if (genoType_early == "plink") {
            if (!config["plinkFile"]) {
                throw std::runtime_error("Config missing required key: plinkFile (needed for genoType=plink)");
            }
        } else if (genoType_early == "vcf") {
            if (!config["vcfFile"]) {
                throw std::runtime_error("Config missing required key: vcfFile (needed for genoType=vcf)");
            }
        } else if (genoType_early == "bgen") {
            if (!config["bgenFile"]) {
                throw std::runtime_error("Config missing required key: bgenFile (needed for genoType=bgen)");
            }
            if (!config["bgenSampleFile"]) {
                throw std::runtime_error("Config missing required key: bgenSampleFile (needed for genoType=bgen)");
            }
        } else if (genoType_early == "pgen") {
            if (!config["pgenFile"]) {
                throw std::runtime_error("Config missing required key: pgenFile (needed for genoType=pgen)");
            }
            if (!config["pvarFile"]) {
                throw std::runtime_error("Config missing required key: pvarFile (needed for genoType=pgen)");
            }
            if (!config["psamFile"]) {
                throw std::runtime_error("Config missing required key: psamFile (needed for genoType=pgen)");
            }
        } else {
            throw std::runtime_error("Unsupported genoType: " + genoType_early +
                                     ". Supported types: plink, vcf, bgen, pgen");
        }

        std::string modelFile = config["modelFile"].as<std::string>();
        std::string varianceRatioFile = config["varianceRatioFile"].as<std::string>();
        std::string plinkPrefix = config["plinkFile"] ? config["plinkFile"].as<std::string>() : "";
        std::string vcfFile = config["vcfFile"] ? config["vcfFile"].as<std::string>() : "";
        std::string vcfField = config["vcfField"] ? config["vcfField"].as<std::string>() : "GT";
        std::string bgenFile = config["bgenFile"] ? config["bgenFile"].as<std::string>() : "";
        std::string bgenSampleFile = config["bgenSampleFile"] ? config["bgenSampleFile"].as<std::string>() : "";
        std::string pgenFile = config["pgenFile"] ? config["pgenFile"].as<std::string>() : "";
        std::string pvarFile = config["pvarFile"] ? config["pvarFile"].as<std::string>() : "";
        std::string psamFile = config["psamFile"] ? config["psamFile"].as<std::string>() : "";
        std::string outputFile = config["outputFile"].as<std::string>();

        // Optional keys with defaults
        double minMAF = config["minMAF"] ? config["minMAF"].as<double>() : 0.0;
        double minMAC = config["minMAC"] ? config["minMAC"].as<double>() : 0.5;
        double maxMissRate = config["maxMissRate"] ? config["maxMissRate"].as<double>() : 0.15;
        double minINFO = config["minINFO"] ? config["minINFO"].as<double>() : 0.0;
        std::string alleleOrder = config["AlleleOrder"] ? config["AlleleOrder"].as<std::string>() : "alt-first";
        // P3 fix (2026-05-09): align defaults with R's step2_SPAtests.R.
        // R defaults: dosage_zerod_cutoff=0.2, dosage_zerod_MAC_cutoff=10
        // cpp previously defaulted both to 0, which kept ~12k extra low-info
        // imputed markers in the output (where R zeros out the noise dosages
        // on MAC<=10 markers, dropping them post re-filter). Aligning defaults
        // closes the cpp-vs-R row-count gap on UKB chr6.
        double dosage_zerod_cutoff = config["dosage_zerod_cutoff"] ? config["dosage_zerod_cutoff"].as<double>() : 0.2;
        double dosage_zerod_MAC_cutoff = config["dosage_zerod_MAC_cutoff"] ? config["dosage_zerod_MAC_cutoff"].as<double>() : 10.0;
        std::string genoType = config["genoType"] ? config["genoType"].as<std::string>() : "plink";
        bool isImputation = config["isImputation"] ? config["isImputation"].as<bool>() : false;
        bool isMoreOutput = config["isMoreOutput"] ? config["isMoreOutput"].as<bool>() : false;
        int marker_chunksize = config["marker_chunksize"] ? config["marker_chunksize"].as<int>() : 10000;
        int nThreads = config["nThreads"] ? config["nThreads"].as<int>() : 1;
        if (nThreads < 1) nThreads = 1;
        g_nThreads = nThreads;
        omp_set_num_threads(g_nThreads);
        // Prevent BLAS oversubscription when OMP threads call BLAS underneath.
        openblas_set_num_threads(1);
        // Phase D (Wave 1.3): BGEN streaming decoder thread count.
        // Only used when genoType=="bgen" in single-variant mode.
        int bgenDecoders = config["bgenDecoders"] ? config["bgenDecoders"].as<int>() : 4;
        g_bgenDecoders = bgenDecoders;
        double MACCutoffforER = config["MACCutoffforER"] ? config["MACCutoffforER"].as<double>() : 4.0;
        bool isFirth = config["isFirth"] ? config["isFirth"].as<bool>() : false;

        // Weights beta parameters (default Beta(1,25))
        arma::vec weights_beta(2);
        if (config["weights_beta"] && config["weights_beta"].IsSequence() &&
            config["weights_beta"].size() == 2) {
            weights_beta(0) = config["weights_beta"][0].as<double>();
            weights_beta(1) = config["weights_beta"][1].as<double>();
        } else {
            weights_beta(0) = 1.0;
            weights_beta(1) = 25.0;
        }

        // Checkpoint configuration
        std::string checkpointDir = config["checkpointDir"] ? config["checkpointDir"].as<std::string>() : "";
        if (!checkpointDir.empty()) {
            g_writeCheckpoints = true;
            g_checkpointDir = checkpointDir;
            std::cout << "Checkpoints enabled, output to: " << checkpointDir << std::endl;
        }

        // ---- Region/gene-based testing config keys ----
        std::string groupFile = config["groupFile"] ? config["groupFile"].as<std::string>() : "";
        bool isRegionTest = !groupFile.empty();

        // Annotation list (e.g., ["lof", "lof;missense", "lof;missense;synonymous"])
        std::vector<std::string> annotationList;
        if (config["annotationList"] && config["annotationList"].IsSequence()) {
            for (size_t i = 0; i < config["annotationList"].size(); i++) {
                annotationList.push_back(config["annotationList"][i].as<std::string>());
            }
        }

        // Max MAF list (e.g., [0.0001, 0.001, 0.01])
        arma::vec maxMAFList;
        if (config["maxMAFList"] && config["maxMAFList"].IsSequence()) {
            maxMAFList.set_size(config["maxMAFList"].size());
            for (size_t i = 0; i < config["maxMAFList"].size(); i++) {
                maxMAFList(i) = config["maxMAFList"][i].as<double>();
            }
        }

        std::string regionTestType = config["regionTestType"] ? config["regionTestType"].as<std::string>() : "SKATO";
        double MACCutoff_to_CollapseUltraRare = config["MACCutoff_to_CollapseUltraRare"] ?
            config["MACCutoff_to_CollapseUltraRare"].as<double>() : 10.0;
        int markers_per_chunk_in_groupTest = config["markers_per_chunk_in_groupTest"] ?
            config["markers_per_chunk_in_groupTest"].as<int>() : 500;
        int groups_per_chunk = config["groups_per_chunk"] ?
            config["groups_per_chunk"].as<int>() : 100;

        // r_corr: 0 means SKAT-O (optimal.adj), 1 means BURDEN
        double r_corr_val = config["r_corr"] ? config["r_corr"].as<double>() : 0.0;
        bool isSingleInGroupTest = config["isSingleInGroupTest"] ?
            config["isSingleInGroupTest"].as<bool>() : true;
        bool isOutputMarkerList = config["isOutputMarkerList"] ?
            config["isOutputMarkerList"].as<bool>() : false;
        unsigned int max_markers_region = config["max_markers_region"] ?
            config["max_markers_region"].as<unsigned int>() : 100000;
        double min_gourpmac_for_burdenonly = config["min_gourpmac_for_burdenonly"] ?
            config["min_gourpmac_for_burdenonly"].as<double>() : 5.0;

        // ---- LD matrix generation config ----
        bool isLDMatrix = config["isLDMatrix"] ? config["isLDMatrix"].as<bool>() : false;
        double ldmat_maxMAF = config["ldmat_maxMAF"] ? config["ldmat_maxMAF"].as<double>() : 0.5;

        // ---- Conditional analysis config ----
        // condition: list of marker IDs to condition on (rsIDs or chr:pos:ref:alt format)
        std::vector<std::string> conditionMarkerIDs;
        if (config["condition"] && config["condition"].IsSequence()) {
            for (size_t i = 0; i < config["condition"].size(); i++) {
                conditionMarkerIDs.push_back(config["condition"][i].as<std::string>());
            }
        } else if (config["condition"] && config["condition"].IsScalar()) {
            // Also support comma-separated string like R: "rs1,rs3,rs5"
            std::string condStr = config["condition"].as<std::string>();
            if (!condStr.empty()) {
                std::stringstream ss(condStr);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    // Trim whitespace
                    size_t start = token.find_first_not_of(" \t");
                    size_t end = token.find_last_not_of(" \t");
                    if (start != std::string::npos) {
                        conditionMarkerIDs.push_back(token.substr(start, end - start + 1));
                    }
                }
            }
        }
        // Condition weights (optional, like R's weights_for_condition)
        arma::vec condition_weights;
        if (config["weights_for_condition"] && config["weights_for_condition"].IsSequence()) {
            condition_weights.set_size(config["weights_for_condition"].size());
            for (size_t i = 0; i < config["weights_for_condition"].size(); i++) {
                condition_weights(i) = config["weights_for_condition"][i].as<double>();
            }
        }

        // Build r_corr vector: if r_corr_val == 0, use SKAT-O optimal.adj rho grid
        // If r_corr_val == 1, use BURDEN (single rho = 1)
        arma::vec r_corr_vec;
        if (r_corr_val == 0.0) {
            // SKAT-O optimal.adj rho grid (matches SKAT:::SKAT_Check_Method("optimal.adj", 0))
            // R returns: {0, 0.01, 0.04, 0.09, 0.25, 0.5, 1}
            r_corr_vec = {0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0};
            regionTestType = "SKATO";
            // For SKAT-O, single-variant results are always output
            isSingleInGroupTest = true;
        } else if (r_corr_val == 1.0) {
            r_corr_vec = {1.0};
            regionTestType = "BURDEN";
        } else {
            throw std::runtime_error("r_corr must be either 0 (SKAT-O) or 1 (BURDEN)");
        }

        // Validate region test config
        if (isRegionTest) {
            if (annotationList.empty()) {
                throw std::runtime_error("Region testing requires annotationList in config");
            }
            if (maxMAFList.n_elem == 0) {
                throw std::runtime_error("Region testing requires maxMAFList in config");
            }
        }

        // Print configuration
        std::cout << std::endl;
        std::cout << "Configuration:" << std::endl;
        std::cout << "  modelFile:         " << modelFile << std::endl;
        std::cout << "  varianceRatioFile: " << varianceRatioFile << std::endl;
        if (genoType == "plink") {
            std::cout << "  plinkFile:         " << plinkPrefix << std::endl;
        } else if (genoType == "vcf") {
            std::cout << "  vcfFile:           " << vcfFile << std::endl;
            std::cout << "  vcfField:          " << vcfField << std::endl;
        } else if (genoType == "bgen") {
            std::cout << "  bgenFile:          " << bgenFile << std::endl;
            std::cout << "  bgenSampleFile:    " << bgenSampleFile << std::endl;
        } else if (genoType == "pgen") {
            std::cout << "  pgenFile:          " << pgenFile << std::endl;
            std::cout << "  pvarFile:          " << pvarFile << std::endl;
            std::cout << "  psamFile:          " << psamFile << std::endl;
        }
        std::cout << "  outputFile:        " << outputFile << std::endl;
        std::cout << "  genoType:          " << genoType << std::endl;
        std::cout << "  alleleOrder:       " << alleleOrder << std::endl;
        std::cout << "  minMAF:            " << minMAF << std::endl;
        std::cout << "  minMAC:            " << minMAC << std::endl;
        std::cout << "  maxMissRate:       " << maxMissRate << std::endl;
        std::cout << "  minINFO:           " << minINFO << std::endl;
        std::cout << "  isImputation:      " << std::boolalpha << isImputation << std::endl;
        std::cout << "  isMoreOutput:      " << std::boolalpha << isMoreOutput << std::endl;
        std::cout << "  isFirth:           " << std::boolalpha << isFirth << std::endl;
        std::cout << "  MACCutoffforER:    " << MACCutoffforER << std::endl;
        std::cout << "  weights_beta:      [" << weights_beta(0) << ", " << weights_beta(1) << "]" << std::endl;
        std::cout << "  marker_chunksize:  " << marker_chunksize << std::endl;
        std::cout << "  nThreads:          " << g_nThreads << std::endl;
        if (isRegionTest) {
            std::cout << std::endl;
            std::cout << "  --- Region testing ---" << std::endl;
            std::cout << "  groupFile:         " << groupFile << std::endl;
            std::cout << "  regionTestType:    " << regionTestType << std::endl;
            std::cout << "  annotationList:    [";
            for (size_t i = 0; i < annotationList.size(); i++) {
                if (i > 0) std::cout << ", ";
                std::cout << annotationList[i];
            }
            std::cout << "]" << std::endl;
            std::cout << "  maxMAFList:        [";
            for (size_t i = 0; i < maxMAFList.n_elem; i++) {
                if (i > 0) std::cout << ", ";
                std::cout << maxMAFList(i);
            }
            std::cout << "]" << std::endl;
            std::cout << "  MACCutoff_to_CollapseUltraRare: " << MACCutoff_to_CollapseUltraRare << std::endl;
            std::cout << "  markers_per_chunk_in_groupTest:  " << markers_per_chunk_in_groupTest << std::endl;
            std::cout << "  groups_per_chunk:  " << groups_per_chunk << std::endl;
            std::cout << "  r_corr:            " << r_corr_val << std::endl;
            std::cout << "  isSingleInGroupTest:   " << std::boolalpha << isSingleInGroupTest << std::endl;
            std::cout << "  isOutputMarkerList:    " << std::boolalpha << isOutputMarkerList << std::endl;
            std::cout << "  max_markers_region:    " << max_markers_region << std::endl;
        }
        if (isLDMatrix) {
            std::cout << std::endl;
            std::cout << "  --- LD Matrix generation ---" << std::endl;
            std::cout << "  isLDMatrix:        " << std::boolalpha << isLDMatrix << std::endl;
            std::cout << "  ldmat_maxMAF:      " << ldmat_maxMAF << std::endl;
        }
        if (!conditionMarkerIDs.empty()) {
            std::cout << std::endl;
            std::cout << "  --- Conditional analysis ---" << std::endl;
            std::cout << "  condition markers: [";
            for (size_t i = 0; i < conditionMarkerIDs.size(); i++) {
                if (i > 0) std::cout << ", ";
                std::cout << conditionMarkerIDs[i];
            }
            std::cout << "]" << std::endl;
            if (condition_weights.n_elem > 0) {
                std::cout << "  weights_for_condition: [";
                for (size_t i = 0; i < condition_weights.n_elem; i++) {
                    if (i > 0) std::cout << ", ";
                    std::cout << condition_weights(i);
                }
                std::cout << "]" << std::endl;
            }
        }
        std::cout << std::endl;

        // ---- 3. Load null model from Step 1 ----
        std::cout << "===== Loading null model =====" << std::endl;
        NullModelData nullModel = loadNullModel(modelFile, varianceRatioFile);

        std::cout << "  Trait type:   " << nullModel.traitType << std::endl;
        std::cout << "  Sample size:  " << nullModel.n << std::endl;
        std::cout << "  Covariates:   " << nullModel.p << std::endl;
        std::cout << "  tau[0]:       " << nullModel.tau0 << std::endl;
        std::cout << "  tau[1]:       " << nullModel.tau1 << std::endl;
        std::cout << "  SPA_Cutoff:   " << nullModel.SPA_Cutoff << std::endl;
        std::cout << "  flagSparseGRM:" << std::boolalpha << nullModel.flagSparseGRM << std::endl;
        std::cout << "  isFastTest:   " << std::boolalpha << nullModel.isFastTest << std::endl;
        std::cout << "  isCondition:  " << std::boolalpha << nullModel.isCondition << std::endl;
        std::cout << "  is_Firth_beta:" << std::boolalpha << nullModel.is_Firth_beta << std::endl;
        std::cout << std::endl;

        // ---- 4. Set global variables ----
        std::cout << "===== Setting global variables =====" << std::endl;
        setAssocTest_GlobalVarsInCPP(
            nullModel.impute_method,
            maxMissRate,
            minMAF,
            minMAC,
            minINFO,
            dosage_zerod_cutoff,
            dosage_zerod_MAC_cutoff,
            weights_beta,
            outputFile,
            MACCutoffforER);

        setMarker_GlobalVarsInCPP(isMoreOutput, marker_chunksize);

        // Override Firth settings from null model, then let YAML config override
        // (mirrors R Step-2's --is_Firth_beta / --pCutoffforFirth CLI behavior).
        g_is_Firth_beta = nullModel.is_Firth_beta;
        g_pCutoffforFirth = nullModel.pCutoffforFirth;
        if (config["is_Firth_beta"]) {
            nullModel.is_Firth_beta = config["is_Firth_beta"].as<bool>();
            g_is_Firth_beta = nullModel.is_Firth_beta;
            std::cout << "  is_Firth_beta overridden from config: "
                      << std::boolalpha << g_is_Firth_beta << std::endl;
        }
        if (config["pCutoffforFirth"]) {
            nullModel.pCutoffforFirth = config["pCutoffforFirth"].as<double>();
            g_pCutoffforFirth = nullModel.pCutoffforFirth;
            std::cout << "  pCutoffforFirth overridden from config: "
                      << g_pCutoffforFirth << std::endl;
        }

        // ---- Override cateVarRatioMinMACVecExclude / cateVarRatioMaxMACVecInclude ----
        // R SAIGE defaults: cateVarRatioMinMACVecExclude = c(10.5, 20.5)
        //                   cateVarRatioMaxMACVecInclude = c(20.5, N)
        // These control which MAC threshold triggers the PCG (sparse GRM) path
        // during re-evaluation for markers with p < pval_cutoff_for_fastTest.
        // If not specified in config, the values from the VR file are used.
        if (config["cateVarRatioMinMACVecExclude"] && config["cateVarRatioMinMACVecExclude"].IsSequence()) {
            size_t n = config["cateVarRatioMinMACVecExclude"].size();
            nullModel.cateVarRatioMinMACVecExclude.set_size(n);
            for (size_t i = 0; i < n; i++) {
                nullModel.cateVarRatioMinMACVecExclude(i) = config["cateVarRatioMinMACVecExclude"][i].as<double>();
            }
            std::cout << "  cateVarRatioMinMACVecExclude overridden from config: [";
            for (size_t i = 0; i < n; i++) {
                if (i > 0) std::cout << ", ";
                std::cout << nullModel.cateVarRatioMinMACVecExclude(i);
            }
            std::cout << "]" << std::endl;
        }
        if (config["cateVarRatioMaxMACVecInclude"] && config["cateVarRatioMaxMACVecInclude"].IsSequence()) {
            size_t n = config["cateVarRatioMaxMACVecInclude"].size();
            nullModel.cateVarRatioMaxMACVecInclude.set_size(n);
            for (size_t i = 0; i < n; i++) {
                nullModel.cateVarRatioMaxMACVecInclude(i) = config["cateVarRatioMaxMACVecInclude"][i].as<double>();
            }
            std::cout << "  cateVarRatioMaxMACVecInclude overridden from config: [";
            for (size_t i = 0; i < n; i++) {
                if (i > 0) std::cout << ", ";
                std::cout << nullModel.cateVarRatioMaxMACVecInclude(i);
            }
            std::cout << "]" << std::endl;
        }

        // ---- Override isnoadjCov from YAML config ----
        // The JSON null model may have isnoadjCov=false, but the YAML config
        // can override it (e.g., for testing the scoreTestFast_noadjCov path).
        if (config["isnoadjCov"]) {
            nullModel.isnoadjCov = config["isnoadjCov"].as<bool>();
            std::cout << "  isnoadjCov overridden from config: " << std::boolalpha << nullModel.isnoadjCov << std::endl;
        }

        // ---- Set isCondition from YAML condition markers ----
        // Note: isCondition is set to true if condition markers are specified.
        // The actual condition_genoIndex will be populated after PLINK setup
        // (since we need the genotype file to look up marker indices).
        if (!conditionMarkerIDs.empty()) {
            nullModel.isCondition = true;
            // Use a dummy condition_genoIndex for now; will be populated after PLINK setup.
            // SAIGEClass just needs to know m_numMarker_cond.
            nullModel.condition_genoIndex.resize(conditionMarkerIDs.size(), 0);
            std::cout << "  Conditional analysis enabled with " << conditionMarkerIDs.size()
                      << " conditioning markers." << std::endl;
        }

        // ---- 5. Construct SAIGEClass from NullModelData ----
        std::cout << "===== Constructing SAIGEClass =====" << std::endl;
        setSAIGEobjInCPP(
            nullModel.XVX,
            nullModel.XXVX_inv,
            nullModel.XV,
            nullModel.XVX_inv_XV,
            nullModel.Sigma_iXXSigma_iX,
            nullModel.X,
            nullModel.S_a,
            nullModel.res,
            nullModel.mu2,
            nullModel.mu,
            nullModel.varRatio_sparse,
            nullModel.varRatio_null,
            nullModel.varRatio_null_noXadj,
            nullModel.cateVarRatioMinMACVecExclude,
            nullModel.cateVarRatioMaxMACVecInclude,
            nullModel.SPA_Cutoff,
            nullModel.tauvec,
            nullModel.traitType,
            nullModel.y,
            nullModel.impute_method,
            nullModel.flagSparseGRM,
            nullModel.isFastTest,
            nullModel.isnoadjCov,
            nullModel.pval_cutoff_for_fastTest,
            nullModel.locationMat,
            nullModel.valueVec,
            nullModel.dimNum,
            nullModel.isCondition,
            nullModel.condition_genoIndex,
            nullModel.is_Firth_beta,
            nullModel.pCutoffforFirth,
            nullModel.offset,
            nullModel.resout);

        std::cout << "  SAIGEClass constructed successfully." << std::endl;
        std::cout << "  n = " << ptr_gSAIGEobj->m_n << ", p = " << ptr_gSAIGEobj->m_p << std::endl;
        std::cout << std::endl;

        // ---- 6. Set up genotype reader (PLINK, VCF, or BGEN) ----
        std::cout << "===== Setting up genotype reader =====" << std::endl;
        uint32_t numMarkers = 0;

        if (genoType == "plink") {
            std::string bimFile = plinkPrefix + ".bim";
            std::string famFile = plinkPrefix + ".fam";
            std::string bedFile = plinkPrefix + ".bed";

            setPLINKobjInCPP(bimFile, famFile, bedFile, nullModel.sampleIDs, alleleOrder);
            numMarkers = ptr_gPLINKobj->getM();
        } else if (genoType == "vcf") {
            setVCFobjInCPP(vcfFile, vcfField, nullModel.sampleIDs);
            // Pre-scan to count markers (needed for index generation)
            numMarkers = ptr_gVCFobj->prescanMarkerCount();
            // Reset to beginning for actual reading
            ptr_gVCFobj->resetFile();
            // Re-set sample mapping after reset
            ptr_gVCFobj->setPosSampleInVcf(nullModel.sampleIDs);
        } else if (genoType == "bgen") {
            // Read sample IDs from .sample file
            std::vector<std::string> sampleInBgen;
            std::ifstream sampleFile(bgenSampleFile);
            if (!sampleFile.is_open()) {
                throw std::runtime_error("Cannot open BGEN sample file: " + bgenSampleFile);
            }
            std::string line;
            int lineNum = 0;
            while (std::getline(sampleFile, line)) {
                lineNum++;
                // Skip first two header lines of .sample file
                if (lineNum <= 2) continue;
                // First column is FID, second is IID
                std::istringstream iss(line);
                std::string fid, iid;
                if (iss >> fid >> iid) {
                    sampleInBgen.push_back(iid);
                }
            }
            sampleFile.close();
            std::cout << "Read " << sampleInBgen.size() << " sample IDs from " << bgenSampleFile << std::endl;

            setBGENobjInCPP(bgenFile, sampleInBgen, nullModel.sampleIDs, alleleOrder);
            numMarkers = ptr_gBGENobj->getM0();
        } else if (genoType == "pgen") {
            setPGENobjInCPP(pgenFile, psamFile, pvarFile, nullModel.sampleIDs);
            numMarkers = ptr_gPGENobj->getM();
        }

        uint32_t numSamplesGeno = Unified_getSampleSizeinGeno(genoType);
        uint32_t numSamplesAnalysis = Unified_getSampleSizeinAnalysis(genoType);

        std::cout << "  Total markers in genotype file: " << numMarkers << std::endl;
        std::cout << "  Total samples in genotype file: " << numSamplesGeno << std::endl;
        std::cout << "  Samples in analysis:            " << numSamplesAnalysis << std::endl;
        std::cout << std::endl;

        // Verify sample size consistency
        if ((int)numSamplesAnalysis != nullModel.n) {
            std::cerr << "WARNING: Sample size mismatch. Null model n=" << nullModel.n
                      << " but genotype file analysis n=" << numSamplesAnalysis << std::endl;
        }

        // ---- 6b. Set up conditional analysis (if condition markers specified) ----
        if (!conditionMarkerIDs.empty()) {
            std::cout << "===== Setting up conditional analysis =====" << std::endl;
            std::cout << "  Looking up " << conditionMarkerIDs.size()
                      << " conditioning marker(s) in genotype file..." << std::endl;

            // Build lookup maps: rsID -> index and chr:pos:ref:alt -> index
            // Note: For VCF, conditional analysis requires a pre-scan of all markers first.
            // Currently conditional analysis is only fully supported with PLINK input.
            if (genoType == "vcf" || genoType == "bgen") {
                throw std::runtime_error(
                    "Conditional analysis is not yet supported with " + genoType + " input. "
                    "Please use PLINK format for conditional analysis.");
            }
            std::unordered_map<std::string, uint32_t> markerNameMap = ptr_gPLINKobj->getMarkerNameToIndex();
            std::unordered_map<std::string, uint32_t> markerIDMap = ptr_gPLINKobj->getMarkerIDToIndex();

            std::vector<uint32_t> condGenoIndices;
            for (size_t i = 0; i < conditionMarkerIDs.size(); i++) {
                const std::string& markerID = conditionMarkerIDs[i];
                uint32_t idx = 0;
                bool found = false;

                // Try rsID lookup first (m_MarkerInPlink, column 2 of .bim)
                auto it1 = markerNameMap.find(markerID);
                if (it1 != markerNameMap.end()) {
                    idx = it1->second;
                    found = true;
                }

                // Try chr:pos:ref:alt lookup if rsID not found
                if (!found) {
                    auto it2 = markerIDMap.find(markerID);
                    if (it2 != markerIDMap.end()) {
                        idx = it2->second;
                        found = true;
                    }
                }

                if (!found) {
                    throw std::runtime_error(
                        "Conditioning marker '" + markerID +
                        "' not found in genotype file. "
                        "Use rsID or chr:pos:ref:alt format.");
                }

                condGenoIndices.push_back(idx);
                std::cout << "    " << markerID << " -> index " << idx << std::endl;
            }

            // Update SAIGEClass with actual condition indices
            ptr_gSAIGEobj->m_condition_genoIndex = condGenoIndices;
            ptr_gSAIGEobj->m_numMarker_cond = condGenoIndices.size();

            // Set up condition weights
            arma::vec cond_weights;
            if (condition_weights.n_elem > 0) {
                if (condition_weights.n_elem != condGenoIndices.size()) {
                    throw std::runtime_error(
                        "weights_for_condition length (" +
                        std::to_string(condition_weights.n_elem) +
                        ") does not match number of conditioning markers (" +
                        std::to_string(condGenoIndices.size()) + ")");
                }
                cond_weights = condition_weights;
            } else {
                // Default: zero weights (will use Beta(MAF, 1, 25) in assign_conditionMarkers_factors)
                cond_weights.zeros(condGenoIndices.size());
            }

            // Read conditioning marker genotypes and compute factors
            std::cout << "  Reading conditioning marker genotypes..." << std::endl;
            assign_conditionMarkers_factors(genoType, condGenoIndices,
                                            numSamplesAnalysis, cond_weights);
            std::cout << std::endl;
        }

        // ================================================================
        // Trifurcate: single-variant / LD matrix / region testing
        // ================================================================

        if (isLDMatrix && isRegionTest) {
            // ============================================================
            // LD MATRIX GENERATION PATH
            // Computes G^T * G for each region in the group file.
            // Ported from SAIGE/R/SAIGE_SPATest_Region_LDMat.R + LDmat.cpp
            // ============================================================

            std::cout << "===== LD Matrix generation mode =====" << std::endl;
            std::cout << std::endl;

            // Set LDmat global variables
            setGlobalVarsInCPP_LDmat(
                nullModel.impute_method,
                dosage_zerod_cutoff,
                dosage_zerod_MAC_cutoff,
                maxMissRate,
                ldmat_maxMAF,
                minMAF,
                minMAC,
                minINFO,
                max_markers_region,
                outputFile);

            // Build marker ID to index map
            std::cout << "===== Building marker ID to index map =====" << std::endl;
            std::unordered_map<std::string, uint32_t> markerIDToIndex = Unified_getMarkerIDToIndex(genoType);
            std::cout << "  Built map with " << markerIDToIndex.size() << " entries." << std::endl;
            std::cout << std::endl;

            // Check group file
            std::cout << "===== Checking group file =====" << std::endl;
            GroupFileInfo gfInfo = checkGroupFile(groupFile);
            int nRegions = gfInfo.nRegions;
            bool is_weight_included = gfInfo.is_weight_included;
            int nline_per_gene = is_weight_included ? 3 : 2;

            std::cout << "  Group file: " << groupFile << std::endl;
            std::cout << "  Number of regions: " << nRegions << std::endl;
            std::cout << "  Weights included: " << std::boolalpha << is_weight_included << std::endl;
            std::cout << "  Lines per gene: " << nline_per_gene << std::endl;
            std::cout << std::endl;

            // Open LDmat output files
            bool isOpenMarkerInfo = openOutfile_single_LDmat(false);
            if (!isOpenMarkerInfo) {
                throw std::runtime_error("Cannot open marker info output file: " + outputFile + ".marker_info.txt");
            }
            std::cout << "  Marker info file opened: " << outputFile << ".marker_info.txt" << std::endl;

            bool isOpenLDmat = openOutfile_LDmat(false);
            if (!isOpenLDmat) {
                throw std::runtime_error("Cannot open LDmat output file: " + outputFile + ".LDmat.txt");
            }
            std::cout << "  LDmat file opened: " << outputFile << ".LDmat.txt" << std::endl;

            bool isOpenIndex = openOutfile_index_LDmat(false);
            if (!isOpenIndex) {
                throw std::runtime_error("Cannot open index output file: " + outputFile + ".index.txt");
            }
            std::cout << "  Index file opened: " << outputFile << ".index.txt" << std::endl;
            std::cout << std::endl;

            // Loop over regions
            std::cout << "===== Starting LD matrix computation =====" << std::endl;
            arma::vec timeStart = getTime();

            std::ifstream gf(groupFile);
            if (!gf.is_open()) {
                throw std::runtime_error("Cannot open group file: " + groupFile);
            }

            unsigned int t_n = (unsigned int)numSamplesAnalysis;
            int regionsProcessed = 0;
            int regionsSkipped = 0;

            while (regionsProcessed + regionsSkipped < nRegions) {
                int remaining = nRegions - regionsProcessed - regionsSkipped;
                int nregions_to_read = std::min(groups_per_chunk, remaining);

                std::vector<RegionData> regionChunk = readRegionChunk(
                    gf, nregions_to_read, nline_per_gene, annotationList, markerIDToIndex);

                for (int r = 0; r < (int)regionChunk.size(); r++) {
                    RegionData& region = regionChunk[r];
                    int totalIdx = regionsProcessed + regionsSkipped + 1;

                    if (region.variantIDs.empty() || region.annoVec.empty()) {
                        std::cout << "  Skipping region " << region.regionName
                                  << " (" << totalIdx << "/" << nRegions
                                  << "): no matching variants." << std::endl;
                        regionsSkipped++;
                        continue;
                    }

                    std::cout << "  Computing LD matrix for region " << region.regionName
                              << " (" << totalIdx << "/" << nRegions
                              << "), " << region.variantIDs.size() << " variants."
                              << std::endl;

                    // Build annoIndicatorMat as arma::mat (LDmatRegionInCPP expects arma::mat)
                    arma::mat annoIndicatorMat = arma::conv_to<arma::mat>::from(region.annoIndicatorMat);

                    LDmatRegionInCPP(
                        genoType,
                        region.genoIndex_prev,
                        region.genoIndex,
                        annoIndicatorMat,
                        outputFile,
                        t_n,
                        isImputation,
                        region.annoVec,
                        region.regionName);

                    regionsProcessed++;

                    if (regionsProcessed % 100 == 0) {
                        std::cout << "    Processed " << regionsProcessed << " regions ("
                                  << regionsSkipped << " skipped)." << std::endl;
                    }
                }
            }

            gf.close();

            arma::vec timeEnd = getTime();
            printTime(timeStart, timeEnd, "complete LD matrix computation");
            std::cout << std::endl;

            std::cout << "  Total regions processed: " << regionsProcessed << std::endl;
            std::cout << "  Total regions skipped:   " << regionsSkipped << std::endl;
            std::cout << std::endl;

            // Close LDmat output files
            closeOutfile_single_LDmat();
            closeOutfile_LDmat();
            closeOutfile_index_LDmat();

            std::cout << "  LD matrix output:    " << outputFile << ".LDmat.txt" << std::endl;
            std::cout << "  Marker info output:  " << outputFile << ".marker_info.txt" << std::endl;
            std::cout << "  Index output:        " << outputFile << ".index.txt" << std::endl;

        } else if (!isRegionTest) {
            // ============================================================
            // SINGLE-VARIANT TESTING PATH
            // ============================================================

            // ---- 7a. Build genotype index vectors ----
            // In the R pipeline, these are passed as vectors of marker indices.
            // For the standalone version, we test ALL markers in the file sequentially.
            //
            // For PLINK/VCF: genoIndex = 0, 1, 2, ... (marker index for seeking)
            // For BGEN: genoIndex = 0 for all (sequential read, no seeking needed;
            //           BGEN uses byte positions for seeking, but we read sequentially)
            std::cout << "===== Building marker index vectors =====" << std::endl;
            std::vector<std::string> genoIndex(numMarkers);
            std::vector<std::string> genoIndex_prev(numMarkers);
            if (genoType == "bgen") {
                // BGEN: all zeros = "just read next variant, no seeking"
                for (uint32_t i = 0; i < numMarkers; i++) {
                    genoIndex[i] = "0";
                    genoIndex_prev[i] = "0";
                }
            } else {
                for (uint32_t i = 0; i < numMarkers; i++) {
                    genoIndex[i] = std::to_string(i);
                    if (i > 0) {
                        genoIndex_prev[i] = std::to_string(i - 1);
                    } else {
                        genoIndex_prev[i] = "0";
                    }
                }
            }

            std::cout << "  Will test " << numMarkers << " markers." << std::endl;
            std::cout << std::endl;

            // ---- 8a. Open output file ----
            std::cout << "===== Opening output file =====" << std::endl;
            bool isopen = openOutfile_single(nullModel.traitType, isImputation, false, isMoreOutput);
            if (!isopen) {
                throw std::runtime_error("Cannot open output file: " + outputFile);
            }
            std::cout << "  Output file opened: " << g_outputFilePrefixSingle << std::endl;
            std::cout << std::endl;

            // ---- 9a. Run single-variant testing ----
            std::cout << "===== Starting single-variant testing =====" << std::endl;
            arma::vec timeStart = getTime();

            mainMarkerInCPP(
                genoType,
                nullModel.traitType,
                genoIndex_prev,
                genoIndex,
                isMoreOutput,
                isImputation,
                isFirth);

            arma::vec timeEnd = getTime();
            printTime(timeStart, timeEnd, "complete single-variant testing");
            std::cout << std::endl;

            // ---- 10a. Close output file ----
            OutFile_single.close();

        } else {
            // ============================================================
            // REGION/GENE-BASED TESTING PATH
            // ============================================================

            std::cout << "===== Region/gene-based testing mode =====" << std::endl;
            std::cout << std::endl;

            // ---- 7b. Set region global variables ----
            setRegion_GlobalVarsInCPP(
                maxMAFList,
                max_markers_region,
                MACCutoff_to_CollapseUltraRare,
                min_gourpmac_for_burdenonly);

            // ---- 8b. Build marker ID to index map ----
            std::cout << "===== Building marker ID to index map =====" << std::endl;
            std::unordered_map<std::string, uint32_t> markerIDToIndex = Unified_getMarkerIDToIndex(genoType);
            std::cout << "  Built map with " << markerIDToIndex.size() << " entries." << std::endl;
            std::cout << std::endl;

            // ---- 9b. Check group file ----
            std::cout << "===== Checking group file =====" << std::endl;
            GroupFileInfo gfInfo = checkGroupFile(groupFile);
            int nRegions = gfInfo.nRegions;
            bool is_weight_included = gfInfo.is_weight_included;
            int nline_per_gene = is_weight_included ? 3 : 2;

            std::cout << "  Group file: " << groupFile << std::endl;
            std::cout << "  Number of regions: " << nRegions << std::endl;
            std::cout << "  Weights included: " << std::boolalpha << is_weight_included << std::endl;
            std::cout << "  Lines per gene: " << nline_per_gene << std::endl;
            std::cout << std::endl;

            // ---- 10b. Determine test type and open output files ----
            // Mirrors R logic: r.corr == 0 -> SKAT-O, r.corr == 1 -> BURDEN
            if (regionTestType == "BURDEN") {
                std::cout << "BURDEN test will be performed." << std::endl;
                bool isOpenOutFile = openOutfile(nullModel.traitType, false);
                if (!isOpenOutFile) {
                    throw std::runtime_error("Cannot open region output file: " + g_outputFilePrefixGroup);
                }
                std::cout << "  Region output file opened: " << g_outputFilePrefixGroup << std::endl;
            } else {
                // SKAT-O or SKAT
                std::cout << "SKAT-O test will be performed. P-values for BURDEN and SKAT will also be output." << std::endl;
                bool isOpenOutFile = openOutfile_SKATO(nullModel.traitType, false);
                if (!isOpenOutFile) {
                    throw std::runtime_error("Cannot open region output file: " + g_outputFilePrefixGroup);
                }
                std::cout << "  Region output file opened: " << g_outputFilePrefixGroup << std::endl;
            }

            if (isSingleInGroupTest) {
                std::cout << "  Single-variant results within groups will be output." << std::endl;
                bool isOpenSingle = openOutfile_singleinGroup(
                    nullModel.traitType, isImputation, false, isMoreOutput);
                if (!isOpenSingle) {
                    throw std::runtime_error("Cannot open single-in-group output file: " + g_outputFilePrefixSingleInGroup);
                }
                std::cout << "  Single-in-group output file opened: " << g_outputFilePrefixSingleInGroup << std::endl;
            } else {
                std::cout << "  Single-variant results within groups will NOT be output." << std::endl;
            }
            std::cout << std::endl;

            // ---- 11b. Allocate P1Mat and P2Mat ----
            unsigned int t_n = (unsigned int)nullModel.n;
            arma::mat P1Mat, P2Mat;
            if (regionTestType != "BURDEN") {
                P1Mat.set_size(markers_per_chunk_in_groupTest, t_n);
                P2Mat.set_size(t_n, markers_per_chunk_in_groupTest);
                P1Mat.zeros();
                P2Mat.zeros();
                std::cout << "  P1Mat size: " << P1Mat.n_rows << " x " << P1Mat.n_cols << std::endl;
                std::cout << "  P2Mat size: " << P2Mat.n_rows << " x " << P2Mat.n_cols << std::endl;
            } else {
                // BURDEN doesn't need P1Mat/P2Mat (1x1 placeholder)
                P1Mat.set_size(1, 1);
                P2Mat.set_size(1, 1);
                P1Mat.zeros();
                P2Mat.zeros();
            }

            // ---- 12b. Loop over regions ----
            std::cout << "===== Starting region-based testing =====" << std::endl;
            arma::vec timeStart = getTime();

            std::ifstream gf(groupFile);
            if (!gf.is_open()) {
                throw std::runtime_error("Cannot open group file: " + groupFile);
            }

            int regionsProcessed = 0;
            int regionsSkipped = 0;

            while (regionsProcessed + regionsSkipped < nRegions) {
                // Determine how many regions to read in this chunk
                int remaining = nRegions - regionsProcessed - regionsSkipped;
                int nregions_to_read = std::min(groups_per_chunk, remaining);

                // Read a chunk of regions from the group file
                std::vector<RegionData> regionChunk = readRegionChunk(
                    gf, nregions_to_read, nline_per_gene, annotationList, markerIDToIndex);

                for (int r = 0; r < (int)regionChunk.size(); r++) {
                    RegionData& region = regionChunk[r];
                    int totalIdx = regionsProcessed + regionsSkipped + 1;

                    // Skip regions with no matching variants
                    if (region.variantIDs.empty() || region.annoVec.empty()) {
                        std::cout << "  Skipping region " << region.regionName
                                  << " (" << totalIdx << "/" << nRegions
                                  << "): no matching variants." << std::endl;
                        regionsSkipped++;
                        continue;
                    }

                    std::cout << "  Analyzing region " << region.regionName
                              << " (" << totalIdx << "/" << nRegions
                              << "), " << region.variantIDs.size() << " variants."
                              << std::endl;

                    // Set sparseGRM flag for this region
                    if (!nullModel.isFastTest) {
                        ptr_gSAIGEobj->set_flagSparseGRM_cur(nullModel.flagSparseGRM);
                    } else {
                        ptr_gSAIGEobj->set_flagSparseGRM_cur(false);
                    }

                    // Build weight vector
                    arma::vec weightVec;
                    if (!region.weights.empty()) {
                        weightVec.set_size(region.weights.size());
                        for (size_t w = 0; w < region.weights.size(); w++) {
                            weightVec(w) = region.weights[w];
                        }
                    } else {
                        // No weights provided in group file: use zeros to signal
                        // that default Beta(MAF, 1, 25) weights should be used.
                        // This matches R SAIGE behavior where WEIGHT=c(0) when
                        // the group file has no weight line (nline_per_gene != 3).
                        weightVec = arma::zeros<arma::vec>(1);
                    }

                    // Call mainRegionInCPP for this region
                    mainRegionInCPP(
                        genoType,
                        nullModel.traitType,
                        region,
                        maxMAFList,
                        outputFile,
                        t_n,
                        P1Mat,
                        P2Mat,
                        regionTestType,
                        isImputation,
                        weightVec,
                        isSingleInGroupTest,
                        isMoreOutput,
                        r_corr_vec,
                        nullModel.mu);

                    regionsProcessed++;

                    // Progress report
                    if (regionsProcessed % 100 == 0) {
                        std::cout << "    Processed " << regionsProcessed << " regions ("
                                  << regionsSkipped << " skipped)." << std::endl;
                    }
                }
            }

            gf.close();

            arma::vec timeEnd = getTime();
            printTime(timeStart, timeEnd, "complete region-based testing");
            std::cout << std::endl;

            std::cout << "  Total regions processed: " << regionsProcessed << std::endl;
            std::cout << "  Total regions skipped:   " << regionsSkipped << std::endl;
            std::cout << std::endl;

            // ---- 13b. Close output files ----
            OutFile.close();
            if (isSingleInGroupTest) {
                OutFile_singleInGroup.close();
            }

        } // end of region testing path

        // ---- Close genotype file ----
        closeGenoFile(genoType);

        // ---- Cleanup ----
        if (ptr_gSAIGEobj) {
            delete ptr_gSAIGEobj;
            ptr_gSAIGEobj = NULL;
        }

        // Memory usage report
        double vm_usage, resident_set;
        process_mem_usage(vm_usage, resident_set);
        std::cout << "===== Done =====" << std::endl;
        std::cout << "  Virtual memory:  " << vm_usage / (1024.0 * 1024.0) << " MB" << std::endl;
        std::cout << "  Resident memory: " << resident_set / (1024.0 * 1024.0) << " MB" << std::endl;
        if (isLDMatrix && isRegionTest) {
            std::cout << "  LD matrix output:    " << outputFile << ".LDmat.txt" << std::endl;
            std::cout << "  Marker info output:  " << outputFile << ".marker_info.txt" << std::endl;
            std::cout << "  Index output:        " << outputFile << ".index.txt" << std::endl;
        } else if (!isRegionTest) {
            std::cout << "  Output written to: " << g_outputFilePrefixSingle << std::endl;
        } else {
            std::cout << "  Region output: " << g_outputFilePrefixGroup << std::endl;
            if (isSingleInGroupTest) {
                std::cout << "  Single-in-group output: " << g_outputFilePrefixSingleInGroup << std::endl;
            }
        }
        std::cout << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "ERROR: Unknown exception occurred." << std::endl;
        return 1;
    }
}
