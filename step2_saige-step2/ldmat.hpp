// Standalone port of SAIGE/src/LDmat.hpp
// LD matrix computation for rare-variant meta-analysis
//
// Ported from: SAIGE/src/LDmat.hpp (~53 lines)
//
// Conversions from original SAIGE LDmat.hpp:
//   1. #include <RcppArmadillo.h>  -->  #include <armadillo>
//   2. // [[Rcpp::depends(RcppArmadillo)]]  --> removed
//   3. // [[Rcpp::export]]  --> removed
//   4. Added outputFile parameter to setGlobalVarsInCPP_LDmat signature
//      (was separate in hpp vs cpp in original)
//   5. Added openOutfile_LDmat and openOutfile_index_LDmat declarations
//   6. All algorithm logic preserved exactly

#ifndef LDMAT_HPP
#define LDMAT_HPP

#include <armadillo>
#include <string>
#include <vector>
#include <fstream>

// Forward declaration of RegionData from group_file.hpp
struct RegionData;

// Set global variables for LDmat computation
void setGlobalVarsInCPP_LDmat(std::string t_impute_method,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
                               double t_missing_cutoff,
                               double t_maxMAFLimit,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               unsigned int t_max_markers_region,
                               std::string t_outputFile);

// Compute LD matrix for one region
// Direct port from SAIGE/src/LDmat.cpp LDmatRegionInCPP()
void LDmatRegionInCPP(
    std::string t_genoType,
    std::vector<std::string>& t_genoIndex_prev,
    std::vector<std::string>& t_genoIndex,
    arma::mat& annoIndicatorMat,
    std::string t_outputFile,
    unsigned int t_n,
    bool t_isImputation,
    std::vector<std::string>& annoStringVec,
    std::string regionName);

// Write marker information for one region
void writeOutfile_single_LDmat(
    std::vector<std::string>& chrVec,
    std::vector<std::string>& posVec,
    std::vector<std::string>& refVec,
    std::vector<std::string>& altVec,
    std::vector<std::string>& infoVec,
    std::vector<double>& altCountsVec,
    std::vector<double>& missingRateVec,
    std::vector<uint32_t>& N_Vec,
    std::string regionName);

// Open output files
bool openOutfile_single_LDmat(bool isappend);
bool openOutfile_LDmat(bool isappend);
bool openOutfile_index_LDmat(bool isappend);

// Close output files
void closeOutfile_single_LDmat();
void closeOutfile_LDmat();
void closeOutfile_index_LDmat();

#endif
