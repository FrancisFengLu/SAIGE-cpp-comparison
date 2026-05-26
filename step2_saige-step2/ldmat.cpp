// Standalone port of SAIGE/src/LDmat.cpp
// LD matrix computation for rare-variant meta-analysis
//
// Ported from: SAIGE/src/LDmat.cpp (~580 lines)
//
// Conversions from original SAIGE LDmat.cpp:
//   1. #include <RcppArmadillo.h>  -->  #include <armadillo>
//   2. Rcpp::stop(...)             -->  throw std::runtime_error(...)
//   3. Rcpp::Rcout                 -->  std::cout
//   4. Rcpp::checkUserInterrupt()  -->  removed (no R event loop)
//   5. // [[Rcpp::export]]         -->  removed
//   6. Uses existing Unified_getOneMarker from genotype_reader.hpp
//   7. Uses existing imputeGenoAndFlip from UTIL.hpp
//   8. All algorithm logic preserved EXACTLY as in SAIGE

#include <armadillo>

#include <vector>
#include <thread>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <iostream>
#include <cmath>

#include "ldmat.hpp"
#include "genotype_reader.hpp"
#include "UTIL.hpp"

// ============================================================
// Global output file streams for LDmat (exact match from SAIGE)
// ============================================================
std::ofstream OutFile_single_LDmat;
std::string g_outputFilePrefixSingle_LDmat;
std::ofstream OutFile_index_LDmat;
std::string g_outputFilePrefixIndex_LDmat;
std::ofstream OutFile_LDmat;
std::string g_outputFilePrefix_LDmat;

// Global variables for analysis (exact match from SAIGE)
std::string g_impute_method_LDmat;      // "mean", "minor", or "drop"
double g_dosage_zerod_MAC_cutoff_LDmat;
double g_dosage_zerod_cutoff_LDmat;

double g_missingRate_cutoff_LDmat;
double g_maxMAFLimit_LDmat;
double g_marker_minMAF_cutoff_LDmat;
double g_marker_minMAC_cutoff_LDmat;
double g_marker_minINFO_cutoff_LDmat;

unsigned int g_region_maxMarkers_cutoff_LDmat;
unsigned int g_startpos;
unsigned int g_endpos;


// ============================================================
// setGlobalVarsInCPP_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 64-96
// ============================================================
void setGlobalVarsInCPP_LDmat(std::string t_impute_method,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,
                               double t_missing_cutoff,
                               double t_maxMAFLimit,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               unsigned int t_max_markers_region,
                               std::string t_outputFile)
{
    g_impute_method_LDmat = t_impute_method;
    g_dosage_zerod_cutoff_LDmat = t_dosage_zerod_cutoff;
    g_dosage_zerod_MAC_cutoff_LDmat = t_dosage_zerod_MAC_cutoff;

    g_missingRate_cutoff_LDmat = t_missing_cutoff;
    g_maxMAFLimit_LDmat = t_maxMAFLimit;

    g_marker_minMAF_cutoff_LDmat = t_min_maf_marker;
    g_marker_minMAC_cutoff_LDmat = t_min_mac_marker;
    g_marker_minINFO_cutoff_LDmat = t_min_info_marker;
    g_region_maxMarkers_cutoff_LDmat = t_max_markers_region;

    g_outputFilePrefixSingle_LDmat = t_outputFile + ".marker_info.txt";
    g_outputFilePrefixIndex_LDmat = t_outputFile + ".index.txt";
    g_outputFilePrefix_LDmat = t_outputFile + ".LDmat.txt";
    g_startpos = -1;
    g_endpos = -1;
}


// ============================================================
// LDmatRegionInCPP
// Direct port from SAIGE/src/LDmat.cpp lines 99-493
// Computes G^T * G (LD matrix) for markers in one region,
// using chunked sparse matrix multiplication.
// ============================================================
void LDmatRegionInCPP(
    std::string t_genoType,
    std::vector<std::string>& t_genoIndex_prev,
    std::vector<std::string>& t_genoIndex,
    arma::mat& annoIndicatorMat,
    std::string t_outputFile,
    unsigned int t_n,
    bool t_isImputation,
    std::vector<std::string>& annoStringVec,
    std::string regionName)
{
    unsigned int q0 = t_genoIndex.size();
    unsigned int q = q0;
    std::vector<std::string> markerVec(q);
    std::vector<std::string> chrVec(q);
    std::vector<std::string> posVec(q);
    std::vector<std::string> refVec(q);
    std::vector<std::string> altVec(q);
    std::vector<double> altCountsVec(q, arma::datum::nan);
    std::vector<double> imputationInfoVec(q, arma::datum::nan);
    std::vector<double> missingRateVec(q, arma::datum::nan);
    std::vector<double> altFreqVec(q, arma::datum::nan);
    std::vector<std::string> infoVec(q);
    std::vector<double> MACVec(q, arma::datum::nan);
    std::vector<double> MAFVec(q, arma::datum::nan);
    std::vector<uint32_t> N_Vec(q, 0);

    unsigned int m1 = g_region_maxMarkers_cutoff_LDmat;
    std::vector<unsigned int> mPassCVVec;
    arma::vec GVec(t_n);
    arma::vec GZeroVec(t_n);
    std::vector<uint> indexZeroVec;
    std::vector<uint> indexNonZeroVec;
    unsigned int nchunks = 0;
    unsigned int ichunk = 0;
    unsigned int i1InChunk = 0;
    unsigned int i1 = 0;

    unsigned int jm;
    std::vector<uint32_t> location_m_P1Mat;
    std::vector<uint32_t> location_n_P1Mat;
    std::vector<double> value_P1Mat;

    for (unsigned int i = 0; i < q0; i++)
    {
        // marker-level information
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
            } else if (t_genoType == "plink") {
                t_genoIndex_prev_str = t_genoIndex.at(i - 1);
            }
            gIndex_prev = std::strtoull(t_genoIndex_prev_str.c_str(), &end_prev, 10);
        }

        bool isReadMarker = Unified_getOneMarker(
            t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr,
            altFreq, altCounts, missingRate, imputeInfo,
            isOutputIndexForMissing,
            indexForMissing,
            isOnlyOutputNonZero,
            indexNonZeroVec,
            GVec,
            t_isImputation);

        if (!isReadMarker) {
            std::cout << "ERROR: Reading " << i << "th marker failed." << std::endl;
            break;
        }
        std::string pds = std::to_string(pd);
        std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;

        double MAF = std::min(altFreq, 1 - altFreq);
        double w0;
        double MAC = MAF * 2 * t_n * (1 - missingRate);
        flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing,
                                  g_impute_method_LDmat, g_dosage_zerod_cutoff_LDmat,
                                  g_dosage_zerod_MAC_cutoff_LDmat, MAC,
                                  indexZeroVec, indexNonZeroVec);
        arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
        MAF = std::min(altFreq, 1 - altFreq);
        MAC = std::min(altCounts, t_n * 2.0 - altCounts);
        chrVec.at(i) = chr;
        posVec.at(i) = pds;
        if (flip) {
            refVec.at(i) = alt;
            altVec.at(i) = ref;
        } else {
            refVec.at(i) = ref;
            altVec.at(i) = alt;
        }

        markerVec.at(i) = marker;
        altFreqVec.at(i) = altFreq;
        missingRateVec.at(i) = missingRate;
        MAFVec.at(i) = MAF;
        N_Vec.at(i) = t_n;
        imputationInfoVec.at(i) = imputeInfo;
        infoVec.at(i) = info;

        if ((missingRate > g_missingRate_cutoff_LDmat) ||
            (MAF > g_maxMAFLimit_LDmat) ||
            (MAF < g_marker_minMAF_cutoff_LDmat) ||
            (MAC < g_marker_minMAC_cutoff_LDmat) ||
            (imputeInfo < g_marker_minINFO_cutoff_LDmat)) {
            continue;
        } else {
            MACVec.at(i) = MAC;
            altCountsVec.at(i) = altCounts;
            indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
            location_m_P1Mat.insert(std::end(location_m_P1Mat), indexNonZeroVec.size(), i1InChunk);
            location_n_P1Mat.insert(std::end(location_n_P1Mat), std::begin(indexNonZeroVec), std::end(indexNonZeroVec));
            indexNonZeroVec.clear();
            arma::vec GVecnonZero = GVec(indexNonZeroVec_arma);
            typedef std::vector<double> stdvec;
            stdvec GVecnonZero_stdvec = arma::conv_to<stdvec>::from(GVecnonZero);
            value_P1Mat.insert(std::end(value_P1Mat), std::begin(GVecnonZero_stdvec), std::end(GVecnonZero_stdvec));
            GVecnonZero_stdvec.clear();
            i1 += 1;
            i1InChunk += 1;
        }

        if (i1InChunk == m1) {
            std::cout << "In chunks 0-" << ichunk << ", " << i1 << " markers are included." << std::endl;
            int nonzeroSize = location_m_P1Mat.size();
            arma::uvec location_m_P1Mat_arma = arma::conv_to<arma::uvec>::from(location_m_P1Mat);
            arma::uvec location_n_P1Mat_arma = arma::conv_to<arma::uvec>::from(location_n_P1Mat);

            arma::umat location_P1Mat_arma(2, location_m_P1Mat_arma.n_elem);
            location_P1Mat_arma.row(0) = location_m_P1Mat_arma.t();
            location_P1Mat_arma.row(1) = location_n_P1Mat_arma.t();
            arma::vec value_P1Mat_arma_d = arma::conv_to<arma::vec>::from(value_P1Mat);
            arma::ivec value_P1Mat_arma = arma::conv_to<arma::ivec>::from(value_P1Mat);
            arma::sp_imat P1Mat(location_P1Mat_arma, value_P1Mat_arma, i1InChunk, t_n);

            P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
            value_P1Mat.clear();
            location_m_P1Mat.clear();
            location_n_P1Mat.clear();

            mPassCVVec.push_back(m1);
            ichunk += 1;
            i1InChunk = 0;
            nchunks = nchunks + 1;
        }
        // Rcpp::checkUserInterrupt() removed -- no R event loop
    } // for(unsigned int i = 0; i < q0; i++)


    // the last chunk
    if (i1InChunk != 0) {
        std::cout << "In chunks 0-" << ichunk << ", " << i1 << " markers are included." << std::endl;
        int nonzeroSize = location_m_P1Mat.size();
        arma::uvec location_m_P1Mat_arma = arma::conv_to<arma::uvec>::from(location_m_P1Mat);
        arma::uvec location_n_P1Mat_arma = arma::conv_to<arma::uvec>::from(location_n_P1Mat);
        arma::umat location_P1Mat_arma(2, location_m_P1Mat_arma.n_elem);
        location_P1Mat_arma.row(0) = location_m_P1Mat_arma.t();
        location_P1Mat_arma.row(1) = location_n_P1Mat_arma.t();

        arma::ivec value_P1Mat_arma = arma::conv_to<arma::ivec>::from(value_P1Mat);
        arma::sp_imat P1Mat(location_P1Mat_arma, value_P1Mat_arma, i1InChunk, t_n);
        P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        value_P1Mat.clear();
        location_m_P1Mat.clear();
        location_n_P1Mat.clear();

        ichunk = ichunk + 1;
        mPassCVVec.push_back(i1InChunk);
        nchunks = nchunks + 1;
        i1InChunk = 0;
    }


    int mPassCVVecsize = mPassCVVec.size();
    nchunks = mPassCVVecsize;

    // Compute LD matrix from chunked P1Mat sparse matrices
    std::vector<unsigned int> rowIndices_VarMat_vec;
    std::vector<unsigned int> colIndices_VarMat_vec;
    unsigned int rowind_new;
    unsigned int colind_new;
    std::vector<int> values_VarMat_vec;

    arma::sp_imat P1Mat;
    arma::sp_imat P2Mat;

    if (nchunks >= 1)
    {
        int first_row = 0, first_col = 0, last_row = 0, last_col = 0;

        for (unsigned int index1 = 0; index1 < nchunks; index1++)
        {
            last_row = first_row + mPassCVVec.at(index1) - 1;

            std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";

            P1Mat.load(P1MatFile);

            if (P1Mat.n_cols == 0) continue;

            // off-diagonal sub-matrix
            for (unsigned int index2 = 0; index2 < index1; index2++)
            {
                std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;

                P2Mat.load(t_outputFile + "_P1Mat_Chunk_" + std::to_string(index2) + ".bin");

                if (P2Mat.n_cols == 0) continue;
                arma::sp_imat offVarMat = P1Mat * (P2Mat.t());

                last_col = first_col + mPassCVVec.at(index2) - 1;

                arma::sp_imat::const_iterator it = offVarMat.begin();
                arma::sp_imat::const_iterator it_end = offVarMat.end();

                for (arma::uword iof = 0; it != it_end; ++it, ++iof) {
                    rowind_new = it.row() + first_row;
                    colind_new = it.col() + first_col;
                    if (rowind_new >= colind_new) {
                        rowIndices_VarMat_vec.push_back(it.row() + first_row);
                        colIndices_VarMat_vec.push_back(it.col() + first_col);
                        values_VarMat_vec.push_back(*it);
                    }
                }

                first_col = last_col + 1;
            }

            // diagonal sub-matrix
            last_col = first_col + mPassCVVec.at(index1) - 1;
            std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
            P2Mat.load(t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin");

            arma::sp_imat diagVarMat = P1Mat * (P2Mat.t());

            {
                arma::sp_imat::const_iterator it = diagVarMat.begin();
                arma::sp_imat::const_iterator it_end = diagVarMat.end();

                for (arma::uword iof = 0; it != it_end; ++it, ++iof) {
                    rowind_new = it.row() + first_row;
                    colind_new = it.col() + first_col;
                    if (rowind_new >= colind_new) {
                        rowIndices_VarMat_vec.push_back(it.row() + first_row);
                        colIndices_VarMat_vec.push_back(it.col() + first_col);
                        values_VarMat_vec.push_back(*it);
                    }
                }
            }

            first_row = last_row + 1;
            first_col = 0;
            // Rcpp::checkUserInterrupt() removed
        }

        // Clean up temporary chunk files
        for (unsigned int index1 = 0; index1 < nchunks; index1++)
        {
            std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
            const char* File1 = P1MatFile.c_str();
            std::remove(File1);
        }
    }

    // Write LD matrix to output file (lower-triangular sparse COO format)
    if (OutFile_LDmat.is_open()) {
        for (unsigned int k = 0; k < rowIndices_VarMat_vec.size(); k++) {
            OutFile_LDmat << rowIndices_VarMat_vec.at(k);
            OutFile_LDmat << " ";
            OutFile_LDmat << colIndices_VarMat_vec.at(k);
            OutFile_LDmat << " ";
            OutFile_LDmat << values_VarMat_vec.at(k);
            OutFile_LDmat << "\n";
        }
    } else {
        std::cout << g_outputFilePrefix_LDmat << " is not opened" << std::endl;
    }

    // Write index file entry for this region
    if (OutFile_index_LDmat.is_open()) {
        g_startpos = g_endpos + 1;
        g_endpos = g_startpos + rowIndices_VarMat_vec.size() - 1;
        OutFile_index_LDmat << g_startpos;
        OutFile_index_LDmat << " ";
        OutFile_index_LDmat << g_endpos;
        OutFile_index_LDmat << " ";
        OutFile_index_LDmat << regionName;
        OutFile_index_LDmat << "\n";
    } else {
        std::cout << g_outputFilePrefixIndex_LDmat << " is not opened" << std::endl;
    }

    // Write marker info for this region
    writeOutfile_single_LDmat(chrVec,
                               posVec,
                               refVec,
                               altVec,
                               infoVec,
                               MACVec,
                               missingRateVec,
                               N_Vec,
                               regionName);
}


// ============================================================
// writeOutfile_single_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 496-532
// ============================================================
void writeOutfile_single_LDmat(
    std::vector<std::string>& chrVec,
    std::vector<std::string>& posVec,
    std::vector<std::string>& refVec,
    std::vector<std::string>& altVec,
    std::vector<std::string>& infoVec,
    std::vector<double>& altCountsVec,
    std::vector<double>& missingRateVec,
    std::vector<uint32_t>& N_Vec,
    std::string regionName)
{
    int numtest = 0;
    for (unsigned int k = 0; k < chrVec.size(); k++) {
        if (!std::isnan(altCountsVec.at(k))) {
            OutFile_single_LDmat << chrVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << posVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << refVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << altVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << altCountsVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << N_Vec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << missingRateVec.at(k);
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << regionName;
            OutFile_single_LDmat << "\t";
            OutFile_single_LDmat << numtest;
            OutFile_single_LDmat << "\n";
            numtest = numtest + 1;
        }
    }
}


// ============================================================
// openOutfile_single_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 535-549
// ============================================================
bool openOutfile_single_LDmat(bool isappend)
{
    bool isopen;
    if (!isappend) {
        OutFile_single_LDmat.open(g_outputFilePrefixSingle_LDmat.c_str());
        isopen = OutFile_single_LDmat.is_open();
        if (isopen) {
            OutFile_single_LDmat << "CHR\tPOS\tMajor_Allele\tMinor_Allele\tMAC\tN\tMissing_rate\tSet\tIndex\n";
        }
    } else {
        OutFile_single_LDmat.open(g_outputFilePrefixSingle_LDmat.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_single_LDmat.is_open();
    }
    return isopen;
}


// ============================================================
// openOutfile_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 551-562
// ============================================================
bool openOutfile_LDmat(bool isappend)
{
    bool isopen;
    if (!isappend) {
        OutFile_LDmat.open(g_outputFilePrefix_LDmat.c_str());
        isopen = OutFile_LDmat.is_open();
    } else {
        OutFile_LDmat.open(g_outputFilePrefix_LDmat.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_LDmat.is_open();
    }
    return isopen;
}


// ============================================================
// openOutfile_index_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 564-575
// ============================================================
bool openOutfile_index_LDmat(bool isappend)
{
    bool isopen;
    if (!isappend) {
        OutFile_index_LDmat.open(g_outputFilePrefixIndex_LDmat.c_str());
        isopen = OutFile_index_LDmat.is_open();
    } else {
        OutFile_index_LDmat.open(g_outputFilePrefixIndex_LDmat.c_str(), std::ofstream::out | std::ofstream::app);
        isopen = OutFile_index_LDmat.is_open();
    }
    return isopen;
}


// ============================================================
// closeOutfile_single_LDmat
// Direct port from SAIGE/src/LDmat.cpp lines 577-580
// ============================================================
void closeOutfile_single_LDmat()
{
    OutFile_single_LDmat.close();
}


// ============================================================
// closeOutfile_LDmat
// ============================================================
void closeOutfile_LDmat()
{
    OutFile_LDmat.close();
}


// ============================================================
// closeOutfile_index_LDmat
// ============================================================
void closeOutfile_index_LDmat()
{
    OutFile_index_LDmat.close();
}
