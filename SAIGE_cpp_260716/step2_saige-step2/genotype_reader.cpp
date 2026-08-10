// Standalone port of SAIGE/src/PLINK.cpp, SAIGE/src/VCF.cpp, and SAIGE/src/BGEN.cpp
// PLINK + VCF + BGEN genotype readers -- all Rcpp dependencies removed
//
// Conversions from original SAIGE:
//   1. #include <RcppArmadillo.h>  -->  #include <armadillo>
//   2. Rcpp::stop(...)             -->  throw std::runtime_error(...)
//   3. boost::split(...)           -->  splitLine() using std::istringstream
//   4. boost::replace_all(...)     -->  manual \r removal
//   5. Rcpp::match(...)            -->  std::unordered_map lookup
//   6. Rcpp::CharacterVector       -->  std::vector<std::string>
//   7. Rcpp::IntegerVector         -->  std::vector<uint32_t>
//   8. savvy library (VCF)         -->  htslib (VCF/BCF/VCF.GZ)
//   9. Rcpp BGEN                   -->  standalone BGEN with zstd/zlib

#include <armadillo>
#include <sqlite3.h>
#include <sys/stat.h>
#include "genotype_reader.hpp"
#include "UTIL.hpp"
// locoChromLabelsMatch(): the single canonical lenient chromosome-name
// comparison ("chr1" == "1" == "01"). Reused here so the VCF/BGEN LOCO filters
// agree exactly with the PLINK/PGEN one in main.cpp.
#include "null_model_loader.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>
#include <cstdio>
#include <cstdint>
#include <cmath>
#include <sys/stat.h>
#include <filesystem>
#include <iomanip>

// htslib headers for VCF reading
#include <htslib/hts.h>
#include <htslib/vcf.h>

// zstd and zlib headers for BGEN decompression
#include <zstd.h>
#include <zlib.h>

// Global checkpoint flag for debugging
bool g_writeCheckpoints = false;
std::string g_checkpointDir = "";

// Global pointer to PLINK object (mirrors SAIGE's ptr_gPLINKobj in Main.cpp)
PLINK::PlinkClass* ptr_gPLINKobj = nullptr;

// Global pointer to VCF object (mirrors SAIGE's ptr_gVCFobj in Main.cpp)
VCF::VcfClass* ptr_gVCFobj = nullptr;

// Global pointer to BGEN object (mirrors SAIGE's ptr_gBGENobj in Main.cpp)
BGEN::BgenClass* ptr_gBGENobj = nullptr;

// Global pointer to PGEN object (mirrors SAIGE's ptr_gPGENobj in Main.cpp)
PGEN::PgenClass* ptr_gPGENobj = nullptr;

// Static counter for checkpoint tracking
static uint64_t s_checkpointMarkerCount = 0;

// ============================================================
// Helper: split a line on whitespace (replaces boost::split)
// ============================================================
static std::vector<std::string> splitLine(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

// ============================================================
// Helper: remove trailing \r from string (replaces boost::replace_all for \r)
// ============================================================
static void removeTrailingCR(std::string& s) {
    if (!s.empty() && s.back() == '\r') {
        s.pop_back();
    }
}

namespace PLINK {

// ============================================================
// Constructor
// ============================================================
PlinkClass::PlinkClass(std::string t_bimFile,
                       std::string t_famFile,
                       std::string t_bedFile,
                       std::string t_AlleleOrder)
{
    setPlinkobj(t_bimFile, t_famFile, t_bedFile);
    m_AlleleOrder = t_AlleleOrder;
}

// ============================================================
// setPlinkobj: read PLINK files and open .bed
// ============================================================
void PlinkClass::setPlinkobj(std::string t_bimFile,
                             std::string t_famFile,
                             std::string t_bedFile)
{
    m_bimFile = t_bimFile;
    m_famFile = t_famFile;
    m_bedFile = t_bedFile;

    readBimFile();
    readFamFile();

    m_fin = fopen(t_bedFile.c_str(), "rb");
    if (!m_fin) {
        throw std::runtime_error("Cannot open PLINK .bed file: " + t_bedFile);
    }

    char magicNumber_2[65536];
    fread(magicNumber_2, 2, 1, m_fin);
    char magicNumber3;
    fread(&magicNumber3, 1, 1, m_fin);
    fseek(m_fin, 3, SEEK_SET);

    if (magicNumber3 != 1) {
        throw std::runtime_error(
            "The third magic number of the plink bed file is not 00000001. "
            "Please use SNP-major plink (plink version >= 1.9) files.");
    }
}

// ============================================================
// readBimFile: parse the .bim file
// ============================================================
void PlinkClass::readBimFile()
{
    std::cout << "Reading bim file...." << std::endl;
    std::ifstream bim(m_bimFile);
    if (!bim.is_open()) {
        throw std::runtime_error("Cannot open PLINK .bim file: " + m_bimFile);
    }

    m_M0 = 0;
    std::string line;

    while (getline(bim, line)) {
        m_M0++;
        std::vector<std::string> line_elements = splitLine(line);

        if (line_elements.size() < 6) {
            throw std::runtime_error(
                "PLINK .bim file has fewer than 6 columns at line " +
                std::to_string(m_M0));
        }

        // Remove trailing \r from the last element (Windows line endings)
        removeTrailingCR(line_elements[line_elements.size() - 1]);

        m_chr.push_back(line_elements[0]);
        m_MarkerInPlink.push_back(line_elements[1]);
        m_gd.push_back(std::stof(line_elements[2]));
        m_pd.push_back(std::stoi(line_elements[3]));

        // Convert alleles to uppercase (matching SAIGE behavior)
        std::transform(line_elements[4].begin(), line_elements[4].end(),
                       line_elements[4].begin(), ::toupper);
        std::transform(line_elements[5].begin(), line_elements[5].end(),
                       line_elements[5].begin(), ::toupper);

        m_alt.push_back(line_elements[4]);  // allele 1, usually minor allele, alt-first
        m_ref.push_back(line_elements[5]);  // allele 2, usually major allele
    }
    m_M = m_M0;

    std::cout << "Number of markers in bim file: " << m_M0 << std::endl;
}

// ============================================================
// readFamFile: parse the .fam file
// ============================================================
void PlinkClass::readFamFile()
{
    std::cout << "Reading fam file...." << std::endl;
    std::ifstream fam(m_famFile);
    if (!fam.is_open()) {
        throw std::runtime_error("Cannot open PLINK .fam file: " + m_famFile);
    }

    m_N0 = 0;
    std::string line;

    while (getline(fam, line)) {
        m_N0++;
        std::vector<std::string> line_elements = splitLine(line);

        if (line_elements.size() < 2) {
            throw std::runtime_error(
                "PLINK .fam file has fewer than 2 columns at line " +
                std::to_string(m_N0));
        }

        m_SampleInPlink.push_back(line_elements[1]);  // put IID to m_SampleInPlink
    }

    m_N = m_N0;
    m_numBytesofEachMarker0 = (m_N0 + 3) / 4;
    m_OneMarkerG4.reserve(m_numBytesofEachMarker0);
    m_OneMarkerG4.resize(m_numBytesofEachMarker0);

    std::cout << "Number of samples in fam file: " << m_N0 << std::endl;
}

// ============================================================
// setPosSampleInPlink: find positions of model samples in PLINK file
//
// Replaces Rcpp::match() with std::unordered_map lookup
// ============================================================
void PlinkClass::setPosSampleInPlink(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in PLINK files...." << std::endl;

    // If no sample IDs provided, use all fam samples in order (identity mapping)
    if (t_SampleInModel.empty()) {
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from fam file." << std::endl;
        m_N = m_N0;
        m_numBytesofEachMarker = (m_N + 3) / 4;
        m_posSampleInPlink.resize(m_N);
        for (uint32_t i = 0; i < m_N; i++) {
            m_posSampleInPlink[i] = i;
        }
        std::cout << "Number of samples in analysis: " << m_N << std::endl;
        return;
    }

    m_N = t_SampleInModel.size();
    m_numBytesofEachMarker = (m_N + 3) / 4;

    // Build map of plink sample IID -> index (0-based)
    std::unordered_map<std::string, uint32_t> plinkSampleMap;
    for (uint32_t i = 0; i < m_SampleInPlink.size(); i++) {
        plinkSampleMap[m_SampleInPlink[i]] = i;
    }

    m_posSampleInPlink.resize(m_N);
    for (uint32_t i = 0; i < m_N; i++) {
        auto it = plinkSampleMap.find(t_SampleInModel[i]);
        if (it == plinkSampleMap.end()) {
            throw std::runtime_error(
                "At least one subject requested is not in Plink file. "
                "Sample not found: " + t_SampleInModel[i]);
        }
        m_posSampleInPlink[i] = it->second;  // already 0-based
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// closegenofile: close the .bed file handle
// ============================================================
void PlinkClass::closegenofile()
{
    if (m_fin) {
        fclose(m_fin);
        m_fin = nullptr;
    }
}

// ============================================================
// getOneMarker_ts: thread-safe variant of getOneMarker.
//
// Each OMP thread keeps its own FILE* + scratch buffer in thread_local
// storage. On first use the thread opens its own fd to the .bed file
// (the kernel maintains an independent file position per fd, so multiple
// threads can read concurrently with no contention).
//
// Always performs absolute SEEK_SET (no sequential CUR seek). This is
// correct under OMP parallel-for where iterations execute out of order.
//
// The thread_local FILE* is leaked at thread exit (no destructor for
// raw FILE*). For SAIGE step 2 this is fine: the process exits shortly
// after the parallel region, and the OS reclaims the fd.
// ============================================================
void PlinkClass::getOneMarker_ts(uint64_t t_gIndex,
                                 std::string& t_ref,
                                 std::string& t_alt,
                                 std::string& t_marker,
                                 uint32_t& t_pd,
                                 std::string& t_chr,
                                 double& t_altFreq,
                                 double& t_altCounts,
                                 double& t_missingRate,
                                 double& t_imputeInfo,
                                 bool t_isOutputIndexForMissing,
                                 std::vector<uint>& t_indexForMissing,
                                 bool t_isOnlyOutputNonZero,
                                 std::vector<uint>& t_indexForNonZero,
                                 bool t_isTrueGenotype,
                                 arma::vec& OneMarkerG1)
{
    if (!t_isTrueGenotype) {
        if (t_isOutputIndexForMissing) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker_ts, if t_isTrueGenotype = FALSE, "
                "then t_isOutputIndexForMissing should be FALSE.");
        }
        if (t_isOnlyOutputNonZero) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker_ts, if t_isTrueGenotype = FALSE, "
                "then t_isOnlyOutputNonZero should be FALSE.");
        }
    }

    // Per-thread file handle and scratch buffer. Opened lazily on first
    // call from each thread.
    thread_local FILE* tlsFin = nullptr;
    thread_local std::vector<unsigned char> tlsOneMarkerG4;

    if (tlsFin == nullptr) {
        tlsFin = fopen(m_bedFile.c_str(), "rb");
        if (!tlsFin) {
            throw std::runtime_error(
                "PlinkClass::getOneMarker_ts: cannot open .bed file: " + m_bedFile);
        }
    }
    if (tlsOneMarkerG4.size() < m_numBytesofEachMarker0) {
        tlsOneMarkerG4.resize(m_numBytesofEachMarker0);
    }

    // Always absolute seek (thread-safe; no dependence on prior file pos)
    uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
    fseek(tlsFin, posSeek, SEEK_SET);

    size_t nRead = fread((char*)(tlsOneMarkerG4.data()), 1,
                         m_numBytesofEachMarker0, tlsFin);
    if (nRead != m_numBytesofEachMarker0) {
        throw std::runtime_error(
            "PlinkClass::getOneMarker_ts: short read on marker " +
            std::to_string(t_gIndex));
    }

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    t_marker = m_MarkerInPlink[t_gIndex];
    t_pd     = m_pd[t_gIndex];
    t_chr    = m_chr[t_gIndex];

    std::vector<int8_t> genoMaps;
    if (m_AlleleOrder == "alt-first") {
        t_ref = m_ref[t_gIndex];
        t_alt = m_alt[t_gIndex];
        genoMaps = m_genoMaps_alt_first;
    }
    if (m_AlleleOrder == "ref-first") {
        t_ref = m_alt[t_gIndex];
        t_alt = m_ref[t_gIndex];
        genoMaps = m_genoMaps_ref_first;
    }

    uint j = 0;
    int counts[] = {0, 0, 0, 0};

    for (uint32_t i = 0; i < m_N; i++) {
        auto ind = m_posSampleInPlink[i];
        unsigned char bufferG4 = tlsOneMarkerG4[ind / 4];
        size_t bufferG1;

        getGenotype(bufferG4, ind % 4, bufferG1);
        counts[bufferG1]++;

        if (bufferG1 == MISSING && t_isOutputIndexForMissing) {
            t_indexForMissing.push_back(i);
        }
        if (t_isTrueGenotype) {
            bufferG1 = genoMaps[bufferG1];
        }
        if (bufferG1 > 0) {
            t_indexForNonZero.push_back(i);
        }
        if (t_isOnlyOutputNonZero) {
            if (bufferG1 > 0) {
                OneMarkerG1[j] = bufferG1;
                j = j + 1;
            }
        } else {
            OneMarkerG1[i] = bufferG1;
        }
    }

    int numMissing = counts[MISSING];
    int count = m_N - numMissing;
    t_missingRate = (double)numMissing / (double)m_N;
    t_imputeInfo  = 1;
    t_altCounts   = (double)(counts[HET] + 2 * counts[HOM_ALT]);

    if (count > 0) {
        t_altFreq = t_altCounts / (double)count / 2;
    } else {
        t_altFreq = 0;
    }

    if (m_AlleleOrder == "ref-first") {
        t_altFreq = 1 - t_altFreq;
        t_altCounts = 2 * (double)count * t_altFreq;
    }

    if (t_isOnlyOutputNonZero) {
        OneMarkerG1.resize(j);
    }
}

// ============================================================
// B2: getOneMarker_carriers_ts — decode directly into the post-
// imputeGenoAndFlip carrier (nonzero) representation, no dense N-vector.
// Two passes over the L2-resident packed buffer. Faithful reproduction of
// getOneMarker_ts + imputeGenoAndFlip for hard-call genotypes.
// ============================================================
bool PlinkClass::getOneMarker_carriers_ts(uint64_t t_gIndex,
                                          double t_dosage_zerod_cutoff,
                                          double t_dosage_zerod_MAC_cutoff,
                                          int t_impute_case,
                                          CarrierMarker& out)
{
    thread_local FILE* tlsFinC = nullptr;
    thread_local std::vector<unsigned char> tlsBufC;
    if (tlsFinC == nullptr) {
        tlsFinC = fopen(m_bedFile.c_str(), "rb");
        if (!tlsFinC) throw std::runtime_error(
            "PlinkClass::getOneMarker_carriers_ts: cannot open .bed: " + m_bedFile);
    }
    if (tlsBufC.size() < m_numBytesofEachMarker0) tlsBufC.resize(m_numBytesofEachMarker0);

    uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
    fseek(tlsFinC, posSeek, SEEK_SET);
    if (fread((char*)tlsBufC.data(), 1, m_numBytesofEachMarker0, tlsFinC)
        != m_numBytesofEachMarker0) return false;

    out.idx.clear(); out.val.clear();
    out.marker = m_MarkerInPlink[t_gIndex];
    out.pd     = m_pd[t_gIndex];
    out.chr    = m_chr[t_gIndex];
    const std::vector<int8_t>* genoMaps;
    if (m_AlleleOrder == "ref-first") { out.ref = m_alt[t_gIndex]; out.alt = m_ref[t_gIndex]; genoMaps = &m_genoMaps_ref_first; }
    else                              { out.ref = m_ref[t_gIndex]; out.alt = m_alt[t_gIndex]; genoMaps = &m_genoMaps_alt_first; }

    // ---- pass 1: counts (over the packed buffer) ----
    int counts[4] = {0, 0, 0, 0};
    const unsigned char* buf = tlsBufC.data();
    for (uint32_t i = 0; i < m_N; i++) {
        const uint32_t ind = m_posSampleInPlink[i];
        size_t code;
        getGenotype(buf[ind >> 2], ind & 3u, code);
        counts[code]++;
    }
    const int numMissing = counts[MISSING];
    const int nonMissing = (int)m_N - numMissing;
    out.nMissing     = (uint32_t)numMissing;
    out.missingRate  = (double)numMissing / (double)m_N;
    out.imputeInfo   = 1.0;

    // reader-style pre-impute altFreq (matches getOneMarker_ts)
    double altCounts_pre = (double)(counts[HET] + 2 * counts[HOM_ALT]);
    double altFreq_pre   = (nonMissing > 0) ? altCounts_pre / (double)nonMissing / 2.0 : 0.0;
    if (m_AlleleOrder == "ref-first") altFreq_pre = 1.0 - altFreq_pre;

    // ---- imputeGenoAndFlip replication ----
    out.flip = (altFreq_pre > 0.5);
    double altFreq_flip = out.flip ? (1.0 - altFreq_pre) : altFreq_pre;
    double imputeG = 0.0;
    if      (t_impute_case == 1) imputeG = std::round(2.0 * altFreq_flip);
    else if (t_impute_case == 2) imputeG = 2.0 * altFreq_flip;
    else                         imputeG = 0.0;   // "none" -> 0 (hard-call path)

    // MAC that gates .clean(): QC1 MAC (min-allele) + imputeG*nMissing
    double MAF_pre  = std::min(altFreq_pre, 1.0 - altFreq_pre);
    double MAC_qc1  = MAF_pre * (double)m_N * (1.0 - out.missingRate) * 2.0;
    double MAC_clean = MAC_qc1 + imputeG * (double)numMissing;
    bool   doClean  = (t_dosage_zerod_cutoff > 0.0) && (MAC_clean <= t_dosage_zerod_MAC_cutoff);

    // ---- pass 2: emit carriers with correct polarity + imputed missing ----
    const bool flip = out.flip;
    for (uint32_t i = 0; i < m_N; i++) {
        const uint32_t ind = m_posSampleInPlink[i];
        size_t code;
        getGenotype(buf[ind >> 2], ind & 3u, code);
        double d;
        if (code == MISSING) {
            d = imputeG;
        } else {
            int d0 = (*genoMaps)[code];       // 0/1/2 pre-flip dosage
            d = flip ? (2.0 - (double)d0) : (double)d0;
        }
        if (d == 0.0) continue;
        if (doClean && std::abs(d) < t_dosage_zerod_cutoff) continue;
        out.idx.push_back(i);
        out.val.push_back(d);
    }

    // ---- post-impute stats (== arma::sum over full GVec; non-carriers are 0) ----
    double sumG = 0.0;
    for (double v : out.val) sumG += v;
    double altFreq_post = sumG / (2.0 * (double)m_N);
    if (flip) { out.altFreq = 1.0 - altFreq_post; out.altCount = 2.0 * (double)m_N - sumG; }
    else      { out.altFreq = altFreq_post;       out.altCount = sumG; }
    out.MAC = std::min(out.altCount, 2.0 * (double)m_N - out.altCount);
    return true;
}

// ============================================================
// getOneMarker: read one marker from .bed and decode genotypes
//
// This is a direct port of SAIGE/src/PLINK.cpp::getOneMarker
// with the EXACT same algorithm:
//   - Seek in .bed file based on previous/current index
//   - Read m_numBytesofEachMarker0 bytes
//   - Decode 2-bit genotypes
//   - Apply allele order mapping
//   - Calculate alt freq, alt counts, missing rate
// ============================================================
void PlinkClass::getOneMarker(uint64_t& t_gIndex_prev,
                              uint64_t& t_gIndex,
                              std::string& t_ref,
                              std::string& t_alt,
                              std::string& t_marker,
                              uint32_t& t_pd,
                              std::string& t_chr,
                              double& t_altFreq,
                              double& t_altCounts,
                              double& t_missingRate,
                              double& t_imputeInfo,
                              bool& t_isOutputIndexForMissing,
                              std::vector<uint>& t_indexForMissing,
                              bool& t_isOnlyOutputNonZero,
                              std::vector<uint>& t_indexForNonZero,
                              bool& t_isTrueGenotype,
                              arma::vec& OneMarkerG1)
{
    // t_isTrueGenotype = FALSE is used only when calculating GRM
    if (!t_isTrueGenotype) {
        if (t_isOutputIndexForMissing) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, "
                "then t_isOutputIndexForMissing should be FALSE.");
        }
        if (t_isOnlyOutputNonZero) {
            throw std::runtime_error(
                "Check PlinkClass::getOneMarker, if t_isTrueGenotype = FALSE, "
                "then t_isOnlyOutputNonZero should be FALSE.");
        }
    }

    // Seek to the correct position in the .bed file.
    // Phase E: when t_gIndex_prev == 0 we ALWAYS force an absolute SEEK_SET
    // (even when t_gIndex == 0), so callers under OpenMP (where reads from
    // different regions interleave at a critical section) can request any
    // marker without depending on the prior file position.
    if (t_gIndex_prev == 0) {
        uint64_t posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex;
        fseek(m_fin, posSeek, SEEK_SET);
    } else if (t_gIndex > 0) {
        uint64_t posSeek = m_numBytesofEachMarker0 * (t_gIndex - t_gIndex_prev - 1);
        if (posSeek > 0) {
            fseek(m_fin, posSeek, SEEK_CUR);
        }
    }

    // Read the raw genotype bytes for this marker
    fread((char*)(&m_OneMarkerG4[0]), 1, m_numBytesofEachMarker0, m_fin);

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    t_marker = m_MarkerInPlink[t_gIndex];
    t_pd = m_pd[t_gIndex];
    t_chr = m_chr[t_gIndex];

    // Select genotype mapping based on allele order
    std::vector<int8_t> genoMaps;

    if (m_AlleleOrder == "alt-first") {
        t_ref = m_ref[t_gIndex];
        t_alt = m_alt[t_gIndex];
        genoMaps = m_genoMaps_alt_first;
    }

    if (m_AlleleOrder == "ref-first") {
        t_ref = m_alt[t_gIndex];
        t_alt = m_ref[t_gIndex];
        genoMaps = m_genoMaps_ref_first;
    }

    // Decode genotypes for each sample in the analysis
    uint j = 0;
    int counts[] = {0, 0, 0, 0};

    for (uint32_t i = 0; i < m_N; i++) {
        auto ind = m_posSampleInPlink[i];             // C++ start from 0
        unsigned char bufferG4 = m_OneMarkerG4[ind / 4];  // 1 byte for 4 genotypes
        size_t bufferG1;                                   // 1 genotype (1 sample)

        // https://www.cog-genomics.org/plink/1.9/formats#bed
        getGenotype(bufferG4, ind % 4, bufferG1);     // bufferG4 -> bufferG1

        counts[bufferG1]++;

        if (bufferG1 == MISSING && t_isOutputIndexForMissing) {
            t_indexForMissing.push_back(i);
        }

        if (t_isTrueGenotype) {
            bufferG1 = genoMaps[bufferG1];
        }

        if (bufferG1 > 0) {
            t_indexForNonZero.push_back(i);
        }

        if (t_isOnlyOutputNonZero) {
            if (bufferG1 > 0) {
                OneMarkerG1[j] = bufferG1;
                j = j + 1;
            }
        } else {
            OneMarkerG1[i] = bufferG1;
        }
    }

    int numMissing = counts[MISSING];

    int count = m_N - numMissing;
    t_missingRate = (double)numMissing / (double)m_N;
    t_imputeInfo = 1;
    t_altCounts = (double)(counts[HET] + 2 * counts[HOM_ALT]);

    if (count > 0) {
        t_altFreq = t_altCounts / (double)count / 2;
    } else {
        t_altFreq = 0;
    }

    // updated on 03/14/2021
    if (m_AlleleOrder == "ref-first") {
        t_altFreq = 1 - t_altFreq;
        t_altCounts = 2 * (double)count * t_altFreq;
    }

    if (t_isOnlyOutputNonZero) {
        OneMarkerG1.resize(j);
    }

    // === CHECKPOINT OUTPUT BEGIN ===
    if (g_writeCheckpoints && s_checkpointMarkerCount == 0 && !g_checkpointDir.empty()) {
        // Create checkpoint directory if it does not exist
        std::filesystem::create_directories(g_checkpointDir);

        // Write marker info checkpoint
        std::string infoFile = g_checkpointDir + "/ckpt_07_marker0_raw.txt";
        std::ofstream ofs_info(infoFile);
        if (ofs_info.is_open()) {
            ofs_info << std::setprecision(15);
            ofs_info << "field\tvalue" << std::endl;
            ofs_info << "chr\t" << t_chr << std::endl;
            ofs_info << "pos\t" << t_pd << std::endl;
            ofs_info << "ref\t" << t_ref << std::endl;
            ofs_info << "alt\t" << t_alt << std::endl;
            ofs_info << "marker\t" << t_marker << std::endl;
            ofs_info << "altFreq\t" << t_altFreq << std::endl;
            ofs_info << "altCounts\t" << t_altCounts << std::endl;
            ofs_info << "missingRate\t" << t_missingRate << std::endl;
            ofs_info << "imputeInfo\t" << t_imputeInfo << std::endl;
            ofs_info.close();
            std::cout << "[CHECKPOINT] Wrote ckpt_07_marker0_raw.txt to " << g_checkpointDir << std::endl;
        }

        // Write first 20 genotype values checkpoint
        std::string genoFile = g_checkpointDir + "/ckpt_08_marker0_geno.txt";
        std::ofstream ofs_geno(genoFile);
        if (ofs_geno.is_open()) {
            ofs_geno << std::setprecision(15);
            ofs_geno << "index\tgenotype" << std::endl;
            uint32_t nToWrite = std::min((uint32_t)20, (uint32_t)OneMarkerG1.n_elem);
            for (uint32_t k = 0; k < nToWrite; k++) {
                ofs_geno << (k + 1) << "\t" << OneMarkerG1[k] << std::endl;
            }
            ofs_geno.close();
            std::cout << "[CHECKPOINT] Wrote ckpt_08_marker0_geno.txt to " << g_checkpointDir << std::endl;
        }
    }
    s_checkpointMarkerCount++;
    // === CHECKPOINT OUTPUT END ===
}

} // namespace PLINK


// ============================================================
// VCF namespace: VCF/BCF/VCF.GZ reader using htslib
// Ported from SAIGE/src/VCF.cpp with savvy -> htslib conversion
// ============================================================
namespace VCF {

// ============================================================
// Constructor
// ============================================================
VcfClass::VcfClass(const std::string& t_vcfFileName,
                   const std::string& t_vcfField,
                   std::vector<std::string>& t_SampleInModel)
    : m_htsFile(nullptr), m_hdr(nullptr), m_rec(nullptr),
      m_vcfFileName(t_vcfFileName), m_fmtField(t_vcfField),
      m_N0(0), m_N(0), m_M0(0), m_totalMarkers(0), m_isPreScanned(false),
      m_metaPrepopulated(false)
{
    // Open VCF/BCF file
    m_htsFile = hts_open(t_vcfFileName.c_str(), "r");
    if (!m_htsFile) {
        throw std::runtime_error("Cannot open VCF file: " + t_vcfFileName);
    }

    // Read header
    m_hdr = bcf_hdr_read(m_htsFile);
    if (!m_hdr) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
        throw std::runtime_error("Cannot read VCF header from: " + t_vcfFileName);
    }

    // Allocate record
    m_rec = bcf_init();

    // Get sample IDs
    getSampleIDlist();
    std::cout << "Open VCF done" << std::endl;
    std::cout << "To read the field " << t_vcfField << std::endl;
    std::cout << "Number of samples in the vcf file: " << m_N0 << std::endl;

    // Verify the format field exists
    // Check if the requested format field is in the header
    int fmt_id = bcf_hdr_id2int(m_hdr, BCF_DT_ID, t_vcfField.c_str());
    if (fmt_id < 0 || !bcf_hdr_idinfo_exists(m_hdr, BCF_HL_FMT, fmt_id)) {
        // Try fallback to HDS if DS or GT was requested
        if ((t_vcfField == "DS" || t_vcfField == "GT")) {
            int hds_id = bcf_hdr_id2int(m_hdr, BCF_DT_ID, "HDS");
            if (hds_id >= 0 && bcf_hdr_idinfo_exists(m_hdr, BCF_HL_FMT, hds_id)) {
                m_fmtField = "HDS";
                std::cout << "Format field " << t_vcfField << " not found, using HDS instead" << std::endl;
            } else {
                std::cerr << "WARNING: vcfField (" << t_vcfField << ") not present in genotype file." << std::endl;
            }
        } else {
            std::cerr << "WARNING: vcfField (" << t_vcfField << ") not present in genotype file." << std::endl;
        }
    }

    // Set up sample mapping
    setPosSampleInVcf(t_SampleInModel);
}

// ============================================================
// Destructor
// ============================================================
VcfClass::~VcfClass()
{
    closegenofile();
}

// ============================================================
// getSampleIDlist: extract sample IDs from VCF header
// ============================================================
void VcfClass::getSampleIDlist()
{
    m_N0 = bcf_hdr_nsamples(m_hdr);
    m_SampleInVcf.resize(m_N0);
    for (uint32_t i = 0; i < m_N0; i++) {
        m_SampleInVcf[i] = m_hdr->samples[i];
    }
}

// ============================================================
// setPosSampleInVcf: find positions of model samples in VCF file
// Mirrors SAIGE VCF.cpp::setPosSampleInVcf but replaces Rcpp::match
// with std::unordered_map lookup
// ============================================================
void VcfClass::setPosSampleInVcf(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in VCF files...." << std::endl;

    // If no sample IDs provided, use all VCF samples in order
    if (t_SampleInModel.empty()) {
        m_N = m_N0;
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from VCF file." << std::endl;

        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            m_posSampleInModel[i] = i;  // identity mapping
        }
    } else {
        m_N = t_SampleInModel.size();
        std::cout << "m_N " << m_N << std::endl;

        // Build map of model sample ID -> model index (0-based)
        std::unordered_map<std::string, int32_t> modelSampleMap;
        for (uint32_t i = 0; i < m_N; i++) {
            modelSampleMap[t_SampleInModel[i]] = (int32_t)i;
        }

        // Verify all model samples exist in VCF
        for (uint32_t i = 0; i < m_N; i++) {
            bool found = false;
            for (uint32_t j = 0; j < m_N0; j++) {
                if (m_SampleInVcf[j] == t_SampleInModel[i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "At least one subject requested is not in VCF file. "
                    "Sample not found: " + t_SampleInModel[i]);
            }
        }

        // Build mapping: for each VCF sample, what is its position in the model?
        // m_posSampleInModel[vcf_idx] = model_idx, or -1 if not in model
        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            auto it = modelSampleMap.find(m_SampleInVcf[i]);
            if (it != modelSampleMap.end()) {
                m_posSampleInModel[i] = it->second;
            } else {
                m_posSampleInModel[i] = -1;
            }
        }
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// getOneMarker: read one marker from VCF and extract dosages
//
// This mirrors SAIGE VCF.cpp::getOneMarker but uses htslib
// instead of savvy. Reads sequentially (no random access).
//
// Supports GT (hard calls -> 0/1/2) and DS (dosage) fields.
// ============================================================
bool VcfClass::getOneMarker(
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    arma::vec& dosages,
    bool t_isImputation)
{
    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    // Read next record. Under a LOCO chromosome restriction, skip forward over
    // every record on another chromosome (the reader is strictly sequential;
    // see VcfClass::m_chromRestrict).
    while (true) {
        int ret = bcf_read(m_htsFile, m_hdr, m_rec);
        if (ret < 0) {
            // End of file or error
            std::cout << "Reach the end of the vcf file" << std::endl;
            return false;
        }
        if (m_chromRestrict.empty()) break;
        const char* recChr = bcf_hdr_id2name(m_hdr, m_rec->rid);
        if (locoChromLabelsMatch(recChr ? recChr : "", m_chromRestrict)) break;
    }

    // Unpack the record (we need INFO, FORMAT, and shared fields)
    bcf_unpack(m_rec, BCF_UN_ALL);

    // Extract chromosome
    t_chr = bcf_hdr_id2name(m_hdr, m_rec->rid);

    // Extract position (VCF is 1-based, htslib stores 0-based)
    t_pd = (uint32_t)(m_rec->pos + 1);

    // Extract REF
    t_ref = m_rec->d.allele[0];

    // Extract ALT (first ALT allele only, skip multiallelic)
    if (m_rec->n_allele < 2) {
        t_alt = ".";
    } else {
        if (m_rec->n_allele > 2) {
            std::cerr << "Warning: skipping multiallelic variant at "
                      << t_chr << ":" << t_pd << std::endl;
        }
        t_alt = m_rec->d.allele[1];
    }

    // Extract marker ID
    if (m_rec->d.id && std::string(m_rec->d.id) != ".") {
        t_marker = m_rec->d.id;
    } else {
        t_marker = t_chr + ":" + std::to_string(t_pd) + ":" + t_ref + ":" + t_alt;
    }

    // Store marker info for later lookup.
    //
    // When prescanMarkerCount() has already recorded every record (the normal
    // case), appending here would duplicate the metadata and break the
    // index space getMarkerIDToIndex() hands out. Instead, assert that the
    // record we just streamed IS the one recorded at this position: that is
    // the run-time proof that the prescan index space and the streamed index
    // space are the same. A mismatch is fatal -- a silently shifted index
    // would attribute a variant's genotypes to a different variant.
    if (m_metaPrepopulated) {
        if (m_M0 < m_chr.size()) {
            if (m_chr[m_M0] != t_chr || m_pd[m_M0] != t_pd ||
                m_ref[m_M0] != t_ref || m_alt[m_M0] != t_alt) {
                throw std::runtime_error(
                    "VCF prescan/stream index mismatch at marker index " +
                    std::to_string(m_M0) + ": prescan recorded " +
                    m_chr[m_M0] + ":" + std::to_string(m_pd[m_M0]) + ":" +
                    m_ref[m_M0] + ":" + m_alt[m_M0] + " but the stream returned " +
                    t_chr + ":" + std::to_string(t_pd) + ":" + t_ref + ":" + t_alt +
                    ". Marker-ID lookups would point at the wrong variant.");
            }
        } else {
            throw std::runtime_error(
                "VCF stream returned more markers (" + std::to_string(m_M0 + 1) +
                ") than the prescan recorded (" + std::to_string(m_chr.size()) +
                "). The file changed underneath us; refusing to continue.");
        }
    } else {
        m_chr.push_back(t_chr);
        m_pd.push_back(t_pd);
        m_ref.push_back(t_ref);
        m_alt.push_back(t_alt);
        m_MarkerInVcf.push_back(t_marker);
    }
    m_M0++;

    // Initialize dosage vector
    dosages.set_size(m_N);
    dosages.fill(0.0);

    t_altCounts = 0;
    int missing_cnt = 0;

    // Extract imputation info (R2) if requested
    if (t_isImputation) {
        float* info_val = nullptr;
        int n_info = 0;
        int info_ret = bcf_get_info_float(m_hdr, m_rec, "R2", &info_val, &n_info);
        if (info_ret > 0 && n_info > 0) {
            t_imputeInfo = (double)info_val[0];
        } else {
            t_imputeInfo = 1.0;
        }
        if (info_val) free(info_val);
    } else {
        t_imputeInfo = 1.0;
    }

    // Read genotype data based on format field
    if (m_fmtField == "GT") {
        // Read GT field (hard-call genotypes)
        int32_t* gt_arr = nullptr;
        int n_gt = 0;
        int ngt_ret = bcf_get_genotypes(m_hdr, m_rec, &gt_arr, &n_gt);

        if (ngt_ret <= 0) {
            // No GT data available
            dosages.fill(-1.0);
            missing_cnt = m_N;
        } else {
            int ploidy = n_gt / m_N0;

            for (uint32_t i = 0; i < m_N0; i++) {
                int32_t model_idx = m_posSampleInModel[i];
                if (model_idx < 0) continue;  // sample not in model

                int dose = 0;
                bool is_missing = false;

                for (int p = 0; p < ploidy; p++) {
                    int32_t allele = gt_arr[i * ploidy + p];

                    if (allele == bcf_int32_vector_end) {
                        break;  // fewer ploidy for this sample
                    }
                    if (bcf_gt_is_missing(allele)) {
                        is_missing = true;
                        break;
                    }
                    // bcf_gt_allele extracts the allele index (0 = REF, 1 = ALT)
                    dose += bcf_gt_allele(allele);
                }

                if (is_missing) {
                    dosages[model_idx] = -1.0;
                    missing_cnt++;
                    if (t_isOutputIndexForMissing) {
                        t_indexForMissing.push_back((uint)model_idx);
                    }
                } else {
                    dosages[model_idx] = (double)dose;
                    t_altCounts += dose;
                    if (dose > 0) {
                        t_indexForNonZero.push_back((uint)model_idx);
                    }
                }
            }
        }

        if (gt_arr) free(gt_arr);

    } else {
        // Read DS or HDS field (dosage)
        float* ds_arr = nullptr;
        int n_ds = 0;
        int nds_ret = bcf_get_format_float(m_hdr, m_rec, m_fmtField.c_str(), &ds_arr, &n_ds);

        if (nds_ret <= 0) {
            // No dosage data available
            dosages.fill(-1.0);
            missing_cnt = m_N;
        } else {
            int stride = n_ds / m_N0;

            for (uint32_t i = 0; i < m_N0; i++) {
                int32_t model_idx = m_posSampleInModel[i];
                if (model_idx < 0) continue;  // sample not in model

                float dose_val;
                if (stride == 1) {
                    // DS field: single dosage value per sample
                    dose_val = ds_arr[i];
                } else {
                    // HDS field: sum haplotype dosages (like savvy stride_reduce)
                    dose_val = 0.0f;
                    for (int p = 0; p < stride; p++) {
                        float v = ds_arr[i * stride + p];
                        if (!bcf_float_is_missing(v) && !bcf_float_is_vector_end(v)) {
                            dose_val += v;
                        }
                    }
                }

                if (bcf_float_is_missing(dose_val) || std::isnan(dose_val)) {
                    dosages[model_idx] = -1.0;
                    missing_cnt++;
                    if (t_isOutputIndexForMissing) {
                        t_indexForMissing.push_back((uint)model_idx);
                    }
                } else {
                    dosages[model_idx] = (double)dose_val;
                    t_altCounts += dose_val;
                    if (dose_val > 0) {
                        t_indexForNonZero.push_back((uint)model_idx);
                    }
                }
            }
        }

        if (ds_arr) free(ds_arr);
    }

    // Calculate alt frequency and missing rate (same logic as SAIGE VCF.cpp)
    if (missing_cnt > 0) {
        if (missing_cnt == (int)m_N) {
            t_altFreq = 0;
        } else {
            t_altFreq = t_altCounts / 2.0 / (double)(m_N - missing_cnt);
        }
        t_missingRate = (double)missing_cnt / (double)m_N;
    } else {
        t_altFreq = t_altCounts / 2.0 / (double)m_N;
        t_missingRate = 0;
    }

    // Convert -1 (missing) values in dosages to the proper SAIGE convention
    // SAIGE VCF.cpp stores -1 for missing, but the downstream pipeline uses
    // the same missing convention as PLINK: impute or exclude.
    // We leave -1 in the dosage vector; the imputation step will handle it.
    // Actually, let's match PLINK convention: missing genotypes get value -1
    // (which the downstream imputeGenoAndFlip handles via indexForMissing).

    return true;
}

// ============================================================
// prescanMarkerCount: count total markers in VCF by scanning, AND record
// each retained record's chr/pos/ref/alt/ID.
//
// Why the metadata is captured here: VcfClass::getMarkerIDToIndex() is built
// from m_chr/m_pd/m_ref/m_alt, but those used to be filled only by
// getOneMarker() as markers stream past. SAIGE-GENE+ region tests need the
// map BEFORE any streaming (to resolve the group file's variant IDs), so
// every VCF region run found 0 markers. Same chicken-and-egg BGEN had, which
// is solved by populateFromBgi(); this prescan is the VCF analogue and it is
// already walking every record, so it now keeps what it reads.
//
// INDEX-SPACE INVARIANT: the k-th entry of m_chr/m_pd/m_ref/m_alt must be the
// k-th record getOneMarker() returns. Guaranteed by (a) walking the file in
// the same sequential order with no tabix/seek, (b) applying the SAME
// m_chromRestrict + locoChromLabelsMatch skip predicate as getOneMarker(),
// and (c) deriving chr/pos/ref/alt/ID with the same expressions. It is also
// CHECKED at run time: getOneMarker() compares each record it reads against
// the prescanned entry at that index and throws on any mismatch.
// ============================================================
uint32_t VcfClass::prescanMarkerCount()
{
    if (m_isPreScanned) {
        return m_totalMarkers;
    }

    // Save current position and re-open file for counting
    htsFile* countFile = hts_open(m_vcfFileName.c_str(), "r");
    if (!countFile) {
        throw std::runtime_error("Cannot re-open VCF file for marker counting: " + m_vcfFileName);
    }
    bcf_hdr_t* countHdr = bcf_hdr_read(countFile);
    bcf1_t* countRec = bcf_init();

    m_chr.clear(); m_pd.clear(); m_ref.clear(); m_alt.clear(); m_MarkerInVcf.clear();
    m_offChromIDToChrom.clear();

    m_totalMarkers = 0;
    uint32_t nSeen = 0;
    while (bcf_read(countFile, countHdr, countRec) >= 0) {
        nSeen++;
        // BCF_UN_STR unpacks ID + REF/ALT only (not INFO/FORMAT), which is the
        // cheapest level that exposes rec->d.allele / rec->d.id.
        bcf_unpack(countRec, BCF_UN_STR);

        const char* recChrC = bcf_hdr_id2name(countHdr, countRec->rid);
        std::string chr = recChrC ? recChrC : "";
        uint32_t pd = (uint32_t)(countRec->pos + 1);
        std::string ref = countRec->n_allele > 0 ? countRec->d.allele[0] : ".";
        std::string alt = countRec->n_allele < 2 ? std::string(".")
                                                 : std::string(countRec->d.allele[1]);

        if (!m_chromRestrict.empty() && !locoChromLabelsMatch(chr, m_chromRestrict)) {
            // Filtered out of the index space, but remember its chromosome so
            // the region-level LOCO check can still see that a group-file
            // variant lives on another chromosome (rather than mistaking it
            // for "absent from the genotype file" and silently trimming the
            // gene's burden statistic).
            std::string base = chr + ":" + std::to_string(pd) + ":";
            m_offChromIDToChrom.emplace(base + ref + ":" + alt, chr);
            m_offChromIDToChrom.emplace(base + alt + ":" + ref, chr);
            continue;
        }

        std::string marker;
        if (countRec->d.id && std::string(countRec->d.id) != ".") {
            marker = countRec->d.id;
        } else {
            marker = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;
        }

        m_chr.push_back(chr);
        m_pd.push_back(pd);
        m_ref.push_back(ref);
        m_alt.push_back(alt);
        m_MarkerInVcf.push_back(marker);
        m_totalMarkers++;
    }

    bcf_destroy(countRec);
    bcf_hdr_destroy(countHdr);
    hts_close(countFile);

    m_isPreScanned = true;
    m_metaPrepopulated = true;
    if (!m_chromRestrict.empty()) {
        std::cout << "Number of markers in VCF file on chromosome "
                  << m_chromRestrict << ": " << m_totalMarkers
                  << " (of " << nSeen << ")" << std::endl;
    } else {
        std::cout << "Number of markers in VCF file: " << m_totalMarkers << std::endl;
    }
    return m_totalMarkers;
}

// ============================================================
// setChromRestriction: LOCO chromosome filter (see header)
// ============================================================
void VcfClass::setChromRestriction(const std::string& t_chrom)
{
    m_chromRestrict = t_chrom;
    // A restriction change invalidates any cached marker count AND any
    // prescanned metadata: the retained set (hence the index space) changes.
    m_isPreScanned = false;
    m_totalMarkers = 0;
    m_metaPrepopulated = false;
    m_chr.clear(); m_pd.clear(); m_ref.clear(); m_alt.clear(); m_MarkerInVcf.clear();
    m_offChromIDToChrom.clear();
    m_M0 = 0;
}

// ============================================================
// resetFile: close and re-open the VCF file to restart reading
// ============================================================
void VcfClass::resetFile()
{
    // Close current file
    if (m_rec) {
        bcf_destroy(m_rec);
        m_rec = nullptr;
    }
    if (m_hdr) {
        bcf_hdr_destroy(m_hdr);
        m_hdr = nullptr;
    }
    if (m_htsFile) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
    }

    // Clear accumulated marker info -- but NOT when prescanMarkerCount() has
    // already populated it for the whole (restricted) file. That metadata is
    // what getMarkerIDToIndex() serves to the region path, and resetFile() is
    // called right after the prescan; wiping it would restore the original
    // "0 markers in every VCF region test" bug.
    if (!m_metaPrepopulated) {
        m_chr.clear();
        m_pd.clear();
        m_ref.clear();
        m_alt.clear();
        m_MarkerInVcf.clear();
    }
    m_M0 = 0;

    // Re-open
    m_htsFile = hts_open(m_vcfFileName.c_str(), "r");
    if (!m_htsFile) {
        throw std::runtime_error("Cannot re-open VCF file: " + m_vcfFileName);
    }
    m_hdr = bcf_hdr_read(m_htsFile);
    if (!m_hdr) {
        throw std::runtime_error("Cannot re-read VCF header from: " + m_vcfFileName);
    }
    m_rec = bcf_init();
}

// ============================================================
// getMarkerIDToIndex: build chr:pos:ref:alt -> index map
// Requires markers to have been read (or pre-scanned)
// ============================================================
std::unordered_map<std::string, uint32_t> VcfClass::getMarkerIDToIndex()
{
    std::unordered_map<std::string, uint32_t> idMap;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        // Primary key: chr:pos:ref:alt
        std::string id1 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_ref[i] + ":" + m_alt[i];
        idMap[id1] = i;
        // Secondary key: chr:pos:alt:ref (reversed allele order)
        std::string id2 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_alt[i] + ":" + m_ref[i];
        if (idMap.find(id2) == idMap.end()) {
            idMap[id2] = i;
        }
    }
    return idMap;
}

// ============================================================
// getMarkerNameToIndex: build marker ID -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> VcfClass::getMarkerNameToIndex()
{
    std::unordered_map<std::string, uint32_t> nameMap;
    for (uint32_t i = 0; i < m_MarkerInVcf.size(); i++) {
        nameMap[m_MarkerInVcf[i]] = i;
    }
    return nameMap;
}

// ============================================================
// getMarkerIDToChrom: variant ID -> CHROM, same keys as getMarkerIDToIndex().
// Also includes the records the LOCO restriction filtered out of the index
// space, so the region-level LOCO check can distinguish "this variant is on
// another chromosome" from "this variant is not in the genotype file".
// ============================================================
std::unordered_map<std::string, std::string> VcfClass::getMarkerIDToChrom()
{
    std::unordered_map<std::string, std::string> chrMap = m_offChromIDToChrom;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        std::string base = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":";
        chrMap[base + m_ref[i] + ":" + m_alt[i]] = m_chr[i];
        chrMap.emplace(base + m_alt[i] + ":" + m_ref[i], m_chr[i]);
    }
    return chrMap;
}

// ============================================================
// closegenofile: clean up htslib resources
// ============================================================
void VcfClass::closegenofile()
{
    if (m_rec) {
        bcf_destroy(m_rec);
        m_rec = nullptr;
    }
    if (m_hdr) {
        bcf_hdr_destroy(m_hdr);
        m_hdr = nullptr;
    }
    if (m_htsFile) {
        hts_close(m_htsFile);
        m_htsFile = nullptr;
    }
}

} // namespace VCF


// ============================================================
// BGEN namespace: BGEN v1.2 reader (manual binary parsing)
// Ported from SAIGE/src/BGEN.cpp with all Rcpp dependencies removed
// ============================================================
namespace BGEN {

// ============================================================
// Constructor
// ============================================================
BgenClass::BgenClass(const std::string& t_bgenFileName,
                     const std::vector<std::string>& t_SampleInBgen,
                     std::vector<std::string>& t_SampleInModel,
                     const std::string& t_AlleleOrder)
    : m_fin(nullptr), m_zBufLens(0), m_bufLens(0),
      CompressedSNPBlocks(0), m_N0(0), m_N(0), m_M0(0)
{
    setBgenObj(t_bgenFileName, t_SampleInBgen);
    setPosSampleInBgen(t_SampleInModel);
    m_AlleleOrder = t_AlleleOrder;
    allele0.resize(65536);
    allele1.resize(65536);
}

// ============================================================
// Destructor
// ============================================================
BgenClass::~BgenClass()
{
    closegenofile();
}

// ============================================================
// setBgenObj: read BGEN header, validate format
// Ported from SAIGE/src/BGEN.cpp::setBgenObj
// ============================================================
void BgenClass::setBgenObj(const std::string& t_bgenFileName,
                            const std::vector<std::string>& t_SampleInBgen)
{
    m_bgenFileName = t_bgenFileName;
    m_SampleInBgen = t_SampleInBgen;

    m_fin = fopen(t_bgenFileName.c_str(), "rb");
    if (!m_fin) {
        throw std::runtime_error("Cannot open BGEN file: " + t_bgenFileName);
    }

    // Read BGEN v1.2 header
    uint32_t offset;
    fread(&offset, 4, 1, m_fin);

    uint32_t L_H;
    fread(&L_H, 4, 1, m_fin);

    fread(&m_M0, 4, 1, m_fin);
    std::cout << "snpBlocks (Mbgen): " << m_M0 << std::endl;
    if (m_M0 == 0) {
        throw std::runtime_error("No snpBlocks in BGEN file.");
    }

    fread(&m_N0, 4, 1, m_fin);
    std::cout << "samples (Nbgen): " << m_N0 << std::endl;

    unsigned int m_Nsample = t_SampleInBgen.size();
    if (m_N0 != m_Nsample) {
        throw std::runtime_error(
            "Number of samples in BGEN header (" + std::to_string(m_N0) +
            ") does not match sample file (" + std::to_string(m_Nsample) + ")");
    }

    char magic[5];
    fread(magic, 1, 4, m_fin);
    magic[4] = '\0';

    // Skip free data area
    fseek(m_fin, L_H - 20, SEEK_CUR);

    uint32_t flags;
    fread(&flags, 4, 1, m_fin);
    CompressedSNPBlocks = flags & 3;
    std::cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;

    if (CompressedSNPBlocks == COMPRESSION_ZLIB) {
        std::cout << "Warning: your bgen file uses zlib compression, which is slow." << std::endl;
    } else if (CompressedSNPBlocks != COMPRESSION_ZSTD) {
        throw std::runtime_error("Bgenreader only supports zlib or zstd compression.");
    }

    uint32_t Layout = (flags >> 2) & 0xf;
    std::cout << "Layout: " << Layout << std::endl;
    if (Layout != 1 && Layout != 2) {
        throw std::runtime_error("Bgenreader only supports Layout 1 or 2.");
    }

    // Seek past the header to the data area
    fseek(m_fin, offset + 4, SEEK_SET);
}

// ============================================================
// setPosSampleInBgen: find positions of model samples in BGEN
// Replaces Rcpp::match() with std::unordered_map lookup
// ============================================================
void BgenClass::setPosSampleInBgen(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in Bgen files...." << std::endl;

    // If no sample IDs provided, use all BGEN samples in order
    if (t_SampleInModel.empty()) {
        m_N = m_N0;
        std::cout << "  No sampleIDs in null model; using all " << m_N0
                  << " samples from BGEN file." << std::endl;
        m_posSampleInModel.resize(m_N0);
        for (uint32_t i = 0; i < m_N0; i++) {
            m_posSampleInModel[i] = i;
        }
        std::cout << "Number of samples in analysis: " << m_N << std::endl;
        return;
    }

    m_N = t_SampleInModel.size();

    // Build map of BGEN sample -> index (0-based)
    std::unordered_map<std::string, uint32_t> bgenSampleMap;
    for (uint32_t i = 0; i < m_N0; i++) {
        bgenSampleMap[m_SampleInBgen[i]] = i;
    }

    // Verify all model samples exist in BGEN
    for (uint32_t i = 0; i < m_N; i++) {
        if (bgenSampleMap.find(t_SampleInModel[i]) == bgenSampleMap.end()) {
            throw std::runtime_error(
                "At least one subject requested is not in BGEN file. "
                "Sample not found: " + t_SampleInModel[i]);
        }
    }

    // Build map of model sample -> model index (0-based)
    std::unordered_map<std::string, int32_t> modelSampleMap;
    for (uint32_t i = 0; i < m_N; i++) {
        modelSampleMap[t_SampleInModel[i]] = (int32_t)i;
    }

    // For each BGEN sample, what is its position in the model?
    // m_posSampleInModel[bgen_idx] = model_idx, or -1 if not in model
    m_posSampleInModel.resize(m_N0);
    for (uint32_t i = 0; i < m_N0; i++) {
        auto it = modelSampleMap.find(m_SampleInBgen[i]);
        if (it != modelSampleMap.end()) {
            m_posSampleInModel[i] = it->second;
        } else {
            m_posSampleInModel[i] = -1;
        }
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// Parse2: decompress and parse genotype probabilities
// Direct port from SAIGE/src/BGEN.cpp::Parse2
// ============================================================
void BgenClass::Parse2(unsigned char* buf, uint32_t bufLen,
                        const unsigned char* zBuf, uint32_t zBufLen,
                        std::string& snpName,
                        arma::vec& dosages,
                        double& AC, double& AF,
                        std::vector<uint>& indexforMissing,
                        double& info,
                        std::vector<uint>& indexNonZero,
                        bool isImputation)
{
    // Decompress
    if (CompressedSNPBlocks == COMPRESSION_ZLIB) {
        z_stream strm = {};
        strm.next_in = const_cast<Bytef*>(zBuf);
        strm.avail_in = zBufLen;
        strm.next_out = buf;
        strm.avail_out = bufLen;

        if (inflateInit(&strm) != Z_OK) {
            std::cerr << "inflateInit failed" << std::endl;
            return;
        }

        int ret = inflate(&strm, Z_FINISH);
        if (ret != Z_STREAM_END) {
            std::cerr << "inflate failed with code " << ret << std::endl;
        }

        inflateEnd(&strm);
    }
    if (CompressedSNPBlocks == COMPRESSION_ZSTD) {
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        size_t actual_decompressed_size = ZSTD_decompressDCtx(dctx, buf, bufLen, zBuf, zBufLen);
        if (ZSTD_isError(actual_decompressed_size)) {
            std::cerr << "Decompression failed: " << ZSTD_getErrorName(actual_decompressed_size) << std::endl;
            ZSTD_freeDCtx(dctx);
            return;
        }
        ZSTD_freeDCtx(dctx);
    }

    unsigned char* bufAt = buf;
    uint32_t N = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24);
    bufAt += 4;

    if (N != m_N0) {
        std::cerr << "ERROR: " << snpName << " has N = " << N
                  << " (mismatch with header block)" << std::endl;
        throw std::runtime_error("BGEN sample count mismatch in variant " + snpName);
    }

    uint32_t K = bufAt[0] | (bufAt[1] << 8);
    bufAt += 2;
    if (K != 2U) {
        std::cerr << "ERROR: " << snpName << " has K = " << K
                  << " (non-bi-allelic)" << std::endl;
        throw std::runtime_error("Non-bi-allelic variant in BGEN: " + snpName);
    }

    uint32_t Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
        std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin
                  << " (not 2)" << std::endl;
        throw std::runtime_error("Unsupported minimum ploidy in BGEN: " + snpName);
    }

    uint32_t Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
        std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax
                  << " (not 2)" << std::endl;
        throw std::runtime_error("Unsupported maximum ploidy in BGEN: " + snpName);
    }

    // Read per-sample ploidy/missingness bytes
    const unsigned char* ploidyMissBytes = bufAt;
    for (uint32_t i = 0; i < N; i++) {
        uint32_t ploidyMiss = *bufAt; bufAt++;
        if (ploidyMiss != 2U && ploidyMiss != 130U) {
            std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = "
                      << ploidyMiss << " (not 2 or 130)" << std::endl;
            throw std::runtime_error("Unsupported ploidy/missingness in BGEN: " + snpName);
        }
    }

    uint32_t Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
        std::cerr << "ERROR: " << snpName << " has Phased = " << Phased
                  << " (not 0)" << std::endl;
        throw std::runtime_error("Phased data not supported in BGEN reader: " + snpName);
    }

    uint32_t B = *bufAt; bufAt++;
    if (B != 8U) {
        std::cerr << "ERROR: " << snpName << " has B = " << B
                  << " (not 8)" << std::endl;
        throw std::runtime_error("Unsupported bit depth in BGEN: " + snpName);
    }

    // Build lookup table for 8-bit probabilities
    double lut[256];
    for (int i = 0; i <= 255; i++) {
        lut[i] = i / 255.0;
    }

    // Parse genotype probabilities
    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0;
    double p11, p10, dosage, eij, fij;
    double dosage_new;

    dosages.set_size(m_N);
    dosages.fill(arma::datum::nan);

    std::size_t missing_cnt = 0;

    for (uint32_t i = 0; i < N; i++) {
        if (ploidyMissBytes[i] != 130U) {
            // Not missing
            p11 = lut[*bufAt]; bufAt++;
            p10 = lut[*bufAt]; bufAt++;

            if (m_posSampleInModel[i] >= 0) {
                dosage = 2 * p11 + p10;

                // SAIGE default is ref-first: dosage_new = 2 - dosage
                dosage_new = 2 - dosage;

                eij = dosage;
                fij = 4 * p11 + p10;
                sum_eij += eij;
                sum_fij_minus_eij2 += fij - eij * eij;

                dosages[m_posSampleInModel[i]] = dosage_new;
                if (dosage_new > 0) {
                    indexNonZero.push_back(m_posSampleInModel[i]);
                }
                sum_eij_sub += eij;
            }
        } else {
            // Missing (ploidy = 130)
            bufAt += 2;
            if (m_posSampleInModel[i] >= 0) {
                indexforMissing.push_back(m_posSampleInModel[i]);
                ++missing_cnt;
                dosages[m_posSampleInModel[i]] = -1;
            }
        }
    }

    // Compute AC and AF
    AC = 2 * ((double)(m_N - missing_cnt)) - sum_eij_sub;
    if (m_N == missing_cnt) {
        AF = 0;
    } else {
        AF = AC / 2 / ((double)(m_N - missing_cnt));
    }

    // Compute imputation info
    double thetaHat = sum_eij / (2 * (m_N - missing_cnt));
    if (isImputation) {
        info = (thetaHat == 0 || thetaHat == 1) ? 1 :
            1 - sum_fij_minus_eij2 / (2 * (m_N - missing_cnt) * thetaHat * (1 - thetaHat));
    } else {
        info = 1.0;
    }
}

// ============================================================
// getOneMarker: read one marker from BGEN
// Direct port from SAIGE/src/BGEN.cpp::getOneMarker
// t_gIndex / t_gIndex_prev are byte positions (BGEN uses byte seeking)
// ============================================================
void BgenClass::getOneMarker(uint64_t& t_gIndex_prev,
                              uint64_t& t_gIndex,
                              std::string& t_ref,
                              std::string& t_alt,
                              std::string& t_marker,
                              uint32_t& t_pd,
                              std::string& t_chr,
                              double& t_altFreq,
                              double& t_altCounts,
                              double& t_missingRate,
                              double& t_imputeInfo,
                              bool& t_isOutputIndexForMissing,
                              std::vector<uint>& t_indexForMissing,
                              bool& t_isOnlyOutputNonZero,
                              std::vector<uint>& t_indexForNonZero,
                              bool& t_isBoolRead,
                              arma::vec& dosages,
                              bool t_isImputation)
{
    // Seek to the variant's ABSOLUTE byte offset (file_start_position from .bgi).
    // t_gIndex == 0 => "no seek, read the next variant sequentially" (single-variant
    // streaming sets all indices to 0; byte 0 is never a valid variant-block start
    // because the BGEN header precedes the data). The previous relative SEEK_CUR by
    // (t_gIndex - t_gIndex_prev) overshot by the previous block's size — the file
    // pointer sits at the END of the previous block, not at t_gIndex_prev — which
    // desynced region/by-ID reads and produced garbage variants ("inflate failed").
    (void)t_gIndex_prev;
    if (t_gIndex > 0) {
        fseek(m_fin, t_gIndex, SEEK_SET);
    }

    std::string SNPID, RSID, chromosome, first_allele, second_allele;
    uint32_t position;
    double AC, AF, info;

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    char snpID[65536], rsID[65536], chrStr[65536];
    uint16_t LS;
    size_t numBoolRead = fread(&LS, 2, 1, m_fin);

    if (numBoolRead > 0) {
        t_isBoolRead = true;

        fread(snpID, 1, LS, m_fin); snpID[LS] = '\0';

        uint16_t LR;
        fread(&LR, 2, 1, m_fin);
        fread(rsID, 1, LR, m_fin); rsID[LR] = '\0';
        RSID = std::string(rsID) == "." ? snpID : rsID;

        uint16_t LC;
        fread(&LC, 2, 1, m_fin);
        fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
        chromosome = std::string(chrStr);

        uint32_t physpos;
        fread(&physpos, 4, 1, m_fin);
        position = physpos;

        uint16_t K;
        fread(&K, 2, 1, m_fin);

        uint32_t LA;
        fread(&LA, 4, 1, m_fin);
        if (LA >= allele1.size()) {
            allele1.resize(2 * LA);
        }
        fread(allele1.data(), 1, LA, m_fin); allele1[LA] = '\0';
        first_allele = std::string(allele1.data());

        uint32_t LB;
        fread(&LB, 4, 1, m_fin);
        if (LB >= allele0.size()) {
            allele0.resize(2 * LB);
        }
        fread(allele0.data(), 1, LB, m_fin); allele0[LB] = '\0';
        second_allele = std::string(allele0.data());

        uint32_t C;
        fread(&C, 4, 1, m_fin);
        if (C > m_zBuf.size()) m_zBuf.resize(C - 4);

        uint32_t D;
        fread(&D, 4, 1, m_fin);
        m_zBufLens = C - 4;
        m_bufLens = D;
        fread(&m_zBuf[0], 1, C - 4, m_fin);

        AC = 0;
        AF = 0;
        info = 0;

        if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens);

        Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens,
               RSID, dosages, AC, AF, t_indexForMissing, info,
               t_indexForNonZero, t_isImputation);

        // Default setting is "ref-first" (matching SAIGE 03-14-2021 update)
        t_alt = second_allele;
        t_ref = first_allele;
        t_marker = RSID;
        t_pd = position;
        t_chr = chromosome;
        t_altFreq = AF;
        t_altCounts = AC;
        t_imputeInfo = info;
        t_missingRate = (double)t_indexForMissing.size() / (double)m_N;

        if (m_AlleleOrder == "alt-first") {
            t_alt = first_allele;
            t_ref = second_allele;
            t_altFreq = 1 - t_altFreq;
            t_altCounts = t_altFreq * 2 * ((double)m_N - (double)t_indexForMissing.size());
            for (unsigned int i = 0; i < dosages.n_elem; i++) {
                if (dosages[i] >= 0) {  // don't flip missing values (-1)
                    dosages[i] = 2 - dosages[i];
                }
            }
        }

        // Store marker info for index maps
        m_chr.push_back(t_chr);
        m_pd.push_back(t_pd);
        m_ref.push_back(t_ref);
        m_alt.push_back(t_alt);
        m_MarkerInBgen.push_back(t_marker);

    } else {
        t_isBoolRead = false;
    }
}

// ============================================================
// Phase D: readRawBlock — pure I/O, runs single-threaded on the streamer's
// reader thread. Mutates only m_fin; emits a BgenRawBlock that decoders can
// process independently (they own scratch + zBuf).
// ============================================================
bool BgenClass::readRawBlock(uint64_t t_markerIndex,
                              uint64_t t_byteOffset,
                              BgenRawBlock& out)
{
    out.markerIndex = t_markerIndex;
    out.byteOffset  = t_byteOffset;
    out.eof         = false;

    if (t_byteOffset > 0) {
        fseek(m_fin, t_byteOffset, SEEK_SET);
    }

    char snpID[65536], rsID[65536], chrStr[65536];
    uint16_t LS;
    size_t numBoolRead = fread(&LS, 2, 1, m_fin);
    if (numBoolRead == 0) {
        out.eof = true;
        return false;
    }

    fread(snpID, 1, LS, m_fin); snpID[LS] = '\0';
    out.snpID = std::string(snpID);

    uint16_t LR;
    fread(&LR, 2, 1, m_fin);
    fread(rsID, 1, LR, m_fin); rsID[LR] = '\0';
    out.rsID = (std::string(rsID) == "." ? out.snpID : std::string(rsID));

    uint16_t LC;
    fread(&LC, 2, 1, m_fin);
    fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
    out.chr = std::string(chrStr);

    uint32_t physpos;
    fread(&physpos, 4, 1, m_fin);
    out.physpos = physpos;

    uint16_t K;
    fread(&K, 2, 1, m_fin);
    (void)K;  // validated in parseBlock

    uint32_t LA;
    fread(&LA, 4, 1, m_fin);
    if (LA >= allele1.size()) allele1.resize(2 * LA);
    fread(allele1.data(), 1, LA, m_fin); allele1[LA] = '\0';
    std::string first_allele(allele1.data());

    uint32_t LB;
    fread(&LB, 4, 1, m_fin);
    if (LB >= allele0.size()) allele0.resize(2 * LB);
    fread(allele0.data(), 1, LB, m_fin); allele0[LB] = '\0';
    std::string second_allele(allele0.data());

    out.alleles.clear();
    out.alleles.push_back(first_allele);
    out.alleles.push_back(second_allele);

    uint32_t C;
    fread(&C, 4, 1, m_fin);
    out.C = C;

    uint32_t D;
    fread(&D, 4, 1, m_fin);
    out.D = D;

    out.zBuf.resize(C - 4);
    fread(out.zBuf.data(), 1, C - 4, m_fin);

    return true;
}

// ============================================================
// Phase D: parseBlock — pure decode, callable from any decoder thread.
// Uses caller-supplied scratch buffer for the decompressed bytes.
// Mirrors Parse2 + the post-Parse2 allele-order-swap logic from getOneMarker.
// ============================================================
void BgenClass::parseBlock(const BgenClass::ParseCtx& ctx,
                            const BgenRawBlock& in,
                            BgenDecodedMarker& out,
                            std::vector<unsigned char>& scratchBuf)
{
    out.markerIndex = in.markerIndex;
    out.snpID       = in.snpID;
    out.rsID        = in.rsID;
    out.chr         = in.chr;
    out.physpos     = in.physpos;
    out.alleles     = in.alleles;       // [first_allele, second_allele] for now
    out.indexForMissing.clear();
    out.indexForNonZero.clear();
    out.valid       = false;

    if (in.eof) return;

    // Grow scratch as needed (per-thread, reused across markers)
    if (in.D > scratchBuf.size()) scratchBuf.resize(in.D);

    // Decompress into scratchBuf
    if (ctx.compressedSNPBlocks == COMPRESSION_ZLIB) {
        z_stream strm = {};
        strm.next_in   = const_cast<Bytef*>(in.zBuf.data());
        strm.avail_in  = in.zBuf.size();
        strm.next_out  = scratchBuf.data();
        strm.avail_out = in.D;
        if (inflateInit(&strm) != Z_OK) {
            std::cerr << "inflateInit failed" << std::endl;
            return;
        }
        int ret = inflate(&strm, Z_FINISH);
        if (ret != Z_STREAM_END) {
            std::cerr << "inflate failed with code " << ret << std::endl;
        }
        inflateEnd(&strm);
    } else if (ctx.compressedSNPBlocks == COMPRESSION_ZSTD) {
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        size_t actual = ZSTD_decompressDCtx(dctx, scratchBuf.data(), in.D,
                                             in.zBuf.data(), in.zBuf.size());
        if (ZSTD_isError(actual)) {
            std::cerr << "Decompression failed: " << ZSTD_getErrorName(actual) << std::endl;
            ZSTD_freeDCtx(dctx);
            return;
        }
        ZSTD_freeDCtx(dctx);
    }

    unsigned char* bufAt = scratchBuf.data();
    uint32_t N = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24);
    bufAt += 4;
    if (N != ctx.N0) {
        std::cerr << "ERROR: " << in.rsID << " has N = " << N
                  << " (mismatch with header block)" << std::endl;
        throw std::runtime_error("BGEN sample count mismatch in variant " + in.rsID);
    }

    uint32_t K = bufAt[0] | (bufAt[1] << 8); bufAt += 2;
    if (K != 2U) {
        throw std::runtime_error("Non-bi-allelic variant in BGEN: " + in.rsID);
    }
    uint32_t Pmin = *bufAt++; if (Pmin != 2U) {
        throw std::runtime_error("Unsupported minimum ploidy in BGEN: " + in.rsID);
    }
    uint32_t Pmax = *bufAt++; if (Pmax != 2U) {
        throw std::runtime_error("Unsupported maximum ploidy in BGEN: " + in.rsID);
    }

    const unsigned char* ploidyMissBytes = bufAt;
    for (uint32_t i = 0; i < N; i++) {
        uint32_t ploidyMiss = *bufAt++;
        if (ploidyMiss != 2U && ploidyMiss != 130U) {
            throw std::runtime_error("Unsupported ploidy/missingness in BGEN: " + in.rsID);
        }
    }

    uint32_t Phased = *bufAt++;
    if (Phased != 0U) {
        throw std::runtime_error("Phased data not supported in BGEN reader: " + in.rsID);
    }
    uint32_t B = *bufAt++;
    if (B != 8U) {
        throw std::runtime_error("Unsupported bit depth in BGEN: " + in.rsID);
    }

    double lut[256];
    for (int i = 0; i <= 255; i++) lut[i] = i / 255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0;
    double p11, p10, dosage, eij, fij, dosage_new;

    out.dosages.set_size(ctx.N);
    out.dosages.fill(arma::datum::nan);
    std::size_t missing_cnt = 0;

    const auto& posMap = *ctx.posSampleInModel;
    for (uint32_t i = 0; i < N; i++) {
        if (ploidyMissBytes[i] != 130U) {
            p11 = lut[*bufAt++];
            p10 = lut[*bufAt++];
            if (posMap[i] >= 0) {
                dosage = 2 * p11 + p10;
                dosage_new = 2 - dosage;          // ref-first default
                eij = dosage;
                fij = 4 * p11 + p10;
                sum_eij += eij;
                sum_fij_minus_eij2 += fij - eij * eij;
                out.dosages[posMap[i]] = dosage_new;
                if (dosage_new > 0) out.indexForNonZero.push_back(posMap[i]);
                sum_eij_sub += eij;
            }
        } else {
            bufAt += 2;
            if (posMap[i] >= 0) {
                out.indexForMissing.push_back(posMap[i]);
                ++missing_cnt;
                out.dosages[posMap[i]] = -1;
            }
        }
    }

    double AC = 2 * ((double)(ctx.N - missing_cnt)) - sum_eij_sub;
    double AF = (ctx.N == missing_cnt) ? 0.0 : AC / 2 / ((double)(ctx.N - missing_cnt));

    double thetaHat = sum_eij / (2 * (ctx.N - missing_cnt));
    double info;
    if (ctx.isImputation) {
        info = (thetaHat == 0 || thetaHat == 1) ? 1.0 :
            1.0 - sum_fij_minus_eij2 / (2 * (ctx.N - missing_cnt) * thetaHat * (1 - thetaHat));
    } else {
        info = 1.0;
    }

    // Default ref-first: ref=second_allele, alt=first_allele inverted? No —
    // mirror existing getOneMarker: ref-first means t_ref=first, t_alt=second
    // (see getOneMarker lines ~1378-1380).
    std::string ref = in.alleles[0];      // first_allele
    std::string alt = in.alleles[1];      // second_allele
    double altFreq = AF;
    double altCounts = AC;

    if (ctx.alleleOrder == "alt-first") {
        std::swap(ref, alt);
        altFreq = 1 - altFreq;
        altCounts = altFreq * 2 * ((double)ctx.N - (double)out.indexForMissing.size());
        for (unsigned int i = 0; i < out.dosages.n_elem; i++) {
            if (out.dosages[i] >= 0) {
                out.dosages[i] = 2 - out.dosages[i];
            }
        }
    }

    out.alleles.clear();
    out.alleles.push_back(ref);
    out.alleles.push_back(alt);
    out.altFreq     = altFreq;
    out.altCounts   = altCounts;
    out.info        = info;
    out.missingRate = (double)out.indexForMissing.size() / (double)ctx.N;
    out.valid       = true;
}

// ============================================================
// Phase D: BgenStreamer — 1 reader + N decoders + bounded queue,
// order-preserving consumer via decodedMap[markerIndex].
// ============================================================
BgenStreamer::BgenStreamer(BgenClass* bgen,
                           const std::vector<std::string>& genoIndex,
                           bool isImputation,
                           int nDecoders,
                           size_t queueCapacity)
    : m_bgen(bgen), m_queueCap(queueCapacity), m_nDecoders(nDecoders)
{
    m_byteOffsets.reserve(genoIndex.size());
    for (const auto& s : genoIndex) {
        char* end;
        m_byteOffsets.push_back(std::strtoull(s.c_str(), &end, 10));
    }

    m_ctx.N0                  = bgen->getN0();
    m_ctx.N                   = bgen->getN();
    m_ctx.compressedSNPBlocks = bgen->getCompressedSNPBlocks();
    m_ctx.posSampleInModel    = &bgen->getPosSampleInModel();
    m_ctx.alleleOrder         = bgen->getAlleleOrder();
    m_ctx.isImputation        = isImputation;

    m_readerThread = std::thread(&BgenStreamer::readerLoop, this);
    m_decoderThreads.reserve(nDecoders);
    for (int i = 0; i < nDecoders; i++) {
        m_decoderThreads.emplace_back(&BgenStreamer::decoderLoop, this);
    }
}

BgenStreamer::~BgenStreamer()
{
    m_stop.store(true);
    // Wake everyone so they can observe m_stop and exit
    m_rawNotEmpty.notify_all();
    m_rawNotFull.notify_all();
    m_decReady.notify_all();
    m_decNotFull.notify_all();
    if (m_readerThread.joinable()) m_readerThread.join();
    for (auto& t : m_decoderThreads) {
        if (t.joinable()) t.join();
    }
}

void BgenStreamer::readerLoop()
{
    const uint64_t M = m_byteOffsets.size();
    for (uint64_t i = 0; i < M; i++) {
        if (m_stop.load()) break;

        // Memory backpressure: don't read so far ahead that m_decoded grows
        // unbounded. If the consumer is slow / out-of-order, pause here until
        // it has drained some entries.
        {
            std::unique_lock<std::mutex> dlk(m_decMtx);
            m_decNotFull.wait(dlk, [&]() {
                return m_decoded.size() < m_decodedCap || m_stop.load();
            });
            if (m_stop.load()) break;
        }

        BgenRawBlock blk;
        bool ok = m_bgen->readRawBlock(i, m_byteOffsets[i], blk);
        if (!ok) {
            // EOF mid-stream: record where decoders should stop and bail out
            std::unique_lock<std::mutex> dlk(m_decMtx);
            m_decodersFinishedAt = i;
            m_decReady.notify_all();
            break;
        }

        std::unique_lock<std::mutex> lk(m_rawMtx);
        m_rawNotFull.wait(lk, [&]() {
            return m_rawQueue.size() < m_queueCap || m_stop.load();
        });
        if (m_stop.load()) break;
        m_rawQueue.push(std::move(blk));
        m_rawNotEmpty.notify_one();
    }

    {
        std::unique_lock<std::mutex> lk(m_rawMtx);
        m_readerDone = true;
        m_rawNotEmpty.notify_all();
    }
    // Wake any consumer waiting on a marker that will never arrive (e.g. if
    // m_stop fired or we exited early without setting m_decodersFinishedAt).
    {
        std::unique_lock<std::mutex> dlk(m_decMtx);
        m_decReady.notify_all();
    }
}

void BgenStreamer::decoderLoop()
{
    std::vector<unsigned char> scratch;     // per-thread, reused across markers
    while (true) {
        BgenRawBlock raw;
        {
            std::unique_lock<std::mutex> lk(m_rawMtx);
            m_rawNotEmpty.wait(lk, [&]() {
                return !m_rawQueue.empty() || m_readerDone || m_stop.load();
            });
            if (m_stop.load()) return;
            if (m_rawQueue.empty()) {
                if (m_readerDone) return;
                continue;
            }
            raw = std::move(m_rawQueue.front());
            m_rawQueue.pop();
            m_rawNotFull.notify_one();
        }

        BgenDecodedMarker dm;
        try {
            BgenClass::parseBlock(m_ctx, raw, dm, scratch);
        } catch (const std::exception& e) {
            std::cerr << "BgenStreamer decoder error at idx "
                      << raw.markerIndex << ": " << e.what() << std::endl;
            dm.markerIndex = raw.markerIndex;
            dm.valid = false;
        }

        {
            std::unique_lock<std::mutex> lk(m_decMtx);
            m_decoded[raw.markerIndex] = std::move(dm);
            m_decReady.notify_all();
        }
    }
}

bool BgenStreamer::getMarker(uint64_t markerIndex, BgenDecodedMarker& out)
{
    // Block until the specific marker we asked for is in m_decoded, or until
    // it is provably never going to arrive (out-of-range, truncated stream,
    // or stop signal). Each markerIndex is requested at most once by the
    // OpenMP consumer, so we evict from the map after retrieval to bound
    // memory.
    std::unique_lock<std::mutex> lk(m_decMtx);
    const uint64_t total = m_byteOffsets.size();

    if (markerIndex >= total) return false;

    m_decReady.wait(lk, [&]() {
        if (m_stop.load()) return true;
        if (markerIndex >= m_decodersFinishedAt) return true;  // truncated past us
        return m_decoded.find(markerIndex) != m_decoded.end();
    });

    if (m_stop.load()) return false;
    if (markerIndex >= m_decodersFinishedAt) return false;

    auto it = m_decoded.find(markerIndex);
    if (it == m_decoded.end()) return false;
    out = std::move(it->second);
    m_decoded.erase(it);
    // Tell reader thread the decoded-map now has room (backpressure release).
    m_decNotFull.notify_one();
    return out.valid;
}

// Kept as a no-op (helper present only to satisfy the header declaration; the
// real cross-queue check would require holding m_rawMtx while m_decMtx is
// held, which we deliberately avoid).
bool BgenStreamer::m_rawQueueIsEmptyHelper() const { return true; }

// ============================================================
// populateFromBgi: pre-populate marker metadata from BGEN .bgi SQLite index
// so by-ID lookups (e.g. SAIGE-GENE+ region tests) work without first
// streaming every marker. Schema (bgenix v1.1.x):
//   Variant(chromosome TEXT, position INT, rsid TEXT, number_of_alleles INT,
//           allele1 TEXT, allele2 TEXT, file_start_position INT,
//           size_in_bytes INT, PRIMARY KEY(...))
// We read (chromosome, position, rsid, allele1, allele2) for every row in
// file_start_position order so the resulting indices match the order the
// rest of BgenClass expects (sequential marker stream).
// ============================================================
uint32_t BgenClass::populateFromBgi(const std::string& t_bgiFileName)
{
    sqlite3* db = nullptr;
    if (sqlite3_open_v2(t_bgiFileName.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
        if (db) sqlite3_close(db);
        throw std::runtime_error("populateFromBgi: cannot open .bgi file: " + t_bgiFileName);
    }
    const char* sql =
        "SELECT chromosome, position, rsid, allele1, allele2, file_start_position "
        "FROM Variant ORDER BY file_start_position";
    sqlite3_stmt* stmt = nullptr;
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        std::string err = sqlite3_errmsg(db);
        sqlite3_close(db);
        throw std::runtime_error("populateFromBgi: prepare failed: " + err);
    }
    // Reserve up to the BGEN-header-declared marker count to avoid repeated growth.
    if (m_M0 > 0) {
        m_chr.reserve(m_M0); m_pd.reserve(m_M0); m_ref.reserve(m_M0);
        m_alt.reserve(m_M0); m_MarkerInBgen.reserve(m_M0);
    }
    m_chr.clear(); m_pd.clear(); m_ref.clear(); m_alt.clear(); m_MarkerInBgen.clear();
    m_byteOffset.clear();
    if (m_M0 > 0) m_byteOffset.reserve(m_M0);
    uint32_t n = 0;
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        const unsigned char* chr_  = sqlite3_column_text(stmt, 0);
        int                  pos   = sqlite3_column_int (stmt, 1);
        const unsigned char* rsid  = sqlite3_column_text(stmt, 2);
        const unsigned char* a1    = sqlite3_column_text(stmt, 3);
        const unsigned char* a2    = sqlite3_column_text(stmt, 4);
        sqlite3_int64        fpos  = sqlite3_column_int64(stmt, 5);  // byte offset of variant block
        m_chr.emplace_back(chr_ ? reinterpret_cast<const char*>(chr_) : "");
        m_pd.push_back(static_cast<uint32_t>(pos));
        m_ref.emplace_back(a1 ? reinterpret_cast<const char*>(a1) : "");
        m_alt.emplace_back(a2 ? reinterpret_cast<const char*>(a2) : "");
        m_MarkerInBgen.emplace_back(rsid ? reinterpret_cast<const char*>(rsid) : "");
        m_byteOffset.push_back(static_cast<uint64_t>(fpos));
        ++n;
    }
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return n;
}

// ============================================================
// getMarkerIDToIndex: build chr:pos:ref:alt -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> BgenClass::getMarkerIDToIndex()
{
    // The region/by-ID read path seeks via getOneMarker(t_gIndex = BYTE OFFSET).
    // So the map must return each variant's .bgi file_start_position, NOT the
    // sequential row index i. Returning i (the old bug) made getOneMarker fseek to
    // byte == row index -> garbage block -> "inflate failed" / N=0 / segfault.
    const bool haveOffsets = (m_byteOffset.size() == m_chr.size());
    std::unordered_map<std::string, uint32_t> idMap;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        uint64_t v = haveOffsets ? m_byteOffset[i] : static_cast<uint64_t>(i);
        if (v > 0xFFFFFFFFULL) {
            throw std::runtime_error("BGEN variant byte offset exceeds 32 bits; the "
                "marker-ID->index map is uint32_t. Widen it to uint64_t for bgens >4GB.");
        }
        uint32_t off = static_cast<uint32_t>(v);
        std::string id1 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_ref[i] + ":" + m_alt[i];
        idMap[id1] = off;
        std::string id2 = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":"
                         + m_alt[i] + ":" + m_ref[i];
        if (idMap.find(id2) == idMap.end()) {
            idMap[id2] = off;
        }
    }
    return idMap;
}

// ============================================================
// getMarkerNameToIndex: build marker ID -> index map
// ============================================================
std::unordered_map<std::string, uint32_t> BgenClass::getMarkerNameToIndex()
{
    // Same as getMarkerIDToIndex: return the .bgi byte offset (what getOneMarker
    // seeks to), not the sequential row index.
    const bool haveOffsets = (m_byteOffset.size() == m_MarkerInBgen.size());
    std::unordered_map<std::string, uint32_t> nameMap;
    for (uint32_t i = 0; i < m_MarkerInBgen.size(); i++) {
        uint64_t v = haveOffsets ? m_byteOffset[i] : static_cast<uint64_t>(i);
        if (v > 0xFFFFFFFFULL) {
            throw std::runtime_error("BGEN variant byte offset exceeds 32 bits; the "
                "marker-name->index map is uint32_t. Widen it to uint64_t for bgens >4GB.");
        }
        nameMap[m_MarkerInBgen[i]] = static_cast<uint32_t>(v);
    }
    return nameMap;
}

// ============================================================
// getMarkerIDToChrom: variant ID -> CHROM, same keys as getMarkerIDToIndex().
// Populated from the .bgi index (populateFromBgi); empty if none was loaded.
// ============================================================
std::unordered_map<std::string, std::string> BgenClass::getMarkerIDToChrom()
{
    std::unordered_map<std::string, std::string> chrMap;
    for (uint32_t i = 0; i < m_chr.size(); i++) {
        std::string base = m_chr[i] + ":" + std::to_string(m_pd[i]) + ":";
        chrMap[base + m_ref[i] + ":" + m_alt[i]] = m_chr[i];
        chrMap.emplace(base + m_alt[i] + ":" + m_ref[i], m_chr[i]);
    }
    return chrMap;
}

// ============================================================
// closegenofile: close the BGEN file handle
// ============================================================
void BgenClass::closegenofile()
{
    if (m_fin) {
        fclose(m_fin);
        m_fin = nullptr;
    }
}

} // namespace BGEN


// ============================================================
// PGEN namespace: PGEN v2 reader (mode 0x02 basic variant-major)
// Standalone implementation without pgenlib dependency
// ============================================================
namespace PGEN {

// ============================================================
// Constructor
// ============================================================
PgenClass::PgenClass(const std::string& t_pgenFile,
                     const std::string& t_psamFile,
                     const std::string& t_pvarFile,
                     std::vector<std::string>& t_SampleInModel)
    : m_pgenFile(t_pgenFile), m_pvarFile(t_pvarFile), m_psamFile(t_psamFile),
      m_fin(nullptr), m_M0(0), m_N0(0), m_mode(0), m_dataOffset(0),
      m_bytesPerVariant(0), m_M(0), m_N(0)
{
    std::cout << "pgenFile: " << t_pgenFile << std::endl;
    std::cout << "psamFile: " << t_psamFile << std::endl;
    std::cout << "pvarFile: " << t_pvarFile << std::endl;

    readPvarFile();
    readPsamFile();
    readPgenHeader();
    setPosSampleInPgen(t_SampleInModel);

    // Allocate raw buffer for one variant
    m_bytesPerVariant = (m_N0 + 3) / 4;
    m_OneMarkerRaw.resize(m_bytesPerVariant);
}

// ============================================================
// Destructor
// ============================================================
PgenClass::~PgenClass()
{
    closegenofile();
}

// ============================================================
// readPgenHeader: read and validate .pgen file header
//
// PGEN format (mode 0x02):
//   Bytes 0-1: magic 0x6c 0x1b
//   Byte 2:    mode (0x02 = basic variant-major)
//   Bytes 3-6: variant count (4 bytes, little-endian)
//   Bytes 7-10: sample count (4 bytes, little-endian)
//   Byte 11:   header control byte (0x00 for mode 0x02)
//   Data starts at byte 12
// ============================================================
void PgenClass::readPgenHeader()
{
    m_fin = fopen(m_pgenFile.c_str(), "rb");
    if (!m_fin) {
        throw std::runtime_error("Cannot open PGEN file: " + m_pgenFile);
    }

    // Read magic bytes
    uint8_t magic[2];
    if (fread(magic, 1, 2, m_fin) != 2) {
        throw std::runtime_error("Cannot read magic bytes from PGEN file: " + m_pgenFile);
    }
    if (magic[0] != 0x6c || magic[1] != 0x1b) {
        throw std::runtime_error("Invalid PGEN magic bytes. Expected 0x6c 0x1b, got 0x"
            + std::to_string(magic[0]) + " 0x" + std::to_string(magic[1]));
    }

    // Read mode byte
    if (fread(&m_mode, 1, 1, m_fin) != 1) {
        throw std::runtime_error("Cannot read mode byte from PGEN file: " + m_pgenFile);
    }

    if (m_mode == 0x01) {
        // PLINK 1 variant-major mode (essentially .bed format in .pgen container)
        // Header is just 3 bytes (magic + mode), data starts at byte 3
        // Use fam/bim sample/variant counts from .psam/.pvar
        m_N0 = m_SampleInPgen.size();
        m_M0 = m_M;
        m_dataOffset = 3;
        m_bytesPerVariant = (m_N0 + 3) / 4;
        std::cout << "PGEN mode 0x01 (PLINK 1 variant-major)" << std::endl;
        std::cout << "  Variants: " << m_M0 << " (from .pvar)" << std::endl;
        std::cout << "  Samples:  " << m_N0 << " (from .psam)" << std::endl;
    } else if (m_mode == 0x02) {
        // PLINK 2 basic variant-major mode
        // Read variant count
        if (fread(&m_M0, 4, 1, m_fin) != 1) {
            throw std::runtime_error("Cannot read variant count from PGEN header");
        }
        // Read sample count
        if (fread(&m_N0, 4, 1, m_fin) != 1) {
            throw std::runtime_error("Cannot read sample count from PGEN header");
        }
        // Read header control byte
        uint8_t headerCtrl;
        if (fread(&headerCtrl, 1, 1, m_fin) != 1) {
            throw std::runtime_error("Cannot read header control byte from PGEN header");
        }

        // For mode 0x02, header control bits 0-5 should be zero
        // Bits 6-7 encode nonref flags: 00=unstored, 01=all ref/alt, 10=never, 11=explicit
        uint8_t nonrefFlags = (headerCtrl >> 6) & 0x03;
        m_dataOffset = 12;  // 2(magic) + 1(mode) + 4(M) + 4(N) + 1(ctrl)

        // If nonref flags are explicitly stored (bits 6-7 = 11), skip the bitarray
        if (nonrefFlags == 3) {
            uint64_t nonrefBytes = (m_M0 + 7) / 8;
            m_dataOffset += nonrefBytes;
        }

        std::cout << "PGEN mode 0x02 (basic variant-major)" << std::endl;
        std::cout << "  Variants in header: " << m_M0 << std::endl;
        std::cout << "  Samples in header:  " << m_N0 << std::endl;
        std::cout << "  Header control:     0x" << std::hex << (int)headerCtrl << std::dec << std::endl;
        std::cout << "  Data offset:        " << m_dataOffset << std::endl;
    } else {
        throw std::runtime_error(
            "Unsupported PGEN mode 0x" + std::to_string((int)m_mode) +
            ". This standalone reader only supports mode 0x01 (PLINK 1) and 0x02 (basic variant-major). "
            "For mode 0x10/0x11 (variable-type records with LD compression or dosage), "
            "the full pgenlib library would be required.");
    }

    // Validate consistency with .pvar
    if (m_M0 != m_M) {
        std::cerr << "WARNING: PGEN header variant count (" << m_M0
                  << ") differs from .pvar count (" << m_M << "). Using .pvar count." << std::endl;
        // Trust .pvar count (pgen may have extra variants)
    }

    // Validate consistency with .psam
    uint32_t psamN = m_SampleInPgen.size();
    if (m_N0 != psamN) {
        std::cerr << "WARNING: PGEN header sample count (" << m_N0
                  << ") differs from .psam count (" << psamN << ")." << std::endl;
    }

    m_bytesPerVariant = (m_N0 + 3) / 4;
    m_OneMarkerRaw.resize(m_bytesPerVariant);

    // Seek to start of variant data
    fseek(m_fin, m_dataOffset, SEEK_SET);
}

// ============================================================
// readPvarFile: parse the .pvar file
// Supports both:
//   - Header with #CHROM line (standard .pvar)
//   - No header, 6-column .bim-like format (fallback)
// ============================================================
void PgenClass::readPvarFile()
{
    std::cout << "Reading pvar file...." << std::endl;
    std::ifstream pvar(m_pvarFile);
    if (!pvar.is_open()) {
        throw std::runtime_error("Cannot open .pvar file: " + m_pvarFile);
    }

    std::string line;
    int col_chr = -1, col_pos = -1, col_id = -1, col_ref = -1, col_alt = -1;
    bool hasHeader = false;

    while (std::getline(pvar, line)) {
        if (line.empty()) continue;

        // Remove trailing \r
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Parse header line
        if (line[0] == '#') {
            if (line.rfind("#CHROM", 0) == 0) {
                hasHeader = true;
                // Parse column positions from header
                std::istringstream iss(line);
                std::string token;
                int colIdx = 0;
                while (iss >> token) {
                    if (token == "#CHROM" || token == "CHROM") col_chr = colIdx;
                    else if (token == "POS") col_pos = colIdx;
                    else if (token == "ID") col_id = colIdx;
                    else if (token == "REF") col_ref = colIdx;
                    else if (token == "ALT") col_alt = colIdx;
                    colIdx++;
                }
                if (col_pos == -1 || col_id == -1 || col_ref == -1 || col_alt == -1) {
                    throw std::runtime_error(
                        ".pvar file header missing required columns (POS/ID/REF/ALT)");
                }
            }
            continue;  // skip all comment/header lines
        }

        // If no header found, assume .bim-like 6-column format:
        // CHR ID CM POS ALT REF
        if (!hasHeader) {
            col_chr = 0;
            col_id = 1;
            // col 2 = CM (ignored)
            col_pos = 3;
            col_alt = 4;
            col_ref = 5;
            hasHeader = true;  // only set defaults once
        }

        // Parse data line
        std::vector<std::string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }

        int maxCol = std::max({col_chr, col_pos, col_id, col_ref, col_alt});
        if ((int)tokens.size() <= maxCol) {
            throw std::runtime_error(
                ".pvar file has too few columns at variant " +
                std::to_string(m_chr.size() + 1));
        }

        m_chr.push_back(tokens[col_chr]);
        m_position.push_back(std::stoul(tokens[col_pos]));
        m_variantId.push_back(tokens[col_id]);

        // Convert alleles to uppercase (matching SAIGE behavior)
        std::string refAllele = tokens[col_ref];
        std::string altAllele = tokens[col_alt];
        std::transform(refAllele.begin(), refAllele.end(), refAllele.begin(), ::toupper);
        std::transform(altAllele.begin(), altAllele.end(), altAllele.begin(), ::toupper);

        m_ref.push_back(refAllele);
        m_alt.push_back(altAllele);
    }

    m_M = m_chr.size();
    std::cout << "Number of markers in pvar file: " << m_M << std::endl;
}

// ============================================================
// readPsamFile: parse the .psam file
// Supports both:
//   - Header with #FID or #IID line (standard .psam)
//   - No header, .fam-like format (uses first 2 columns as FID IID)
// ============================================================
void PgenClass::readPsamFile()
{
    std::cout << "Reading psam file (" << m_psamFile << ")...." << std::endl;
    std::ifstream psam(m_psamFile);
    if (!psam.is_open()) {
        throw std::runtime_error("Cannot open .psam file: " + m_psamFile);
    }

    std::string line;
    int col_iid = -1;
    bool hasHeader = false;

    while (std::getline(psam, line)) {
        if (line.empty()) continue;

        // Remove trailing \r
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        // Parse header line
        if (line[0] == '#') {
            // Look for #FID or #IID in header
            std::istringstream iss(line);
            std::string token;
            int colIdx = 0;
            bool hasFID = false;
            while (iss >> token) {
                if (token == "#IID" || token == "IID") col_iid = colIdx;
                if (token == "#FID" || token == "FID") hasFID = true;
                colIdx++;
            }
            // If #FID is present but #IID not found as separate column,
            // IID is typically the second column
            if (col_iid == -1 && hasFID) {
                col_iid = 1;
            }
            hasHeader = true;
            continue;
        }

        // If no header found, assume .fam-like format (FID IID ...)
        if (!hasHeader) {
            col_iid = 1;  // second column = IID
            // But if file has only IID column (plink2 --make-psam with --no-fid)
            // we check if there's only one column
            std::istringstream testIss(line);
            std::string testTok;
            int ncols = 0;
            while (testIss >> testTok) ncols++;
            if (ncols == 1) {
                col_iid = 0;
            }
            hasHeader = true;  // only set defaults once
        }

        // Parse data line
        std::istringstream iss(line);
        std::string token;
        int colIdx = 0;
        std::string iid;
        while (iss >> token) {
            if (colIdx == col_iid) {
                iid = token;
                break;
            }
            colIdx++;
        }

        if (!iid.empty()) {
            m_SampleInPgen.push_back(iid);
        }
    }

    std::cout << "Number of samples in psam file: " << m_SampleInPgen.size() << std::endl;
}

// ============================================================
// setPosSampleInPgen: find positions of model samples in PGEN file
// ============================================================
void PgenClass::setPosSampleInPgen(std::vector<std::string>& t_SampleInModel)
{
    std::cout << "Setting position of samples in PGEN file...." << std::endl;

    // If no sample IDs provided, use all .psam samples in order
    if (t_SampleInModel.empty()) {
        m_N = m_SampleInPgen.size();
        std::cout << "  No sampleIDs in null model; using all " << m_N
                  << " samples from .psam file." << std::endl;
        m_posSampleInPgen.resize(m_N);
        for (uint32_t i = 0; i < m_N; i++) {
            m_posSampleInPgen[i] = i;
        }
        std::cout << "Number of samples in analysis: " << m_N << std::endl;
        return;
    }

    m_N = t_SampleInModel.size();

    // Build map of .psam sample IID -> index (0-based)
    std::unordered_map<std::string, uint32_t> pgenSampleMap;
    for (uint32_t i = 0; i < m_SampleInPgen.size(); i++) {
        pgenSampleMap[m_SampleInPgen[i]] = i;
    }

    m_posSampleInPgen.resize(m_N);
    for (uint32_t i = 0; i < m_N; i++) {
        auto it = pgenSampleMap.find(t_SampleInModel[i]);
        if (it == pgenSampleMap.end()) {
            throw std::runtime_error(
                "At least one subject requested is not in .psam file. "
                "Sample not found: " + t_SampleInModel[i]);
        }
        m_posSampleInPgen[i] = it->second;
    }

    std::cout << "Number of samples in analysis: " << m_N << std::endl;
}

// ============================================================
// getOneMarker: read one marker from .pgen and decode genotypes
//
// For PGEN mode 0x02, the genotype encoding is:
//   00 = hom ref (genotype 0)
//   01 = het     (genotype 1)
//   10 = hom alt (genotype 2)
//   11 = missing
//
// Note: Unlike PLINK .bed (which encodes 00=hom_alt, 01=missing,
//   10=het, 11=hom_ref), PGEN mode 0x02 uses a cleaner encoding.
//   The SAIGE PGEN reader in PGEN.cpp calls Read() which calls
//   PgrGet1D -> Dosage16ToDoubles -> kGenoRDoublePairs which maps
//   to {0.0, 1.0, 2.0, NA_REAL}, i.e. 00->0, 01->1, 10->2, 11->NA.
// ============================================================
void PgenClass::getOneMarker(
    uint64_t& t_gIndex,
    std::string& t_ref,
    std::string& t_alt,
    std::string& t_marker,
    uint32_t& t_pd,
    std::string& t_chr,
    double& t_altFreq,
    double& t_altCounts,
    double& t_missingRate,
    double& t_imputeInfo,
    bool& t_isOutputIndexForMissing,
    std::vector<uint>& t_indexForMissing,
    bool& t_isOnlyOutputNonZero,
    std::vector<uint>& t_indexForNonZero,
    arma::vec& OneMarkerG1)
{
    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    // Variant metadata from .pvar
    t_chr = m_chr[t_gIndex];
    t_pd = m_position[t_gIndex];
    t_ref = m_ref[t_gIndex];
    t_alt = m_alt[t_gIndex];
    t_marker = m_variantId[t_gIndex];

    // Seek to the variant record in .pgen
    uint64_t filePos = m_dataOffset + m_bytesPerVariant * t_gIndex;
    fseek(m_fin, filePos, SEEK_SET);

    // Read raw genotype bytes
    if (fread(m_OneMarkerRaw.data(), 1, m_bytesPerVariant, m_fin) != m_bytesPerVariant) {
        throw std::runtime_error(
            "Failed to read variant " + std::to_string(t_gIndex) + " from PGEN file");
    }

    // Decode genotypes for each sample in the analysis
    // PGEN mode 0x02 encoding: 00=0(ref), 01=1(het), 10=2(alt), 11=missing
    uint32_t numMissing = 0;
    t_altCounts = 0;
    uint j = 0;

    for (uint32_t i = 0; i < m_N; i++) {
        uint32_t sampleIdx = m_posSampleInPgen[i];
        uint8_t geno = getGenotype(sampleIdx);

        double genoVal;
        if (geno == 0x03) {
            // Missing
            genoVal = std::numeric_limits<double>::quiet_NaN();
            numMissing++;
            if (t_isOutputIndexForMissing) {
                t_indexForMissing.push_back(i);
            }
        } else {
            // 00->0, 01->1, 10->2
            genoVal = (double)geno;
            t_altCounts += genoVal;
        }

        if (geno > 0 && geno != 0x03) {
            t_indexForNonZero.push_back(i);
        }

        if (t_isOnlyOutputNonZero) {
            if (geno > 0 && geno != 0x03) {
                OneMarkerG1[j] = genoVal;
                j++;
            }
        } else {
            OneMarkerG1[i] = genoVal;
        }
    }

    // Compute frequency and missing rate
    uint32_t count = m_N - numMissing;
    t_missingRate = (double)numMissing / (double)m_N;
    t_imputeInfo = 1.0;

    if (count > 0) {
        t_altFreq = t_altCounts / (double)count / 2.0;
    } else {
        t_altFreq = 0;
    }

    if (t_isOnlyOutputNonZero) {
        OneMarkerG1.resize(j);
    }
}

// ============================================================
// closegenofile: close the .pgen file handle
// ============================================================
void PgenClass::closegenofile()
{
    if (m_fin) {
        fclose(m_fin);
        m_fin = nullptr;
    }
}

// ============================================================
// getOneMarker_ts: thread-safe variant of getOneMarker.
//
// Each OMP thread keeps its own FILE* + scratch buffer in thread_local
// storage. On first use the thread opens its own fd to the same .pgen
// file. PGEN mode 0x01/0x02 records are fixed-size (m_bytesPerVariant
// bytes per variant after m_dataOffset), so a thread can seek directly
// by variant index without coordinating with other threads.
//
// The thread_local FILE* is leaked at thread exit (cleaned up on
// process exit by the OS). Acceptable for SAIGE step 2's lifecycle.
// ============================================================
void PgenClass::getOneMarker_ts(uint64_t t_gIndex,
                                std::string& t_ref,
                                std::string& t_alt,
                                std::string& t_marker,
                                uint32_t& t_pd,
                                std::string& t_chr,
                                double& t_altFreq,
                                double& t_altCounts,
                                double& t_missingRate,
                                double& t_imputeInfo,
                                bool t_isOutputIndexForMissing,
                                std::vector<uint>& t_indexForMissing,
                                bool t_isOnlyOutputNonZero,
                                std::vector<uint>& t_indexForNonZero,
                                arma::vec& OneMarkerG1)
{
    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    t_chr    = m_chr[t_gIndex];
    t_pd     = m_position[t_gIndex];
    t_ref    = m_ref[t_gIndex];
    t_alt    = m_alt[t_gIndex];
    t_marker = m_variantId[t_gIndex];

    thread_local FILE* tlsFin = nullptr;
    thread_local std::vector<unsigned char> tlsOneMarkerRaw;

    if (tlsFin == nullptr) {
        tlsFin = fopen(m_pgenFile.c_str(), "rb");
        if (!tlsFin) {
            throw std::runtime_error(
                "PgenClass::getOneMarker_ts: cannot open .pgen file: " + m_pgenFile);
        }
    }
    if (tlsOneMarkerRaw.size() < m_bytesPerVariant) {
        tlsOneMarkerRaw.resize(m_bytesPerVariant);
    }

    uint64_t filePos = m_dataOffset + m_bytesPerVariant * t_gIndex;
    fseek(tlsFin, filePos, SEEK_SET);

    if (fread(tlsOneMarkerRaw.data(), 1, m_bytesPerVariant, tlsFin) != m_bytesPerVariant) {
        throw std::runtime_error(
            "PgenClass::getOneMarker_ts: failed to read variant " +
            std::to_string(t_gIndex) + " from PGEN file");
    }

    // Decode genotypes (PGEN mode 0x02 encoding: 00=0, 01=1, 10=2, 11=missing)
    // We inline the 2-bit extraction because getGenotype() reads m_OneMarkerRaw.
    uint32_t numMissing = 0;
    t_altCounts = 0;
    uint j = 0;

    for (uint32_t i = 0; i < m_N; i++) {
        uint32_t sampleIdx = m_posSampleInPgen[i];
        uint32_t byteIdx  = sampleIdx / 4;
        uint32_t bitShift = (sampleIdx % 4) * 2;
        uint8_t  geno     = (tlsOneMarkerRaw[byteIdx] >> bitShift) & 0x03;

        double genoVal;
        if (geno == 0x03) {
            genoVal = std::numeric_limits<double>::quiet_NaN();
            numMissing++;
            if (t_isOutputIndexForMissing) {
                t_indexForMissing.push_back(i);
            }
        } else {
            genoVal = (double)geno;
            t_altCounts += genoVal;
        }

        if (geno > 0 && geno != 0x03) {
            t_indexForNonZero.push_back(i);
        }

        if (t_isOnlyOutputNonZero) {
            if (geno > 0 && geno != 0x03) {
                OneMarkerG1[j] = genoVal;
                j++;
            }
        } else {
            OneMarkerG1[i] = genoVal;
        }
    }

    uint32_t count = m_N - numMissing;
    t_missingRate = (double)numMissing / (double)m_N;
    t_imputeInfo  = 1.0;

    if (count > 0) {
        t_altFreq = t_altCounts / (double)count / 2.0;
    } else {
        t_altFreq = 0;
    }

    if (t_isOnlyOutputNonZero) {
        OneMarkerG1.resize(j);
    }
}

} // namespace PGEN


// ============================================================
// whichCPP: C++ version of which(). Note: start from 0, not 1
// (Outside PLINK namespace, as declared in header)
// ============================================================
std::vector<unsigned int> whichCPP(std::vector<std::string>& strVec,
                                   std::string strValue)
{
    std::vector<unsigned int> indexVec;
    for (unsigned int i = 0; i < strVec.size(); i++) {
        if (strVec.at(i) == strValue)
            indexVec.push_back(i);
    }
    return indexVec;
}


// ============================================================
// setPLINKobjInCPP: create and configure the global PLINK object
// Mirrors SAIGE Main.cpp::setPLINKobjInCPP
// ============================================================
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string>& t_SampleInModel,
                      std::string t_AlleleOrder)
{
    // Clean up any existing object
    if (ptr_gPLINKobj != nullptr) {
        ptr_gPLINKobj->closegenofile();
        delete ptr_gPLINKobj;
        ptr_gPLINKobj = nullptr;
    }

    ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                          t_famFile,
                                          t_bedFile,
                                          t_AlleleOrder);

    ptr_gPLINKobj->setPosSampleInPlink(t_SampleInModel);

    int n = ptr_gPLINKobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// setVCFobjInCPP: create and configure the global VCF object
// Mirrors SAIGE Main.cpp::setVCFobjInCPP
// ============================================================
void setVCFobjInCPP(const std::string& t_vcfFileName,
                     const std::string& t_vcfField,
                     std::vector<std::string>& t_SampleInModel)
{
    // Clean up any existing object
    if (ptr_gVCFobj != nullptr) {
        ptr_gVCFobj->closegenofile();
        delete ptr_gVCFobj;
        ptr_gVCFobj = nullptr;
    }

    ptr_gVCFobj = new VCF::VcfClass(t_vcfFileName, t_vcfField, t_SampleInModel);

    int n = ptr_gVCFobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// setBGENobjInCPP: create and configure the global BGEN object
// Mirrors SAIGE Main.cpp::setBGENobjInCPP
// ============================================================
void setBGENobjInCPP(const std::string& t_bgenFileName,
                      const std::vector<std::string>& t_SampleInBgen,
                      std::vector<std::string>& t_SampleInModel,
                      const std::string& t_AlleleOrder)
{
    // Clean up any existing object
    if (ptr_gBGENobj != nullptr) {
        ptr_gBGENobj->closegenofile();
        delete ptr_gBGENobj;
        ptr_gBGENobj = nullptr;
    }

    std::cout << "t_SampleInBgen " << t_SampleInBgen.size() << std::endl;
    ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                        t_SampleInBgen,
                                        t_SampleInModel,
                                        t_AlleleOrder);

    int n = ptr_gBGENobj->getN();
    std::cout << "n:\t" << n << std::endl;

    // Pre-populate marker metadata from the BGEN .bgi index if present.
    // Without this, getMarkerIDToIndex() returns an empty map until the
    // marker stream has been walked once, which breaks SAIGE-GENE+ region
    // tests (group-file variant lookup happens before any streaming).
    std::string bgiPath = t_bgenFileName + ".bgi";
    struct stat sb;
    if (stat(bgiPath.c_str(), &sb) == 0) {
        try {
            uint32_t nVar = ptr_gBGENobj->populateFromBgi(bgiPath);
            std::cout << "  Pre-loaded " << nVar
                      << " variant metadata rows from " << bgiPath << std::endl;
        } catch (const std::exception& e) {
            std::cout << "  Warning: .bgi pre-load failed (" << e.what()
                      << "); falling back to stream-on-demand metadata." << std::endl;
        }
    } else {
        std::cout << "  .bgi not found at " << bgiPath
                  << " — marker-ID lookups won't work until stream completes." << std::endl;
    }
}


// ============================================================
// setPGENobjInCPP: create and configure the global PGEN object
// Mirrors SAIGE Main.cpp::setPGENobjInCPP
// ============================================================
void setPGENobjInCPP(const std::string& t_pgenFile,
                      const std::string& t_psamFile,
                      const std::string& t_pvarFile,
                      std::vector<std::string>& t_SampleInModel)
{
    // Clean up any existing object
    if (ptr_gPGENobj != nullptr) {
        ptr_gPGENobj->closegenofile();
        delete ptr_gPGENobj;
        ptr_gPGENobj = nullptr;
    }

    ptr_gPGENobj = new PGEN::PgenClass(t_pgenFile, t_psamFile, t_pvarFile, t_SampleInModel);

    int n = ptr_gPGENobj->getN();
    std::cout << "n:\t" << n << std::endl;
}


// ============================================================
// Unified_getOneMarker: unified dispatcher for genotype reading
//
// Ported from SAIGE Main.cpp. Supports "plink", "vcf", "bgen", and "pgen".
// ============================================================
bool Unified_getOneMarker(std::string& t_genoType,
                          uint64_t& t_gIndex_prev,
                          uint64_t& t_gIndex,
                          std::string& t_ref,
                          std::string& t_alt,
                          std::string& t_marker,
                          uint32_t& t_pd,
                          std::string& t_chr,
                          double& t_altFreq,
                          double& t_altCounts,
                          double& t_missingRate,
                          double& t_imputeInfo,
                          bool& t_isOutputIndexForMissing,
                          std::vector<uint>& t_indexForMissing,
                          bool& t_isOnlyOutputNonZero,
                          std::vector<uint>& t_indexForNonZero,
                          arma::vec& t_GVec,
                          bool t_isImputation)
{
    bool isBoolRead = true;

    if (t_genoType == "plink") {
        bool isTrueGenotype = true;
        // t_gIndex_prev is after reading the last marker

        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: PLINK object not initialized. "
                "Call setPLINKobjInCPP first.");
        }

        ptr_gPLINKobj->getOneMarker(t_gIndex_prev, t_gIndex,
                                    t_ref, t_alt, t_marker, t_pd, t_chr,
                                    t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                    t_isOutputIndexForMissing, t_indexForMissing,
                                    t_isOnlyOutputNonZero, t_indexForNonZero,
                                    isTrueGenotype, t_GVec);
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: VCF object not initialized. "
                "Call setVCFobjInCPP first.");
        }

        isBoolRead = ptr_gVCFobj->getOneMarker(
            t_ref, t_alt, t_marker, t_pd, t_chr,
            t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
            t_isOutputIndexForMissing, t_indexForMissing,
            t_isOnlyOutputNonZero, t_indexForNonZero,
            t_GVec, t_isImputation);
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: BGEN object not initialized. "
                "Call setBGENobjInCPP first.");
        }

        ptr_gBGENobj->getOneMarker(t_gIndex_prev, t_gIndex,
                                    t_ref, t_alt, t_marker, t_pd, t_chr,
                                    t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                    t_isOutputIndexForMissing, t_indexForMissing,
                                    t_isOnlyOutputNonZero, t_indexForNonZero,
                                    isBoolRead, t_GVec, t_isImputation);
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker: PGEN object not initialized. "
                "Call setPGENobjInCPP first.");
        }

        ptr_gPGENobj->getOneMarker(t_gIndex,
                                    t_ref, t_alt, t_marker, t_pd, t_chr,
                                    t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
                                    t_isOutputIndexForMissing, t_indexForMissing,
                                    t_isOnlyOutputNonZero, t_indexForNonZero,
                                    t_GVec);
    } else {
        throw std::runtime_error(
            "Unified_getOneMarker: Unknown genotype type '" + t_genoType + "'. "
            "Supported types: plink, vcf, bgen");
    }

    return isBoolRead;
}


// ============================================================
// Unified_getOneMarker_ts: thread-safe unified dispatcher.
//
// For PLINK / PGEN: dispatches to per-format _ts variants that use
//   thread_local FILE* + scratch buffers. No critical section required.
// For VCF: falls back to the locked getOneMarker — htslib bcf_read is
//   sequential and not thread-safe in our usage. The CALLER is responsible
//   for wrapping the VCF call in a critical section.
// For BGEN: should not be called via this dispatcher — the BgenStreamer
//   already provides thread-safe access. Throws if invoked.
// ============================================================
bool Unified_getOneMarker_ts(std::string& t_genoType,
                             uint64_t t_gIndex,
                             std::string& t_ref,
                             std::string& t_alt,
                             std::string& t_marker,
                             uint32_t& t_pd,
                             std::string& t_chr,
                             double& t_altFreq,
                             double& t_altCounts,
                             double& t_missingRate,
                             double& t_imputeInfo,
                             bool t_isOutputIndexForMissing,
                             std::vector<uint>& t_indexForMissing,
                             bool t_isOnlyOutputNonZero,
                             std::vector<uint>& t_indexForNonZero,
                             arma::vec& t_GVec,
                             bool t_isImputation)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker_ts: PLINK object not initialized.");
        }
        bool isTrueGenotype = true;
        ptr_gPLINKobj->getOneMarker_ts(
            t_gIndex,
            t_ref, t_alt, t_marker, t_pd, t_chr,
            t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
            t_isOutputIndexForMissing, t_indexForMissing,
            t_isOnlyOutputNonZero,    t_indexForNonZero,
            isTrueGenotype, t_GVec);
        return true;
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker_ts: PGEN object not initialized.");
        }
        ptr_gPGENobj->getOneMarker_ts(
            t_gIndex,
            t_ref, t_alt, t_marker, t_pd, t_chr,
            t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
            t_isOutputIndexForMissing, t_indexForMissing,
            t_isOnlyOutputNonZero,    t_indexForNonZero,
            t_GVec);
        return true;
    } else if (t_genoType == "vcf") {
        // VCF cannot be parallelized here — caller must serialize.
        // Reuse the (locked-by-caller) getOneMarker.
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getOneMarker_ts: VCF object not initialized.");
        }
        return ptr_gVCFobj->getOneMarker(
            t_ref, t_alt, t_marker, t_pd, t_chr,
            t_altFreq, t_altCounts, t_missingRate, t_imputeInfo,
            t_isOutputIndexForMissing, t_indexForMissing,
            t_isOnlyOutputNonZero,    t_indexForNonZero,
            t_GVec, t_isImputation);
    } else if (t_genoType == "bgen") {
        throw std::runtime_error(
            "Unified_getOneMarker_ts: BGEN should be read via BgenStreamer, "
            "not this dispatcher.");
    } else {
        throw std::runtime_error(
            "Unified_getOneMarker_ts: Unknown genotype type '" + t_genoType + "'.");
    }
}


// ============================================================
// Unified_getMarkerCount: get total marker count from genotype file
// ============================================================
uint32_t Unified_getMarkerCount(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: PLINK object not initialized.");
        }
        return ptr_gPLINKobj->getM();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: VCF object not initialized.");
        }
        return ptr_gVCFobj->getM0();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: BGEN object not initialized.");
        }
        return ptr_gBGENobj->getM0();
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerCount: PGEN object not initialized.");
        }
        return ptr_gPGENobj->getM();
    } else {
        throw std::runtime_error(
            "Unified_getMarkerCount: Unsupported type: " + t_genoType);
    }
}


// ============================================================
// Unified_getSampleSizeinGeno: get total sample size in genotype file
// ============================================================
uint32_t Unified_getSampleSizeinGeno(std::string& t_genoType)
{
    uint32_t N0;
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: PLINK object not initialized.");
        }
        N0 = ptr_gPLINKobj->getN0();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: VCF object not initialized.");
        }
        N0 = ptr_gVCFobj->getN0();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: BGEN object not initialized.");
        }
        N0 = ptr_gBGENobj->getN0();
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinGeno: PGEN object not initialized.");
        }
        N0 = ptr_gPGENobj->getN0();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinGeno: Unsupported type: " + t_genoType);
    }
    return N0;
}


// ============================================================
// Unified_getSampleSizeinAnalysis: get sample size in analysis (after subsetting)
// ============================================================
uint32_t Unified_getSampleSizeinAnalysis(std::string& t_genoType)
{
    uint32_t N;
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: PLINK object not initialized.");
        }
        N = ptr_gPLINKobj->getN();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: VCF object not initialized.");
        }
        N = ptr_gVCFobj->getN();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: BGEN object not initialized.");
        }
        N = ptr_gBGENobj->getN();
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error(
                "Unified_getSampleSizeinAnalysis: PGEN object not initialized.");
        }
        N = ptr_gPGENobj->getN();
    } else {
        throw std::runtime_error(
            "Unified_getSampleSizeinAnalysis: Unsupported type: " + t_genoType);
    }
    return N;
}


// ============================================================
// Unified_getMarkerIDToIndex: get marker ID to index map
// ============================================================
std::unordered_map<std::string, uint32_t> Unified_getMarkerIDToIndex(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: PLINK object not initialized.");
        }
        return ptr_gPLINKobj->getMarkerIDToIndex();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: VCF object not initialized.");
        }
        return ptr_gVCFobj->getMarkerIDToIndex();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: BGEN object not initialized.");
        }
        return ptr_gBGENobj->getMarkerIDToIndex();
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToIndex: PGEN object not initialized.");
        }
        return ptr_gPGENobj->getMarkerIDToIndex();
    } else {
        throw std::runtime_error(
            "Unified_getMarkerIDToIndex: Unsupported type: " + t_genoType);
    }
}


// ============================================================
// Unified_getMarkerIDToChrom: get marker ID -> CHROM map
// ============================================================
std::unordered_map<std::string, std::string> Unified_getMarkerIDToChrom(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToChrom: PLINK object not initialized.");
        }
        return ptr_gPLINKobj->getMarkerIDToChrom();
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToChrom: VCF object not initialized.");
        }
        return ptr_gVCFobj->getMarkerIDToChrom();
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToChrom: BGEN object not initialized.");
        }
        return ptr_gBGENobj->getMarkerIDToChrom();
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj == nullptr) {
            throw std::runtime_error("Unified_getMarkerIDToChrom: PGEN object not initialized.");
        }
        return ptr_gPGENobj->getMarkerIDToChrom();
    } else {
        throw std::runtime_error(
            "Unified_getMarkerIDToChrom: Unsupported type: " + t_genoType);
    }
}


// ============================================================
// closeGenoFile: close the genotype file
// ============================================================
void closeGenoFile(std::string& t_genoType)
{
    if (t_genoType == "plink") {
        if (ptr_gPLINKobj != nullptr) {
            ptr_gPLINKobj->closegenofile();
        }
    } else if (t_genoType == "vcf") {
        if (ptr_gVCFobj != nullptr) {
            ptr_gVCFobj->closegenofile();
        }
    } else if (t_genoType == "bgen") {
        if (ptr_gBGENobj != nullptr) {
            ptr_gBGENobj->closegenofile();
        }
    } else if (t_genoType == "pgen") {
        if (ptr_gPGENobj != nullptr) {
            ptr_gPGENobj->closegenofile();
        }
    } else {
        throw std::runtime_error(
            "closeGenoFile: Unsupported type: " + t_genoType);
    }
}
