// Standalone port of SAIGE/src/PLINK.hpp, SAIGE/src/VCF.hpp, and SAIGE/src/BGEN.hpp
// PLINK + VCF + BGEN genotype readers -- all Rcpp dependencies removed
// Ported from: /SAIGE/src/PLINK.hpp, /SAIGE/src/PLINK.cpp,
//              /SAIGE/src/VCF.hpp, /SAIGE/src/VCF.cpp,
//              /SAIGE/src/BGEN.hpp, /SAIGE/src/BGEN.cpp

#ifndef GENOTYPE_READER_HPP
#define GENOTYPE_READER_HPP

#include <armadillo>
#include <string>
#include <vector>
#include <cstdio>
#include <unordered_map>

// Forward-declare htslib types to avoid including htslib headers in the header
struct htsFile;
struct bcf_hdr_t;
struct bcf1_t;

namespace PLINK {

class PlinkClass {
private:

    // added on 03/14/2021
    std::string m_AlleleOrder;  // "alt-first" or "ref-first"

    // information from bim file
    uint32_t m_M0, m_M;  // total markers, markers in analysis
    std::vector<std::string> m_chr;               // Chromosome code
    std::vector<std::string> m_MarkerInPlink;     // Variant identifier
    std::vector<float> m_gd;                      // Position in morgans or centimorgans
    std::vector<uint32_t> m_pd;                   // Base-pair coordinate (1-based)
    std::vector<std::string> m_alt;               // Allele 1 (clear bits in .bed; usually minor)
    std::vector<std::string> m_ref;               // Allele 2 (set bits in .bed; usually major)

    // information from fam file
    std::vector<std::string> m_SampleInPlink;
    uint32_t m_N0, m_N;
    unsigned long long int m_numBytesofEachMarker0, m_numBytesofEachMarker;

    // input file stream of .bed file
    FILE* m_fin;

    // PLINK files
    std::string m_bimFile, m_famFile, m_bedFile;
    std::vector<uint32_t> m_posSampleInPlink;

    // https://www.cog-genomics.org/plink/1.9/formats#bed
    // PLINK format
    // The two-bit genotype codes have the following meanings:
    // 00  Homozygous for first allele in .bim file
    // 01  Missing genotype
    // 10  Heterozygous
    // 11  Homozygous for second allele in .bim file
    static const unsigned char HOM_REF = 0x3;   // 0b11
    static const unsigned char HET = 0x2;       // 0b10
    static const unsigned char HOM_ALT = 0x0;   // 0b00
    static const unsigned char MISSING = 0x1;   // 0b01

    // Genotype mapping vectors
    std::vector<int8_t> m_genoMaps_alt_first = {2, -1, 1, 0};
    std::vector<int8_t> m_genoMaps_ref_first = {0, -1, 1, 2};

    // pipeline: OneMarkerG4 --> bufferG4 --> bufferG1 --> OneMarkerG1
    std::vector<unsigned char> m_OneMarkerG4;

    void readBimFile();
    void readFamFile();

    // extract geno (0,1,2,3) at specific pos (0,1,2,3) of address c (1 byte)
    inline void getGenotype(const unsigned char c, const uint32_t pos, size_t& geno) {
        geno = (c >> (pos << 1)) & 0x3;  // 0b11 = 0x3
    }

public:

    PlinkClass(std::string t_bimFile,
               std::string t_famFile,
               std::string t_bedFile,
               std::string t_AlleleOrder);

    // setup PlinkClass
    void setPlinkobj(std::string t_bimFile,
                     std::string t_famFile,
                     std::string t_bedFile);

    void setPosSampleInPlink(std::vector<std::string>& t_SampleInModel);

    void getOneMarker(uint64_t& t_gIndex_prev,
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
                      arma::vec& OneMarkerG1);

    // Convenience overload: simplified getOneMarker (like SAIGE PLINK.hpp inline overloads)
    void getOneMarker(uint64_t t_gIndex_prev,
                      uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::string& t_chr,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForMissing, indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        bool isTrueGenotype = false;
        getOneMarker(t_gIndex_prev, t_gIndex, ref, alt, marker, pd, t_chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero,
                     isTrueGenotype, OneMarkerG1);
    }

    // Convenience overload: with indexForMissing output
    void getOneMarker(uint64_t t_gIndex_prev,
                      uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::vector<uint>& t_indexForMissing,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker, chr;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        bool isTrueGenotype = true;
        getOneMarker(t_gIndex_prev, t_gIndex, ref, alt, marker, pd, chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, t_indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero,
                     isTrueGenotype, OneMarkerG1);
    }

    uint32_t getN0() { return m_N0; }
    uint32_t getN() { return m_N; }
    uint32_t getM0() { return m_M0; }
    uint32_t getM() { return m_M; }
    uint32_t getnumBytesofEachMarker0() { return m_numBytesofEachMarker0; }
    uint32_t getnumBytesofEachMarker() { return m_numBytesofEachMarker; }

    std::vector<std::string> getChrVec() { return m_chr; }

    // Build a map from "chr:pos:ref:alt" (and "chr:pos:alt:ref") to marker index
    // Used by group file parser to look up variant genotype indices
    std::unordered_map<std::string, uint32_t> getMarkerIDToIndex() {
        std::unordered_map<std::string, uint32_t> idMap;
        for (uint32_t i = 0; i < m_M0; i++) {
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

    // Build a map from rsID (column 2 of .bim, m_MarkerInPlink) to marker index
    // Used by conditional analysis to look up conditioning markers by name
    std::unordered_map<std::string, uint32_t> getMarkerNameToIndex() {
        std::unordered_map<std::string, uint32_t> nameMap;
        for (uint32_t i = 0; i < m_M0; i++) {
            nameMap[m_MarkerInPlink[i]] = i;
        }
        return nameMap;
    }

    void closegenofile();
};

} // namespace PLINK


// ============================================================
// VCF namespace: VCF/BCF/VCF.GZ reader using htslib
// Ported from SAIGE/src/VCF.hpp and VCF.cpp
// ============================================================
namespace VCF {

class VcfClass {
private:
    // htslib file handles
    htsFile*   m_htsFile;
    bcf_hdr_t* m_hdr;
    bcf1_t*    m_rec;

    // VCF file path
    std::string m_vcfFileName;

    // Format field to read: "GT", "DS", or "HDS"
    std::string m_fmtField;

    // Sample information
    std::vector<std::string> m_SampleInVcf;  // sample IDs from VCF header
    uint32_t m_N0;  // total samples in VCF
    uint32_t m_N;   // samples in analysis

    // Mapping from VCF sample index -> model sample index
    // m_posSampleInModel[vcf_idx] = model_idx, or -1 if not in model
    std::vector<int32_t> m_posSampleInModel;

    // Marker information (accumulated as markers are read)
    uint32_t m_M0;          // total markers read so far
    std::vector<std::string> m_chr;
    std::vector<uint32_t> m_pd;
    std::vector<std::string> m_ref;
    std::vector<std::string> m_alt;
    std::vector<std::string> m_MarkerInVcf;

    // Pre-scanned marker count (total markers in file)
    uint32_t m_totalMarkers;
    bool m_isPreScanned;

    // Internal: read sample IDs from VCF header
    void getSampleIDlist();

public:
    VcfClass(const std::string& t_vcfFileName,
             const std::string& t_vcfField,
             std::vector<std::string>& t_SampleInModel);

    ~VcfClass();

    // Set up sample position mapping (model samples -> VCF samples)
    void setPosSampleInVcf(std::vector<std::string>& t_SampleInModel);

    // Read one marker from VCF (sequential read, ignores gIndex for seeking)
    // Returns false if end of file reached
    bool getOneMarker(
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
        bool t_isImputation);

    uint32_t getN0() { return m_N0; }
    uint32_t getN()  { return m_N; }
    uint32_t getM0() { return m_totalMarkers; }

    // Pre-scan VCF to count total markers (needed for index generation)
    uint32_t prescanMarkerCount();

    // Reset file to beginning for re-reading
    void resetFile();

    // Build marker ID to index map (chr:pos:ref:alt -> index)
    // Note: requires a full prescan first
    std::unordered_map<std::string, uint32_t> getMarkerIDToIndex();

    // Build marker name to index map (ID field -> index)
    std::unordered_map<std::string, uint32_t> getMarkerNameToIndex();

    std::vector<std::string> getChrVec() { return m_chr; }

    void closegenofile();
};

} // namespace VCF


// ============================================================
// BGEN namespace: BGEN v1.2 reader (manual binary parsing)
// Ported from SAIGE/src/BGEN.hpp and BGEN.cpp
// Supports zstd and zlib decompression
// ============================================================
namespace BGEN {

static constexpr uint32_t COMPRESSION_ZLIB = 1;
static constexpr uint32_t COMPRESSION_ZSTD = 2;

class BgenClass {
private:
    // Allele order: "alt-first" or "ref-first"
    std::string m_AlleleOrder;

    // BGEN file path
    std::string m_bgenFileName;

    // File handle for binary reading
    FILE* m_fin;

    // Decompression buffers
    std::vector<unsigned char> m_buf;
    std::vector<unsigned char> m_zBuf;
    uint32_t m_zBufLens;
    uint32_t m_bufLens;
    uint32_t CompressedSNPBlocks;  // 1=zlib, 2=zstd

    // Sample mapping: m_posSampleInModel[bgen_idx] = model_idx, or -1 if not in model
    std::vector<int32_t> m_posSampleInModel;

    // Sample information
    std::vector<std::string> m_SampleInBgen;
    uint32_t m_N0;   // total samples in BGEN
    uint32_t m_N;    // samples in analysis
    uint32_t m_M0;   // total markers in BGEN header

    // Marker metadata (accumulated during sequential reading for index maps)
    std::vector<std::string> m_chr;
    std::vector<uint32_t> m_pd;
    std::vector<std::string> m_ref;
    std::vector<std::string> m_alt;
    std::vector<std::string> m_MarkerInBgen;
    // Byte offset (file_start_position) of each variant block, parallel to the
    // metadata vectors. Populated from the .bgi index in populateFromBgi.
    // getMarkerIDToIndex()/getMarkerNameToIndex() return THIS (the byte offset
    // that getOneMarker() seeks to), NOT the sequential row index.
    std::vector<uint64_t> m_byteOffset;

    // Pre-allocated allele buffers (to avoid repeated allocation)
    std::vector<char> allele0, allele1;

    // Internal: decompress and parse genotype probabilities
    void Parse2(unsigned char* buf, uint32_t bufLen,
                const unsigned char* zBuf, uint32_t zBufLen,
                std::string& snpName,
                arma::vec& dosages,
                double& AC, double& AF,
                std::vector<uint>& indexforMissing,
                double& info,
                std::vector<uint>& indexNonZero,
                bool isImputation);

public:
    BgenClass(const std::string& t_bgenFileName,
              const std::vector<std::string>& t_SampleInBgen,
              std::vector<std::string>& t_SampleInModel,
              const std::string& t_AlleleOrder);

    ~BgenClass();

    // Set up BGEN file: read header, validate
    void setBgenObj(const std::string& t_bgenFileName,
                    const std::vector<std::string>& t_SampleInBgen);

    // Pre-populate marker metadata vectors (m_chr/m_pd/m_ref/m_alt/m_MarkerInBgen)
    // from the BGEN .bgi SQLite index, so getMarkerIDToIndex() and other
    // by-ID lookups work BEFORE any markers have been streamed. Required for
    // SAIGE-GENE+ region tests where the group-file marker IDs must be resolved
    // up-front. Returns the number of variants loaded.
    uint32_t populateFromBgi(const std::string& t_bgiFileName);

    // Set up sample position mapping
    void setPosSampleInBgen(std::vector<std::string>& t_SampleInModel);

    // Read one marker from BGEN (sequential read using byte-position seeking)
    // t_gIndex / t_gIndex_prev are byte positions (for BGEN), not marker indices
    // Returns false via t_isBoolRead if EOF
    void getOneMarker(uint64_t& t_gIndex_prev,
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
                      bool t_isImputation);

    uint32_t getN0() { return m_N0; }
    uint32_t getN()  { return m_N; }
    uint32_t getM0() { return m_M0; }

    // Build marker ID to index map (chr:pos:ref:alt -> marker read order index)
    // Note: only contains markers that have been read so far
    std::unordered_map<std::string, uint32_t> getMarkerIDToIndex();

    // Build marker name to index map (RSID -> marker read order index)
    std::unordered_map<std::string, uint32_t> getMarkerNameToIndex();

    std::vector<std::string> getChrVec() { return m_chr; }

    void closegenofile();
};

} // namespace BGEN


// ============================================================
// PGEN namespace: PGEN v2 reader (mode 0x02 basic variant-major)
// Standalone implementation that reads .pgen + .pvar + .psam files
// Ported from SAIGE/src/PGEN.hpp and PGEN.cpp without pgenlib dependency
//
// Supports mode 0x02 (basic variant-major, hard-calls only):
//   00 = hom ref, 01 = het, 10 = hom alt, 11 = missing
// This is the most common output of `plink2 --make-pgen` for hard-call data.
// For dosage data (mode 0x03/0x04/0x10/0x11), the full pgenlib would be needed.
// ============================================================
namespace PGEN {

class PgenClass {
private:
    // PGEN file paths
    std::string m_pgenFile, m_pvarFile, m_psamFile;

    // File handle for .pgen binary
    FILE* m_fin;

    // PGEN header info
    uint32_t m_M0;       // total variants in .pgen header
    uint32_t m_N0;       // total samples in .pgen header
    uint8_t  m_mode;     // storage mode (must be 0x02 for our reader)

    // Offset to first variant record in .pgen file
    uint64_t m_dataOffset;

    // Bytes per variant record = ceil(N0 / 4)
    uint64_t m_bytesPerVariant;

    // Information from .pvar file
    std::vector<std::string> m_chr;
    std::vector<std::string> m_variantId;
    std::vector<uint32_t>    m_position;
    std::vector<std::string> m_ref;
    std::vector<std::string> m_alt;
    uint32_t m_M;   // marker count from pvar (should match m_M0)

    // Information from .psam file
    std::vector<std::string> m_SampleInPgen;   // IID list
    uint32_t m_N;   // samples in analysis (after subsetting)

    // Sample mapping: model sample index -> pgen sample index
    std::vector<uint32_t> m_posSampleInPgen;

    // Raw genotype buffer for one variant
    std::vector<unsigned char> m_OneMarkerRaw;

    // Internal parsers
    void readPvarFile();
    void readPsamFile();
    void readPgenHeader();

    // Extract 2-bit genotype at sample index from raw buffer
    inline uint8_t getGenotype(uint32_t sampleIdx) const {
        uint32_t byteIdx = sampleIdx / 4;
        uint32_t bitShift = (sampleIdx % 4) * 2;
        return (m_OneMarkerRaw[byteIdx] >> bitShift) & 0x03;
    }

public:
    PgenClass(const std::string& t_pgenFile,
              const std::string& t_psamFile,
              const std::string& t_pvarFile,
              std::vector<std::string>& t_SampleInModel);

    ~PgenClass();

    void setPosSampleInPgen(std::vector<std::string>& t_SampleInModel);

    void getOneMarker(
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
        arma::vec& OneMarkerG1);

    // Convenience overload: simplified getOneMarker
    void getOneMarker(uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::string& t_chr,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForMissing, indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        getOneMarker(t_gIndex, ref, alt, marker, pd, t_chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero, OneMarkerG1);
    }

    // Convenience overload: with indexForMissing output
    void getOneMarker(uint64_t t_gIndex,
                      double& t_altFreq,
                      double& t_missingRate,
                      std::vector<uint>& t_indexForMissing,
                      arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker, chr;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        getOneMarker(t_gIndex, ref, alt, marker, pd, chr,
                     t_altFreq, altCounts, t_missingRate, imputeInfo,
                     isOutputIndexForMissing, t_indexForMissing,
                     isOnlyOutputNonZero, indexForNonZero, OneMarkerG1);
    }

    uint32_t getN0() { return m_N0; }
    uint32_t getN()  { return m_N; }
    uint32_t getM()  { return m_M; }

    std::vector<std::string> getChrVec() { return m_chr; }

    // Build a map from "chr:pos:ref:alt" (and "chr:pos:alt:ref") to marker index
    std::unordered_map<std::string, uint32_t> getMarkerIDToIndex() {
        std::unordered_map<std::string, uint32_t> idMap;
        for (uint32_t i = 0; i < m_M; i++) {
            std::string id1 = m_chr[i] + ":" + std::to_string(m_position[i]) + ":"
                             + m_ref[i] + ":" + m_alt[i];
            idMap[id1] = i;
            std::string id2 = m_chr[i] + ":" + std::to_string(m_position[i]) + ":"
                             + m_alt[i] + ":" + m_ref[i];
            if (idMap.find(id2) == idMap.end()) {
                idMap[id2] = i;
            }
        }
        return idMap;
    }

    // Build a map from variant ID to marker index
    std::unordered_map<std::string, uint32_t> getMarkerNameToIndex() {
        std::unordered_map<std::string, uint32_t> nameMap;
        for (uint32_t i = 0; i < m_M; i++) {
            nameMap[m_variantId[i]] = i;
        }
        return nameMap;
    }

    void closegenofile();
};

} // namespace PGEN


// Global checkpoint flag for debugging
extern bool g_writeCheckpoints;
extern std::string g_checkpointDir;

// Global pointer to PLINK object (mirrors SAIGE's ptr_gPLINKobj in Main.cpp)
extern PLINK::PlinkClass* ptr_gPLINKobj;

// Global pointer to VCF object (mirrors SAIGE's ptr_gVCFobj in Main.cpp)
extern VCF::VcfClass* ptr_gVCFobj;

// Global pointer to BGEN object (mirrors SAIGE's ptr_gBGENobj in Main.cpp)
extern BGEN::BgenClass* ptr_gBGENobj;

// Global pointer to PGEN object (mirrors SAIGE's ptr_gPGENobj in Main.cpp)
extern PGEN::PgenClass* ptr_gPGENobj;

// Unified dispatcher (like SAIGE Main.cpp)
// Supports "plink", "vcf", "bgen", and "pgen" types.
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
                          bool t_isImputation);

// Helper: set up PLINK object (mirrors SAIGE Main.cpp::setPLINKobjInCPP)
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string>& t_SampleInModel,
                      std::string t_AlleleOrder);

// Helper: set up VCF object (mirrors SAIGE Main.cpp::setVCFobjInCPP)
void setVCFobjInCPP(const std::string& t_vcfFileName,
                     const std::string& t_vcfField,
                     std::vector<std::string>& t_SampleInModel);

// Helper: set up BGEN object (mirrors SAIGE Main.cpp::setBGENobjInCPP)
void setBGENobjInCPP(const std::string& t_bgenFileName,
                      const std::vector<std::string>& t_SampleInBgen,
                      std::vector<std::string>& t_SampleInModel,
                      const std::string& t_AlleleOrder);

// Helper: set up PGEN object (mirrors SAIGE Main.cpp::setPGENobjInCPP)
void setPGENobjInCPP(const std::string& t_pgenFile,
                      const std::string& t_psamFile,
                      const std::string& t_pvarFile,
                      std::vector<std::string>& t_SampleInModel);

// Helper: get total marker count from genotype file
uint32_t Unified_getMarkerCount(std::string& t_genoType);

// Helper: get sample size from genotype file
uint32_t Unified_getSampleSizeinGeno(std::string& t_genoType);

// Helper: get sample size in analysis (after subsetting)
uint32_t Unified_getSampleSizeinAnalysis(std::string& t_genoType);

// Helper: close genotype file
void closeGenoFile(std::string& t_genoType);

// Helper: get marker ID to index map (for region testing)
std::unordered_map<std::string, uint32_t> Unified_getMarkerIDToIndex(std::string& t_genoType);

// C++ version of which(). Note: start from 0, not 1
std::vector<unsigned int> whichCPP(std::vector<std::string>& strVec,
                                   std::string strValue);

#endif
