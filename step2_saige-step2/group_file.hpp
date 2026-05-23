// Standalone port of SAIGE group file parsing
// Ported from: SAIGE/R/SAIGE_SPATest_Region_Func.R (ReadGroupFile, checkGroupFile)
//              SAIGE/R/SAIGE_SPATest_Region.R (SAIGE.getRegionList_new)
//
// Conversions from original SAIGE R code:
//   1. R list / data.frame      -->  C++ structs + vectors
//   2. strsplit(..., "[\ \t]+") -->  std::istringstream (whitespace splitting)
//   3. stop(...)                -->  throw std::runtime_error(...)
//   4. R is.null / NA checks    -->  empty vector checks
//   5. data.table merge         -->  std::unordered_map lookup

#ifndef GROUP_FILE_HPP
#define GROUP_FILE_HPP

#include <armadillo>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

// Data for one region/gene after parsing and filtering
struct RegionData {
    std::string regionName;
    std::vector<std::string> variantIDs;    // chr:pos:ref:alt identifiers
    std::vector<std::string> annotations;   // per-variant annotation label
    std::vector<double> weights;            // per-variant weights (empty if not provided)
    arma::imat annoIndicatorMat;            // (num_variants x num_anno_groups) binary indicator
    std::vector<std::string> annoVec;       // annotation group names that have markers
    std::vector<std::string> genoIndex;     // genotype file index per variant (as string)
    std::vector<std::string> genoIndex_prev; // previous genotype index (as string)
};

// Summary information about a group file's format
struct GroupFileInfo {
    int nRegions;
    bool is_weight_included;
    int nline_per_gene; // 2 or 3
};

// Check group file format and count regions.
// Direct port from SAIGE R: checkGroupFile()
GroupFileInfo checkGroupFile(const std::string& groupFile);

// Read a chunk of regions from an already-open file stream.
// This is a C++ port of SAIGE.getRegionList_new() from SAIGE/R/SAIGE_SPATest_Region.R.
//
// Parameters:
//   gf                - open ifstream positioned at the start of the next chunk
//   nregions_to_read  - number of genes/regions to read in this chunk
//   nline_per_gene    - 2 (var+anno) or 3 (var+anno+weight)
//   annoVec           - annotation groups, e.g. {"lof", "lof;missense", "lof;missense;synonymous"}
//                        For "lof;missense", a variant with annotation "lof" OR "missense" matches.
//   markerIDToIndex   - maps "chr:pos:ref:alt" -> genotype file marker index
//
// Returns vector of RegionData, one per gene/region.
// Regions with no matching variants are returned with empty variantIDs.
std::vector<RegionData> readRegionChunk(
    std::ifstream& gf,
    int nregions_to_read,
    int nline_per_gene,
    const std::vector<std::string>& annoVec,
    const std::unordered_map<std::string, uint32_t>& markerIDToIndex);

#endif
