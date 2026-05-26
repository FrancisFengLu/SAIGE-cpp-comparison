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
//   6. R matrix indexing        -->  arma::imat with row/col access

#include "group_file.hpp"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <cmath>


// ============================================================
// Helper: split a line on whitespace (same pattern as genotype_reader.cpp)
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
// Helper: split a string by a delimiter character
// Used to split annotation groups like "lof;missense" by ";"
// ============================================================
static std::vector<std::string> splitByChar(const std::string& s, char delim) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string token;
    while (std::getline(iss, token, delim)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}


// ============================================================
// Helper: strip trailing \r from a line (for Windows line endings)
// ============================================================
static void stripCR(std::string& line) {
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
}


// ============================================================
// checkGroupFile
// Direct port from SAIGE/R/SAIGE_SPATest_Region_Func.R: checkGroupFile()
//
// Opens the group file, detects whether weights are included
// (3rd line's type field == "weight"), counts total regions,
// and validates the var/anno/weight structure.
// ============================================================
GroupFileInfo checkGroupFile(const std::string& groupFile) {
    std::cout << "Start extracting marker-level information from 'groupFile' of "
              << groupFile << " ...." << std::endl;

    // --- Pass 1: Detect if weight line is present ---
    std::ifstream gf(groupFile);
    if (!gf.is_open()) {
        throw std::runtime_error("Cannot open group file: " + groupFile);
    }

    // Read first 2 lines (first gene's var + anno)
    std::string line1, line2;
    if (!std::getline(gf, line1) || !std::getline(gf, line2)) {
        throw std::runtime_error("Group file has fewer than 2 lines: " + groupFile);
    }

    bool is_weight_included = false;
    int a = 2; // lines per gene

    // Read 3rd line to check if it's a weight line
    std::string line3;
    if (std::getline(gf, line3)) {
        stripCR(line3);
        if (!line3.empty()) {
            std::vector<std::string> tokens = splitLine(line3);
            if (tokens.size() < 3) {
                throw std::runtime_error(
                    "Error, group file line 3: Each line should have a region name "
                    "and a category (var, anno or weight)");
            }
            std::string type = tokens[1];
            if (type == "weight") {
                is_weight_included = true;
                std::cout << "weights are included for markers" << std::endl;
                a = 3;
            }
        }
    }
    gf.close();

    // --- Pass 2: Count regions and validate structure ---
    std::ifstream gf2(groupFile);
    if (!gf2.is_open()) {
        throw std::runtime_error("Cannot reopen group file: " + groupFile);
    }

    int nregion = 0;
    int lineNum = 0;

    while (true) {
        // Read 'a' lines at a time (one gene block)
        std::vector<std::string> lines(a);
        int linesRead = 0;
        for (int k = 0; k < a; k++) {
            if (std::getline(gf2, lines[k])) {
                stripCR(lines[k]);
                linesRead++;
            } else {
                break;
            }
        }

        if (linesRead == 0) {
            break; // EOF
        }

        if (linesRead < a) {
            // Incomplete gene block
            std::vector<std::string> tokens = splitLine(lines[0]);
            std::string geneID = tokens.size() >= 1 ? tokens[0] : "UNKNOWN";
            throw std::runtime_error(
                "Group file is incomplete for " + geneID + ".");
        }

        // Validate var line
        std::vector<std::string> tokens_var = splitLine(lines[0]);
        lineNum++;
        if (tokens_var.size() < 3) {
            throw std::runtime_error(
                "Error, group file line " + std::to_string(lineNum) +
                " is incomplete.");
        }
        std::string geneID0 = tokens_var[0];
        std::string type_var = tokens_var[1];
        if (type_var != "var") {
            throw std::runtime_error(
                "Error! No var is specified for " + geneID0);
        }
        int numMarkers = static_cast<int>(tokens_var.size()) - 2;

        // Validate anno line
        std::vector<std::string> tokens_anno = splitLine(lines[1]);
        lineNum++;
        if (tokens_anno.size() < 3) {
            throw std::runtime_error(
                "Error, group file line " + std::to_string(lineNum) +
                " is incomplete.");
        }
        std::string geneID_anno = tokens_anno[0];
        std::string type_anno = tokens_anno[1];
        if (type_anno != "anno") {
            throw std::runtime_error(
                "Error! No anno is specified for " + geneID_anno);
        }
        if (geneID_anno != geneID0) {
            throw std::runtime_error(
                "anno for " + geneID0 + " is missing.");
        }
        int numAnnos = static_cast<int>(tokens_anno.size()) - 2;
        if (numAnnos != numMarkers) {
            throw std::runtime_error(
                "The length of annotations for markers in region " + geneID0 +
                " is not equal to the length of marker IDs");
        }

        // Validate weight line (if applicable)
        if (is_weight_included) {
            std::vector<std::string> tokens_wt = splitLine(lines[2]);
            lineNum++;
            if (tokens_wt.size() < 3) {
                throw std::runtime_error(
                    "Error, group file line " + std::to_string(lineNum) +
                    " is incomplete.");
            }
            std::string geneID_wt = tokens_wt[0];
            std::string type_wt = tokens_wt[1];
            if (type_wt != "weight") {
                throw std::runtime_error(
                    "Error! No weight is specified for " + geneID_wt);
            }
            if (geneID_wt != geneID0) {
                throw std::runtime_error(
                    "weight for " + geneID0 + " is missing.");
            }
            int numWeights = static_cast<int>(tokens_wt.size()) - 2;
            if (numWeights != numMarkers) {
                throw std::runtime_error(
                    "The length of weights for markers in region " + geneID0 +
                    " is not equal to the length of marker IDs");
            }
        }

        nregion++;
    }
    gf2.close();

    GroupFileInfo info;
    info.nRegions = nregion;
    info.is_weight_included = is_weight_included;
    info.nline_per_gene = a;
    return info;
}


// ============================================================
// readRegionChunk
// Direct port from SAIGE/R/SAIGE_SPATest_Region.R: SAIGE.getRegionList_new()
//
// Reads nregions_to_read * nline_per_gene lines from the file stream.
// For each gene: parses var line, anno line, optional weight line.
// Looks up each variant ID in markerIDToIndex to get genotype indices.
// Removes variants not found in the genotype file.
// Builds annoIndicatorMat: for each annotation group (e.g. "lof;missense"),
//   splits by ";", marks variants whose annotation matches any part.
// Removes variants that don't match any annotation group.
// Filters out annotation groups that have no variants.
// ============================================================
std::vector<RegionData> readRegionChunk(
    std::ifstream& gf,
    int nregions_to_read,
    int nline_per_gene,
    const std::vector<std::string>& annoVec,
    const std::unordered_map<std::string, uint32_t>& markerIDToIndex)
{
    if (annoVec.empty()) {
        throw std::runtime_error("At least one annotation is required");
    }

    // --- Step 1: Read all lines for this chunk ---
    int nlinetoread = nregions_to_read * nline_per_gene;
    std::vector<std::string> allLines;
    allLines.reserve(nlinetoread);

    for (int k = 0; k < nlinetoread; k++) {
        std::string line;
        if (std::getline(gf, line)) {
            stripCR(line);
            allLines.push_back(line);
        } else {
            break; // EOF
        }
    }

    int ngroup = static_cast<int>(allLines.size()) / nline_per_gene;

    // --- Step 2: Pre-split each annotation group by ";" for matching ---
    // annoVec contains groups like {"lof", "lof;missense", "lof;missense;synonymous"}
    // For "lof;missense", a variant with annotation "lof" OR "missense" should match.
    std::vector<std::unordered_set<std::string>> annoGroupSets(annoVec.size());
    for (size_t q = 0; q < annoVec.size(); q++) {
        std::vector<std::string> parts = splitByChar(annoVec[q], ';');
        for (const auto& p : parts) {
            annoGroupSets[q].insert(p);
        }
    }

    // --- Step 3: Parse each gene block ---
    std::vector<RegionData> results;
    results.reserve(ngroup);

    std::unordered_set<std::string> genesSeen;

    for (int i = 0; i < ngroup; i++) {
        RegionData region;

        // Parse var line
        int varLineIdx = i * nline_per_gene;
        std::vector<std::string> tokensVar = splitLine(allLines[varLineIdx]);
        if (tokensVar.size() < 3) {
            throw std::runtime_error(
                "Error, group file: var line too short for gene block " +
                std::to_string(i + 1));
        }
        std::string gene = tokensVar[0];
        // Check for duplicate genes
        if (genesSeen.count(gene)) {
            throw std::runtime_error(gene + " is duplicated in the group file");
        }
        genesSeen.insert(gene);
        region.regionName = gene;

        int nvar = static_cast<int>(tokensVar.size()) - 2;
        std::vector<std::string> varIDs(tokensVar.begin() + 2, tokensVar.end());

        // Parse anno line
        int annoLineIdx = i * nline_per_gene + 1;
        std::vector<std::string> tokensAnno = splitLine(allLines[annoLineIdx]);
        if (tokensAnno.size() < 3) {
            throw std::runtime_error(
                "Error, group file: anno line too short for gene " + gene);
        }
        int nanno = static_cast<int>(tokensAnno.size()) - 2;
        if (nanno != nvar) {
            throw std::runtime_error(
                "The length of annotations for markers in region " + gene +
                " is not equal to the length of marker IDs");
        }
        std::vector<std::string> annoLabels(tokensAnno.begin() + 2, tokensAnno.end());

        // Parse weight line (if present)
        std::vector<double> weightVals;
        if (nline_per_gene == 3) {
            int weightLineIdx = i * nline_per_gene + 2;
            std::vector<std::string> tokensWeight = splitLine(allLines[weightLineIdx]);
            if (tokensWeight.size() < 3) {
                throw std::runtime_error(
                    "Error, group file: weight line too short for gene " + gene);
            }
            int nwt = static_cast<int>(tokensWeight.size()) - 2;
            if (nwt != nvar) {
                throw std::runtime_error(
                    "The length of weights for markers in region " + gene +
                    " is not equal to the length of marker IDs");
            }
            weightVals.resize(nwt);
            for (int w = 0; w < nwt; w++) {
                try {
                    weightVals[w] = std::stod(tokensWeight[w + 2]);
                } catch (const std::exception& e) {
                    throw std::runtime_error(
                        "Invalid weight value '" + tokensWeight[w + 2] +
                        "' for gene " + gene);
                }
            }
        }

        // --- Step 4: Look up variant IDs in markerIDToIndex ---
        // Remove variants not found in the genotype file.
        // Mirrors the R code: merge(RegionData, markerInfo, ...) + removing NA genoIndex
        std::vector<std::string> filteredVar;
        std::vector<std::string> filteredAnno;
        std::vector<double> filteredWeight;
        std::vector<std::string> filteredGenoIndex;
        std::vector<std::string> filteredGenoIndexPrev;

        int notFoundCount = 0;
        for (int v = 0; v < nvar; v++) {
            auto it = markerIDToIndex.find(varIDs[v]);
            if (it != markerIDToIndex.end()) {
                filteredVar.push_back(varIDs[v]);
                filteredAnno.push_back(annoLabels[v]);
                if (!weightVals.empty()) {
                    filteredWeight.push_back(weightVals[v]);
                }
                uint32_t idx = it->second;
                filteredGenoIndex.push_back(std::to_string(idx));
                // For PLINK, genoIndex_prev is the previous marker's index.
                // We use the same index here; the caller/mainRegionInCPP handles
                // the sequential access pattern.
                filteredGenoIndexPrev.push_back(std::to_string(idx));
            } else {
                notFoundCount++;
            }
        }

        if (notFoundCount > 0) {
            std::cout << notFoundCount
                      << " markers in region " << gene
                      << " from 'RegionFile' are not in 'GenoFile'."
                      << std::endl;
        }

        int nFiltered = static_cast<int>(filteredVar.size());

        if (nFiltered == 0) {
            // No variants found -- return an empty region
            region.variantIDs.clear();
            region.annotations.clear();
            region.weights.clear();
            region.annoIndicatorMat.reset();
            region.annoVec.clear();
            region.genoIndex.clear();
            region.genoIndex_prev.clear();
            results.push_back(std::move(region));
            continue;
        }

        // --- Step 5: Build annoIndicatorMat ---
        // Matrix is (nFiltered x annoVec.size()), value 1 if variant's annotation
        // matches any sub-annotation in the group.
        //
        // Mirrors R code:
        //   for(q in 1:length(annoVec)) {
        //     indiceVec = which(RegionData$ANNO[posSNP] %in% RegionAnnoHeaderList[[q]])
        //     annoIndicatorMat[indiceVec, q] = 1
        //   }
        arma::imat annoIndMat(nFiltered, static_cast<int>(annoVec.size()), arma::fill::zeros);

        // Track which annotation groups actually have at least one variant
        std::vector<bool> annoGroupHasVariants(annoVec.size(), false);

        for (size_t q = 0; q < annoVec.size(); q++) {
            for (int v = 0; v < nFiltered; v++) {
                if (annoGroupSets[q].count(filteredAnno[v])) {
                    annoIndMat(v, static_cast<int>(q)) = 1;
                    annoGroupHasVariants[q] = true;
                }
            }
        }

        // --- Step 6: Filter out annotation groups that have no variants ---
        // Mirrors R code:
        //   if(length(annoVecNew) < length(annoVec)) { rebuild annoIndicatorMat }
        std::vector<std::string> annoVecNew;
        std::vector<size_t> keptAnnoIndices;
        for (size_t q = 0; q < annoVec.size(); q++) {
            if (annoGroupHasVariants[q]) {
                annoVecNew.push_back(annoVec[q]);
                keptAnnoIndices.push_back(q);
            }
        }

        if (annoVecNew.empty()) {
            // No annotation groups matched any variant
            std::cerr << "WARNING: No markers are found for at least one annotation, "
                      << "so region " << gene << " is skipped" << std::endl;
            region.variantIDs = filteredVar;
            region.annotations = filteredAnno;
            region.weights = filteredWeight;
            region.annoIndicatorMat.reset();
            region.annoVec.clear();
            region.genoIndex = filteredGenoIndex;
            region.genoIndex_prev = filteredGenoIndexPrev;
            results.push_back(std::move(region));
            continue;
        }

        // If some groups were removed, rebuild the indicator matrix with only kept columns
        // Also re-split the kept groups for the rebuild (R does this in the
        // "if(length(annoVecNew) < length(annoVec))" block)
        arma::imat annoIndMatFinal;
        if (keptAnnoIndices.size() < annoVec.size()) {
            annoIndMatFinal.set_size(nFiltered, static_cast<int>(keptAnnoIndices.size()));
            annoIndMatFinal.zeros();
            for (size_t qi = 0; qi < keptAnnoIndices.size(); qi++) {
                size_t origQ = keptAnnoIndices[qi];
                // Re-split and re-match (exactly as R does in the rebuild loop)
                std::vector<std::string> parts = splitByChar(annoVecNew[qi], ';');
                std::unordered_set<std::string> partSet(parts.begin(), parts.end());
                for (int v = 0; v < nFiltered; v++) {
                    if (partSet.count(filteredAnno[v])) {
                        annoIndMatFinal(v, static_cast<int>(qi)) = 1;
                    }
                }
            }
        } else {
            annoIndMatFinal = annoIndMat;
        }

        // --- Step 7: Remove variants that don't match ANY annotation group ---
        // Mirrors R code:
        //   annoIndicatorMat_rmind = which(rowSums(annoIndicatorMat) == 0)
        //   if(length(annoIndicatorMat_rmind) > 0) { remove those rows }
        std::vector<int> keepRows;
        for (int v = 0; v < nFiltered; v++) {
            int rowSum = 0;
            for (int c = 0; c < static_cast<int>(annoIndMatFinal.n_cols); c++) {
                rowSum += annoIndMatFinal(v, c);
            }
            if (rowSum > 0) {
                keepRows.push_back(v);
            }
        }

        if (static_cast<int>(keepRows.size()) < nFiltered) {
            // Subset everything to keep only matching rows
            std::vector<std::string> keptVar, keptAnno, keptGenoIdx, keptGenoIdxPrev;
            std::vector<double> keptWeight;
            int nKept = static_cast<int>(keepRows.size());
            arma::imat keptIndMat(nKept, static_cast<int>(annoIndMatFinal.n_cols),
                                  arma::fill::zeros);

            for (int ki = 0; ki < nKept; ki++) {
                int origRow = keepRows[ki];
                keptVar.push_back(filteredVar[origRow]);
                keptAnno.push_back(filteredAnno[origRow]);
                keptGenoIdx.push_back(filteredGenoIndex[origRow]);
                keptGenoIdxPrev.push_back(filteredGenoIndexPrev[origRow]);
                if (!filteredWeight.empty()) {
                    keptWeight.push_back(filteredWeight[origRow]);
                }
                for (int c = 0; c < static_cast<int>(annoIndMatFinal.n_cols); c++) {
                    keptIndMat(ki, c) = annoIndMatFinal(origRow, c);
                }
            }

            region.variantIDs = std::move(keptVar);
            region.annotations = std::move(keptAnno);
            region.weights = std::move(keptWeight);
            region.annoIndicatorMat = std::move(keptIndMat);
            region.annoVec = std::move(annoVecNew);
            region.genoIndex = std::move(keptGenoIdx);
            region.genoIndex_prev = std::move(keptGenoIdxPrev);
        } else {
            // All variants match at least one group -- no removal needed
            region.variantIDs = std::move(filteredVar);
            region.annotations = std::move(filteredAnno);
            region.weights = std::move(filteredWeight);
            region.annoIndicatorMat = std::move(annoIndMatFinal);
            region.annoVec = std::move(annoVecNew);
            region.genoIndex = std::move(filteredGenoIndex);
            region.genoIndex_prev = std::move(filteredGenoIndexPrev);
        }

        results.push_back(std::move(region));
    }

    return results;
}
