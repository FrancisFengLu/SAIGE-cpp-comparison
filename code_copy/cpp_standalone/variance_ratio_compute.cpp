// variance_ratio_compute.cpp
// ------------------------------------------------------------------
// Real variance ratio computation matching R's extractVarianceRatio()
// in SAIGE/R/SAIGE_fitGLMM_fast.R lines 2569-2956.
//
// Algorithm:
//   For each MAC category k:
//     For each marker i (default 30, adaptive):
//       1. G0 = raw genotype from PLINK
//       2. Flip to minor allele if AF > 0.5
//       3. AC = sum(G0), skip if AC < 2
//       4. G = G0 - XXVX_inv * (XV * G0)   (covariate-adjusted)
//       5. g = G / sqrt(AC)                 (normalized)
//       6. Sigma_iG = PCG solve Σ^{-1} G
//       7. var1 = (G' Sigma_iG - G' Sigma_iX (X' Sigma_iX)^{-1} X' Sigma_iG) / AC
//       8. var2 = innerProduct(mu*(1-mu), g*g) for binary
//              or innerProduct(g, g)           for quantitative
//       9. ratio_i = var1 / var2
//     Average ratios; if CV > threshold, add 10 more markers and repeat.
//   Write output: tab-separated, no header, 3 columns: value \t type \t category
// ------------------------------------------------------------------

#include "variance_ratio_compute.hpp"
#include "SAIGE_step1_fast.hpp"
#include "score.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <algorithm>
#include <cmath>
#include <filesystem>

namespace saige {

// ---- helpers ----

// innerProduct: sum(a .* b)
static double innerProduct(const arma::fvec& a, const arma::fvec& b) {
    return arma::dot(a, b);
}

void compute_variance_ratio(const Paths& paths,
                            const FitNullConfig& cfg,
                            const LocoRanges& /*chr*/,
                            const FitNullResult& fit,
                            const Design& design,
                            std::string& out_vr_path,
                            std::string& out_marker_results_path)
{
    std::cout << "\n===== Variance Ratio Computation =====" << std::endl;

    // VR marker bypass: when true, read marker indices from R's output
    // Set to true for testing exact match with R, then set back to false
    bool use_r_vr_bypass = true;

    const int n = design.n;
    const int p = design.p;
    const int numMarkers_default = cfg.num_markers_for_vr;  // default 30
    const float tolPCG = static_cast<float>(cfg.tolPCG);
    const int maxiterPCG = cfg.maxiterPCG;
    const double ratioCVcutoff = cfg.ratio_cv_cutoff;
    const bool is_binary = (cfg.trait == "binary");
    // const bool is_quant = (cfg.trait == "quantitative");

    // --- Extract tau and build W vector ---
    arma::fvec tauVec(fit.theta.size());
    for (size_t i = 0; i < fit.theta.size(); ++i)
        tauVec(i) = static_cast<float>(fit.theta[i]);

    std::cout << "[VR] tau = [";
    for (size_t i = 0; i < fit.theta.size(); ++i)
        std::cout << fit.theta[i] << (i+1 < fit.theta.size() ? ", " : "");
    std::cout << "]\n";

    // --- Build mu, y vectors from the GLMM fit ---
    // We need to reconstruct mu from the null model.
    // The GLMM solver should have stored mu in ScoreNullPack (obj_noK.mu).
    // However, obj_noK might not be populated yet at this point.
    // Instead, we reconstruct from X*alpha + offset through the link function.

    // Build X matrix (arma::fmat) from design
    arma::fmat X(n, p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < p; ++j)
            X(i, j) = static_cast<float>(design.X[i * p + j]);  // row-major

    arma::fvec y_vec(n);
    for (int i = 0; i < n; ++i)
        y_vec(i) = static_cast<float>(design.y[i]);

    // Compute eta = X * alpha + offset
    arma::fvec alpha_f(p);
    for (int i = 0; i < p; ++i)
        alpha_f(i) = static_cast<float>(fit.alpha[i]);

    arma::fvec eta = X * alpha_f;
    if (!fit.offset.empty()) {
        for (int i = 0; i < n; ++i)
            eta(i) += static_cast<float>(fit.offset[i]);
    }

    // Compute mu from eta
    arma::fvec mu(n);
    if (is_binary) {
        for (int i = 0; i < n; ++i) {
            float e = std::exp(eta(i));
            mu(i) = e / (1.0f + e);  // logistic
        }
    } else {
        // quantitative: identity link
        mu = eta;
    }

    // W = working weights
    arma::fvec W(n);
    if (is_binary) {
        W = mu % (1.0f - mu);  // mu*(1-mu)
    } else {
        // quantitative: W = 1.0 (Gaussian IRLS working weight)
        W.fill(1.0f);
    }

    std::cout << "[VR] mu[0:5]: ";
    for (int i = 0; i < std::min(5, n); ++i) std::cout << mu(i) << " ";
    std::cout << "\n";
    std::cout << "[VR] W[0:5]: ";
    for (int i = 0; i < std::min(5, n); ++i) std::cout << W(i) << " ";
    std::cout << "\n";

    // --- Compute Sigma_iX (global, no LOCO) ---
    // getSigma_X solves Σ^{-1} X column by column via PCG
    arma::fvec W_copy = W;
    arma::fvec tau_copy = tauVec;
    arma::fmat X_copy = X;
    arma::fmat Sigma_iX = getSigma_X(W_copy, tau_copy, X_copy, maxiterPCG, tolPCG);

    std::cout << "[VR] Sigma_iX computed: " << Sigma_iX.n_rows << " x " << Sigma_iX.n_cols << "\n";

    // Precompute (X' Sigma_iX)^{-1}
    arma::fmat XtSiX = X.t() * Sigma_iX;  // p x p
    arma::fmat XtSiX_inv = arma::inv_sympd(arma::symmatu(XtSiX));

    // Precompute XXVX_inv and XV for covariate adjustment of genotype
    // G_adj = G0 - XXVX_inv * (XV * G0)
    // where XXVX_inv = X * (X'VX)^{-1}, XV = (X ⊙ V)'
    // But for VR, we use the score-null versions.
    // Build ScoreNull for covariate adjustment
    ScoreNull sn;
    if (is_binary) {
        sn = build_score_null_binary(X, y_vec, mu);
    } else {
        sn = build_score_null_quant(X, y_vec, mu, 1.0f / tauVec(0));
    }

    std::cout << "[VR] ScoreNull built. XV: " << sn.XV.n_rows << "x" << sn.XV.n_cols
              << ", XXVX_inv: " << sn.XXVX_inv.n_rows << "x" << sn.XXVX_inv.n_cols << "\n";

    // --- Get markers for VR ---
    bool isVarRatioGeno = getIsVarRatioGeno();
    std::cout << "[VR] isVarRatioGeno: " << isVarRatioGeno << "\n";

    // Marker indices and MAC values
    arma::ivec macVec, indexVec;
    int numAvailMarkers;

    if (isVarRatioGeno) {
        // Separate VR genotype data is loaded
        macVec = getMACVec_forVarRatio();
        indexVec = getIndexVec_forVarRatio();
        numAvailMarkers = macVec.n_elem;
    } else {
        // Use main genotype data
        macVec = getMACVec();
        numAvailMarkers = macVec.n_elem;
        indexVec.set_size(numAvailMarkers);
        for (int i = 0; i < numAvailMarkers; ++i)
            indexVec(i) = i;
    }

    std::cout << "[VR] Available markers for VR: " << numAvailMarkers << "\n";

    // --- VR marker bypass: read marker indices from R's output ---
    // When use_r_vr_bypass is true, we read the exact marker order from R
    // instead of using our own random shuffle. This lets us verify that
    // the per-marker var1/var2/ratio computations match R exactly.
    struct BypassMarker {
        int snp_index;
        int geno_ind;
        int orig_plink_index;  // 0-based index into main plink file
    };
    std::vector<BypassMarker> bypass_markers;
    bool bypass_active = false;

    if (use_r_vr_bypass) {
        std::string bypass_path = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/bypass/vr_marker_indices.csv";
        std::ifstream bypass_ifs(bypass_path);
        if (bypass_ifs.is_open()) {
            std::string header_line;
            std::getline(bypass_ifs, header_line);  // skip header: snp_index,geno_ind,mac,var1,var2null,ratio,orig_plink_index

            std::string line;
            while (std::getline(bypass_ifs, line)) {
                if (line.empty()) continue;
                std::istringstream iss(line);
                std::string token;

                // Parse snp_index (column 0)
                if (!std::getline(iss, token, ',')) continue;
                int snp_index = std::stoi(token);

                // Parse geno_ind (column 1)
                if (!std::getline(iss, token, ',')) continue;
                int geno_ind = std::stoi(token);

                // Skip mac (column 2), var1 (3), var2null (4), ratio (5)
                for (int skip = 0; skip < 4; ++skip) {
                    if (!std::getline(iss, token, ',')) break;
                }
                // Parse orig_plink_index (column 6)
                int orig_plink_idx = -1;
                if (std::getline(iss, token, ',')) {
                    orig_plink_idx = std::stoi(token);
                }

                bypass_markers.push_back({snp_index, geno_ind, orig_plink_idx});
            }
            bypass_ifs.close();

            std::cout << "[VR] Bypass: loaded " << bypass_markers.size()
                      << " marker indices from R" << std::endl;
            bypass_active = true;
        } else {
            std::cout << "[VR] WARNING: Bypass file not found: " << bypass_path << std::endl;
            std::cout << "[VR] Falling back to random marker selection." << std::endl;
        }
    }

    // --- Marker ordering ---
    // If bypass is active, we use the bypass_markers vector directly.
    // Otherwise, we shuffle marker indices as before.
    std::vector<int> markerOrder;
    if (!bypass_active) {
        markerOrder.resize(numAvailMarkers);
        std::iota(markerOrder.begin(), markerOrder.end(), 0);

        // Use a deterministic seed for reproducibility (R uses set.seed internally)
        std::mt19937 rng(200);
        std::shuffle(markerOrder.begin(), markerOrder.end(), rng);
    }

    // --- VR computation loop ---
    int numMarkers0 = numMarkers_default;
    float ratioCV = 1.0f;  // start above threshold
    int numTestedMarker = 0;
    int indexInMarkerList = 0;
    int totalAvailable = bypass_active ? static_cast<int>(bypass_markers.size()) : numAvailMarkers;

    arma::fvec varRatio_NULL_vec;
    arma::fvec varRatio_NULL_noXadj_vec;

    int Nnomissing = getNnomissingOut();
    std::cout << "[VR] Nnomissing: " << Nnomissing << "\n";

    // Marker results for output
    struct MarkerResult {
        int snpIdx;
        double mac;
        double af;
        double var1;
        double var2;
        double ratio;
    };
    std::vector<MarkerResult> markerResults;

    while (ratioCV > ratioCVcutoff) {
        while (numTestedMarker < numMarkers0 && indexInMarkerList < totalAvailable) {
            int snpIdx;
            bool genoInd;

            if (bypass_active) {
                // Use the original plink marker index from R's output
                // R's snp_index is a VR-array index (useless without the same VR array).
                // orig_plink_index is the 0-based index into the main plink file.
                int orig_plink_idx = bypass_markers[indexInMarkerList].orig_plink_index;
                if (orig_plink_idx >= 0) {
                    int mainIdx = findMainArrayIdx(orig_plink_idx);
                    if (mainIdx < 0) {
                        std::cout << "[VR] WARNING: plink index " << orig_plink_idx
                                  << " not found in main geno array, skipping" << std::endl;
                        indexInMarkerList++;
                        continue;
                    }
                    snpIdx = mainIdx;
                    genoInd = false;  // Always use main plink Get_OneSNP_Geno
                } else {
                    // Fallback: try VR index (will fail if VR geno not loaded)
                    snpIdx = bypass_markers[indexInMarkerList].snp_index;
                    genoInd = (bypass_markers[indexInMarkerList].geno_ind != 0);
                }
            } else {
                // Original random shuffle path
                int macdata_i = markerOrder[indexInMarkerList];
                snpIdx = indexVec(macdata_i);
                genoInd = isVarRatioGeno;
            }
            indexInMarkerList++;

            // Get raw genotype
            arma::ivec G0;
            if (!genoInd) {
                G0 = Get_OneSNP_Geno(snpIdx);
            } else {
                G0 = Get_OneSNP_Geno_forVarRatio(snpIdx);
            }

            // Flip to minor allele if AF > 0.5
            int sumG0 = arma::sum(G0);
            if (static_cast<double>(sumG0) / (2.0 * Nnomissing) > 0.5) {
                G0 = 2 - G0;
                sumG0 = arma::sum(G0);
            }

            double AC = static_cast<double>(sumG0);
            double AF = AC / (2.0 * Nnomissing);

            // Skip markers with very low AC
            if (AC < 2.0) continue;

            // Convert to float for computation
            arma::fvec G0f = arma::conv_to<arma::fvec>::from(G0);

            // Covariate-adjusted genotype: G = G0 - XXVX_inv * (XV * G0)
            arma::fvec G = G0f - sn.XXVX_inv * (sn.XV * G0f);

            // Normalized genotype
            arma::fvec g = G / std::sqrt(static_cast<float>(AC));

            // Also compute non-X-adjusted version for comparison
            arma::fvec G_noXadj = G0f - arma::mean(G0f);
            arma::fvec g_noXadj = G_noXadj / std::sqrt(static_cast<float>(AC));

            // --- var1 (exact): Sigma^{-1} based ---
            // Sigma_iG = PCG solve for Σ^{-1} G
            arma::fvec G_for_pcg = G0f;  // R uses G (not G0) — but actually R uses G (the raw genotype, not covariate-adjusted)
            // Actually, looking at R code more carefully:
            // Line 2850: Sigma_iG = getSigma_G(W, tauVecNew, G, maxiterPCG, tolPCG)
            // where G was already covariate-adjusted at line 2833.
            // Wait - R's line 2833: G = G0 - obj.noK$XXVX_inv %*% (obj.noK$XV %*% G0)
            // And line 2850: Sigma_iG = getSigma_G(W, tauVecNew, G, ...)
            // So it uses the covariate-adjusted G for PCG solve.
            // But then line 2861: var1a = t(G)%*%Sigma_iG - t(G)%*%Sigma_iX%*%(solve(t(X)%*%Sigma_iX))%*%t(X)%*%Sigma_iG
            // This seems redundant (adjusting both G and subtracting projection), but we must match R exactly.

            arma::fvec Sigma_iG = getPCG1ofSigmaAndVector(W_copy, tau_copy, G, maxiterPCG, tolPCG);

            // var1a = G' * Sigma_iG - G' * Sigma_iX * (X' * Sigma_iX)^{-1} * X' * Sigma_iG
            float GtSiG = arma::dot(G, Sigma_iG);
            arma::fvec XtSiG = Sigma_iX.t() * G;     // p x 1 -- wait, this should be X' * Sigma_iG
            // Actually R: t(X) %*% Sigma_iG  and  t(G) %*% Sigma_iX
            // Let me re-read:
            //   var1a = t(G) %*% Sigma_iG
            //         - t(G) %*% Sigma_iX %*% solve(t(X) %*% Sigma_iX) %*% t(X) %*% Sigma_iG
            // where Sigma_iX = Σ^{-1} X
            // So: term2 = G' * (Σ^{-1} X) * (X' Σ^{-1} X)^{-1} * X' * (Σ^{-1} G)

            arma::fvec GtSiX = Sigma_iX.t() * G;     // p x 1: (Σ^{-1}X)' G = X' Σ^{-1} G... no.
            // Sigma_iX is n x p = Σ^{-1} X
            // G' Sigma_iX = G' (Σ^{-1} X) which is 1 x p
            arma::fvec GtSiX_vec(p);
            for (int j = 0; j < p; ++j)
                GtSiX_vec(j) = arma::dot(G, Sigma_iX.col(j));  // G' * col_j of Sigma_iX

            arma::fvec XtSiG_vec(p);
            for (int j = 0; j < p; ++j)
                XtSiG_vec(j) = arma::dot(X.col(j), Sigma_iG);  // X(:,j)' * Sigma_iG

            // term2 = GtSiX' * XtSiX_inv * XtSiG
            float term2 = arma::dot(GtSiX_vec, XtSiX_inv * XtSiG_vec);

            double var1 = (GtSiG - term2) / AC;

            // --- var2 (approximate): null-model based ---
            double var2;
            if (is_binary) {
                // var2 = sum(mu*(1-mu) * g^2)
                var2 = innerProduct(mu % (1.0f - mu), g % g);
            } else {
                // var2 = sum(g^2) for quantitative
                var2 = innerProduct(g, g);
            }

            // Also compute noXadj version
            double var2_noXadj;
            if (is_binary) {
                var2_noXadj = innerProduct(mu % (1.0f - mu), g_noXadj % g_noXadj);
            } else {
                var2_noXadj = innerProduct(g_noXadj, g_noXadj);
            }

            double ratio = var1 / var2;
            double ratio_noXadj = var1 / var2_noXadj;

            varRatio_NULL_vec.resize(varRatio_NULL_vec.n_elem + 1);
            varRatio_NULL_vec(varRatio_NULL_vec.n_elem - 1) = static_cast<float>(ratio);
            varRatio_NULL_noXadj_vec.resize(varRatio_NULL_noXadj_vec.n_elem + 1);
            varRatio_NULL_noXadj_vec(varRatio_NULL_noXadj_vec.n_elem - 1) = static_cast<float>(ratio_noXadj);

            numTestedMarker++;

            markerResults.push_back({snpIdx, AC, AF, var1, var2, ratio});

            if (numTestedMarker % 10 == 0 || numTestedMarker == numMarkers0) {
                std::cout << "[VR] Marker " << numTestedMarker << "/" << numMarkers0
                          << ": AC=" << AC << " AF=" << AF
                          << " var1=" << var1 << " var2=" << var2
                          << " ratio=" << ratio << "\n";
            }
        }

        // Check CV
        if (numTestedMarker > 0) {
            ratioCV = calCV(varRatio_NULL_vec);
        }

        if (ratioCV > ratioCVcutoff) {
            std::cout << "[VR] CV for variance ratio estimate using " << numMarkers0
                      << " markers is " << ratioCV << " > " << ratioCVcutoff << "\n";
            numMarkers0 += 10;
            std::cout << "[VR] Trying " << numMarkers0 << " markers\n";

            // Check if we've run out of markers
            if (indexInMarkerList >= totalAvailable) {
                std::cout << "[VR] No more markers available. Stopping.\n";
                break;
            }
        } else {
            std::cout << "[VR] CV for variance ratio estimate using " << numMarkers0
                      << " markers is " << ratioCV << " <= " << ratioCVcutoff << "\n";
        }
    }

    // --- Average the ratios ---
    double varRatio_null = (varRatio_NULL_vec.n_elem > 0)
                             ? arma::mean(varRatio_NULL_vec)
                             : 1.0;
    double varRatio_null_noXadj = (varRatio_NULL_noXadj_vec.n_elem > 0)
                                    ? arma::mean(varRatio_NULL_noXadj_vec)
                                    : 1.0;

    std::cout << "[VR] Final varRatio_null: " << varRatio_null
              << " (using " << numTestedMarker << " markers)\n";
    std::cout << "[VR] Final varRatio_null_noXadj: " << varRatio_null_noXadj << "\n";

    // --- Write output ---
    // Match R format: write.table(varRatioTable, file, quote=F, col.names=F, row.names=F)
    // Each row: value \t type \t category
    std::string vr_out = paths.out_prefix_vr.empty()
                           ? (paths.out_prefix + ".varianceRatio.txt")
                           : (paths.out_prefix_vr + ".varianceRatio.txt");

    std::filesystem::create_directories(std::filesystem::path(vr_out).parent_path());

    {
        std::ofstream ofs(vr_out);
        if (!ofs) throw std::runtime_error("Cannot write VR file: " + vr_out);

        // Category k=1 (only one category for non-categorical mode)
        int k = 1;
        ofs << varRatio_null << "\tnull\t" << k << "\n";
        ofs << varRatio_null_noXadj << "\tnull_noXadj\t" << k << "\n";
    }

    std::cout << "[VR] Wrote variance ratios to: " << vr_out << "\n";

    // --- Write marker results (optional diagnostic file) ---
    std::string marker_out = paths.out_prefix_vr.empty()
                               ? (paths.out_prefix + "." + std::to_string(numTestedMarker) + "markers.SAIGE.results.txt")
                               : (paths.out_prefix_vr + "." + std::to_string(numTestedMarker) + "markers.SAIGE.results.txt");

    {
        std::ofstream ofs(marker_out);
        if (!ofs) throw std::runtime_error("Cannot write marker results: " + marker_out);

        ofs << "SNPIdx\tMAC\tAF\tvar1\tvar2\tratio\n";
        for (const auto& mr : markerResults) {
            ofs << mr.snpIdx << "\t" << mr.mac << "\t" << mr.af
                << "\t" << mr.var1 << "\t" << mr.var2 << "\t" << mr.ratio << "\n";
        }
    }

    std::cout << "[VR] Wrote marker results to: " << marker_out << "\n";
    std::cout << "===== Variance Ratio Complete =====\n" << std::endl;

    out_vr_path = vr_out;
    out_marker_results_path = marker_out;
}

} // namespace saige
