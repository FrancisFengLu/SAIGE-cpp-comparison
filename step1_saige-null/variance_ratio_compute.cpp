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

// Forward declaration: sparse-GRM linear solve defined in SAIGE_step1_fast.cpp.
// Solves Σ_sparse^{-1} · y where Σ_sparse = tau0·diag(1/w) + tau1·G_sparse.
// Used to compute var2_sparse for the SAIGE-GENE+ "sparse" variance-ratio row.
arma::fvec gen_spsolve_v4(arma::fvec& wVec, arma::fvec& tauVec, arma::fvec& yvec);

namespace saige {

// ---- helpers ----

// innerProduct: sum(a .* b)
static double innerProduct(const arma::fvec& a, const arma::fvec& b) {
    return arma::dot(a, b);
}

// Categorical-VR bin assignment from MAC.
// R: for i in [0..numCate-2): MAC > minExclude[i]  AND  MAC <= maxInclude[i]
//    for i == numCate-1:     MAC > minExclude[numCate-1]  (open upper)
// Returns 0..numBins-1, or -1 if marker is excluded (MAC <= minExclude[0]).
// When isCateVarianceRatio=false: returns 0 (single bin).
static int bin_of_mac(double mac, const FitNullConfig& cfg) {
    if (!cfg.isCateVarianceRatio) return 0;
    const auto& lo = cfg.cateVarRatioMinMACVecExclude;
    const auto& hi = cfg.cateVarRatioMaxMACVecInclude;
    if (lo.empty()) return 0;
    int nBins = static_cast<int>(lo.size());
    for (int b = 0; b < nBins; ++b) {
        if (mac <= lo[b]) continue;            // below this bin's lower bound
        if (b < static_cast<int>(hi.size())) {
            if (mac <= hi[b]) return b;        // bounded bin
        } else {
            return b;                          // open-upper last bin
        }
    }
    return -1;                                 // didn't fit any bin (mac <= lo[0])
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
    // ===== Step 6: VR bypass disabled for production (was true for R-comparison testing) =====
    bool use_r_vr_bypass = false;

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

    // Use the converged mu from the GLMM solver (includes BLUP b̂).
    // For VR on quantitative trait, W=1 and the formula doesn't reach mu, so
    // this is mainly a binary-trait correctness fix; for quant we still prefer
    // the converged mu so any downstream sn/score-null builds on the right values.
    // See glmm.cpp::stash_eta_mu_into and STEP2_PVALUE_FIX_2026-05-06.md.
    arma::fvec alpha_f(p);
    for (int i = 0; i < p; ++i)
        alpha_f(i) = static_cast<float>(fit.alpha[i]);

    arma::fvec mu(n);
    if (fit.mu.size() == static_cast<size_t>(n)) {
        for (int i = 0; i < n; ++i)
            mu(i) = static_cast<float>(fit.mu[i]);
    } else {
        // Fallback (legacy path — keeps the build green if solver didn't populate).
        std::cerr << "[VR] WARNING: fit.mu unset; falling back to X*alpha "
                     "(b̂ dropped — VR on binary trait may drift).\n";
        arma::fvec eta_fallback = X * alpha_f;
        if (!fit.offset.empty())
            for (int i = 0; i < n; ++i) eta_fallback(i) += static_cast<float>(fit.offset[i]);
        if (is_binary) {
            for (int i = 0; i < n; ++i) {
                float e = std::exp(eta_fallback(i));
                mu(i) = e / (1.0f + e);
            }
        } else {
            mu = eta_fallback;
        }
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

    // ===== Step 8: No markers found for VR (R line 3021) =====
    // R: stop("No markers were found for variance ratio estimation...")
    if (numAvailMarkers == 0) {
      throw std::runtime_error(
          "ERROR: No markers were found for variance ratio estimation. "
          "Please make sure there are markers with MAC >= "
          + std::to_string(cfg.vr_min_mac) + " in the plink file.");
    }

    // ===== Step 9: Insufficient markers for VR (R lines 3055-3063) =====
    // R: if(length(listOfMarkersForVarRatio[[k]]) < numMarkers) stop(...)
    if (numAvailMarkers < numMarkers_default) {
      std::cerr << "[warning] Only " << numAvailMarkers
                << " markers available for variance ratio estimation, but "
                << numMarkers_default << " requested. Using all available markers.\n";
    }

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

    // Per-bin VR accumulators (vec-of-vecs).
    // Bin count: numBins = cfg.cateVarRatioMinMACVecExclude.size() when categorical, else 1.
    const int numBins = cfg.isCateVarianceRatio
                          ? std::max<int>(1, static_cast<int>(cfg.cateVarRatioMinMACVecExclude.size()))
                          : 1;
    std::vector<arma::fvec> varRatio_NULL_vec_per_bin(numBins);
    std::vector<arma::fvec> varRatio_NULL_noXadj_vec_per_bin(numBins);
    std::vector<arma::fvec> varRatio_sparse_vec_per_bin(numBins);
    std::vector<int> numTestedMarker_per_bin(numBins, 0);
    std::vector<int> numTarget_per_bin(numBins, numMarkers_default);
    std::vector<float> ratioCV_per_bin(numBins, 1.0f);
    std::vector<bool> binConverged(numBins, false);
    // Compatibility aliases kept as flat views for the diagnostic output later.
    arma::fvec varRatio_NULL_vec;          // global pool (kept for legacy diagnostics)
    arma::fvec varRatio_NULL_noXadj_vec;   // global pool (kept for legacy diagnostics)

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
                // Original random shuffle path.
                // markerOrder[i] is the *VR-array sequential index* (0..numAvailMarkers-1),
                // which is exactly what Get_OneSNP_Geno_forVarRatio expects when
                // genoInd = true. The previously-used `indexVec(macdata_i)` returns the
                // ORIGINAL plink BIM index (~0..340 k) and would access way out-of-bounds
                // inside `genoVecofPointers_forVarRatio` (only ~numAvailMarkers entries).
                snpIdx  = markerOrder[indexInMarkerList];
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

            // --- Sparse-GRM VR: var2_sparse = g' · Σ_sparse^{-1} · g  (only if requested) ---
            // R: SAIGE_fitGLMM_fast.R:3164-3181 — emits the "sparse" VR row when
            // useSparseGRMforVarRatio=TRUE. var1 is the same exact-PCG numerator;
            // var2_sparse uses the sparse-only Σ (no dense GRM term).
            double ratio_sparse = 1.0;
            bool   have_sparse  = false;
            if (cfg.use_sparse_grm_for_vr) {
                // gen_spsolve_v4 needs non-const refs; copy w/tau locally.
                arma::fvec w_for_sp   = W_copy;
                arma::fvec tau_for_sp = tau_copy;
                arma::fvec G_for_sp   = G;  // covariate-adjusted genotype
                arma::fvec Sigma_sp_iG = gen_spsolve_v4(w_for_sp, tau_for_sp, G_for_sp);
                double var2_sparse = arma::dot(G, Sigma_sp_iG) / AC;
                if (var2_sparse > 0.0 && std::isfinite(var2_sparse)) {
                    ratio_sparse = var1 / var2_sparse;
                    have_sparse  = true;
                }
            }

            // --- Per-bin accumulation (categorical VR) ---
            int binId = bin_of_mac(AC, cfg);
            if (binId < 0) {
                // Marker below the lowest bin — skip without incrementing target counts.
                continue;
            }

            // Append to this bin's accumulators (and the legacy flat view).
            auto append_f = [](arma::fvec& v, double x) {
                v.resize(v.n_elem + 1);
                v(v.n_elem - 1) = static_cast<float>(x);
            };
            append_f(varRatio_NULL_vec_per_bin[binId],        ratio);
            append_f(varRatio_NULL_noXadj_vec_per_bin[binId], ratio_noXadj);
            if (have_sparse)
                append_f(varRatio_sparse_vec_per_bin[binId], ratio_sparse);
            append_f(varRatio_NULL_vec, ratio);
            append_f(varRatio_NULL_noXadj_vec, ratio_noXadj);
            numTestedMarker_per_bin[binId]++;
            numTestedMarker++;

            markerResults.push_back({snpIdx, AC, AF, var1, var2, ratio});

            if (numTestedMarker % 10 == 0) {
                std::cout << "[VR] Marker " << numTestedMarker
                          << " (bin " << (binId+1) << "/" << numBins << ")"
                          << ": MAC=" << AC << " AF=" << AF
                          << " var1=" << var1 << " var2=" << var2
                          << " ratio=" << ratio;
                if (have_sparse) std::cout << " ratio_sparse=" << ratio_sparse;
                std::cout << "\n";
            }
        }

        // ---- Per-bin convergence check ----
        bool all_converged = true;
        for (int b = 0; b < numBins; ++b) {
            // Skip-bin via cateVarRatioIndexVec[b]==0 → mark converged with default 1.0.
            bool skip_bin = (!cfg.cateVarRatioIndexVec.empty()
                              && b < static_cast<int>(cfg.cateVarRatioIndexVec.size())
                              && cfg.cateVarRatioIndexVec[b] == 0);
            if (skip_bin) { binConverged[b] = true; continue; }
            if (numTestedMarker_per_bin[b] >= numTarget_per_bin[b] && varRatio_NULL_vec_per_bin[b].n_elem > 0) {
                ratioCV_per_bin[b] = calCV(varRatio_NULL_vec_per_bin[b]);
                if (ratioCV_per_bin[b] <= ratioCVcutoff) {
                    binConverged[b] = true;
                    std::cout << "[VR] Bin " << (b+1) << ": CV=" << ratioCV_per_bin[b]
                              << " <= " << ratioCVcutoff << " using "
                              << numTestedMarker_per_bin[b] << " markers (converged)\n";
                } else {
                    // grow this bin's target
                    numTarget_per_bin[b] += 10;
                    std::cout << "[VR] Bin " << (b+1) << ": CV=" << ratioCV_per_bin[b]
                              << " > " << ratioCVcutoff << "; trying "
                              << numTarget_per_bin[b] << " markers in this bin\n";
                }
            }
            if (!binConverged[b]) all_converged = false;
        }
        // mirror to legacy single-bin ratioCV for log compatibility
        ratioCV = 0.0f;
        for (int b = 0; b < numBins; ++b) ratioCV = std::max(ratioCV, ratioCV_per_bin[b]);

        if (all_converged) {
            std::cout << "[VR] All " << numBins << " bin(s) converged.\n";
            break;
        }
        if (indexInMarkerList >= totalAvailable) {
            std::cout << "[VR] No more markers available. Stopping with "
                      << numTestedMarker << " markers tested across " << numBins << " bin(s).\n";
            break;
        }
        // ensure outer-loop guard fires even if some bin never advances numMarkers0
        numMarkers0 = numTestedMarker + 10;
    }

    // --- Per-bin aggregation + multi-row output (R SAIGE_fitGLMM_fast.R:3240-3255) ---
    // For each bin k (1-based in the output), emit:
    //   <vr>\tnull\t<k>
    //   <vr>\tnull_noXadj\t<k>
    //   <vr>\tsparse\t<k>            (only if cfg.use_sparse_grm_for_vr)
    // Bins flagged skip via cateVarRatioIndexVec[b]==0 emit default 1.0.
    std::vector<double> bin_null(numBins, 1.0), bin_noXadj(numBins, 1.0), bin_sparse(numBins, 1.0);
    for (int b = 0; b < numBins; ++b) {
        bool skip_bin = (!cfg.cateVarRatioIndexVec.empty()
                          && b < static_cast<int>(cfg.cateVarRatioIndexVec.size())
                          && cfg.cateVarRatioIndexVec[b] == 0);
        if (!skip_bin) {
            if (varRatio_NULL_vec_per_bin[b].n_elem        > 0) bin_null  [b] = arma::mean(varRatio_NULL_vec_per_bin[b]);
            if (varRatio_NULL_noXadj_vec_per_bin[b].n_elem > 0) bin_noXadj[b] = arma::mean(varRatio_NULL_noXadj_vec_per_bin[b]);
            if (varRatio_sparse_vec_per_bin[b].n_elem      > 0) bin_sparse[b] = arma::mean(varRatio_sparse_vec_per_bin[b]);
        }
        std::cout << "[VR] Bin " << (b+1) << ": null=" << bin_null[b]
                  << " null_noXadj=" << bin_noXadj[b];
        if (cfg.use_sparse_grm_for_vr) std::cout << " sparse=" << bin_sparse[b];
        std::cout << " (n=" << numTestedMarker_per_bin[b] << ")\n";
    }

    std::string vr_out = paths.out_prefix_vr.empty()
                           ? (paths.out_prefix + ".varianceRatio.txt")
                           : (paths.out_prefix_vr + ".varianceRatio.txt");

    std::filesystem::create_directories(std::filesystem::path(vr_out).parent_path());

    {
        std::ofstream ofs(vr_out);
        if (!ofs) throw std::runtime_error("Cannot write VR file: " + vr_out);
        // 3 rows per bin (sparse row only if use_sparse_grm_for_vr);
        // R order is: null, null_noXadj, [sparse]
        for (int b = 0; b < numBins; ++b) {
            int k = b + 1; // 1-based
            ofs << bin_null[b]   << "\tnull\t"        << k << "\n";
            ofs << bin_noXadj[b] << "\tnull_noXadj\t" << k << "\n";
            if (cfg.use_sparse_grm_for_vr) {
                ofs << bin_sparse[b] << "\tsparse\t" << k << "\n";
            }
        }
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
