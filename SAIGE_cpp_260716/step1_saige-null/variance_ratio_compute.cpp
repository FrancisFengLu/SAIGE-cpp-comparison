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
#include "fused_variance_ratio.hpp"
#include <fstream>
#include <cstdlib>
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
    int numMarkers_default = cfg.num_markers_for_vr;  // default 30 (lowered below when the fused anchor is on)
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

    // ---------------- fused variance ratio: the closed-form anchor ------------
    // fit.fused_variance_ratio. The per-bin ratio is anchored on the EXACT
    // closed form anchor = tr(P Psi)/tr((I-H) Psi) and markers are spent only on
    // testing the per-bin correction delta = mean(ratio)/anchor against 1.
    // Every gate below names a reason the closed form would not describe what
    // the sampled ratio measures; when one fires the flag is ignored out loud
    // and the ordinary sampled path runs unchanged.
    FusedVrAnchor fused;
    bool fused_on = false;
    if (cfg.fused_variance_ratio) {
        std::string why;
        if (is_binary)
            why = "the trait is binary -- var2 carries the working weights "
                  "mu(1-mu), so the closed-form denominator is neither tr(Psi) "
                  "nor tr(I); that case has not been derived or measured";
        else if (cfg.trait != "quantitative")
            why = "fit.trait is '" + cfg.trait + "'; only quantitative is derived";
        else if (!cfg.use_sparse_grm_to_fit)
            why = "fit.use_sparse_grm_to_fit is off -- with a dense GRM tr(P Psi) "
                  "is O(N^2) and Sigma is not block diagonal, so no trace here is exact";
        else if (!get_isUseSparseSigmaforModelFitting())
            why = "the sparse Sigma was not installed for the model fit";
        else if (cfg.use_sparse_grm_for_vr && cfg.fused_vr_markers <= 0)
            why = "fit.use_sparse_grm_for_vr needs the per-marker sparse row, which "
                  "has no closed form here, but fit.fused_vr_markers is 0";
        if (why.empty()) {
            arma::vec Wd(n), taud(2);
            for (int i = 0; i < n; ++i) Wd(i) = static_cast<double>(W(i));
            taud(0) = fit.theta.size() > 0 ? fit.theta[0] : 0.0;
            taud(1) = fit.theta.size() > 1 ? fit.theta[1] : 0.0;
            // The SAME X the var1/var2 loop below uses -- on the default
            // covariate_offset path that is the collapsed intercept-only design
            // (design.p == 1, the covariates having moved into the offset), and
            // the anchor has to describe what the sampled ratio measures, not
            // what step 2 later projects with. Measured cost of the difference
            // on the fair simulation: 1.2e-6 relative, both on the anchor and
            // on the per-marker mean.
            arma::mat Xd = arma::conv_to<arma::mat>::from(X);
            fused = compute_fused_vr_anchor(Wd, taud, Xd, cfg.fused_vr_max_block);
            if (!fused.ok) why = fused.why;
        }
        if (!why.empty()) {
            std::cout << "[fusedVR] fit.fused_variance_ratio requested but NOT used: "
                      << why << ". Falling back to the sampled variance ratio.\n";
        } else {
            fused_on = true;
            numMarkers_default = std::max(0, cfg.fused_vr_markers);
            const std::streamsize oldprec = std::cout.precision(10);
            std::cout << "[fusedVR] ON  blocks=" << fused.nblocks
                      << " (max " << fused.maxblock << ")  diag floor hits="
                      << fused.floor_hits << "  " << fused.seconds << "s\n"
                      << "[fusedVR]   tr(Psi)=" << fused.trPsi
                      << "  tr((I-H)Psi)=" << fused.trPsi_proj
                      << "  tr(P Psi)=" << fused.trPPsi
                      << "  tr(P)=" << fused.trP << "\n"
                      << "[fusedVR]   anchor = tr(P Psi)/tr((I-H)Psi) = " << fused.anchor
                      << "   (raw tr(P Psi)/tr(Psi) = " << fused.anchor_raw
                      << ", tr(P)/N = " << fused.trP_over_N << ")\n"
                      << "[fusedVR]   noXadj anchor = " << fused.anchor_noXadj
                      << "   delta budget = " << numMarkers_default
                      << " marker(s)/bin, ceiling " << cfg.fused_vr_max_markers
                      << ", se target " << cfg.fused_vr_delta_se
                      << ", keep delta at |delta-1| > " << cfg.fused_vr_delta_z << " se\n";
            std::cout.precision(oldprec);
        }
    }

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
        std::string bypass_path = []{ const char* e = std::getenv("SAIGE_BYPASS_DIR"); return (e && *e) ? std::string(e) + "/vr_marker_indices.csv" : std::string(); }();
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

    // Phase-2 wave batching: markers of one wave are selected first (the
    // selection filters — AC>=2, bin assignment — depend only on the
    // genotype, never on the PCG result, so the selection sequence is
    // identical to the historical serial loop), then their Σ⁻¹G solves run
    // as ONE batched block-PCG, then the per-marker ratio computations are
    // replayed in the original order. SAIGE_NO_BLOCKPCG=1 keeps the wave
    // structure but solves each column serially in the original order,
    // which is bit-identical to the historical per-marker loop.
    struct VrCand {
        int snpIdx;
        double AC;
        double AF;
        int binId;
        arma::fvec G;     // covariate-adjusted genotype
        arma::fvec G0f;   // raw genotype (minor-allele coded), float
    };

    // With the fused anchor and a zero marker budget there is nothing to
    // sample: the closed form IS the answer and no genotype is touched.
    const bool skip_marker_loop = (fused_on && numMarkers_default <= 0);
    if (skip_marker_loop)
        std::cout << "[fusedVR] marker budget 0 -- no PCG solve, no genotype read; "
                     "every bin takes the closed-form anchor.\n";

    while (!skip_marker_loop && ratioCV > ratioCVcutoff) {
        // ---- Phase A: select this wave's markers (no PCG) ----
        std::vector<VrCand> wave;
        while (numTestedMarker + (int)wave.size() < numMarkers0
               && indexInMarkerList < totalAvailable) {
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

            // Bin assignment depends only on MAC — check now so the wave
            // only carries markers that will actually be counted (the
            // historical loop solved PCG for binId<0 markers and then
            // discarded the result; skipping the solve changes nothing).
            int binId = bin_of_mac(AC, cfg);
            if (binId < 0) continue;

            wave.push_back({snpIdx, AC, AF, binId, std::move(G),
                            std::move(G0f)});
        }

        // ---- Phase B: Σ⁻¹G for the whole wave ----
        const int nw = (int)wave.size();
        arma::fmat Sigma_iG_mat(n, std::max(nw, 1));
        if (nw > 0) {
            if (!isBlockPCGdisabled()) {
                arma::fmat Gmat(n, nw);
                for (int c = 0; c < nw; ++c) Gmat.col(c) = wave[c].G;
                Sigma_iG_mat = getPCGofSigmaAndMatrix(W_copy, tau_copy, Gmat,
                                                      maxiterPCG, tolPCG);
            } else {
                // serial fallback: identical per-column solves in the
                // original marker order
                for (int c = 0; c < nw; ++c)
                    Sigma_iG_mat.col(c) = getPCG1ofSigmaAndVector(
                        W_copy, tau_copy, wave[c].G, maxiterPCG, tolPCG);
            }
        }

        // ---- Phase C: per-marker ratios, replayed in selection order ----
        for (int c = 0; c < nw; ++c) {
            const int    snpIdx = wave[c].snpIdx;
            const double AC     = wave[c].AC;
            const double AF     = wave[c].AF;
            const int    binId  = wave[c].binId;
            const arma::fvec& G = wave[c].G;

            // Normalized genotype
            arma::fvec g = G / std::sqrt(static_cast<float>(AC));

            // Also compute non-X-adjusted version for comparison
            arma::fvec G_noXadj = wave[c].G0f - arma::mean(wave[c].G0f);
            arma::fvec g_noXadj = G_noXadj / std::sqrt(static_cast<float>(AC));

            // --- var1 (exact): Sigma^{-1} based ---
            // R: Sigma_iG = getSigma_G(W, tauVecNew, G, ...) with the
            // covariate-adjusted G (SAIGE_fitGLMM_fast.R:2850); solved above
            // for the whole wave.
            arma::fvec Sigma_iG = Sigma_iG_mat.col(c);

            // var1a = G' * Sigma_iG - G' * Sigma_iX * (X' * Sigma_iX)^{-1} * X' * Sigma_iG
            float GtSiG = arma::dot(G, Sigma_iG);
            //   var1a = t(G) %*% Sigma_iG
            //         - t(G) %*% Sigma_iX %*% solve(t(X) %*% Sigma_iX) %*% t(X) %*% Sigma_iG
            // where Sigma_iX = Σ^{-1} X
            // So: term2 = G' * (Σ^{-1} X) * (X' Σ^{-1} X)^{-1} * X' * (Σ^{-1} G)
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
            // binId was assigned in Phase A (selection); binId<0 markers
            // never entered the wave.

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
            if (fused_on) {
                // The stopping rule has to change with the estimand. delta is the
                // bin mean divided by a CONSTANT, so sd/mean -- the CV rule used
                // below -- is numerically identical for delta and for the ratio and
                // would buy nothing. What the anchor buys is a looser accuracy
                // target, so the rule is on the relative standard error of delta,
                // which does fall as 1/sqrt(n).
                const arma::fvec& v = varRatio_NULL_vec_per_bin[b];
                if (numTestedMarker_per_bin[b] >= numTarget_per_bin[b] && v.n_elem > 0) {
                    const double m  = arma::mean(v);
                    const double sd = (v.n_elem > 1) ? arma::stddev(v) : 0.0;
                    const double rse = (m != 0.0 && v.n_elem > 1)
                                         ? sd / std::sqrt((double)v.n_elem) / std::abs(m) : 0.0;
                    ratioCV_per_bin[b] = static_cast<float>(rse);
                    if (rse <= cfg.fused_vr_delta_se
                        || numTestedMarker_per_bin[b] >= cfg.fused_vr_max_markers) {
                        binConverged[b] = true;
                        std::cout << "[fusedVR] Bin " << (b+1) << ": se(delta)/delta=" << rse
                                  << (rse <= cfg.fused_vr_delta_se ? " <= " : " > ")
                                  << cfg.fused_vr_delta_se << " using "
                                  << numTestedMarker_per_bin[b] << " markers"
                                  << (rse <= cfg.fused_vr_delta_se ? " (converged)\n"
                                                                   : " (ceiling reached)\n");
                    } else {
                        numTarget_per_bin[b] += 10;
                        std::cout << "[fusedVR] Bin " << (b+1) << ": se(delta)/delta=" << rse
                                  << " > " << cfg.fused_vr_delta_se << "; trying "
                                  << numTarget_per_bin[b] << " markers in this bin\n";
                    }
                }
                if (!binConverged[b]) all_converged = false;
                continue;
            }
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
        if (fused_on && !skip_bin) {
            // anchor + tested correction. delta is kept only when the markers
            // actually resolve it away from 1; otherwise the exact closed form
            // stands, which is both cheaper and (on the fair simulation) a
            // factor of 3 more accurate than a 30-marker sampled mean.
            const arma::fvec& v = varRatio_NULL_vec_per_bin[b];
            double sampled = bin_null[b], se = 0.0, delta = 1.0, dse = 0.0;
            if (v.n_elem > 0) {
                sampled = arma::mean(v);
                se = (v.n_elem > 1) ? arma::stddev(v) / std::sqrt((double)v.n_elem) : 0.0;
                delta = sampled / fused.anchor;
                dse   = se / fused.anchor;
            }
            const bool keep = (v.n_elem > 1) && (std::abs(delta - 1.0) > cfg.fused_vr_delta_z * dse);
            bin_null[b] = keep ? sampled : fused.anchor;
            if (varRatio_NULL_noXadj_vec_per_bin[b].n_elem > 1) {
                const arma::fvec& vn = varRatio_NULL_noXadj_vec_per_bin[b];
                const double sn = arma::mean(vn);
                const double sen = arma::stddev(vn) / std::sqrt((double)vn.n_elem);
                bin_noXadj[b] = (std::abs(sn - fused.anchor_noXadj) > cfg.fused_vr_delta_z * sen)
                                  ? sn : fused.anchor_noXadj;
            } else {
                bin_noXadj[b] = fused.anchor_noXadj;
            }
            const std::streamsize oldprec2 = std::cout.precision(10);
            std::cout << "[fusedVR] Bin " << (b+1) << ": n=" << v.n_elem
                      << " sampled=" << sampled << " +- " << se
                      << "  anchor=" << fused.anchor
                      << "  delta=" << delta << " +- " << dse
                      << "  -> " << (keep ? "keep delta (bin differs from the anchor)"
                                          : "delta set to 1 (anchor stands)")
                      << ", null=" << bin_null[b] << "\n";
            std::cout.precision(oldprec2);
            if (cfg.use_sparse_grm_for_vr && varRatio_sparse_vec_per_bin[b].n_elem == 0)
                std::cout << "[fusedVR] Bin " << (b+1) << ": WARNING no marker for the "
                             "'sparse' row; it stays at the 1.0 default.\n";
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
