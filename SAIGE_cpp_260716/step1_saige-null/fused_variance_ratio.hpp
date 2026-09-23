// fused_variance_ratio.hpp
// ------------------------------------------------------------------
// Closed-form anchor for the step-1 variance ratio on the
// quantitative + sparse-GRM path.
//
//   VR = E[G' P G] / E[G' G]      G = (I-H) G0 ,  Cov(G0) proportional to K
//
// P X = 0 and H = X(X'X)^-1 X' give (I-H) P (I-H) = P, so with K taken as the
// sparse kinship Psi the ratio is
//
//   anchor = tr(P Psi) / tr((I-H) Psi)
//
// which is the TorchGWAS2 factor 1/c1 = tr(P Psi)/tr(Psi) with the denominator
// corrected for the p covariates the tested genotype is projected against.
// See optimization/FUSED_VARIANCE_RATIO.md and
// optimization/torchgwas2/FAIR_SIM_RELATEDNESS.md.
//
// Every trace here is EXACT, not Hutchinson: on this path
// Sigma = tau0*diag(1/W) + tau1*Psi is block diagonal by the connected
// components of Psi, and biobank sparse GRMs give components of a handful of
// people, so each block is inverted densely.  The 30-probe Hutchinson traces
// the AI-REML step already computes carry ~0.1% noise, an order of magnitude
// worse than the 30-marker sampled ratio they would be replacing.
// ------------------------------------------------------------------
#pragma once
#include <string>
#include <armadillo>

namespace saige {

struct FusedVrAnchor {
    bool        ok{false};
    std::string why;            // why the anchor could not be built (ok == false)

    double trPsi{0};            // tr(Psi)
    double trPsi_proj{0};       // tr((I-H) Psi)          <- the denominator used
    double trPPsi{0};           // tr(P Psi)
    double trSigmaInvPsi{0};    // tr(Sigma^-1 Psi)       (before the covariate term)
    double trP{0};              // tr(P)                  (diagnostic only)
    double trSigmaInv{0};       // tr(Sigma^-1)           (diagnostic only)

    // Binary denominator. var2 for a binary trait is g~' W g~ with the
    // W-weighted covariate projection M = I - X(X'WX)^-1 X'W, and because
    // X'WM = 0 the expansion of tr(W M Psi M') collapses:
    //     M'WM = WM = W - W X (X'WX)^-1 X' W
    // so two of the three covariate terms cancel and one is left.
    //     tr(W Psi) - tr((X'WX)^-1 X'W Psi W X)
    // With W == 1 this is literally tr((I-H) Psi), i.e. the quantitative
    // denominator above -- one formula, not two.
    double trWPsi{0};           // tr(W Psi)
    double trPsi_Wproj{0};      // tr(W Psi) - tr((X'WX)^-1 X'W Psi W X)
    double anchor_binary{0};    // tr(P Psi) / tr(W Psi - ...)

    double trPsi_mean{0};       // tr((I - 11'/n) Psi)    denominator of the noXadj row
    double anchor{0};           // tr(P Psi) / tr((I-H) Psi)     <- what we use
    double anchor_noXadj{0};    // tr(P Psi) / tr((I - 11'/n) Psi)
    double anchor_raw{0};       // tr(P Psi) / tr(Psi)           = TorchGWAS2's 1/c1
    double trP_over_N{0};       // tr(P)/N, the identity-weighted closed form

    int    n{0};
    int    nblocks{0};
    int    maxblock{0};
    int    floor_hits{0};       // diagonal entries clamped at 1e-4 by gen_sp_Sigma
    double seconds{0};
};

// Builds the anchor from the sparse GRM currently installed in
// SAIGE_step1_fast.cpp's globals (locationMat / valueVec / dimNum).
//   W    working weights (quantitative: all ones)
//   tau  [tau0, tau1]
//   X    n x p design matrix
//   max_block  refuse if any connected component is larger than this
FusedVrAnchor compute_fused_vr_anchor(const arma::vec& W,
                                      const arma::vec& tau,
                                      const arma::mat& X,
                                      int max_block);

} // namespace saige
