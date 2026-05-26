#ifndef NULL_MODEL_LOADER_HPP
#define NULL_MODEL_LOADER_HPP

#include <armadillo>
#include <string>
#include <vector>

// Struct holding all Step 1 output needed by Step 2
struct NullModelData {
    // --- Scalars (from nullmodel.json) ---
    double tau0;          // tau[0] -- variance component (=1 for binary)
    double tau1;          // tau[1] -- variance component (genetic)
    std::string traitType;   // "binary" or "quantitative"
    int n;                // sample size
    int p;                // number of covariates (including intercept)
    double SPA_Cutoff;    // default 2.0
    std::string impute_method; // "mean", "best_guess", "minor"
    bool flagSparseGRM;   // whether sparse GRM was used
    bool isFastTest;
    bool isnoadjCov;
    double pval_cutoff_for_fastTest;
    bool isCondition;
    bool is_Firth_beta;
    double pCutoffforFirth;
    std::vector<std::string> sampleIDs;

    // --- Vectors (from .arma binary files) ---
    arma::vec mu;         // [N] fitted values
    arma::vec res;        // [N] residuals (y - mu) / tau[0]
    arma::vec y;          // [N] phenotype
    arma::vec V;          // [N] variance weights: mu*(1-mu) for binary
    arma::vec mu2;        // [N] = V for binary, 1/tau[0] for quantitative
    arma::vec S_a;        // [p] score vector
    arma::vec tauvec;     // [2] tau vector
    arma::vec offset;     // [N] offset vector (for Firth)
    arma::vec resout;     // [N] residual output vector

    // --- Matrices (from .arma binary files) ---
    arma::mat X;          // [N x p] design matrix
    arma::mat XVX;        // [p x p] weighted covariate information
    arma::mat XVX_inv;    // [p x p] inverse of XVX
    arma::mat XXVX_inv;   // [N x p] = X * (X'VX)^{-1}
    arma::mat XV;         // [p x N] weighted covariates transposed
    arma::mat XVX_inv_XV; // [N x p] = (X'VX)^{-1} * X' * V
    arma::mat Sigma_iXXSigma_iX; // [p x p] or [1 x 1] dummy

    // --- Variance ratios ---
    arma::vec varRatio_sparse;
    arma::vec varRatio_null;
    arma::vec varRatio_null_noXadj;
    arma::vec cateVarRatioMinMACVecExclude;
    arma::vec cateVarRatioMaxMACVecInclude;

    // --- Sparse GRM (optional) ---
    arma::umat locationMat;  // sparse GRM location indices
    arma::vec valueVec;      // sparse GRM values
    int dimNum;              // dimension of sparse GRM (0 if not used)

    // --- Condition (optional) ---
    std::vector<uint32_t> condition_genoIndex;
};


// Load null model from Step 1 output directory
// model_dir should contain nullmodel.json and .arma files
// varianceRatio_file is the path to varianceRatio.txt
NullModelData loadNullModel(const std::string & model_dir,
                            const std::string & varianceRatio_file);

// Load a single armadillo vector from binary file
arma::vec loadArmaVec(const std::string & filepath);

// Load a single armadillo matrix from binary file
arma::mat loadArmaMat(const std::string & filepath);

// Load variance ratios from text file
// Format: each line has VR value (possibly with MAC category info)
void loadVarianceRatios(const std::string & filepath,
                        arma::vec & varRatio_null,
                        arma::vec & varRatio_sparse,
                        arma::vec & varRatio_null_noXadj,
                        arma::vec & cateVarRatioMinMACVecExclude,
                        arma::vec & cateVarRatioMaxMACVecInclude);

#endif
