#pragma once

#include <RcppArmadillo.h>
#include <string>
#include <vector>


class genoClass ;
extern genoClass geno;

// Optional: an init function so you can set it up from main
void init_global_geno(const std::string& bed,
                      const std::string& bim,
                      const std::string& fam,
                      std::vector<int> & subSampleInGeno,
                      std::vector<bool> & indicatorGenoSamplesWithPheno,
                      bool setKinDiagtoOne,
                      double minMAFforGRM,
                      double maxMissRateforGRM);

void output_grm_diagonal(const std::string& out_path);

// Forward declaration of genoClass
// class genoClass;

// Global instance declaration
// extern genoClass geno;

// Core genotype data initialization
// void setgeno(std::string bedfile, std::string bimfile, std::string famfile, 
//              std::vector<int> & subSampleInGeno, 
//              std::vector<bool> & indicatorGenoSamplesWithPheno, 
//              float memoryChunk, bool isDiagofKinSetAsOne);

void closeGenoFile_plink();

// Genotype data access functions
arma::ivec Get_OneSNP_Geno(int SNPIdx);
arma::ivec Get_OneSNP_Geno_forVarRatio(int SNPIdx);
arma::fvec Get_OneSNP_StdGeno(int SNPIdx);
arma::fvec Get_OneSNP_StdGeno_forVarRatio(int SNPIdx);

// Marker and sample information
int gettotalMarker();
int getSubMarkerNum();
int getNnomissingOut();
arma::fvec getAlleleFreqVec();
int getMsub_MAFge_minMAFtoConstructGRM();
int getMsub_MAFge_minMAFtoConstructGRM_singleChr();
int getNColStdGenoMultiMarkersMat();
int getNRowStdGenoMultiMarkersMat();

// Variance ratio functions
bool getIsVarRatioGeno();
void setminMAC_VarianceRatio(arma::fvec t_cateVarRatioMinMACVecExclude,
                             arma::fvec t_cateVarRatioMaxMACVecInclude);
void setminMAC_VarianceRatio(float t_minMACVarRatio, float t_maxMACVarRatio,
                             bool t_isVarianceRatioinGeno);

// MAC and marker index vectors for variance ratio
arma::ivec getMACVec();
int findMainArrayIdx(int origPlinkIdx);
arma::ivec getMACVec_forVarRatio();
arma::ivec getIndexVec_forVarRatio();

// Sigma solvers (PCG-based)
arma::fmat getSigma_X(arma::fvec& wVec, arma::fvec& tauVec,
                      arma::fmat& Xmat, int maxiterPCG, float tolPCG);
arma::fvec getSigma_G(arma::fvec& wVec, arma::fvec& tauVec,
                      arma::fvec& Gvec, int maxiterPCG, float tolPCG);

// Multi-marker genotype functions
void Get_MultiMarkersBySample_StdGeno_Mat();
void Get_MultiMarkersBySample_StdGeno(arma::fvec& markerIndexVec, 
                                      std::vector<float> &stdGenoMultiMarkers);

// Sparse GRM functions
void setupSparseGRM(int r, arma::umat & locationMatinR, arma::vec & valueVecinR);
void parallelcalsparseGRM(arma::fvec &GRMvec);
void findIndiceRelatedSample();
void setRelatednessCutoff(float a);

// Kinship matrix operations
void initKinValueVecFinal(int ni);
void muliplyMailman(arma::fvec & bvec, arma::fvec & Gbvec, arma::fvec & kinbvec);
void muliplyMailman_NbyM(arma::fvec & bvec, arma::fvec & tGbvec);

// Matrix manipulation and PCG solver functions
void extractUvecforkthTime(unsigned int kthtime, arma::fvec & RvecIndex, 
                           arma::fvec& NVec, arma::fvec & sqrtDVec, arma::fvec & kthVec);
void extractVecfornthSample(unsigned int nthsample, unsigned int k_uniqTime, 
                            arma::fvec & RvecIndex, arma::fvec & sqrtWinvNVec, 
                            arma::fvec & nthVec);
void extractVecfornthSample_double(unsigned int nthsample, unsigned int k_uniqTime, 
                                   arma::vec & RvecIndex, arma::vec & sqrtWinvNVec, 
                                   arma::vec & nthVec);

// Statistical computation functions
float calCV(arma::fvec& xVec);
float GetTrace(arma::fmat Sigma_iX, arma::fmat& Xmat, arma::fvec& wVec, 
               arma::fvec& tauVec, arma::fmat& cov1, int nrun, int maxiterPCG, 
               float tolPCG, float traceCVcutoff);
float parallelInnerProduct(std::vector<float> &x, std::vector<float> &y);
float calGRMValueforSamplePair(arma::ivec &sampleidsVec);

// Configuration functions
void setminMAFforGRM(float minMAFforGRM);
void setmaxMissingRateforGRM(float maxMissingforGRM);
void setisUsePrecondM(bool isUseSparseSigmaforPCG);
void setisUseSparseSigmaforInitTau(bool isUseSparseSigmaforInitTau0);
void setisUseSparseSigmaforNullModelFitting(bool isUseSparseSigmaforModelFitting0);
void setisUsePCGwithSparseSigma(bool isUsePCGwithSparseSigma0);

// LOCO (Leave One Chromosome Out) functions
void setChromosomeIndicesforLOCO(std::vector<int> chromosomeStartIndexVec, 
                                 std::vector<int> chromosomeEndIndexVec, 
                                 std::vector<int> chromosomeVecVec);
void setStartEndIndex(int startIndex, int endIndex, int chromIndex);
void setStartEndIndexVec(arma::ivec & startIndex_vec, arma::ivec & endIndex_vec);
void set_Diagof_StdGeno_LOCO();

// Utility functions
void freqOverStd(arma::fcolvec& freqOverStdVec);
void setSubMarkerIndex(arma::ivec &subMarkerIndexRandom);
void getstdgenoVectorScalorProduct(int jth, float y, arma::fvec & prodVec);
int computePindex(arma::ivec &ithGeno);

// Matrix-matrix operations for large datasets
void sumPz(arma::fvec & Pbvec, arma::fvec & Ubvec, unsigned int mmchunksize);
void mmGetPb_MbyN(unsigned int cthchunk, unsigned int mmchunksize, 
                  arma::fvec & bvec, arma::fvec & Pbvec, arma::fvec & kinbvec);
void mmGetPb_NbyM(unsigned int cthchunk, unsigned int mmchunksize, 
                  arma::fvec & bvec, arma::fvec & Pbvec);

// Parallel processing utilities
void parallelsumTwoVec(arma::fvec &x);
void printComb(int N);

// ===== Quick hooks for sparse GRM integration (pure C++) =====
// Build sparse GRM into internal globals (locationMat/valueVec/dimNum),
// using current geno state and thresholds.
void build_sparse_grm_in_place(double relatedness_cutoff,
                               double min_maf, double max_miss);

// Export the assembled sparse GRM as COO for writing to disk.
arma::umat export_sparse_grm_locations(); // 2 x nnz (0-based)
arma::vec  export_sparse_grm_values();    // nnz
int        export_sparse_grm_dim();       // n
int        export_sparse_grm_nnz();       // nnz convenience

// Sparse switches (for verification / logging)
void setisUseSparseSigmaforInitTau(bool isUseSparseSigmaforInitTau0);
void setisUseSparseSigmaforNullModelFitting(bool isUseSparseSigmaforModelFitting0);
bool get_isUseSparseSigmaforModelFitting();

arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w,
                                   const arma::fvec& tau,
                                   const arma::fvec& v,
                                   int maxiterPCG, float tolPCG);

arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& w,
                                        const arma::fvec& tau,
                                        const arma::fvec& v,
                                        int maxiterPCG, float tolPCG);

// ---- Legacy symbols expected by old call sites (non-const refs) ----
// arma::fvec getPCG1ofSigmaAndVector(arma::fvec& w,
//                                    arma::fvec& tau,
//                                    arma::fvec& v,
//                                    int maxiterPCG, float tolPCG)
// {
//     // bind to const refs so overload resolution picks the const version
//     const arma::fvec& wc = w;
//     const arma::fvec& tc = tau;
//     const arma::fvec& vc = v;
//     return getPCG1ofSigmaAndVector_mother(wc, tc, vc, maxiterPCG, tolPCG);
// }

// arma::fvec getPCG1ofSigmaAndVector_LOCO(arma::fvec& w,
//                                         arma::fvec& tau,
//                                         arma::fvec& v,
//                                         int maxiterPCG, float tolPCG)
// {
//     const arma::fvec& wc = w;
//     const arma::fvec& tc = tau;
//     const arma::fvec& vc = v;
//     return getPCG1ofSigmaAndVector_LOCO_mother(wc, tc, vc, maxiterPCG, tolPCG);
// }