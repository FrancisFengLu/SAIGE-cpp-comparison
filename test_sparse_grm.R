# Test script for R SAIGE with sparse GRM
# Run with: ~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/pixi.toml Rscript test_sparse_grm.R

library(SAIGE)

# Input paths (same as C++ config)
plinkFile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/extdata/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"
phenoFile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/input/pheno_1000samples.txt"
sparseGRMFile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
sparseGRMSampleIDFile <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/extdata/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
outputPrefix <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/R_sparse_test"

cat("=== Running R SAIGE with Sparse GRM ===\n")
cat("Sparse GRM file:", sparseGRMFile, "\n")
cat("Sample IDs file:", sparseGRMSampleIDFile, "\n\n")

# Run SAIGE with sparse GRM
fitNULLGLMM(
  plinkFile = plinkFile,
  phenoFile = phenoFile,
  phenoCol = "y",
  covarColList = c("x1", "x2"),
  sampleIDColinphenoFile = "IID",
  traitType = "binary",
  outputPrefix = outputPrefix,
  useSparseGRMtoFitNULL = TRUE,
  usePCGwithSparseGRM = TRUE,  # Use PCG like C++ (instead of direct solve)
  sparseGRMFile = sparseGRMFile,
  sparseGRMSampleIDFile = sparseGRMSampleIDFile,
  LOCO = FALSE,
  nThreads = 1,
  tol = 0.02,
  tolPCG = 1e-5,
  maxiterPCG = 500,
  maxiter = 20,
  IsOverwriteVarianceRatioFile = TRUE
)

cat("\n=== R SAIGE Sparse GRM Test Complete ===\n")
