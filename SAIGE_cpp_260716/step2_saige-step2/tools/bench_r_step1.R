#!/usr/bin/env Rscript
# R SAIGE step 1 (fitNULLGLMM), no instrumentation, for timing.
#   bench_r_step1.R <plinkPrefix> <phenoFile> <phenoCol> <traitType> <outPrefix> <nThreads>
# Argument values mirror test_data/loco_reference/run_reference.R except that
# LOCO is off (single-shot fit) so the R and C++ runs do identical work.
args <- commandArgs(trailingOnly = TRUE)
PLINK <- args[1]; PHENO <- args[2]; PCOL <- args[3]
TRAIT <- args[4]; PREFIX <- args[5]; NTHR <- as.integer(args[6])
suppressPackageStartupMessages(library(SAIGE))
t0 <- proc.time()
fitNULLGLMM(
  plinkFile                    = PLINK,
  phenoFile                    = PHENO,
  phenoCol                     = PCOL,
  covarColList                 = c("x1", "x2"),
  sampleIDColinphenoFile       = "IID",
  traitType                    = TRAIT,
  invNormalize                 = FALSE,
  outputPrefix                 = PREFIX,
  nThreads                     = NTHR,
  LOCO                         = FALSE,
  isLowMemLOCO                 = FALSE,
  minMAFforGRM                 = 0.01,
  tol                          = 0.02,
  tolPCG                       = 1e-5,
  maxiter                      = 20,
  maxiterPCG                   = 500,
  nrun                         = 30,
  numMarkersForVarRatio        = 30,
  skipVarianceRatioEstimation  = FALSE,
  IsOverwriteVarianceRatioFile = TRUE
)
el <- proc.time() - t0
load(paste0(PREFIX, ".rda"))
cat("BENCH_R_STEP1_ELAPSED_S ", el[["elapsed"]], "\n", sep = "")
cat("BENCH_R_STEP1_TAU ", paste(format(modglmm$theta, digits = 15), collapse = " "),
    "\n", sep = "")
cat("BENCH_R_STEP1_N ", length(modglmm$sampleID), "\n", sep = "")
