# Reference LOCO run with the R SAIGE package (code_copy/SAIGE_isolated build).
#
# Run this BEFORE the C++ LOCO run from the same working directory: R writes the
# random-vector bypass file that the C++ GetTrace() reads, so both sides use the
# same Hutchinson probes and therefore land on the same tau. Without that, tau
# differs by ~5% and per-chromosome mu cannot be compared meaningfully.
#
# Usage:
#   Rscript r_reference_loco.R <plinkPrefix> <phenoFile> <outPrefix> <phenoCol> <traitType>
#
# Writes <outPrefix>_LOCO_mu.csv with one column per autosome plus a "full"
# column (the non-LOCO fitted values).

suppressPackageStartupMessages(library(SAIGE))

args <- commandArgs(trailingOnly = TRUE)
plinkFile  <- args[1]
phenoFile  <- args[2]
outPrefix  <- args[3]
phenoCol   <- if (length(args) >= 4) args[4] else "y_binary"
traitType  <- if (length(args) >= 5) args[5] else "binary"

fitNULLGLMM(
  plinkFile               = plinkFile,
  phenoFile               = phenoFile,
  phenoCol                = phenoCol,
  covarColList            = c("x1", "x2"),
  sampleIDColinphenoFile  = "IID",
  traitType               = traitType,
  outputPrefix            = outPrefix,
  useSparseGRMtoFitNULL   = FALSE,
  LOCO                    = TRUE,
  nThreads                = 1,
  tol                     = 0.02,
  tolPCG                  = 1e-5,
  maxiterPCG              = 500,
  maxiter                 = 20,
  isCovariateOffset       = TRUE,
  IsOverwriteVarianceRatioFile = TRUE
)

obj <- readRDS(paste0(outPrefix, ".rda"))
cat("R theta:", obj$theta, "\n")
cat("R LOCO:", obj$LOCO, "\n")

out <- data.frame(full = as.numeric(obj$fitted.values))
if (!is.null(obj$LOCOResult)) {
  for (j in 1:22) {
    r <- obj$LOCOResult[[j]]
    if (!is.null(r) && isTRUE(r$isLOCO)) {
      out[[paste0("chr", j)]] <- as.numeric(r$fitted.values)
    }
  }
}
write.csv(out, paste0(outPrefix, "_LOCO_mu.csv"), row.names = FALSE)
cat("wrote ", paste0(outPrefix, "_LOCO_mu.csv"), "\n")
