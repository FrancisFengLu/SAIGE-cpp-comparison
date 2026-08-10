#!/usr/bin/env Rscript
# One R SAIGE SPAGMMATtest run, for timing.  Fresh process per run (calling
# SPAGMMATtest twice in one session empties the second output -- see
# test_data/loco_reference/run_step2.R).
#
#   bench_r_step2.R single <nullPrefix> <bedPrefix> <out>
#   bench_r_step2.R region <nullPrefix> <bedPrefix> <out> <groupFile>
args   <- commandArgs(trailingOnly = TRUE)
MODE   <- args[1]; PREFIX <- args[2]; BED <- args[3]; OUT <- args[4]
suppressPackageStartupMessages(library(SAIGE))

common <- list(
  bedFile               = paste0(BED, ".bed"),
  bimFile               = paste0(BED, ".bim"),
  famFile               = paste0(BED, ".fam"),
  AlleleOrder           = "alt-first",
  min_MAC               = 0.5,
  min_MAF               = 0,
  max_missing           = 0.15,
  GMMATmodelFile        = paste0(PREFIX, ".rda"),
  varianceRatioFile     = paste0(PREFIX, ".varianceRatio.txt"),
  LOCO                  = FALSE,
  is_output_moreDetails = TRUE,
  SAIGEOutputFile       = OUT
)

if (MODE == "region") {
  common <- c(common, list(
    groupFile                      = args[5],
    annotation_in_groupTest        = c("lof", "missense;lof", "missense;lof;synonymous"),
    maxMAF_in_groupTest            = c(0.0001, 0.001, 0.01),
    r.corr                         = 0,
    MACCutoff_to_CollapseUltraRare = 10,
    markers_per_chunk_in_groupTest = 500,
    is_single_in_groupTest         = FALSE,
    is_output_markerList_in_groupTest = FALSE,
    weights.beta                   = c(1, 25)))
}

t0 <- proc.time()
do.call(SPAGMMATtest, common)
el <- proc.time() - t0
cat("BENCH_R_STEP2_ELAPSED_S ", el[["elapsed"]], "\n", sep = "")
