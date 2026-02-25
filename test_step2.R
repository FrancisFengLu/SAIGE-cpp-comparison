library(SAIGE)

cat("=== SAIGE Step 2 Test Run ===\n\n")

# Base path
base <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/extdata"

# Step 1 outputs (model + variance ratio)
model_file <- file.path(base, "output", "example_binary_positive_signal.rda")
vr_file    <- file.path(base, "output", "example_binary_positive_signal.varianceRatio.txt")

# Step 2 genotype input (PLINK format)
bed_file <- file.path(base, "input", "genotype_100markers.bed")
bim_file <- file.path(base, "input", "genotype_100markers.bim")
fam_file <- file.path(base, "input", "genotype_100markers.fam")

# Verify all input files exist
cat("Checking input files:\n")
for (f in c(model_file, vr_file, bed_file, bim_file, fam_file)) {
  exists <- file.exists(f)
  cat(sprintf("  [%s] %s\n", ifelse(exists, "OK", "MISSING"), basename(f)))
  if (!exists) stop(paste("Missing file:", f))
}

# Output file
out_file <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/test_step2_results.txt"

cat("\nRunning Step 2 (PLINK format, chr 1)...\n")

SPAGMMATtest(
  bedFile         = bed_file,
  bimFile         = bim_file,
  famFile         = fam_file,
  AlleleOrder     = "alt-first",
  chrom           = "1",
  GMMATmodelFile  = model_file,
  varianceRatioFile = vr_file,
  SAIGEOutputFile = out_file,
  min_MAF         = 0,
  min_MAC         = 0.5,
  LOCO            = FALSE,
  is_Firth_beta   = TRUE,
  pCutoffforFirth = 0.05,
  is_output_moreDetails = TRUE,
  is_overwrite_output = TRUE
)

cat("\n=== Step 2 Complete ===\n")
cat("Output file:", out_file, "\n\n")

# Read and display results
results <- read.table(out_file, header = TRUE, sep = "\t")
cat(sprintf("Total variants tested: %d\n", nrow(results)))
cat(sprintf("Significant (p < 0.05): %d\n", sum(results$p.value < 0.05, na.rm = TRUE)))
cat(sprintf("Significant (p < 5e-8): %d\n", sum(results$p.value < 5e-8, na.rm = TRUE)))
cat("\nTop 10 results by p-value:\n")
top <- head(results[order(results$p.value), ], 10)
print(top[, c("CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value")])

# Compare with reference output
ref_file <- file.path(base, "output", "example_binary_positive_signal.assoc.step2.txt")
if (file.exists(ref_file)) {
  cat("\n=== Comparing with reference output ===\n")
  ref <- read.table(ref_file, header = TRUE, sep = "\t")
  cat(sprintf("Reference has %d variants\n", nrow(ref)))
  # Find overlapping markers
  common <- merge(ref, results, by = "MarkerID", suffixes = c(".ref", ".new"))
  if (nrow(common) > 0) {
    cat(sprintf("Overlapping markers: %d\n", nrow(common)))
    for (i in 1:nrow(common)) {
      cat(sprintf("  %s: p.ref=%.6e, p.new=%.6e, BETA.ref=%.6f, BETA.new=%.6f\n",
                  common$MarkerID[i],
                  common$p.value.ref[i], common$p.value.new[i],
                  common$BETA.ref[i], common$BETA.new[i]))
    }
  }
}
