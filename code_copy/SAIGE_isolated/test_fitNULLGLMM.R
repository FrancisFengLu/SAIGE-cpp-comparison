# Test script to run fitNULLGLMM and generate checkpoint outputs
# Run with: ~/.pixi/bin/pixi run --manifest-path=pixi.toml Rscript test_fitNULLGLMM.R

cat("=== SAIGE fitNULLGLMM Test with Checkpoints ===\n\n")

# Load SAIGE
library(SAIGE)
cat("SAIGE loaded successfully from:", find.package("SAIGE"), "\n\n")

# Set paths - use the source extdata directory, not the installed package's
extdata_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/code_copy/SAIGE_isolated/extdata/input"
output_dir <- "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)

cat("Input data directory:", extdata_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# Use the larger dataset with 10k markers and 1000 samples
plink_prefix <- file.path(extdata_dir, "plinkforGRM_1000samples_10kMarkers")
pheno_file <- file.path(extdata_dir, "pheno_1000samples.txt")

cat("Checking input files:\n")
cat("  PLINK prefix:", plink_prefix, "\n")
cat("  BED file exists:", file.exists(paste0(plink_prefix, ".bed")), "\n")
cat("  BIM file exists:", file.exists(paste0(plink_prefix, ".bim")), "\n")
cat("  FAM file exists:", file.exists(paste0(plink_prefix, ".fam")), "\n")
cat("  Phenotype file:", pheno_file, "exists:", file.exists(pheno_file), "\n\n")

# Run fitNULLGLMM
cat("Running fitNULLGLMM...\n")
cat("This will generate checkpoint files in:", file.path(output_dir, "checkpoints"), "\n\n")

fitNULLGLMM(
  plinkFile = plink_prefix,
  phenoFile = pheno_file,
  phenoCol = "y",
  covarColList = c("x1", "x2"),
  sampleIDColinphenoFile = "IID",
  traitType = "binary",
  outputPrefix = file.path(output_dir, "test_nullmodel"),
  nThreads = 1,
  IsOverwriteVarianceRatioFile = TRUE
)

cat("\n=== fitNULLGLMM Complete ===\n\n")

# Check checkpoint files
checkpoint_dir <- file.path(output_dir, "checkpoints")
cat("Checkpoint files created:\n")
checkpoint_files <- list.files(checkpoint_dir, pattern = "R_CP.*\\.rds", full.names = TRUE)
if (length(checkpoint_files) > 0) {
  for (f in checkpoint_files) {
    cat("  ", basename(f), "- size:", file.info(f)$size, "bytes\n")
  }
} else {
  cat("  No checkpoint files found!\n")
}

cat("\n=== Test Complete ===\n")
