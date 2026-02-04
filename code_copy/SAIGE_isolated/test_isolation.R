# Test script to verify isolated R environment
# Run with: ~/.pixi/bin/pixi run --manifest-path=pixi.toml Rscript test_isolation.R

cat("=== R Environment Isolation Test ===\n\n")

# Print R location
cat("R.home():", R.home(), "\n")
cat("R executable:", file.path(R.home("bin"), "R"), "\n")
cat(".libPaths():\n")
print(.libPaths())

# Print working directory
cat("\ngetwd():", getwd(), "\n")

# Check if SAIGE is installed
cat("\n=== Checking SAIGE installation ===\n")
if ("SAIGE" %in% installed.packages()[,"Package"]) {
    cat("SAIGE is installed at:", find.package("SAIGE"), "\n")
    library(SAIGE)
    cat("SAIGE loaded successfully!\n")
} else {
    cat("SAIGE is NOT installed yet.\n")
    cat("Please install with: R CMD INSTALL --preclean .\n")
}

cat("\n=== Test Complete ===\n")
