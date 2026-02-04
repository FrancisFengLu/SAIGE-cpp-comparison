#!/usr/bin/env Rscript

# SAIGE R package dependencies installation script
# Based on DESCRIPTION file and pixi.toml

cat("Installing SAIGE R package dependencies...\n")

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# BiocManager for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Core dependencies from DESCRIPTION
packages_cran <- c(
  "Rcpp",
  "RcppParallel",
  "Matrix",
  "data.table",
  "RcppArmadillo",
  "RcppEigen",
  "methods",
  "BH",
  "optparse",
  "SKAT",
  "MetaSKAT",
  "qlcMatrix",
  "RhpcBLASctl",
  "RSQLite",
  "dplyr",
  "lintools",
  "devtools"
)

# Install CRAN packages
for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installing ", pkg, "...\n"))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(paste0(pkg, " is already installed.\n"))
  }
}

# Install SPAtest from GitHub (specific version 3.1.2)
cat("Installing SPAtest 3.1.2 from GitHub...\n")
if (!requireNamespace("SPAtest", quietly = TRUE) || packageVersion("SPAtest") != "3.1.2") {
  devtools::install_github("leeshawn/SPAtest", ref = "3.1.2", upgrade = "never")
} else {
  cat("SPAtest 3.1.2 is already installed.\n")
}

cat("\nAll R package dependencies have been installed successfully!\n")
cat("You can now proceed with: R CMD INSTALL SAIGE\n")
