#!/usr/bin/env Rscript
# Usage: Rscript covariate_offset_compare.R <cpp_out_prefix> [binary|quantitative]
#
# Reads <prefix>_cov_input.csv (written by saige::fit_covariate_offset), runs
# R's glm.fit on the same (X, y), dumps <prefix>_rglm_{beta,fit}.csv and prints
# a diff vs C++'s <prefix>_cov_{beta,fit}.csv.
#
# This lets us reproduce the C++ initial-GLM fit on multiple configs and
# measure exactly where C++ diverges from R's glm.fit.

suppressPackageStartupMessages({
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: covariate_offset_compare.R <cpp_out_prefix> [binary|quantitative]")
}
prefix <- args[1]
family_str <- if (length(args) >= 2) args[2] else "auto"

input_csv <- paste0(prefix, "_cov_input.csv")
beta_csv  <- paste0(prefix, "_cov_beta.csv")
fit_csv   <- paste0(prefix, "_cov_fit.csv")

if (!file.exists(input_csv)) stop("missing input CSV: ", input_csv)

dat <- read.csv(input_csv, check.names = FALSE)
p   <- sum(startsWith(names(dat), "X"))
y   <- dat$y
X   <- as.matrix(dat[, paste0("X", 0:(p-1))])

# Auto-detect family if not given
if (family_str == "auto") {
  uy <- unique(y)
  family_str <- if (length(uy) <= 2 && all(uy %in% c(0,1))) "binary" else "quantitative"
}
fam <- if (family_str == "binary") binomial() else gaussian()

cat("=== covariate_offset_compare.R ===\n")
cat("prefix : ", prefix, "\n")
cat("n      : ", nrow(X), "  p: ", ncol(X), "\n")
cat("family : ", family_str, "  (", fam$family, "/", fam$link, ")\n", sep = "")

# Run glm.fit on the same (X, y)
fit <- glm.fit(X, y, family = fam, intercept = FALSE,
               control = list(epsilon = 1e-10, maxit = 100, trace = FALSE))
cat("glm.fit iters:", fit$iter, "  converged:", fit$converged, "\n")

# Compute offset = X[, -1] %*% beta[-1]
beta_r <- unname(fit$coefficients)
eta_r  <- as.vector(X %*% beta_r)
mu_r   <- if (family_str == "binary") 1/(1+exp(-eta_r)) else eta_r
off_r  <- if (p > 1) as.vector(X[, -1, drop = FALSE] %*% beta_r[-1]) else rep(0, nrow(X))

# Write matching CSVs
write.csv(data.frame(col = 0:(p-1), beta = beta_r),
          paste0(prefix, "_rglm_beta.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(sample = 0:(nrow(X)-1), eta = eta_r, mu = mu_r, offset = off_r),
          paste0(prefix, "_rglm_fit.csv"),  row.names = FALSE, quote = FALSE)

# Diff against C++ if those CSVs exist
diff_summary <- function(label, vr, vc) {
  d <- vr - vc
  cat(sprintf("  %-8s | max|d|=%9.3e  mean|d|=%9.3e  ||d||=%9.3e  cor=%.6f\n",
              label, max(abs(d)), mean(abs(d)), sqrt(sum(d^2)), cor(vr, vc)))
}

cat("\n--- R glm.fit vs C++ covariate_offset ---\n")

if (file.exists(beta_csv)) {
  bc <- read.csv(beta_csv)
  if (nrow(bc) == p) {
    cat("beta:\n")
    for (i in seq_len(p)) {
      cat(sprintf("  beta[%d]  R=%+.6f  C++=%+.6f  diff=%+.3e\n",
                  i-1, beta_r[i], bc$beta[i], beta_r[i] - bc$beta[i]))
    }
  }
} else {
  cat("(C++ beta CSV not found: ", beta_csv, ")\n", sep = "")
}

if (file.exists(fit_csv)) {
  fc <- read.csv(fit_csv)
  cat("per-sample:\n")
  diff_summary("eta",    eta_r,  fc$eta)
  diff_summary("mu",     mu_r,   fc$mu)
  diff_summary("offset", off_r,  fc$offset)
}

cat("\nWrote:\n  ", prefix, "_rglm_beta.csv\n  ", prefix, "_rglm_fit.csv\n", sep = "")
