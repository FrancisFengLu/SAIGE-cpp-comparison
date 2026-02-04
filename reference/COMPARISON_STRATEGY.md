# SAIGE R vs C++ Comparison Strategy

## Overview

This document outlines the step-by-step strategy for comparing SAIGE R code outputs with C++ code outputs to ensure identical behavior during the R-to-C conversion.

## How to Run R Code

Based on `/Users/francis/Desktop/Saige_temp文件/SAIGE_sparse_GRM_handoff_Jan28.txt`:

### 1. Compile R Package

```bash
cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE
rm -rf /Users/francis/Library/R/arm64/4.4/library/00LOCK-SAIGE
~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/pixi.toml R CMD INSTALL --preclean .
```

### 2. Run R Tests

```bash
# Run sparse GRM creation
~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/pixi.toml \
  Rscript /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_28_sparse_grm/R_testing/create_sparse_grm_10k.R

# Run with covariates (x1)
~/.pixi/bin/pixi run --manifest-path=/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/SAIGE/pixi.toml \
  Rscript /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_28_sparse_grm/covar_testing/run_R_128k_x1.R
```

### 3. Run C++ Code

```bash
cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_13_work/iteration_test/cpp_code
make clean && make

# Create sparse GRM only
./saige-null -c /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_28_sparse_grm/C_testing/config_create_sparse_grm_10k.yaml

# Run with covariates
./saige-null -c /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_28_sparse_grm/covar_testing/config_128k_x1.yaml
```

---

## Comparison Checkpoints

### Checkpoint 1: After Genotype Loading

**R code location:** After `setgeno()` call (SAIGE_fitGLMM_fast.R:274)

**Variables to compare:**
- `M` (number of markers)
- `N` (number of samples in genotype file)
- `Nnomissing` (number of samples with non-missing phenotype)
- `alleleFreqVec[1:10]` (first 10 allele frequencies)
- `invstdvVec[1:10]` (first 10 inverse standard deviations)

**R code to add for debugging:**
```r
# Add after line 274 in glmmkin.ai_PCG_Rcpp_Binary
M_total <- gettotalMarker()
alleleFreq <- getAlleleFreqVec()
cat("Checkpoint 1: M =", M_total, "\n")
cat("Checkpoint 1: alleleFreq[1:10] =", alleleFreq[1:10], "\n")
write.csv(data.frame(alleleFreq=alleleFreq[1:100]),
          "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint1_alleleFreq.csv")
```

**C++ code to add for debugging:**
```cpp
// Add after init_global_geno() in main.cpp
std::cout << "Checkpoint 1: M = " << geno.M << std::endl;
std::cout << "Checkpoint 1: N = " << geno.N << std::endl;
std::cout << "Checkpoint 1: Nnomissing = " << geno.Nnomissing << std::endl;
// Save to file
std::ofstream ofs("cpp_checkpoint1_alleleFreq.csv");
for(size_t i = 0; i < std::min(geno.alleleFreqVec.n_elem, (arma::uword)100); i++) {
    ofs << geno.alleleFreqVec(i) << "\n";
}
ofs.close();
```

---

### Checkpoint 2: Initial Values Before GLMM Loop

**R code location:** Lines 339-341 in glmmkin.ai_PCG_Rcpp_Binary

**Variables to compare:**
- `tau` (initial variance components)
- `alpha0` (initial fixed effects)
- `eta0` (initial linear predictor)
- `Y` (initial working response)
- `W` (initial weights)

**R code to add:**
```r
# After line 341
cat("Checkpoint 2: tau =", tau, "\n")
cat("Checkpoint 2: alpha0 =", alpha0, "\n")
cat("Checkpoint 2: sum(Y) =", sum(re.coef$Y), "\n")
cat("Checkpoint 2: sum(W) =", sum(re.coef$W), "\n")
checkpoint2 <- list(tau=tau, alpha=alpha0, Y=re.coef$Y, W=re.coef$W)
saveRDS(checkpoint2, "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint2.rds")
```

---

### Checkpoint 3: After First getCoefficients()

**R code location:** After line 341 (first `Get_Coef()` call)

**Variables to compare:**
- `Sigma_iY` (Sigma inverse times Y)
- `Sigma_iX` (Sigma inverse times X)
- `alpha` (fixed effect coefficients)
- `eta` (linear predictor)
- `cov` (covariance matrix)

**R code to add:**
```r
# After line 341
cat("Checkpoint 3: Sigma_iY[1:5] =", re.coef$Sigma_iY[1:5], "\n")
cat("Checkpoint 3: alpha =", re.coef$alpha, "\n")
cat("Checkpoint 3: cov diag =", diag(re.coef$cov), "\n")
checkpoint3 <- list(Sigma_iY=re.coef$Sigma_iY, Sigma_iX=re.coef$Sigma_iX,
                    alpha=re.coef$alpha, cov=re.coef$cov)
saveRDS(checkpoint3, "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint3.rds")
```

---

### Checkpoint 4: After getAIScore()

**R code location:** After line 342

**Variables to compare:**
- `PY` (projection of Y)
- `APY` (K times PY)
- `YPAPY` (quadratic form)
- `Trace` (trace of P*K)
- `Score` (gradient)
- `AI` (average information/Hessian)

**R code to add:**
```r
# After line 342
cat("Checkpoint 4: YPAPY =", re$YPAPY, "\n")
cat("Checkpoint 4: Trace =", re$Trace, "\n")
cat("Checkpoint 4: Score =", re$YPAPY - re$Trace, "\n")
checkpoint4 <- list(YPAPY=re$YPAPY, Trace=re$Trace)
saveRDS(checkpoint4, "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint4.rds")
```

---

### Checkpoint 5: After Iteration 1

**R code location:** After line 367 (first iteration of main loop)

**Variables to compare:**
- `tau` (updated variance components)
- `Score` / `AI` ratio
- Convergence criterion value

**R code to add:**
```r
# Inside the for loop, after fitglmmaiRPCG
cat("Checkpoint 5 (iter ", i, "): tau =", tau, "\n")
if(i == 1) {
  checkpoint5 <- list(tau=tau, iter=i)
  saveRDS(checkpoint5, "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint5.rds")
}
```

---

### Checkpoint 6: Final Results

**R code location:** After the GLMM loop (around line 420)

**Variables to compare:**
- Final `tau`
- Final `alpha`
- Final `eta`
- Convergence status
- Number of iterations

**R code to add:**
```r
# After the loop
cat("Checkpoint 6: Final tau =", tau, "\n")
cat("Checkpoint 6: Final alpha =", alpha, "\n")
cat("Checkpoint 6: Converged =", converged, "\n")
checkpoint6 <- list(tau=tau, alpha=alpha, eta=eta, converged=converged)
saveRDS(checkpoint6, "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/R_checkpoint6.rds")
```

---

### Checkpoint 7: Sparse GRM

**Variables to compare:**
- Number of non-zero elements (nnz)
- Diagonal values (first 5)
- Off-diagonal values (first 10 pairs)
- Matrix dimensions

**R code:**
```r
# After sparse GRM creation
library(Matrix)
sparseGRM <- readMM("output.sparseGRM.mtx")
cat("Checkpoint 7: dim =", dim(sparseGRM), "\n")
cat("Checkpoint 7: nnz =", nnz(sparseGRM), "\n")
cat("Checkpoint 7: diag[1:5] =", diag(sparseGRM)[1:5], "\n")
```

---

## Key Normalization Differences

### getCrossprodMatAndKin()

**R code (SAIGE_fitGLMM_fast.cpp:1733):**
```cpp
return CorssProd.m_bout/(CorssProd.Msub_mafge1perc);
```

This normalizes by the number of markers with MAF >= minMAF.

### Trace Estimation

The trace is estimated using Monte Carlo with random vectors. Due to randomness, small differences are expected.

---

## Expected Precision Differences

| Variable | Expected Difference | Reason |
|----------|---------------------|--------|
| Allele frequencies | Exact (0) | Deterministic calculation |
| Genotype values | Exact (0) | Deterministic |
| Sigma_iY, Sigma_iX | ~1e-6 | PCG convergence tolerance |
| alpha, eta | ~1e-6 | PCG-based |
| Trace | ~1e-4 | Monte Carlo randomness |
| tau | ~1e-4 | Cumulative from Trace |
| Sparse GRM kinship | ~1e-5 | Float precision |

---

## Comparison Script Template

```r
# compare_checkpoints.R
library(dplyr)

compare_vectors <- function(r_vec, cpp_vec, name, tol=1e-5) {
  if(length(r_vec) != length(cpp_vec)) {
    cat(sprintf("ERROR: %s length mismatch: R=%d, C++=%d\n", name, length(r_vec), length(cpp_vec)))
    return(FALSE)
  }
  max_diff <- max(abs(r_vec - cpp_vec))
  rel_diff <- max(abs(r_vec - cpp_vec) / (abs(r_vec) + 1e-10))

  cat(sprintf("%s: max_abs_diff=%.2e, max_rel_diff=%.2e\n", name, max_diff, rel_diff))

  if(max_diff > tol) {
    cat(sprintf("  WARNING: Difference exceeds tolerance %.2e\n", tol))
    idx <- which.max(abs(r_vec - cpp_vec))
    cat(sprintf("  Max diff at index %d: R=%.6f, C++=%.6f\n", idx, r_vec[idx], cpp_vec[idx]))
    return(FALSE)
  }
  return(TRUE)
}

# Load checkpoints
r_cp1 <- read.csv("R_checkpoint1_alleleFreq.csv")
cpp_cp1 <- read.csv("cpp_checkpoint1_alleleFreq.csv", header=FALSE)

compare_vectors(r_cp1$alleleFreq, cpp_cp1$V1, "alleleFreq")
```

---

## Output Directory

All comparison outputs should be saved to:
```
/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/
```

## Next Steps

1. Add checkpoint code to R functions
2. Add corresponding checkpoint code to C++ functions
3. Run both with same test data
4. Compare outputs at each checkpoint
5. Fix any discrepancies found
6. Document the resolution of each discrepancy
