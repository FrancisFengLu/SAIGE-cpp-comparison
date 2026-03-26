# Corner Case & Input Validation Coverage
## R (`SAIGE_fitGLMM_fast.R`) vs C++ Standalone

---

## 1. Trait Type Validation

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `eventTimeCol != ""` requires `traitType == "survival"` | 1235 | **Partial** — trait type validated, survival not implemented | `null_model_engine.cpp:493` throws on unsupported trait type |
| `eventTime` must not be NULL (survival) | 1743–1744 | **No** | Survival not implemented in C++ |
| `eventTime` must be >= 0 | 1746–1748 | **No** | Survival not implemented |
| Binary phenotype must be only 0 or 1 | 1754–1757 | **Partial** | No explicit 0/1 check; R handles this; C++ reads whatever is in the design CSV |

---

## 2. File Existence Validation

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| PLINK BED file exists | 1344–1345 | **Yes** | `preprocess_engine.cpp:67` (FAM open check); `SAIGE_step1_fast.cpp` BED open |
| PLINK BIM file exists | 1347–1348 | **Yes** | `main.cpp:521`, `preprocess_engine.cpp:226` |
| PLINK FAM file exists | 1364–1365 | **Yes** | `preprocess_engine.cpp:67` |
| Phenotype/design file exists | 1391–1392 | **Yes** | `main.cpp:148`, `null_model_engine.cpp:84` |
| Sparse GRM file exists | 1314–1315 | **Yes** | `main.cpp:433` Matrix Market open check |
| Sparse GRM sample ID file exists | 1317–1319 | **Yes** | `preprocess_engine.cpp:82` |
| Sample ID include file exists | 1433–1434 | **Yes** | `main.cpp:1010` |
| Variance ratio file overwrite check | 1284–1295 | **No** | C++ always overwrites output; no overwrite guard |
| `skipModelFitting` but model file missing | 1253–1254 | **No** | C++ has no `skipModelFitting` mode |

---

## 3. Column Existence in Data

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `phenoCol` exists in phenotype file | 1404–1407 | **Yes** | `main.cpp:671` throws on missing y_col |
| `covarColList` columns exist in file | 1404–1407 | **Yes** | `main.cpp:685` throws on missing covariate column |
| `sampleIDColinphenoFile` exists | 1404–1407 | **Yes** | `main.cpp:674` throws on missing iid_col |
| `eventTimeCol` exists when specified | 1412–1413 | **No** | Survival not implemented |
| Categorical covariate columns all in `covarColList` | 1446–1454 | **No** | C++ has no `qCovarCol` (categorical covariate) concept |

---

## 4. Sex-Specific Filtering

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Both `FemaleOnly` and `MaleOnly` cannot be TRUE | 1329–1330 | **No** | C++ has no sex-specific filtering |
| FemaleOnly=TRUE but no female-coded samples | 1466–1469 | **No** | Not implemented |
| MaleOnly=TRUE but no male-coded samples | 1475–1477 | **No** | Not implemented |

---

## 5. Sample Size & Completeness

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Remove rows with any missing value (`complete.cases`) | 1430 | **Yes** | `main.cpp:718` skips rows where phenotype is empty/NA/NaN; `main.cpp:763` NaN imputation for covariates |
| No samples left after filtering | 1540–1541 | **Yes** | `preprocess_engine.cpp:109` throws if `keep_pos` is empty |
| Duplicate sample IDs removed | 1437 | **No** | C++ assumes IDs are unique; no dedup step |
| Samples in sparse GRM must match model samples | 2785–2802 | **Yes** | `preprocess_engine.cpp:109` intersection check |

---

## 6. Variance Component (Tau) Checks

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Phenotype variance < 0.1 (quantitative) | 879–880 | **No** | No pre-check on phenotype variance; solver will proceed |
| `tau[2] == 0` → break binary loop | 639 | **Partial** | `glmm.cpp:1073–1074` clamps tau ≥ 0 but doesn't break on zero |
| `tau[1] <= 0 \|\| tau[2] <= 0` → break quantitative | 941 | **Partial** | `glmm.cpp:1073–1074` clamps; no explicit break-on-zero |
| Tau convergence check | 944, 641 | **Yes** | `glmm.cpp:722` `rc_tau < tol_coef` → break |
| `max(tau) > tol^(-2)` → warning + force break | 643–646, 947–950 | **No** | No upper bound on tau; solver will keep running |
| Step-halving when tau_new < 0 (binary) | (implicit) | **Yes** | `glmm.cpp:689–694` step-halving while loop |
| Step-halving when tau_new < 0 (quantitative) | (implicit) | **Yes** | `glmm.cpp:1073–1074` clamp + step-halving |

---

## 7. Fixed-Effect Alpha (Coefficient) Convergence

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `max(abs(α - α0)/(abs(α)+abs(α0)+tol)) < tol` | 167–168 | **Yes** | `glmm.cpp:501` `rc_alpha_inner < tol_coef` → break inner IRLS |
| Same check in LOCO variant | 225 | **Partial** | LOCO not implemented as R separate path; shared solver used |

---

## 8. Numerical Stability / Float Checks

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `mu` clamped away from 0/1 (logistic) | (implicit in IRLS) | **Yes** | `null_model_engine.cpp:199` clamp [1e-12, 1-1e-12]; `SAIGE_step1_fast.cpp:8207,8238` clamp [1e-8, 1-1e-8] |
| `eta` overflow prevention (sigmoid) | (implicit) | **Yes** | `glmm.cpp:202` `arma::clamp(eta, -40, 40)` |
| `mu_eta` denominator floor | (implicit) | **Yes** | `glmm.cpp:179` add 1e-20f floor |
| Non-finite values in trace computation | (implicit) | **Yes** | `SAIGE_step1_fast.cpp:5543` throws on non-finite inputs |
| Non-finite entries in trace result vectors | (implicit) | **Yes** | `SAIGE_step1_fast.cpp:5624,5628,5632` |
| Column all-ones detection for intercept | (implicit) | **Yes** | `main.cpp:794` `!std::isfinite(v)` check + all-one flag |

---

## 9. PCG Solver Checks

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| PCG convergence (while loop condition) | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4087` `sumr2 > tolPCG && iter < maxiterPCG` |
| PCG exceeded max iterations warning | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4113–4116` prints warning |
| Input vector `bVec` non-empty | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4045` throws |
| `wVec` length matches `bVec` | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4046` throws |
| `tauVec` has at least 2 elements | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4047` throws |
| Solve result has correct length | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:4078,4103` throws |
| Symmetric positive definite inversion fallback | (Rcpp C++) | **Yes** | `saige_ai.cpp:9–15` `inv_sympd` → `pinv` fallback |

---

## 10. GetTrace / Random Vector Checks

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Adaptive trace: add vectors if CV > cutoff | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:5665` CV check |
| Matrix dimensions match for trace computation | (Rcpp C++) | **Yes** | `SAIGE_step1_fast.cpp:5536–5540` throws |
| Preloaded bypass vector pool exhausted | N/A | **Yes** | `SAIGE_step1_fast.cpp:5472–5475` throws |
| Bypass file unavailable → fallback RNG | N/A | **Yes** | `SAIGE_step1_fast.cpp:5430` cerr + fallback |

---

## 11. Variance Ratio Computation

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| No markers found for VR estimation | 3021 | **No** | C++ will error from empty indices rather than a clear message |
| Insufficient markers per MAC category | 3055–3057 | **No** | No explicit per-category minimum check |
| Categorical VR: insufficient markers in highest MAC | 3062–3063 | **No** | Not checked |
| VR value must be > 0 | (implicit) | **Yes** | `variance_ratio_engine.cpp:109` throws if `!(vr > 0.0)` |
| VR file has no data rows | N/A | **Yes** | `variance_ratio_engine.cpp:112` throws |
| VR file missing expected MAC categories | N/A | **Yes** | `variance_ratio_engine.cpp:116` throws |
| Allele frequency > 0.5 → flip genotype | 3117–3118 | **Yes** | `SAIGE_step1_fast.cpp:488–489` MAF = min(af, 1-af) |
| `use_r_vr_bypass` must be FALSE in production | N/A | **⚠️ TODO** | `variance_ratio_compute.cpp:~55` currently `true` |
| Marker index not found in main array | N/A | **Yes** | `variance_ratio_compute.cpp:288` warning + skip |

---

## 12. Sparse GRM Consistency

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| GRM dimensions match sample count | 2785–2786 | **Partial** | `preprocess_engine.cpp:109` checks sample overlap; no explicit square-matrix dimension assert |
| Sample mismatch between model and GRM | 2799–2802 | **Yes** | `preprocess_engine.cpp:109` throws if intersection is empty |
| Sparse GRM + no plink → skip VR | 1271–1272 | **Yes** | `main.cpp` config logic: `skipVarianceRatioEstimation` path |

---

## 13. Thread & Computational Configuration

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `useSparseGRMtoFitNULL=TRUE` forces `nThreads=1` | 1305 | **Yes** | `main.cpp:929` cerr note; forces nthreads=1 |
| nThreads > 1 → configure parallel | 1324–1326 | **Yes** | Armadillo + OpenBLAS thread config at startup |

---

## 14. Dimension / Shape Guards on Core Matrices

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Design matrix shape consistent with n, p | (implicit) | **Yes** | `glmm.cpp:404–408` `check_dims()` throws on X/y/offset/beta_init shape mismatch |
| Coefficient output shape consistent | (implicit) | **Yes** | `glmm.cpp:521–525` `check_coef()` throws |
| Y, W shape before AI score | (implicit) | **Yes** | `glmm.cpp:574–575` throws |
| Design buffer sizes in null model engine | (implicit) | **Yes** | `null_model_engine.cpp:41,49,328–331`, `preprocess_engine.cpp:18–23` |

---

## 15. Missing Data / NA Imputation for Genotypes

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Missing genotype imputation (mean / bestguess / none) | (Rcpp C++) | **Yes** | `UTIL.cpp:46,110` imputation strategies; MAC cutoff cleaning |

---

## 16. Convergence Diagnostics

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| Set `converged` flag after iteration loop | 668 | **No** | C++ prints convergence messages to stdout but does not expose a `converged` field in output |
| `STOP_AFTER_FIRST_TAU_UPDATE` debug early exit | 585–619 | **No** | No equivalent debug flag |

---

## 17. Initialization Guards

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| `tauInit[fixtau==0] == 0` → default init | 352 | **Yes** | `glmm.cpp:63` zero-fills beta_init if empty |
| `is.null(offset)` → zero offset | 850–851 | **Yes** | `glmm.cpp:49` fills zero offset if empty |
| Conservative tau update on first iteration only | 639–641 | **Yes** | `glmm.cpp:676–685` `it==0` branch |
| `sum(tauInit[fixtau==0]) == 0` quantitative init | 875–883 | **Partial** | C++ initializes tau=[1,0] for binary, [1,0] for quant; no variance-based dynamic init |

---

## 18. Log-Odds / Probability Bounds

| R Check | R Line | C++ Has It? | C++ Location |
|---------|--------|-------------|--------------|
| p ≤ 0 or p ≥ 1 → ±Inf for logit | (implicit) | **Yes** | `main.cpp:487–488` returns ±INFINITY |

---

## Summary

| Category | R Checks | C++ Has | C++ Missing |
|----------|----------|---------|-------------|
| Trait type / phenotype value | 4 | 1 | 3 (survival trait, binary 0/1 check) |
| File existence | 9 | 7 | 2 (overwrite guard, skipModelFitting) |
| Column existence | 5 | 3 | 2 (eventTimeCol, qCovarCol) |
| Sex-specific filtering | 3 | 0 | 3 (feature not implemented) |
| Sample size / completeness | 4 | 3 | 1 (duplicate IDs) |
| Tau / variance component | 7 | 4 | 3 (phenotype variance < 0.1, upper bound warning, break-on-zero) |
| Alpha convergence | 2 | 1 | 1 (LOCO path) |
| Numerical stability | 6 | 6 | 0 |
| PCG solver | 7 | 7 | 0 |
| GetTrace / random vectors | 4 | 4 | 0 |
| Variance ratio | 6 | 3 | 3 (marker count minimums; bypass TODO) |
| Sparse GRM consistency | 3 | 2 | 1 (exact dimension assert) |
| Thread config | 2 | 2 | 0 |
| Core matrix dimensions | 4 | 4 | 0 |
| Genotype missing data | 1 | 1 | 0 |
| Convergence diagnostics | 2 | 0 | 2 |
| Initialization guards | 4 | 2 | 2 |
| Log-odds bounds | 1 | 1 | 0 |
| **Total** | **74** | **51** | **23** |

---

## Notable Gaps (Risk-Ranked)

| Priority | Gap | Risk | Notes |
|----------|-----|------|-------|
| **HIGH** | `use_r_vr_bypass = true` left on | Correctness | `variance_ratio_compute.cpp:~55` — must flip to `false` for production |
| **HIGH** | No phenotype variance < 0.1 check (quantitative) | Silent bad model | R line 879–880 stops with clear message; C++ will proceed with degenerate input |
| **HIGH** | No binary phenotype 0/1 validation | Silent wrong result | C++ reads whatever floats are in the CSV column |
| **MEDIUM** | No `max(tau) > tol^-2` warning | Non-convergence loops forever | R breaks + warns; C++ iterates to maxiter silently |
| **MEDIUM** | No duplicate sample ID removal | Undefined behavior on index lookup | C++ assumes unique IDs |
| **MEDIUM** | No minimum marker count check for VR MAC categories | Runtime crash or empty VR | R gives a clear stop() message |
| **LOW** | No `converged` flag in output | Downstream callers cannot detect non-convergence | R returns `converged=TRUE/FALSE` |
| **LOW** | Sparse GRM no exact dimension assert | Subtle sample mismatch | C++ checks overlap but not exact square dimensions |
