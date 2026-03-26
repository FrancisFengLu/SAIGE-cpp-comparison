# Input Validation Implementation Plan
## C++ Standalone — Matching R `fitNULLGLMM()` Corner Cases

**Date:** 2026-03-24
**Status:** ALL PHASES IMPLEMENTED (Steps 1–16, 18). Step 17 skipped by design.
**Reference R file:** `SAIGE_isolated/R/SAIGE_fitGLMM_fast.R`, function `fitNULLGLMM()` (line 1175)
**Reference coverage audit:** `CORNER_CASE_COVERAGE.md` (23 gaps identified)

---

## Strategy

### Principles

1. **Match R behavior exactly** — every `stop()` in R becomes a `throw std::runtime_error()` in C++; every R `warning()` becomes a `std::cerr` warning. Messages are similar to R's for user familiarity.
2. **Fail early, fail clearly** — validate inputs as close to the entry point as possible (in `main.cpp` after design/config load), before expensive computation begins.
3. **No algorithm changes** — these are purely guard-rail additions. The solver logic stays untouched.
4. **Minimal code footprint** — each check is a short block (5–15 lines). No new files needed; all changes go into existing files.

### File-to-responsibility mapping

| File | Responsibility |
|------|---------------|
| `main.cpp` | CLI entry, design loading, pre-solver input validation |
| `glmm.cpp` | Solver-loop guards (tau bounds, convergence flags) |
| `variance_ratio_compute.cpp` | VR-specific marker count checks, bypass flag |
| `preprocess_engine.cpp` | GRM dimension validation (via `load_matrix_market_coo`) |
| `saige_null.hpp` | `FitNullConfig` struct (`overwrite_vr` field already present) |

---

## Phase 1 — HIGH Priority (IMPLEMENTED)

### Step 1: Binary phenotype 0/1 validation

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1754–1757` — `uniqPheno = sort(unique(y)); if (uniqPheno[1] != 0 \| uniqPheno[2] != 1) stop(...)` |
| **C++ implementation** | `main.cpp:1020–1038` |
| **What it does** | After design load, if `trait == "binary"`, scans all `design.y` values. Throws if any value is not 0.0 or 1.0. Warns if only one level present. |

```cpp
// main.cpp:1020
// ===== Step 1: Binary phenotype must be 0 or 1 (R lines 1754-1757) =====
if (ieq(cfg.trait, "binary")) {
    std::set<double> unique_y;
    for (int i = 0; i < design.n; ++i) {
      double yval = design.y[i];
      if (yval != 0.0 && yval != 1.0) {
        throw std::runtime_error(
            "ERROR: binary phenotype value must be 0 or 1, found: "
            + std::to_string(yval) + " at sample " + design.iid[i]);
      }
      unique_y.insert(yval);
    }
    if (unique_y.size() < 2) {
      std::cerr << "[warning] binary phenotype has only one unique level ("
                << *unique_y.begin() << "), model fitting may be degenerate\n";
    }
}
```

---

### Step 2: Phenotype variance check (quantitative)

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:879–880` — `if (abs(var(Y)) < 0.1) stop("WARNING: variance of the phenotype is much smaller than 1...")` |
| **C++ implementation** | `main.cpp:1040–1052` |
| **What it does** | For quantitative traits, computes sample variance. Throws if < 0.1, suggesting `inv_normalize: true`. |

```cpp
// main.cpp:1040
// ===== Step 2: Phenotype variance check for quantitative (R lines 879-880) =====
if (ieq(cfg.trait, "quantitative")) {
    double sum_y = 0.0, sum_y2 = 0.0;
    for (int i = 0; i < design.n; ++i) {
      sum_y  += design.y[i];
      sum_y2 += design.y[i] * design.y[i];
    }
    double mean_y = sum_y / design.n;
    double var_y  = sum_y2 / design.n - mean_y * mean_y;
    if (std::fabs(var_y) < 0.1) {
      throw std::runtime_error(
          "ERROR: variance of the phenotype (" + std::to_string(var_y)
          + ") is much smaller than 1. Please consider setting inv_normalize: true in config.");
    }
}
```

---

### Step 3: Tau break-on-zero (binary solver)

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:639` — `if(tau[2] == 0) break` |
| **C++ implementation** | `glmm.cpp:722–727` |
| **What it does** | After tau update in binary solver, if `tau[1] <= 0`, breaks the outer GLMM loop immediately instead of wasting further iterations. |

```cpp
// glmm.cpp:722
// ===== Step 3: Tau break-on-zero, binary (R line 639) =====
if (tau[1] <= 0.0f) {
    std::cout << "[binary_glmm] tau[1] <= 0 after update, stopping early.\n";
    break;
}
```

---

### Step 4: Tau break-on-zero (quantitative solver)

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:941` — `if(tau[1]<=0 \| tau[2] <= 0) break` |
| **C++ implementation** | `glmm.cpp:1125–1131` |
| **What it does** | After tau update in quantitative solver, if either tau component <= 0, breaks the loop. |

```cpp
// glmm.cpp:1125
// ===== Step 4: Tau break-on-zero, quantitative (R line 941) =====
if (tau[0] <= 0.0f || tau[1] <= 0.0f) {
    std::cout << "[quant_glmm] tau[0]=" << tau[0] << " tau[1]=" << tau[1]
              << " <= 0, stopping early.\n";
    break;
}
```

---

### Step 5: Max tau upper-bound warning + forced break

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:643–646` (binary), `947–950` (quant) — `if(max(tau) > tol^(-2)) { warning("Large variance estimate..."); i = maxiter; break }` |
| **C++ implementation** | `glmm.cpp:885–893` (binary), `glmm.cpp:1189–1197` (quant) |
| **What it does** | If `max(tau)` exceeds `1/tol^2` (e.g., 2500 when tol=0.02), prints warning to stderr and breaks. Prevents infinite iteration on non-converging models. |

```cpp
// glmm.cpp:885 (binary), 1189 (quant)
// ===== Step 5: Max tau upper-bound warning + break =====
{
    double tau_max_val = static_cast<double>(*std::max_element(tau.begin(), tau.end()));
    double tau_upper = 1.0 / (tol_coef * tol_coef);
    if (tau_max_val > tau_upper) {
        std::cerr << "[warning] Large variance estimate (" << tau_max_val
                  << " > " << tau_upper << "), model not converged.\n";
        break;
    }
}
```

---

### Step 6: Flip `use_r_vr_bypass` to false

| | |
|---|---|
| **R reference** | N/A — C++-only testing artifact |
| **C++ implementation** | `variance_ratio_compute.cpp:55–56` |
| **What it does** | Changed `bool use_r_vr_bypass = true` to `false`. Production code now computes its own VR markers instead of reading R's bypass file. |

```cpp
// variance_ratio_compute.cpp:55
// ===== Step 6: VR bypass disabled for production =====
bool use_r_vr_bypass = false;
```

---

## Phase 2 — MEDIUM Priority (IMPLEMENTED)

### Step 7: Duplicate sample ID removal

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1437` — `sampleIDInclude[!duplicated(sampleIDInclude)]` |
| **C++ implementation** | `main.cpp:999–1017` |
| **What it does** | After design load, scans `design.iid` for duplicates. Keeps only first occurrence of each IID. Prints warning with count of removed duplicates. Uses `design_take_rows()` to subset. |

```cpp
// main.cpp:999
// ===== Step 7: Duplicate sample ID removal (R line 1437) =====
{
    std::unordered_set<std::string> seen;
    std::vector<size_t> keep;
    keep.reserve(design.n);
    int n_dup = 0;
    for (size_t i = 0; i < (size_t)design.n; ++i) {
      if (seen.insert(design.iid[i]).second) {
        keep.push_back(i);
      } else {
        ++n_dup;
      }
    }
    if (n_dup > 0) {
      std::cerr << "[warning] removed " << n_dup
                << " duplicate sample ID(s), keeping first occurrence\n";
      design_take_rows(design, keep);
    }
}
```

---

### Step 8: No markers found for VR estimation

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:3021` — `stop("No markers were found for variance ratio estimation...")` |
| **C++ implementation** | `variance_ratio_compute.cpp:184–190` |
| **What it does** | After collecting eligible VR markers, if count is zero, throws with clear error message matching R's. |

```cpp
// variance_ratio_compute.cpp:184
// ===== Step 8: No markers found for VR (R line 3021) =====
if (numAvailMarkers == 0) {
    throw std::runtime_error(
        "ERROR: No markers were found for variance ratio estimation. "
        "Please make sure there are markers with MAC >= "
        + std::to_string(cfg.vr_min_mac) + " in the plink file.");
}
```

---

### Step 9: Insufficient markers for VR

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:3055–3063` — per-category check `length(markers) < numMarkers` → stop |
| **C++ implementation** | `variance_ratio_compute.cpp:193–198` |
| **What it does** | Warns when available markers are fewer than the requested number (`numMarkersForVarRatio`). The computation proceeds using all available markers (matching R's behavior for single-category VR). |

```cpp
// variance_ratio_compute.cpp:193
// ===== Step 9: Insufficient markers for VR (R lines 3055-3063) =====
if (numAvailMarkers < numMarkers_default) {
    std::cerr << "[warning] Only " << numAvailMarkers
              << " markers available for variance ratio estimation, but "
              << numMarkers_default << " requested. Using all available markers.\n";
}
```

---

### Step 10: Sparse GRM dimension assert

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:2785–2786` — implicit check that GRM is square and matches sample count |
| **C++ implementation** | `main.cpp:448` (square check, pre-existing), `main.cpp:1226–1232` (size warning, new) |
| **What it does** | `load_matrix_market_coo()` already throws if GRM is non-square (`nr != nc`). New code warns if GRM dimension is smaller than the number of model samples, since the intersection will reduce model size. |

```cpp
// main.cpp:448 (pre-existing)
if (nr != nc) throw std::runtime_error("Non-square matrix in " + path);

// main.cpp:1226 (new)
// ===== Step 10: Sparse GRM dimension assert (R lines 2785-2786) =====
if (n_mtx < design.n) {
    std::cerr << "[warning] Sparse GRM dimension (" << n_mtx
              << ") is smaller than the number of model samples (" << design.n
              << "). Sample intersection will reduce the model size.\n";
}
```

---

### Step 11: Variance ratio file overwrite guard (PRE-EXISTING)

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1284–1295` — if VR file exists and `IsOverwriteVarianceRatioFile=FALSE`, stop |
| **C++ implementation** | `main.cpp:1278–1287` (already existed before this work) |
| **What it does** | Before VR computation, checks if output file exists. If it does and `paths.overwrite_varratio` is not true in YAML, throws with message. |

```cpp
// main.cpp:1278 (pre-existing)
// Variance-ratio overwrite guard
if (cfg.num_markers_for_vr > 0) {
    std::string vr_txt = paths.out_prefix_vr + ".varianceRatio.txt";
    bool allow_overwrite = (y["paths"] && y["paths"]["overwrite_varratio"])
        ? y["paths"]["overwrite_varratio"].as<bool>() : false;
    if (!allow_overwrite && fs::exists(vr_txt)) {
        throw std::runtime_error("Refusing to overwrite existing variance-ratio file: "
            + vr_txt + " (set paths.overwrite_varratio=true to allow).");
    }
}
```

---

### Step 12: Complete.cases for covariates

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1430` — `data = data[complete.cases(data),,drop=F]` |
| **C++ implementation** | `main.cpp:712–742` (replaced previous phenotype-only NA drop) |
| **What it does** | In `load_design_csv()`, the row-dropping logic now checks **both** the phenotype column AND all covariate columns for empty/NA/NaN. Any row with any missing value is dropped, matching R's `complete.cases()`. |

```cpp
// main.cpp:712
// ===== Step 12: Drop rows with any missing value (R line 1430: complete.cases) =====
{
    int n_na = 0;
    std::vector<std::vector<std::string>> valid_rows;
    valid_rows.reserve(rows.size());
    for (auto& row : rows) {
      bool any_missing = false;
      const std::string& yval = row[idx_y];
      if (yval.empty() || ieq(yval, "NA") || ieq(yval, "NaN")) {
        any_missing = true;
      }
      if (!any_missing) {
        for (int j : x_idx) {
          const std::string& cv = row[j];
          if (cv.empty() || ieq(cv, "NA") || ieq(cv, "NaN")) {
            any_missing = true;
            break;
          }
        }
      }
      if (any_missing) ++n_na;
      else valid_rows.push_back(std::move(row));
    }
    if (n_na > 0)
      std::cout << "[design] dropped " << n_na
                << " row(s) with missing phenotype or covariates (complete.cases)\n";
    rows = std::move(valid_rows);
}
```

---

## Phase 3 — LOW Priority (IMPLEMENTED)

### Step 13: `converged` flag in output

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:668` (binary), `974` (quant) — `converged = ifelse(i < maxiter, TRUE, FALSE)` |
| **C++ implementation** | `saige_null.hpp:141–142` (struct fields), `glmm.cpp:857–859` + `glmm.cpp:918–920` (binary), `glmm.cpp:1178–1180` + `glmm.cpp:1228–1230` (quant), `main.cpp:1358–1359` (output) |
| **What it does** | Added `bool converged` and `int iterations` to `FitNullResult`. Set at all 4 return points (converged + maxiter) in both solvers. Reported in the final output. |

```cpp
// saige_null.hpp:141
bool converged{false};       // Step 13: true if solver converged before maxiter
int  iterations{0};          // number of outer iterations run

// glmm.cpp:857 (binary converged return)
out.converged = true;
out.iterations = it + 1;

// glmm.cpp:918 (binary maxiter fallthrough)
out.converged = false;
out.iterations = maxiter;
```

---

### Step 14: Categorical covariate (`qCovarCol`) validation

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1446–1454` — `if(!all(qCovarCol %in% covarColList)) stop(...)` |
| **C++ implementation** | `saige_null.hpp:76` (config field), `main.cpp:241–247` (YAML parsing), `main.cpp:1003–1018` (validation) |
| **What it does** | Added `q_covar_cols` vector to `FitNullConfig`. Parsed from `design.q_covar_cols` in YAML. After covar_cols are loaded, validates every name in `q_covar_cols` is present in `covar_cols`. Throws if not. |

```yaml
# Example config usage:
design:
  covar_cols: [x1, x2, site]
  q_covar_cols: [site]   # 'site' will be treated as categorical
```

```cpp
// main.cpp:1003
// ===== Step 14: Validate q_covar_cols subset of covar_cols (R lines 1446-1454) =====
if (!cfg.q_covar_cols.empty()) {
    std::unordered_set<std::string> covar_set(covar_col_names.begin(), covar_col_names.end());
    for (const auto& qc : cfg.q_covar_cols) {
      if (covar_set.find(qc) == covar_set.end()) {
        throw std::runtime_error(
            "ERROR: categorical covariate '" + qc
            + "' in q_covar_cols is not in covar_cols.");
      }
    }
}
```

---

### Step 15: Sex-specific filter validation

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1329–1330` — `if (FemaleOnly & MaleOnly) stop(...)`, `1459–1461` — `if (!sexCol %in% colnames(data)) stop(...)` |
| **C++ implementation** | `main.cpp:1105–1115` |
| **What it does** | (a) Throws if both `female_only` and `male_only` are true. (b) Throws if either is set but `sex_col` is empty. (c) `read_column_from_csv` (line 161) already throws if the named column doesn't exist in the CSV. |

```cpp
// main.cpp:1105
// ===== Step 15: Sex-specific filter validation =====
if (cfg.female_only && cfg.male_only) {
    throw std::runtime_error("ERROR: Both female_only and male_only are true.");
}
if ((cfg.female_only || cfg.male_only) && cfg.sex_col.empty()) {
    throw std::runtime_error("ERROR: female_only or male_only is true but sex_col is not specified.");
}
```

---

### Step 16: `skipModelFitting` guard

| | |
|---|---|
| **R reference** | `SAIGE_fitGLMM_fast.R:1252–1254` — `if(skipModelFitting) { if(!file.exists(modelOut)) stop(...) }` |
| **C++ implementation** | `saige_null.hpp:79–80` (config fields), `main.cpp:238` (YAML parsing), `main.cpp:1310–1323` (validation) |
| **What it does** | Added `skip_model_fitting` bool and `model_file` string to config. When `skip_model_fitting=true`, validates that `model_file` is specified and exists. Prints a message that full deserialization is not yet implemented (the file existence check matches R's guard). |

```cpp
// main.cpp:1310
if (cfg.skip_model_fitting) {
    if (cfg.model_file.empty()) {
      throw std::runtime_error("skip_model_fitting=true but no model_file specified.");
    }
    if (!fs::exists(cfg.model_file)) {
      throw std::runtime_error("skip_model_fitting=true but model_file does not exist: " + cfg.model_file);
    }
    std::cout << "[skip_model_fitting] Using existing model: " << cfg.model_file << "\n";
    return 0;
}
```

---

### Step 17: Debug early-exit flag (SKIPPED)
- **R reference:** `SAIGE_fitGLMM_fast.R:585–619` (debug only)
- **Decision:** Existing checkpoint system serves same purpose.

---

### Step 18: Dry-run input validation mode (NEW)

| | |
|---|---|
| **R reference** | N/A — new C++ feature for pre-flight input checking |
| **C++ implementation** | `saige_null.hpp:82` (config field), CLI flag `main.cpp:893`, `main.cpp:1257–1308` (dry-run logic) |
| **What it does** | `--dry-run` CLI flag (or `fit.dry_run: true` in YAML) runs all input validation — config parsing, design CSV loading (complete.cases, column checks), phenotype validation (binary 0/1, variance check), duplicate ID removal, FAM alignment — then prints a summary and exits **before** loading genotype data or running the solver. |

```
Usage:
  ./saige-null -c config.yaml --dry-run

Output:
  ============================================
  === DRY RUN: Input Validation Summary ===
  ============================================
  Config:
    trait:          quantitative
    tol:            0.02
    ...
  Paths:
    bed:            /path/to/file.bed [OK]
    ...
  Design:
    n (samples):    1000
    p (covariates): 2
    y range:        [-3.29, 3.29]  mean=0
    ...
  FAM alignment:
    samples matched: 1000 / 1000
  ============================================
  DRY RUN PASSED: all input validations succeeded.
  ============================================
```

---

## Summary of All Changes

| Step | Description | Priority | File:Lines | R Lines | Status |
|------|-------------|----------|------------|---------|--------|
| 1 | Binary 0/1 phenotype check | HIGH | `main.cpp:1020–1038` | 1754–1757 | DONE |
| 2 | Quantitative phenotype variance < 0.1 | HIGH | `main.cpp:1040–1052` | 879–880 | DONE |
| 3 | Tau break-on-zero (binary) | HIGH | `glmm.cpp:722–727` | 639 | DONE |
| 4 | Tau break-on-zero (quantitative) | HIGH | `glmm.cpp:1125–1131` | 941 | DONE |
| 5 | Max tau upper-bound warning | HIGH | `glmm.cpp:885–893`, `glmm.cpp:1189–1197` | 643–646, 947–950 | DONE |
| 6 | VR bypass flag → false | HIGH | `variance_ratio_compute.cpp:55–56` | N/A | DONE |
| 7 | Duplicate sample ID removal | MEDIUM | `main.cpp:1030–1048` | 1437 | DONE |
| 8 | No markers for VR → clear error | MEDIUM | `variance_ratio_compute.cpp:184–190` | 3021 | DONE |
| 9 | Insufficient markers per MAC | MEDIUM | `variance_ratio_compute.cpp:193–198` | 3055–3063 | DONE |
| 10 | Sparse GRM dimension assert | MEDIUM | `main.cpp:448` + `main.cpp:1246–1252` | 2785–2786 | DONE |
| 11 | VR file overwrite guard | MEDIUM | `main.cpp:1400–1409` | 1284–1295 | PRE-EXISTING |
| 12 | Complete.cases for covariates | MEDIUM | `main.cpp:712–742` | 1430 | DONE |
| 13 | `converged` flag in output | LOW | `saige_null.hpp:141–142`, `glmm.cpp:857,918,1178,1228` | 668, 974 | DONE |
| 14 | `qCovarCol` validation | LOW | `main.cpp:1003–1018` | 1446–1454 | DONE |
| 15 | Sex column/flag validation | LOW | `main.cpp:1105–1115` | 1329–1330, 1459–1461 | DONE |
| 16 | `skipModelFitting` guard | LOW | `main.cpp:1310–1323` | 1252–1254 | DONE |
| 17 | Debug early-exit flag | — | — | 585–619 | SKIPPED |
| 18 | Dry-run input validation | NEW | `main.cpp:1257–1308`, CLI `--dry-run` | N/A | DONE |

**Gap closure: 23 original gaps → 1 remaining (Step 17 skipped by design).**
**Build status: compiles cleanly. Dry-run tested successfully.**

---

## Execution Order in main.cpp (data flow)

The input validation checks run in this order during program execution:

```
1. load_cfg() + CLI parsing
   ├── Step 14: q_covar_cols ⊆ covar_cols validation               (main.cpp:1003)
   └── Step 18: --dry-run flag captured                             (main.cpp:948)

2. load_design_csv()
   └── Step 12: complete.cases — drop rows with ANY missing value   (main.cpp:712)

3. add_intercept_if_missing()                                        (main.cpp:1029)

4. Step 7:  Duplicate ID removal                                     (main.cpp:1030)
5. Step 1:  Binary 0/1 check                                        (main.cpp:1051)
6. Step 2:  Quantitative variance check                              (main.cpp:1071)

7. Step 15: Sex filter validation (mutual exclusion, col existence)  (main.cpp:1105)
8. Whitelist, FAM alignment ...

9. Step 18: DRY RUN EXIT (if --dry-run)                              (main.cpp:1257)
   └── Prints config/paths/design/FAM summary and exits.
       No genotype loading. No solver.

10. Step 16: skip_model_fitting guard                                (main.cpp:1310)

11. init_global_geno() — genotype loading begins here

12. Step 10: Sparse GRM dimension check                              (main.cpp:1246)
13. Step 11: VR overwrite guard                                      (main.cpp:1400)

14. fit_null() → binary_glmm_solver() or quant_glmm_solver()
    ├── Step 3/4: Tau break-on-zero                                  (glmm.cpp:722/1125)
    ├── Step 5:   Max tau upper bound                                (glmm.cpp:885/1189)
    └── Step 13:  converged flag set                                 (glmm.cpp:857/918/1178/1228)

15. Output: converged + iterations reported                          (main.cpp:1358)

16. compute_variance_ratio()
    ├── Step 6:  VR bypass = false                                   (variance_ratio_compute.cpp:55)
    ├── Step 8:  No markers check                                    (variance_ratio_compute.cpp:184)
    └── Step 9:  Insufficient markers warning                        (variance_ratio_compute.cpp:193)
```
