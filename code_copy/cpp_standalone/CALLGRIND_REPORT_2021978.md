# Callgrind Profiling Report — PID 2021978
# C++ Standalone vs R SAIGE: Performance & Algorithm Comparison

**Profile file:** `callgrind.out.2021978`
**Total instructions (Ir):** 30,380,167,932

---

## 1. Top 20 Hotspots by Instruction Count

| Rank | Ir (abs) | Ir (%) | Function |
|------|----------|--------|----------|
| 1 | 9,063,104,460 | **29.83%** | `genoClass::Get_OneSNP_Geno()` — SAIGE_step1_fast.cpp |
| 2 | 7,275,713,485 | **23.95%** | `genoClass::Get_OneSNP_Geno_atBeginning()` — SAIGE_step1_fast.cpp |
| 3 | 3,779,519,037 | **12.44%** | `0x189480` libc (memcpy/memmove inner loop) |
| 4 | 2,191,786,944 | **7.21%** | `stl_bvector.h` — std::vector<bool> bitwise access inside `Get_OneSNP_Geno_atBeginning` |
| 5 | 2,038,439,770 | **6.71%** | `genoClass::Get_OneSNP_StdGeno()` — SAIGE_step1_fast.cpp |
| 6 | 696,819,380 | **2.29%** | `stl_vector.h` — vector size() inside `Get_OneSNP_Geno` |
| 7 | 618,945,377 | **2.04%** | `blas_thread_server` — OpenBLAS threading overhead |
| 8 | 561,017,520 | **1.85%** | `Mat_meat.hpp` — Armadillo matrix body inside `Get_OneSNP_Geno` |
| 9 | 391,793,600 | **1.29%** | `eglue_core_meat.hpp` — `arma::eglue_plus` elementwise add |
| 10 | 284,153,690 | **0.94%** | `arma::Mat::zeros()` — zero-init of genotype vectors |
| 11 | 281,065,325 | **0.93%** | `arrayops_meat.hpp` — Armadillo array ops inside `Get_OneSNP_Geno` |
| 12 | 242,363,756 | **0.80%** | `Mat_meat.hpp` — Armadillo init inside `Get_OneSNP_Geno_atBeginning` |
| 13 | 228,318,017 | **0.75%** | `blas_memory_alloc` — repeated BLAS thread workspace allocation |
| 14 | 225,058,710 | **0.74%** | `Mat_meat.hpp` inside `Get_OneSNP_StdGeno` |
| 15 | 223,277,830 | **0.73%** | `debug.hpp` — Armadillo bounds checks inside `Get_OneSNP_StdGeno` |
| 16 | 200,928,584 | **0.66%** | `blas_memory_free` — repeated BLAS thread workspace deallocation |
| 17 | 189,096,791 | **0.62%** | `colamd` — SuperLU column ordering for sparse factorization |
| 18 | 167,068,805 | **0.55%** | `memory.hpp` — Armadillo allocator inside `eglue_plus` |
| 19 | 139,242,555 | **0.46%** | `stl_vector.h` inside `Get_OneSNP_StdGeno` |
| 20 | 58,679,900 | **0.19%** | `std::__introsort_loop` — sort inside SpMat construction |

**Genotype decoding alone: `Get_OneSNP_Geno` + `Get_OneSNP_Geno_atBeginning` = 53.78% of all instructions.**

---

## 2. Call Graph Summary

```
main / fit_null_model
├── genoClass::setGenoObj()                         37.15% inclusive
│   ├── Get_OneSNP_Geno_atBeginning()               36.48% inclusive  ← genotype LOADING hotspot
│   │   ├── std::istream::seekg/read                12.74%            ← 128,868 BED file seeks
│   │   ├── std::vector<unsigned char>::reserve     12.74%
│   │   └── stl_bvector bitwise access               7.21%
│   └── get_GRMdiagVec()                            13.15% inclusive
│       └── Get_OneSNP_StdGeno()                     6.71%
│
├── [GLMM fitting loop]
│   ├── getCrossprodMatAndKin()          0.53%  ×150 calls
│   │   └── parallelCrossProd()
│   │       └── Get_OneSNP_StdGeno()         ← inner loop dominant
│   │
│   ├── gen_spsolve_v4()                 4.43%  ×227 calls  ← PCG solve
│   │   ├── gen_sp_Sigma()               0.66%              ← sparse Sigma rebuilt each call
│   │   │   └── SpMat::init_batch_std() + arma_sort_index   0.87%
│   │   └── arma::spsolve() → SuperLU
│   │       ├── colamd                   0.62%
│   │       ├── dgstrf (LU factor)       0.08%
│   │       └── dgstrs (LU solve)        0.06%
│   │
│   └── GetTrace()
│       └── getPCG1ofSigmaAndVector()    2.69%  ×150 calls
│           └── gen_spsolve_v4()
│
└── Get_OneSNP_Geno()                   29.83%  ×556,525 calls  ← VR + GRM diagonal
```

**Call counts:**
- `Get_OneSNP_Geno`: **556,525** calls (VR computation + GRM diagonal)
- `Get_OneSNP_Geno_atBeginning`: **128,868** calls (initial BED loading, one per marker)
- `gen_spsolve_v4`: **227** calls (GLMM outer iterations × PCG solves + GetTrace × nrun)
- `getPCG1ofSigmaAndVector`: **150** calls

---

## 3. Algorithm Comparison: R vs C++

### 3a. GLMM Fitting Loop

**R (`glmmkin.ai_PCG_Rcpp_Binary` in SAIGE_fitGLMM_fast.R):**
```
1. setgeno()           → calls Rcpp setGenoObj() (same BED loading C++ code)
2. Initial Get_Coef()  → inner IRLS for α convergence
3. Initial getAIScore()→ YPAPY, Trace (30 random vectors), AI
4. Conservative τ update (iter 0 only)
5. Outer loop (maxiter=20):
   a. Get_Coef()       → converge α, η, W, Y
   b. fitglmmaiRPCG()  → one AI-REML step: getAIScore + τ update
   c. Check τ convergence: max(|τ - τ0| / (|τ| + |τ0| + tol)) < tol
```

**C++ (`binary_glmm_solver` in glmm.cpp):**
```
Outer loop (maxiter):
  a. Inner IRLS loop  → irls_binary_build() + getCoefficients_cpp() + α convergence
  b. Final IRLS build → recompute W, Y from converged η
  c. getAIScore_cpp() → YPAPY, GetTrace, AI
  d. Conservative τ on iter==0, standard AI-REML after
  e. Step-halving if τ < 0
  f. τ convergence check (same R formula)
```

**Structural differences:**
- R has a pre-loop call to `getAIScore` for the initial conservative τ update; C++ fuses this into iteration 0.
- R separates `Get_Coef` and `fitglmmaiRPCG` as distinct function calls; C++ has one unified outer loop.
- Otherwise algorithmically identical.

### 3b. PCG Solver (`getPCG1ofSigmaAndVector`)

**Both R and C++ use the same Rcpp C++ implementation** (the standalone ported it directly from `src/SAIGE_step1_fast.cpp`).
- `usePCGwithSparseGRM=FALSE` (default): calls `gen_spsolve_v4()` → SuperLU direct sparse solve.
- `usePCGwithSparseGRM=TRUE`: conjugate gradient iterative loop.
- **Critical:** Both implementations rebuild `gen_sp_Sigma` on every call; SuperLU re-runs symbolic + numeric factorization (`colamd` + `dgstrf`) every time. No caching between calls within a τ iteration.

### 3c. Trace Estimation (`GetTrace`)

| Aspect | R | C++ Standalone |
|--------|---|----------------|
| Random vectors | R RNG (Mt19937 with R seed) | C++ RNG, or bypass file from R |
| Vector count | 30, adaptive (+10 if CV > cutoff) | Same |
| Algorithm | u'PKu for each vector | Identical |
| Bypass | N/A | Reads `output/bypass/random_vectors_seed10.csv` when present |

### 3d. Variance Ratio Computation

Both follow the same algorithm:
```
For each selected marker:
  G_adj = G - XXVX_inv * (XV * G)     # covariate-adjusted genotype
  g = G_adj / sqrt(AC)
  Sigma_iG = PCG solve
  var1 = G' Sigma_iG / AC  (minus covariate correction term)
  var2 = mu*(1-mu) * g'g   (binary)
  ratio = var1 / var2
```

**Key difference:** R uses `set.seed` for random marker selection; C++ uses its own RNG.
- Bypass (`use_r_vr_bypass = true`) reads R's marker indices from `output/bypass/vr_marker_indices.csv`.
- **NOTE: `use_r_vr_bypass` must be set to `false` in production** (currently `true` in variance_ratio_compute.cpp ~line 55).

### 3e. Matrix Operations

Both R and C++ use the **same Armadillo library** (via RcppArmadillo in R, or direct Armadillo headers in C++). Both call the same OpenBLAS and SuperLU backends. The Armadillo header paths in the callgrind output (`RcppArmadillo/include/armadillo_bits/`) confirm both use the same headers.

---

## 4. Bottlenecks Unique to C++ vs R

### 4a. Extensive Debug stdout Output (MAJOR wall-clock impact)

The C++ standalone emits hundreds of `std::cout` lines per run from production code paths:
- `getCoefficients_cpp`: ~15 lines per call (called 150+ times)
- `binary_glmm_solver`: ~20 lines per outer + inner IRLS iteration
- `GetTrace`: comparison data for each of 30+ random vectors

**This is not counted in Ir** (I/O is not CPU instructions) but adds significant wall-clock time. R has no equivalent output.
**Fix:** Guard all debug output with `#ifdef DEBUG` or a runtime `--verbose` flag.

### 4b. Sparse Sigma Rebuilt Every Solve (no R-specific advantage here)

Both R and C++ rebuild `gen_sp_Sigma` and re-run SuperLU `colamd`+`dgstrf` on every `getPCG1ofSigmaAndVector` call. Within a single outer GLMM iteration, W is fixed, so the diagonal of Σ changes only with τ values. The symbolic sparsity pattern never changes.

**Optimization:** Cache the `colamd` column permutation and `sp_coletree` across all solves with the same sparsity pattern.

### 4c. BLAS Threading for Small Matrices (~2.79% Ir)

`blas_thread_server` (2.04%) + `blas_memory_alloc/free` (0.75%) show OpenBLAS spawning threads per BLAS call. For N=1000, p=3, single-threaded BLAS outperforms multi-threaded BLAS.
**Fix:** Set `OMP_NUM_THREADS=1` or `OPENBLAS_NUM_THREADS=1` at launch.

### 4d. Redundant `zeros()` per `Get_OneSNP_Geno` Call (~0.94% Ir)

`arma::Mat::zeros()` (memset) is called at the start of every `Get_OneSNP_Geno` call, zero-filling the 1000-element genotype buffer. Since the loop always overwrites every element for non-missing samples, this zero-fill is redundant when there are no missing genotypes.
**Fix:** Remove or conditionalize the zero-fill.

### 4e. Armadillo Bounds Checks (~0.73% Ir)

`debug.hpp` bounds-check calls inside `Get_OneSNP_StdGeno` hot path.
**Fix:** Compile with `-DARMA_NO_DEBUG` (already in many production builds).

---

## 5. Branch Mispredictions

- Total branch mispredictions (Bcm): **117,082,027** — of which **88.65%** occur in `Get_OneSNP_Geno_atBeginning`.
- Caused by: irregular per-marker loop lengths (varying Nnomissing), data-dependent 2-bit decoding branches, and `std::istream::seekg` state changes.
- R has identical misprediction profile in the same Rcpp C++ code.

---

## 6. Ranked Optimization Opportunities

| Priority | Opportunity | Est. Ir Savings | Notes |
|----------|------------|-----------------|-------|
| 1 | SIMD-vectorize PLINK BED decoding (AVX2: 128 genotypes/cycle) | ~30% Ir | `Get_OneSNP_Geno` + `Get_OneSNP_Geno_atBeginning` |
| 2 | Remove debug stdout from production paths | ~large wall-clock | `#ifdef DEBUG` guards |
| 3 | Cache SuperLU symbolic factorization per outer GLMM iter | ~4-5% Ir, 227 calls × colamd+dgstrf | Reuse permutation, redo only numeric factor when W changes |
| 4 | Set `OPENBLAS_NUM_THREADS=1` for N≤5000 | ~2.79% Ir | Threading overhead > parallelism benefit at N=1000 |
| 5 | Remove `m_OneSNP_Geno.zeros()` when no missing | ~0.94% Ir | Conditional on missingness mask |
| 6 | Compile with `-DARMA_NO_DEBUG` | ~0.73% Ir | Standard Armadillo production flag |
| 7 | Disable `use_r_vr_bypass` in production | Correctness | variance_ratio_compute.cpp ~line 55 |

---

## 7. Key Source Locations

| File | Function | Line (approx) |
|------|----------|---------------|
| `SAIGE_step1_fast.cpp` | `Get_OneSNP_Geno` | 216 |
| `SAIGE_step1_fast.cpp` | `Get_OneSNP_Geno_atBeginning` | 319 |
| `SAIGE_step1_fast.cpp` | `Get_OneSNP_StdGeno` | 580 |
| `SAIGE_step1_fast.cpp` | `getCrossprodMatAndKin` | 2040 |
| `SAIGE_step1_fast.cpp` | `GetTrace` | 5490 |
| `SAIGE_step1_fast.cpp` | `getPCG1ofSigmaAndVector` | 4035 |
| `SAIGE_step1_fast.cpp` | `gen_spsolve_v4` | 3972 |
| `glmm.cpp` | `binary_glmm_solver` | 459 |
| `saige_ai.cpp` | `getCoefficients_cpp` | 22 |
| `saige_ai.cpp` | `getAIScore_cpp` | 69 |
| `variance_ratio_compute.cpp` | `compute_variance_ratio` | 43 |
| `variance_ratio_compute.cpp` | `use_r_vr_bypass` flag | ~55 |
| `SAIGE_fitGLMM_fast.R` | `glmmkin.ai_PCG_Rcpp_Binary` | 241 |
| `SAIGE_fitGLMM_fast.R` | `Get_Coef` | 117 |
| `SAIGE_fitGLMM_fast.R` | outer GLMM loop | 414 |
| `SAIGE_fitGLMM_fast.R` | `extractVarianceRatio` | 2858 |
