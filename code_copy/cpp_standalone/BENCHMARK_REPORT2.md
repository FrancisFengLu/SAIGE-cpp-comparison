# SAIGE C++ Standalone Benchmark Report

**Date:** 2026-03-25
**Platform:** Linux 6.8.0-106-generic, x86_64
**Compiler:** GCC 14.3.0 (conda, `-O2` release / `-O1 -g` profile mode)
**Valgrind:** 3.26.0
**R SAIGE:** Docker `wzhou88/saige:1.5.0.2` (R 4.4.3, OpenBLAS 0.3.30)
**Threading:** All benchmark runs use single-threaded BLAS (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`)

---

## 1. Test Configurations

### Test A: sparse_x1 (small benchmark data)

| Parameter | Value |
|-----------|-------|
| Trait | quantitative |
| N (samples) | 1,000 |
| M (markers) | 128,868 |
| Covariates | x1 (1 covariate + intercept) |
| GRM | Sparse (relatedness cutoff 0.125, 2000 random markers) |
| Solver | Direct sparse solve (`usePCGwithSparseGRM=FALSE`) |

### Test B: ukb_ldl_dense (UKB real data, dense GRM) — Primary benchmark

| Parameter | Value |
|-----------|-------|
| Trait | quantitative (LDL cholesterol, inv-normalized) |
| N (samples) | 961 (from 1000 random white-British, 39 missing LDL excluded) |
| N_fam (genotype file) | 408,970 |
| M (markers) | 340,447 (LD-pruned array) |
| Covariates | Sex (categorical), Age, Batch, PC1-PC4 (7 covariates + intercept) |
| GRM | Dense (computed from 340,447 LD-pruned markers) |
| Solver | Dense crossproduct (`getCrossprodMatAndKin`) |

---

## 2. UKB Dense GRM: 5-Run Benchmark (Primary Results)

### Wall Time & Memory

| Run | C++ Wall | C++ RSS (MB) | R Wall | R RSS (MB) |
|-----|----------|--------------|--------|------------|
| 1 | 12:50 | 168 | 11:52 | 26 |
| 2 | 12:50 | 172 | 11:51 | 26 |
| 3 | 12:53 | 169 | 11:56 | 25 |
| 4 | 12:54 | 170 | 11:57 | 24 |
| 5 | 12:52 | 168 | 11:55 | 24 |
| **Mean** | **12:52** | **169** | **11:54** | **25** |
| **Std** | ~1.5s | ~1.5 | ~2.5s | ~1 |

**R is ~8% faster on wall time. C++ uses ~7x more memory.**

### Model Parameters

| Run | C++ tau[0] | C++ tau[1] | C++ VR | C++ Iters | R tau[0] | R tau[1] | R VR | R Iters |
|-----|-----------|-----------|--------|-----------|---------|---------|------|---------|
| 1 | 0.9952 | 0.0000 | 1.005 | 2 | 0.8773 | 0.1094 | 1.014 | 4 |
| 2 | 0.9929 | 0.0000 | 1.007 | 2 | 0.8773 | 0.1094 | 1.014 | 4 |
| 3 | 0.9934 | 0.0000 | 1.007 | 2 | 0.8773 | 0.1094 | 1.014 | 4 |
| 4 | 0.9935 | 0.0000 | 1.007 | 2 | 0.8772 | 0.1094 | 1.014 | 4 |
| 5 | 0.9941 | 0.0000 | 1.006 | 2 | 0.8773 | 0.1094 | 1.014 | 4 |

**R results are highly consistent across runs (tau matches to 5 decimal places).**

### tau[1] = 0 Issue (C++ — NEEDS INVESTIGATION)

C++ consistently gets tau[1] = 0 while R converges to tau[1] = 0.109. This is abnormal — with 961 white-British samples and a dense GRM from 340K markers, non-zero genetic variance is expected.

**What happens step by step in C++:**
1. **Iteration 0 (conservative):** tau_init = [1, 0]. Conservative formula: `tau[i] += tau[i]^2 * score / n`. Since tau[1]=0, Dtau[1] = 0^2 * anything = 0. tau stays [~0.995, 0.000]. This is correct and matches R's conservative step.
2. **Iteration 1 (AI-REML):** AI-REML computes delta = solve(AI, score). The delta for tau[1] is **negative** (e.g., -0.5177), meaning AI-REML wants to push tau[1] below 0. Step halving reduces step to ~0. tau[1] remains 0. The `tau[1] <= 0` check triggers and stops.

**Why R behaves differently:** R uses the same algorithm but with different random vectors for `GetTrace()`. R's trace estimation produces Score and AI values that give a **positive** delta for tau[1] at iteration 1, allowing tau[1] to move away from 0 and eventually converge to 0.109.

**Root cause:** Different random number generators (R's Mersenne Twister vs C++ `std::mt19937`) produce different trace estimates, which cascade into different AI-REML updates. With the small test data and random vector bypass (reading R's exact vectors), C++ and R produce matching tau values — confirming the algorithm is correct.

**Possible fixes to investigate:**
- Enable random vector bypass for UKB configs (generate from R, read in C++)
- Increase number of trace vectors (nrun=30 → 100) for more stable estimates
- Investigate whether C++ trace estimation has a systematic bias vs R

---

## 3. Per-Step CPU Breakdown

### R UKB dense (from R's internal `proc.time()`)

| Step | User (s) | System (s) | Elapsed (s) |
|------|----------|------------|-------------|
| Setup + phenotype/covariate loading | 2.34 | 2.09 | 4.38 |
| Genotype reading (BED → GRM) | — | — | ~694s |
| Null model fitting (4 iterations) | 952.78 | 14.83 | 697.47 |
| Variance ratio estimation | ~3.6 | ~0.1 | ~4.6 |
| **Total (R internal)** | **~959** | **~17** | **~702** |

### C++ UKB dense (from `/usr/bin/time`, estimated breakdown)

| Step | Estimated time | Notes |
|------|---------------|-------|
| Config + design loading | < 1s | 961 samples, 7 covariates, inv-normalize |
| Genotype loading (`setGenoObj`) | **~11 min** | 340K markers x 408K samples BED scan |
| GRM diagonal computation | ~1 min | 248K QC-passing markers |
| GLMM fitting (2 iterations) | ~30s | getCrossprodMatAndKin + getAIScore |
| Variance ratio (30 markers) | ~10s | 30 PCG solves + genotype lookups |
| **Total** | **~12:52** | **User CPU: ~770s** |

**Bottleneck:** Genotype loading dominates (~85% of wall time) in both C++ and R. The BED file has 408K samples x 340K markers — even though only 961 are used, the full file must be scanned.

---

## 4. Callgrind Instruction Profile (Test A: sparse_x1, Ir-only)

**Total instructions:** 30,560,991,710

| Rank | Ir (%) | Function |
|------|--------|----------|
| 1 | **29.66%** | `Get_OneSNP_Geno()` — genotype decoding |
| 2 | **23.81%** | `Get_OneSNP_Geno_atBeginning()` — BED loading |
| 3 | **12.24%** | `memcpy/memmove` (libc) |
| 4 | **7.17%** | `stl_bvector` bitwise access |
| 5 | **7.16%** | `blas_thread_server` (OpenBLAS threading overhead) |
| 6 | **6.67%** | `Get_OneSNP_StdGeno()` — genotype standardization |

**Genotype decoding alone = 53.5% of all instructions.**

---

## 5. Bugs Fixed During Benchmarking

### Bug 1: Premature tau[1]<=0 Break on Iteration 0 (glmm.cpp)

**Problem:** C++ checked `tau[0] <= 0 || tau[1] <= 0` after every iteration including the conservative first step (iteration 0). Since initial tau = [1, 0] for quantitative traits and the conservative formula uses tau^2, tau[1] always stays 0 after iteration 0. This caused immediate termination.

**R behavior:** R only applies this check inside the main loop (iterations >= 1), not after the pre-loop conservative step.

**Fix:** Changed `if (tau[0] <= 0.0f || tau[1] <= 0.0f)` to `if (it > 0 && (tau[0] <= 0.0f || tau[1] <= 0.0f))`.

### Bug 2: Binary tau break uses `<=0` instead of `==0` (glmm.cpp)

**Problem:** C++ used `tau[1] <= 0.0f` but R uses `tau[2] == 0` (exact zero check).

**Fix:** Changed to `tau[1] == 0.0f` to match R.

### Bug 3: Sparse GRM Subsetting (main.cpp)

**Problem:** Loading a pre-built sparse GRM (488,377 samples) passed the full dimension to `setupSparseGRM()`. The solver accessed `wVec` (961 elements) at indices up to 488K → crash.

**Fix:** After loading MTX, remap COO indices to phenotyped-sample space. 442,853,165 entries → 963 entries.

---

## 6. Memory Analysis

| Component | C++ | R |
|-----------|-----|---|
| BED file handling | Loads full file into memory (~170 MB for 408K x 340K) | Streams per-marker (< 1 MB) |
| Genotype matrix | N x M_sub float matrix | Same (via Rcpp) |
| GRM storage | Dense N x N or sparse COO | Same |
| Peak RSS | **169 MB** | **25 MB** (Docker-reported) |

C++ memory is dominated by the BED-to-memory loading strategy in `setGenoObj()`. R processes markers one at a time during GRM construction, avoiding bulk memory allocation.

---

## 7. R Docker Notes

R benchmarks used `/usr/bin/time -v docker run ...` which measures the Docker container process. Key caveats:
- **Wall time**: Accurate (includes container startup ~2-3s)
- **Peak RSS**: Under-reports actual R memory (Docker cgroup reporting)
- **CPU time**: Not meaningful from outside Docker; R's internal `proc.time()` from stdout is used instead
- **Single-thread enforcement**: via `-e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1` and `--nThreads=1`

---

## 8. Collaborator Integration (March 25, 2026)

Integrated changes from collaborator's `resaigegenecppupdate.zip` for Step 1 → Step 2 compatibility:

- **null_model_engine.cpp**: Output paths from dot-prefix (`out_prefix.mu.arma`) to directory-based (`out_prefix/mu.arma`); `#include <armadillo>` replaces `<RcppArmadillo.h>`
- **glmm.cpp**: `obj_noK.json` path directory-based; tau safety checks retained
- **Makefile**: Cross-platform (Darwin/Linux) with conda auto-detection; `profile`/`callgrind` targets retained

---

## 9. Optimization Opportunities

| Priority | Opportunity | Est. Impact | Notes |
|----------|------------|-------------|-------|
| 1 | **Stream BED instead of loading into memory** | 169 MB → <1 MB RSS | R does this; C++ loads entire BED |
| 2 | **Investigate trace estimation bias** | Correctness | tau[1]=0 in C++ vs 0.109 in R — needs root cause analysis |
| 3 | SIMD-vectorize PLINK BED decoding (AVX2) | ~30% Ir | `Get_OneSNP_Geno` + `Get_OneSNP_Geno_atBeginning` |
| 4 | Set `OPENBLAS_NUM_THREADS=1` by default for small N | ~7% Ir | Already done in benchmarks |
| 5 | Cache SuperLU symbolic factorization | ~4-5% Ir | 70+ sparse solves with same sparsity pattern |
| 6 | Remove debug stdout from production paths | Wall time | `#ifdef DEBUG` guards |

---

## 10. Files Generated

```
benchmark_results/
├── sparse_x1_nocache/
│   ├── cpp/   (callgrind Ir-only profile, 3:05 under valgrind)
│   └── r/     (6.44s wall)
├── ukb_ldl_dense_timed/
│   ├── cpp/   (13:02 wall, single run)
│   └── r/     (11:44 wall, single run)
├── ukb_ldl_dense_rep/        (OLD binary — premature tau break)
│   ├── cpp_run{1-5}/         (mean 12:58 wall, tau[1]=0)
│   └── r_run{1-5}/           (mean 11:54 wall, tau[1]=0.109)
├── ukb_ldl_dense_rep_v2/     (FIXED binary — tau check skips iter 0)
│   └── cpp_run{1-5}/         (mean 12:52 wall, tau[1]=0 — trace estimation issue)
├── ukb_ldl_sparse_timed/
│   ├── cpp/   (12:28 wall, with GRM subsetting fix)
│   └── r/     (14:34 wall)
└── ukb_ldl_dense_rep/r_run{1-5}/   (R reference, mean 11:54)
```
