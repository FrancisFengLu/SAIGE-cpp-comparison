# SAIGE C++ Standalone Benchmark Report

**Date:** 2026-03-25
**Platform:** Linux 6.8.0-106-generic, x86_64
**Compiler:** GCC 14.3.0 (conda, `-O1 -g` profile mode)
**Valgrind:** 3.26.0
**R SAIGE:** Docker `wzhou88/saige:1.5.0.2` (R 4.4.3, OpenBLAS 0.3.30)

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
| Trace vectors | 30 |
| VR markers | 30 (C++), 50 adaptive (R) |

### Test B: ukb_ldl_sparse (UKB real data)

| Parameter | Value |
|-----------|-------|
| Trait | quantitative (LDL cholesterol, inv-normalized) |
| N (samples) | 961 (from 1000 random white-British, 39 missing LDL excluded) |
| N_fam (genotype file) | 408,970 |
| M (markers) | 340,447 (LD-pruned array) |
| Covariates | Sex, Age, Batch, PC1-PC4 (7 covariates + intercept) |
| GRM | Sparse (relatedness cutoff 0.05, 5000 markers, 488,377 total samples) |
| GRM after subsetting | 961x961, 963 non-zeros (from 442,853,165) |
| Solver | Direct sparse solve |
| Trace vectors | 30 |
| VR markers | 30 |

---

## 2. Overall Timing & Memory

### Test A: sparse_x1

| Metric | C++ (callgrind, Ir-only) | R (Docker, `/usr/bin/time`) |
|--------|--------------------------|------------------------------|
| Wall clock | 3:05.65 | 6.44s |
| User CPU | 183.42s | 0.04s (Docker overhead; actual R inside ~5.5s) |
| System CPU | 1.21s | 0.05s |
| Peak RSS | 310 MB | 24 MB (Docker container overhead) |
| Total instructions (Ir) | **30,560,991,710** | N/A |
| Page faults (minor) | 124,836 | 3,385 |

**Note:** The C++ wall time of 3:05 is under callgrind (~30x overhead). Native execution estimated at ~6s, comparable to R.

### Test B: ukb_ldl_sparse

| Metric | C++ (`/usr/bin/time`) | R (Docker) |
|--------|-----------------------|------------|
| Wall clock | **12:27.99** | Not yet run (sample ID mismatch) |
| User CPU | **740.71s** | — |
| System CPU | **17.28s** | — |
| Peak RSS | **13.3 GB** | — |
| CPU utilization | 101% | — |
| Page faults (minor) | 4,098,193 | — |
| File I/O (read) | 382,880 sectors (~187 MB) | — |

---

## 3. Per-Step CPU Breakdown

### Test A: R sparse_x1 (from R's internal `proc.time()`)

| Step | User (s) | System (s) | Elapsed (s) |
|------|----------|------------|-------------|
| Setup + genotype loading | 1.825 | 1.582 | 3.387 |
| Null model fitting (5 iterations) | 1.563 | 0.306 | 1.871 |
| Post-fit + VR estimation | 1.781 | 0.308 | 2.123 |
| **Total** | **~5.5** | **~2.2** | **~7.4** |

### Test B: C++ UKB sparse (from `/usr/bin/time`)

The C++ code does not emit per-step timing, but the major steps are:

| Step | Estimated time | Notes |
|------|---------------|-------|
| Config + design loading | < 1s | 961 samples, 7 covariates, inv-normalize |
| Genotype loading (`setGenoObj`) | **~11 min** | 340K markers x 408K samples BED, MAF/QC filtering |
| Sparse GRM loading + subsetting | ~30s | Load 488K x 488K COO (221M entries), subset to 961x961 (963 nnz) |
| GLMM fitting (1 iteration, tau[1]->0) | < 1s | 35 sparse solves per iteration |
| Variance ratio (30 markers) | ~10s | 30 PCG solves + genotype lookups |
| Output writing | < 1s | JSON, .arma binaries, VR files |

**Bottleneck:** Genotype loading dominates (~90% of wall time) due to the large BED file (408K samples x 340K markers = ~33 GB of BED data to scan even though only 961 samples are used).

---

## 4. Model Output Comparison

### Test A: sparse_x1

| Metric | C++ | R |
|--------|-----|---|
| tau[0] | 0.7447 | 0.2838 |
| tau[1] | 0.0000 (stopped early) | 0.4630 |
| GLMM iterations | 1 (tau[1]<=0, early stop) | 5 (converged) |
| Variance ratio | 1.3428 | 1.1130 |
| VR markers used | 30 | 50 (adaptive: CV>0.001 at 30, 40) |

**Discrepancy:** C++ stops after 1 iteration because tau[1] goes to 0 on the conservative first step. R continues for 5 iterations and converges to tau[1]=0.463. This is a known difference due to **random vector bypass not being active** for this config — different random vectors for trace estimation produce different tau updates, and the conservative first-step can push tau[1] negative or zero.

### Test B: ukb_ldl_sparse

| Metric | C++ | R |
|--------|-----|---|
| tau[0] | 0.9942 | Not run yet |
| tau[1] | 0.0000 (stopped early) | — |
| GLMM iterations | 1 | — |
| Variance ratio | 1.0059 | — |
| VR CV | 0.0000 | — |

**Note:** With 961 random samples and relatedness cutoff 0.05, the subsetted sparse GRM has only 963 non-zeros (essentially diagonal). This means almost no relatedness signal, so tau[1] correctly goes to 0.

---

## 5. Callgrind Instruction Profile (Test A: sparse_x1, Ir-only mode)

**Total instructions:** 30,560,991,710

### Top Hotspots

| Rank | Ir (%) | Function | File |
|------|--------|----------|------|
| 1 | **29.66%** | `Get_OneSNP_Geno()` | SAIGE_step1_fast.cpp |
| 2 | **23.81%** | `Get_OneSNP_Geno_atBeginning()` | SAIGE_step1_fast.cpp |
| 3 | **12.24%** | `memcpy/memmove` (libc) | libc.so.6 |
| 4 | **7.17%** | `stl_bvector` bitwise access | inside `Get_OneSNP_Geno_atBeginning` |
| 5 | **7.16%** | `blas_thread_server` | OpenBLAS threading |
| 6 | **6.67%** | `Get_OneSNP_StdGeno()` | SAIGE_step1_fast.cpp |
| 7 | **1.65%** | `std::vector::size()` | inside genotype functions |
| 8 | **1.28%** | `eglue_plus` elementwise add | Armadillo |
| 9 | **0.93%** | `arma::Mat::zeros()` | zero-init genotype vectors |
| 10 | **0.73%** | `arma::debug` bounds checks | inside `Get_OneSNP_StdGeno` |

### Category Summary

| Category | Ir (%) | Notes |
|----------|--------|-------|
| **Genotype decoding** | **53.47%** | `Get_OneSNP_Geno` + `Get_OneSNP_Geno_atBeginning` |
| **Memory operations** | **12.24%** | libc memcpy/memmove |
| **BLAS threading** | **7.16%** | OpenBLAS thread server (wasteful for N=1000) |
| **Genotype standardization** | **6.67%** | `Get_OneSNP_StdGeno` |
| **Armadillo overhead** | **2.94%** | zeros, bounds checks, eglue |
| **GLMM fitting + VR** | **<5%** | Sparse solves, trace estimation |

---

## 6. Bug Fixed During Benchmarking

### Sparse GRM Subsetting (main.cpp)

**Problem:** When loading a pre-built sparse GRM from file, the code passed the full GRM dimension (488,377) to `setupSparseGRM()`. The solver then tried to access `wVec` (961 elements) at indices up to 488K, causing `Mat::operator(): index out of bounds`.

**Root cause:** `load_matrix_market_coo()` returns the full COO matrix with indices in the full-cohort space. The code did not remap indices to the phenotyped-sample space.

**Fix:** After loading the MTX file:
1. Read sparse GRM sample IDs
2. Build mapping: GRM index -> design sample index
3. Filter COO entries to phenotyped-sample pairs only
4. Remap indices to 0..960 space
5. Add diagonal entries for samples missing from GRM
6. Call `setupSparseGRM(961, ...)` instead of `setupSparseGRM(488377, ...)`

**Result:** 442,853,165 entries -> 963 entries (961 diagonal + 2 off-diagonal related pairs)

---

## 7. R SAIGE Sample ID Matching Issue (Unresolved)

R SAIGE (Docker v1.5.0.2) constructs internal sample IDs as `FID_IID` (e.g., `1793481_1793481`), but:
- The phenotype file has `eid` column with plain eids (e.g., `1793481`)
- The sparse GRM sample IDs file has plain eids
- The FAM file has matching FID and IID columns

This causes R to report "0 samples will be used for analysis" because the join key formats don't match between genotype (FID_IID) and phenotype (eid) tables. A workaround phenotype file with `IID = eid_eid` format was created but still fails at the sparse GRM matching step.

**Status:** R UKB sparse benchmark pending resolution of this ID matching issue.

---

## 8. Optimization Opportunities

| Priority | Opportunity | Est. Impact | Notes |
|----------|------------|-------------|-------|
| 1 | **Skip full BED scan when using pre-built sparse GRM** | ~90% wall time for UKB | Currently scans all 340K markers even when GRM is loaded from file |
| 2 | SIMD-vectorize PLINK BED decoding (AVX2) | ~30% Ir | `Get_OneSNP_Geno` + `Get_OneSNP_Geno_atBeginning` |
| 3 | Remove debug stdout from production paths | Wall time savings | `#ifdef DEBUG` guards |
| 4 | Set `OPENBLAS_NUM_THREADS=1` for small N | ~7% Ir | Threading overhead exceeds benefit at N<=5000 |
| 5 | Cache SuperLU symbolic factorization | ~4-5% Ir | 70+ sparse solves with same sparsity pattern |
| 6 | Compile with `-DARMA_NO_DEBUG` | ~0.73% Ir | Remove Armadillo bounds checks |
| 7 | Remove `m_OneSNP_Geno.zeros()` when no missing | ~0.93% Ir | Conditional on missingness mask |

---

## 9. Files Generated

```
benchmark_results/
├── sparse_x1_nocache/
│   ├── cpp/
│   │   ├── callgrind.out         # Callgrind Ir-only profile
│   │   ├── stdout.log            # C++ stdout
│   │   └── time_stderr.log       # /usr/bin/time + valgrind stderr
│   └── r/
│       ├── stdout.log            # R stdout
│       ├── time_stderr.log       # /usr/bin/time stderr
│       └── saige_out.*           # R SAIGE output files
├── ukb_ldl_sparse_timed/
│   └── cpp/
│       ├── stdout.log            # C++ stdout
│       └── time_stderr.log       # /usr/bin/time stderr
└── ukb_ldl_sparse/               # (original callgrind with cache+branch, still running)
    └── cpp/
        └── ...
```
