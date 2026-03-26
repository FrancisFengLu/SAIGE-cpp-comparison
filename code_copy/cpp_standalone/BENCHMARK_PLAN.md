# Benchmark Plan: R SAIGE vs C++ Standalone
## Memory & CPU Profiling with Valgrind

**Date:** 2026-03-24

---

## Overview

Compare R SAIGE (docker `wzhou88/saige:1.5.0.2`) against C++ standalone on:
- **Wall-clock time** — `/usr/bin/time -v`
- **Peak memory (RSS)** — `/usr/bin/time -v` → "Maximum resident set size"
- **CPU hotspots** — `valgrind --tool=callgrind` (C++ only; R profiled via wall time)

## Script

**`benchmark.sh`** — single script that runs all configs, both R and C++.

```
Usage:
  bash benchmark.sh all            # run all 5 config pairs
  bash benchmark.sh sparse_x1      # run one config only
  bash benchmark.sh dense_x1x2     # etc.
```

Output goes to `benchmark_results/<config>/{cpp,r}/`:
```
benchmark_results/
├── sparse_x1/
│   ├── cpp/
│   │   ├── callgrind.out      ← callgrind_annotate or kcachegrind
│   │   ├── stdout.log
│   │   ├── time_stderr.log    ← full /usr/bin/time -v output
│   │   └── time_summary.txt   ← wall clock + peak RSS extracted
│   └── r/
│       ├── stdout.log
│       ├── time_stderr.log
│       ├── time_summary.txt
│       └── saige_out.*        ← R SAIGE output files
├── sparse_x1x2/
│   └── ...
├── dense_x1/
│   └── ...
├── dense_x1x2/
│   └── ...
└── ukb_ldl/
    └── ...
```

At the end, a summary table is printed:
```
Config          Lang    Wall(s)    Peak RSS(KB)
-------------------------------------------------------
sparse_x1       cpp     0:02.15        125000
sparse_x1       r       0:08.43        340000
...
```

---

## Config Mapping: C++ YAML → R SAIGE CLI

Each C++ config YAML has a matching R `step1_fitNULLGLMM.R` command. Key parameter mappings:

### Common parameters (all configs)

| C++ YAML field | R CLI argument | Value |
|---|---|---|
| `fit.trait` | `--traitType` | `quantitative` |
| `fit.loco` | `--LOCO` | `FALSE` |
| `fit.nthreads` | `--nThreads` | `1` |
| `fit.tol` | `--tol` | `0.02` |
| `fit.tolPCG` | `--tolPCG` | `1e-5` |
| `fit.maxiterPCG` | `--maxiterPCG` | `500` |
| `fit.maxiter` | `--maxiter` | `20` |
| `fit.nrun` | (internal, default 30) | `30` |
| `fit.traceCVcutoff` | `--traceCVcutoff` | `0.0025` |
| `fit.num_markers_for_vr` | `--numMarkersForVarRatio` | `30` |
| `fit.relatedness_cutoff` | `--relatednessCutoff` | `0.125` |
| `paths.overwrite_varratio` | `--IsOverwriteVarianceRatioFile` | `TRUE` |

### Per-config differences

| Config | C++ design CSV | R `--phenoFile` | R `--phenoCol` | R `--covarColList` | Sparse GRM? |
|---|---|---|---|---|---|
| `sparse_x1` | `design_dense_x1.csv` | `pheno_1000samples.txt_withdosages_withBothTraitTypes.txt` | `y_quantitative` | `x1` | YES |
| `sparse_x1x2` | `design_test.csv` | same | `y_quantitative` | `x1,x2` | YES |
| `dense_x1` | `design_dense_x1.csv` | same | `y_quantitative` | `x1` | NO |
| `dense_x1x2` | `design_test.csv` | same | `y_quantitative` | `x1,x2` | NO |
| `ukb_ldl` | UKB 1000-sample TSV | same TSV | `f.30780.0.0` | `Sex,Age,Batch,PC1-4` | NO |

### Sparse GRM configs add:

| C++ YAML field | R CLI argument |
|---|---|
| `fit.use_sparse_grm_to_fit: true` | `--useSparseGRMtoFitNULL=TRUE` |
| `fit.use_pcg_with_sparse_grm: false` | `--usePCGwithSparseGRM=FALSE` |
| `paths.sparse_grm` | `--sparseGRMFile=...` |
| `paths.sparse_grm_ids` | `--sparseGRMSampleIDFile=...` |

### UKB LDL config adds:

| C++ YAML field | R CLI argument |
|---|---|
| `fit.inv_normalize: true` | `--invNormalize=TRUE` |
| `design.q_covar_cols: [Sex]` | `--qCovarCol=Sex` |
| `design.iid_col: eid` | `--sampleIDColinphenoFile=eid` |

---

## Data Summary

| Dataset | Samples | Markers | Plink prefix |
|---|---|---|---|
| Test (small) | 1,000 | 128,868 | `nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr` |
| UKB LDL | ~1,000 (random subset) | 340,447 | `ukb_allchr_v2_newID_...pruned` |

---

## What to Compare

### 1. Wall-clock time
- Extracted from `/usr/bin/time -v` ("Elapsed wall clock time")
- Note: R in docker has container startup overhead (~1-2s). For small datasets this matters; for UKB it's negligible.

### 2. Peak memory
- Extracted from `/usr/bin/time -v` ("Maximum resident set size" in KB)
- R includes R runtime + loaded packages (~150-200 MB baseline)
- C++ has no runtime overhead beyond libc + linked libraries

### 3. CPU hotspots (C++ only)
- `callgrind.out` files can be viewed with:
  ```
  callgrind_annotate benchmark_results/<config>/cpp/callgrind.out
  ```
  or interactively with `kcachegrind`.
- Key functions to watch:
  - `getCoefficients_cpp` — PCG solver (inner IRLS)
  - `GetTrace_q` / `GetTrace` — Monte Carlo trace estimation
  - `fitglmmaiRPCG` — one AI-REML iteration
  - `compute_variance_ratio` — VR marker loop
  - `gen_spsolve_v4` — direct sparse solve

### 4. R profiling (optional, beyond this script)
- For R CPU profiling, you can use `Rprof()` inside a custom R script, or:
  ```
  docker run --rm -v /media/leelabsg-storage0:/media/leelabsg-storage0 \
      wzhou88/saige:1.5.0.2 \
      Rscript -e 'Rprof("prof.out"); source("/usr/local/bin/step1_fitNULLGLMM.R"); Rprof(NULL)'
  ```
  Then analyze with `summaryRprof("prof.out")`.

---

## How to Run

```bash
# Full benchmark (all 5 configs, R + C++)
cd /media/leelabsg-storage0/seokho/SAIGE-cpp-comparison/code_copy/cpp_standalone_update
bash benchmark.sh all

# Single config
bash benchmark.sh sparse_x1

# After completion, view callgrind for a specific config:
callgrind_annotate benchmark_results/sparse_x1/cpp/callgrind.out | head -80

# Or open in kcachegrind (GUI):
kcachegrind benchmark_results/sparse_x1/cpp/callgrind.out
```

---

## Notes

- C++ is built with `make profile` (`-O1 -g`) for meaningful callgrind output. Production builds use `-O2`.
- R runs inside docker; `nThreads=1` to match C++.
- All configs use `LOCO=FALSE` to keep the comparison focused on the null model fitting + VR.
- The UKB config uses a 1000-sample random subset, not the full 400K cohort.
