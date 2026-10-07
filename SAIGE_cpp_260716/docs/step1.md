---
title: "Step 1: fit the null model"
nav_order: 2
---

# Step 1: fit the null model

`saige-null` fits one null model per trait (no genotype effect) and estimates
the variance ratio that step 2 needs. All settings are in a YAML config.

```bash
saige-null -c step1.yaml                       # run
saige-null -c step1.yaml --dry-run             # check inputs only, fit nothing
saige-null -c step1.yaml -t 8                  # threads (= fit.nthreads)
saige-null -c step1.yaml -o fit.loco=false     # override any config key (repeatable)
```

`R_HOME` must point to the R installation of the build environment
(`export R_HOME=$CONDA_PREFIX/lib/R`).

## Command-line options

| Option | Meaning |
|---|---|
| `-c, --config FILE` | YAML config (required) |
| `-d, --design FILE` | phenotype file; overrides `design.csv` |
| `-o, --override KEY=VALUE` | set a config key, dotted path, e.g. `-o paths.out_prefix=run2/models`; repeat for several |
| `-t, --threads N` | `fit.nthreads`; 0 = all cores |
| `--gpu` | same as `fit.use_gpu: true` |
| `--dry-run` | read and check the config, phenotype file and `.fam`, print a summary, exit |
| `-h, --help` | option list |

## Important inputs

The config has three sections: `paths`, `design`, `fit`. Relative paths under
`paths:` are taken relative to the config file's directory.

**`paths:`**

- `plinkFile` — PLINK prefix (`.bed/.bim/.fam`). Alternatively `bed`, `bim`, `fam` separately. Required (also with a sparse GRM: the variance-ratio markers are read from it).
- `out_prefix` — where the model goes. One trait: the model directory itself. Several traits: a directory holding one model directory per trait (`<out_prefix>/<trait>/`).
- `out_prefix_vr` — file prefix of the variance-ratio file. One trait: `<out_prefix_vr>.varianceRatio.txt`; several traits: `<out_prefix_vr>_<trait>.varianceRatio.txt`. Default: `out_prefix`.
- `overwrite_varratio` — `true` to allow overwriting an existing variance-ratio file (default `false`: the run stops).
- `sparse_grm`, `sparse_grm_ids` — sparse GRM (`.mtx`) and its sample-ID list. Read if both files exist, written if they do not and the run builds a sparse GRM (see [Sparse GRM](#sparse-grm)).

**`design:`**

- `csv` — phenotype/covariate file (format below).
- `iid_col` — sample-ID column (default `IID`).
- `y_col` — the trait column, for one trait (default `y`).
- `y_cols` — list of trait columns, for several traits in one run. Mutually exclusive with `y_col`.
- `covar_cols` — list of covariate columns (default: none). An intercept is always added.
- `q_covar_cols` — the categorical covariates among `covar_cols` (one-hot coded, first level dropped).
- `whitelist_ids` — file with one sample ID per line; only these samples are used.
- `intersect_samples` — several traits: fit every trait on the samples non-missing in all of them (default `false`; see [Several traits](#several-traits-in-one-run)).
- `sex_col`, `female_only`, `male_only`, `female_code` (default `1`), `male_code` (default `0`) — sex-specific fit.

**`fit:`**

- `trait` — `binary` or `quantitative`. Default `binary`.
- `loco` — leave-one-chromosome-out models. Default **`true`**; switched off automatically (with a log line) when the `.bim` has fewer than 2 autosomes or the fit uses a sparse GRM.
- `nthreads` — threads (default 1).
- `use_gpu` — run the GRM products on the GPU (default `false`; [GPU](gpu.md)).
- `inv_normalize` — rank-based inverse-normal transform of a quantitative trait (default `false`).
- `min_maf_grm` — minimum MAF of the markers used for the full GRM (default 0.01); `max_miss_grm` — maximum missing rate (default 0.15).
- `num_markers_for_vr` — markers drawn for the variance ratio (default 30).
- `use_sparse_grm_to_fit` — fit on the sparse GRM instead of the full GRM (default `false`). Forces `nthreads: 1` and `loco: false`.
- `use_sparse_grm_for_vr` — use the sparse GRM for the variance ratio (default `false`).
- `make_sparse_grm_only` — build the sparse GRM, write it, stop (default `false`).
- `relatedness_cutoff` — sparse GRM entries below this are dropped when building it (default 0.05).
- `tol` (0.02), `maxiter` (20), `tolPCG` (1e-5), `maxiterPCG` (500), `nrun` (30), `trace_seed` — fitting controls, SAIGE's meaning and defaults.

**Stored in the model and used by step 2** (step 2 can override the two Firth keys):

- `spa_cutoff` — saddle-point approximation is used when |test statistic| > this (default 2).
- `fast_test` — fast score test with re-computation of markers with p < 0.05 (default `true`).
- `firth_beta` — Firth-corrected effect sizes for binary traits (default: `true` for binary traits).
- `p_cutoff_for_firth` — Firth is applied when p < this (default 0.01).
- `impute_method` — missing genotypes in step 2: `mean` (default).

## Example phenotype file

Header line, one row per sample. Separator: tab, space or comma (detected from
the header). Missing values: empty, `NA` or `NaN`.

| IID | b1 | b2 | b3 | b4 | q1 | q2 | x1 | x2 |
|---|---|---|---|---|---|---|---|---|
| per0 | 0 | 1 | 0 | 0 | -1.0620 | 0.3673 | -0.4378 | 1 |
| per1 | 0 | 1 | 0 | NA | 0.1620 | 1.3208 | 0.4618 | 0 |

Rules (checked at start-up, the run stops with a message otherwise):

- Binary traits are coded **0 = control, 1 = case**. Any other value (e.g. 1/2) stops the run.
- A quantitative trait with variance < 0.1 stops the run; use `inv_normalize: true` or rescale.
- Every sample in the phenotype file must be in the `.fam` (matched on the `.fam` IID, column 2). A sample that is not stops the run (`IID in design not found in FAM`). To use a phenotype file with extra samples, set `design.whitelist_ids` to the `.fam` IIDs (`cut -f2 geno.fam > ids.txt`).
- `.fam` samples without a phenotype row are left out of the model.
- A row with a missing value in the trait or in any listed covariate is dropped for that trait.
- Duplicate IIDs: the first row is kept, with a warning.

## Example genotype input

PLINK 1 binary files. Hard calls, missing calls allowed (markers above
`max_miss_grm` are not used for the GRM). LOCO needs real chromosome codes
in `.bim` column 1.

| .bim: chr | id | cM | pos | A1 | A2 |
|---|---|---|---|---|---|
| 1 | snp1 | 0 | 1000 | C | A |

## Example config: four binary traits, full GRM, GPU

From [`examples/03_step1_binary.sh`](examples/03_step1_binary.sh)
(`$D` = data directory, `$O` = output directory):

```yaml
paths:
  plinkFile: $D/geno                 # .bed/.bim/.fam prefix
  out_prefix: $O/models              # one model directory per trait: $O/models/<trait>/
  out_prefix_vr: $O/vr               # variance ratios: $O/vr_<trait>.varianceRatio.txt
  overwrite_varratio: true           # allow re-running into the same place
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [b1, b2, b3, b4]
  covar_cols: [x1, x2]
fit:
  trait: binary
  loco: false
  nthreads: 8
  use_gpu: true
  firth_beta: true                   # stored in the model; step 2 can override
  p_cutoff_for_firth: 0.01
  spa_cutoff: 2.0
```

```bash
saige-null -c step1.yaml > step1.log 2>&1
```

## Example config: quantitative traits with LOCO

From [`examples/05_step1_quant_loco.sh`](examples/05_step1_quant_loco.sh):

```yaml
paths:
  plinkFile: $D/geno
  out_prefix: $O/models
  out_prefix_vr: $O/vr
  overwrite_varratio: true
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [q1, q2]
  covar_cols: [x1, x2]
fit:
  trait: quantitative
  loco: true                         # needs >= 2 autosomes in the .bim
  inv_normalize: false               # true: rank-normalise the phenotype first
  nthreads: 8
  use_gpu: true
```

The log confirms LOCO with `LOCO: on  LowMem: no  chroms: 1 2`, and each model
directory gets one `chr<N>/` subdirectory per chromosome. Step 2 is then run once
per chromosome ([Step 2, LOCO](step2.md#loco)).

## One trait

```yaml
paths:
  plinkFile: $D/geno
  out_prefix: $O/q2                  # model directory
design:
  csv: $D/pheno.txt
  y_col: q2
fit:
  trait: quantitative
  loco: false
```

writes `$O/q2/` (model) and `$O/q2.varianceRatio.txt`.

## Several traits in one run

`design.y_cols: [...]` fits each listed trait with the same covariates and
settings (all traits of a run share `fit.trait`; run binary and quantitative
traits as two step-1 runs). Each trait is fitted on its own non-missing
samples. Traits with identical sample sets share one genotype load; the log
shows the grouping:

```
[multi-pheno] 4 phenotypes fall into 2 sample-set groups; each group gets its own genotype load, GRM and GPU upload:
    group 1/2: n=5000  traits(3): b1 b2 b3
    group 2/2: n=4511  traits(1): b4
```

Each trait's result is the same as fitting it alone. `design.intersect_samples: true`
instead fits all traits on the common samples (fewer samples per trait; a warning
says so).

To choose each trait's output paths, use a top-level `models:` list instead of
`y_cols` (mutually exclusive):

```yaml
models:
  - y_col: q1
    out_prefix: $X/s1m/q1_model      # model directory
    out_prefix_vr: $X/s1m/q1         # -> $X/s1m/q1.varianceRatio.txt
```

## Sparse GRM

From [`examples/07_sparse_grm.sh`](examples/07_sparse_grm.sh). First build the
sparse GRM from the genotypes (once per cohort):

```yaml
paths:
  plinkFile: $D/geno
  out_prefix: $O/grm_run
  sparse_grm: $O/grm.mtx             # written
  sparse_grm_ids: $O/grm.ids         # written
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_col: b1
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  make_sparse_grm_only: true
  relatedness_cutoff: 0.05           # entries below this are dropped
  min_maf_grm: 0.01
```

Then fit on it (the files exist, so they are read):

```yaml
paths:
  plinkFile: $D/geno                 # still needed: variance-ratio markers come from here
  out_prefix: $O/models
  out_prefix_vr: $O/vr
  sparse_grm: $O/grm.mtx
  sparse_grm_ids: $O/grm.ids
  overwrite_varratio: true
design:
  csv: $D/pheno.txt
  iid_col: IID
  y_cols: [b1, b2, b3]
  covar_cols: [x1, x2]
fit:
  trait: binary
  use_sparse_grm_to_fit: true
  use_sparse_grm_for_vr: true
  fast_test: true
```

A sparse GRM made elsewhere can be used the same way:

- `.mtx`: MatrixMarket `coordinate real`, `general` (both triangles) or `symmetric`, 1-based, square.
- `.ids`: one sample ID per line, in the matrix's row order, matching the `.fam` IIDs.

```
%%MatrixMarket matrix coordinate real general
%
5000 5000 388682
```

The model is fitted on the samples present in the GRM, the phenotype file and the
`.fam`. With several traits on different sample sets, the sparse GRM files must
already exist (build them once on the whole cohort).

## Outputs

For each trait (paths as configured above):

| File | Contents |
|---|---|
| `<model dir>/nullmodel.json` + `*.arma` | the null model; this directory is step 2's `modelFile` |
| `<model dir>/chr<N>/` | per-chromosome LOCO model (only with LOCO) |
| `<vr prefix>.varianceRatio.txt` | variance ratio; step 2's `varianceRatioFile` |
| `<vr prefix>.30markers.SAIGE.results.txt` | the markers used for the variance ratio |
| `<model dir>.grm_diag.txt`, `<model dir>_cov_*.csv` | diagnostics, not used by step 2 |

The end of the log, per trait:

```
Converged: yes
Model artifact: .../models/b1/nullmodel.json
Variance ratio: .../vr_b1.varianceRatio.txt
```

`nullmodel.json` holds the estimates (`theta` = variance components, `alpha` =
covariate effects) and the step-2 settings (`SPA_Cutoff`, `isFastTest`,
`is_Firth_beta`, `pCutoffforFirth`, `loco`).
