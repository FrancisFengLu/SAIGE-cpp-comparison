# R parity gate (step 2, R SAIGE 1.5.2)

Math-consistency gate, not byte identity: every step-2 path of this port
(CPU single-trait scalar, CPU multi-trait batch, GPU multi-trait) against R
SAIGE 1.5.2 on the SAME null models (R's step 1, converted with
`step2_saige-step2/tools/rda_to_arma.R`), under R's defaults and under
`is_noadjCov=FALSE` + `impute_method=mean`, Firth off and on.

Datasets (simulated; `RDEF`, default `/opt/saige/logs/rdefaults`):

| name | N / M | traits | covariates | what it exercises |
|---|---|---|---|---|
| audit | 3000 / 4000 | 4 binary, prevalence 4.7–50% | 3 | MAC <= 4 (exact test), planted signals |
| bvs | 5000 / 5000 | 4 binary, each with its own missing phenotypes | 4 | different sample sets, 1671 flipped markers, 5% missing genotypes |
| qt12 | 4000 / 5000 | 2 quantitative + 2 binary, two with missing phenotypes | 13 | quantitative traits, p = 13, ~half the markers flipped, best_guess on 1% missing genotypes, MAC <= 4 |

Scripts:

- `gen_qt12.py OUTDIR` — the qt12 data (audit and bvs come from earlier work).
- `r_runs.sh` — R step 1 (qt12) and R step 2 for every dataset x configuration
  (`def` = R defaults, `defF` = + Firth, `adj` = `--is_noadjCov=FALSE --impute_method=mean`,
  `adjF` = adj + Firth). Needs the R environment of `optimization/comp8/env_r.sh`.
- `cpp_runs.sh <saige-step2>` — converts the models and runs the three paths on
  every dataset x configuration (`single`: the legacy one-model config, i.e. the
  P = 1 scalar path; `multi`: all traits in one `models:` run; `gpu`: the same with
  `useGPU: true`).
- `cmp_r.py A B` — two result files row by row: max relative difference of
  p.value / BETA / SE / Tstat / var, rows beyond the printed digits, rows
  crossing 5e-8 / 1e-5 on one side only, Is.SPA and Firth-route mismatches.
- `gate_table.py` — the Markdown tables (C++ vs R, and multi / GPU vs single).

Write-up with the numbers: `SAIGE-work/optimization/torchgwas2/S2_RDEFAULTS.md`.

Own sample lists on the device (`gpuDeviceStats` with `gpuOwnSampleSets`, write-up
`S2_OWNSTATS.md`):

- `gate_bm.sh <new-bin> <base-bin>` — bingpu_test bm (8 binary traits, 4 distinct
  missing-phenotype patterns, N = 50,000, 26,060 markers, 3,640 with missing calls)
  under the four configurations: the new binary's GPU multi-trait run against its own
  CPU single-trait runs and against the base binary's GPU run (host tail), route dumps on.
  `bm_table.py` prints the table. (bm's models come from this port's step 1, so there is
  no R reference; the R comparison for own sample lists is bvs and qt12 above.)
- `timing_ownstats.sh` — g200k x the bm_full models cycled to P = 128 (own sample
  lists), before / after, cold cache, n = 2, TIMING_LOCK; then bt under mean imputation,
  where the missing-cell sums run on the device for every marker with missing calls.

Per-marker preprocessing on the device, quantitative device statistics, the sparse
statistic on the device (branch `gpuprep`, write-up `S2_GPUPREP.md`):

- `cpp_gpu_runs.sh <bin> <OUT>` / `gate_table_gpu.py <OUT>` — the GPU path of one binary on
  the rdefaults datasets x configurations, against R and against the rdefaults single-trait
  outputs, with the per-trait marker counts and the device log lines.
- `gen_fam32.py OUTDIR` — fam32: N = 4,000 in 1,000 families of 4 (Mendelian transmission, so a
  sparse GRM has 4 x 4 blocks), M = 5,000, 1% missing genotypes, 12 covariates, 16 binary + 16
  quantitative traits, EVERY trait with its own missing-phenotype pattern (5–25% NA).
- `r_runs_fam32.sh` — R step 1 (full GRM for the 32; createSparseGRM + sparse step 1 with the
  sparse and categorical variance ratios for b1–b4 / q1–q4) and R step 2 (full: def / adj /
  defF / adjF; sparse: fast test off and on, Firth).
- `cpp_runs_fam32.sh <bin>` — converts the models (`rda_to_arma.R`, with the sparse GRM for the
  sparse ones) and runs the three paths; `gate_table_fam32.py` prints the tables.
- `sparse_gate.sh <bin>` / `sparse_table.py` — the sparse-GRM cases of bingpu_test (bt: 10
  binary, one list; bm: 8 binary own lists; qm: 8 quantitative own lists; mix), fast test off /
  on, Firth, noadjCov both: the GPU path against the CPU block-inverse path of the same binary
  (`blockSparseSigmaSolveMargin: 0.01`, so the CPU never falls back to the 2%-tolerance PCG).
- `subset_model.py` — a model directory with a random subset of its samples dropped, for timing
  runs with 128 distinct sample lists (not a consistent fit; timing only).
- `timing_gpuprep.sh` — g200k P = 128 (same list / 4 lists / 128 lists) and the quantitative
  sparse first pass, before (ownstats) vs after.

Covariate count on the device SPA / Firth (branch `covlimit`, write-up `S2_COVLIMIT.md`): the device
SPA (lib and own, fp64 and fp32), the device Firth (fp64 and fp32) and the own-list device statistics
take up to 64 covariate columns incl. the intercept (were 8 / 8 / 32).

- `gen_covp.py OUTDIR [N] [M] [seed]` -- covp: N = M = 5,000, 40 covariates (sex, age-like, 38 PC-like),
  4 binary traits (prevalence 5-40%), b1 / b2 complete, b3 / b4 with their own missing-phenotype
  patterns; covt is the same at N = M = 20,000 (timing).
- `r_runs_covp.sh [s1|s2|all]` -- R step 1 with the first p - 1 covariates, p = 4 9 13 24 40 (covp
  b1..b4, covt b1 / b2), and R step 2 on covp for every p x def / defF / adj / adjF x trait.
- `cpp_runs_covp.sh <bin> <OUT> [variants]` -- converts the models and runs, per p x configuration:
  `single` (CPU scalar, per trait), `same` (GPU, b1 + b2, one list), `own` (GPU, b1..b4, own lists),
  `sameown` (same with `gpuSpaImpl: own`), `same32` / `own32` (`gpuPrecisionSPA` / `gpuPrecisionFirth`
  fp32). `gate_table_covp.py <OUT>` prints the tables (vs R, vs single, markers tested, the device
  log lines with the SPA / Firth pair counts and us per pair).
- `timing_covp.sh <new> <base> <OUT>` -- covt, defF, device SPA / Firth us per pair: base (gpuprep)
  at p = 4, new at every p (lib fp64 / own fp64 / lib fp32), n = 2, CUDA_MODULE_LOADING=EAGER;
  `timing_table_covp.py <OUT>` prints the table.
