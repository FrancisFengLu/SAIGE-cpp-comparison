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
