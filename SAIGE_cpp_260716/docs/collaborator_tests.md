---
title: Collaborator test protocol
nav_order: 7
---

# Test protocol for collaborators with biobank data

This page is for a collaborator who runs our C++/GPU SAIGE step 2 on real
biobank data that we cannot see. You run a fixed set of tests on one chromosome
and send back **only logs and summary numbers**. Every script is in
[`collab/`](collab/), and every command on this page was run end to end on
simulated data before it was written down (see [What was verified](#what-was-verified)).

**Data compliance.** Individual-level data (genotypes, phenotypes, sample IDs)
and per-variant results never leave your environment. The tar file that
[`collect_results.sh`](#7-collect-and-send-the-results) produces contains logs,
timings, counts, md5 checksums of result files and aggregate comparison
statistics only, and it is scanned for sample and marker IDs before it is
written. Please look at its contents listing before sending it.

Only binary traits are tested.

## Contents

1. [What is tested](#what-is-tested)
2. [Prerequisites and environment record](#1-prerequisites-and-environment-record)
3. [Build](#2-build)
4. [Install self-check](#3-install-self-check)
5. [Prepare the inputs](#4-prepare-the-inputs)
6. [Step 1: the null models](#5-step-1-the-null-models)
7. [Run the test matrix](#6-run-the-test-matrix)
8. [Collect and send the results](#7-collect-and-send-the-results)
9. [Result files we receive](#8-result-files-we-receive)
10. [Run time and what to drop first](#9-run-time-and-what-to-drop-first)

## What is tested

All tests use step-2 genotypes of **one chromosome** (so the CPU runs finish), the
null models from two step-1 runs, and your binary traits. "P" is the number of
traits tested together in one step-2 run (the first P traits of your list).

| Block | Runs | What varies | Question |
|---|---|---|---|
| `main` | 48 | path {CPU, GPU} x P {1, 8, 32, 128} x GRM {full; sparse, fast test off; sparse, fast test on} x Firth {off, on} | speed, memory, CPU vs GPU byte-identity |
| `rare` | 2 | CPU and GPU on the markers with MAC < 20; full GRM, Firth on, P = top level | how many pairs go to the exact test (ER), identity on that path |
| `stage` | 2 | CPU build with stage timers + GPU with `gpuOverlap: false`; full GRM, Firth on, P = top level | where the time goes: read, scan, host tail, SPA, ER, Firth, write, with pair counts |
| `pgen` | 1 | GPU on a hard-call PGEN of the same chromosome; full GRM, Firth on, P = top level (only if you have PGEN) | PGEN input on the GPU |
| `vsr` | 8 + 8 | ours (CPU) and R SAIGE 1.5.2: P {1, 8} x GRM {full; sparse, fast test off} x Firth {off, on} | agreement with R: p.value, BETA, SE |
| `precision` | 8 | GPU, full GRM, Firth on, P = top level: all fp64 (default), scan fp32, scan int8, SPA fp32, ER fp32, Firth fp32, all fp32, scan int8 + SPA / ER / Firth fp32 | speed of the [GPU precision modes](gpu.md#precision) on your card, and their difference from all-fp64 |

P levels at or above your number of binary traits are dropped, and the top level
uses all your binary traits (with 60 traits: P = 1, 8, 32, 60; with 300 traits:
1, 8, 32, 300). With at least 128 binary traits and a PGEN file the matrix is
**77 runs**: 48 + 2 + 2 + 1 + 8 + 8 of ours + 8 of R (an R "run" at P = 8 is 8 R
processes at once, one per trait).

"Fast test" is SAIGE's `isFastTest` for sparse-GRM models: off means every marker
uses the sparse-GRM variance; on means the variance without the GRM is used first
and markers with p < 0.05 are recomputed with the sparse GRM. With a full-GRM
model it has no effect (reported as `na`).

## 1. Prerequisites and environment record

- Linux x86-64, an NVIDIA GPU (compute capability 7.0 or newer), CUDA 12.x, and
  the conda build environment from [Home, Installation](index.md#installation).
- `plink2` (to cut one chromosome and the rare-marker set), GNU `time` at
  `/usr/bin/time`, `python3` (standard library only), `nvidia-smi`.
- For block `vsr`: R with SAIGE **1.5.2** installed, and the path to its
  `extdata/step2_SPAtests.R`.
- Disk: the result files of one run are deleted after they are compared, but the
  largest run (CPU and GPU of P = 128) keeps two sets at once: about
  2 x 128 x (markers on the chromosome) x 200 bytes. Keep 20 GB more than that free;
  the driver stops when less than `MIN_FREE_GB` (20) is left.
- Exclusive use of the node during the runs, if at all possible.

Clone the repository and set the variables used on the rest of this page:

```bash
git clone -b docs-usage https://github.com/FrancisFengLu/SAIGE-cpp-comparison.git
cd SAIGE-cpp-comparison/SAIGE_cpp_260716
export CONDA_ENV=$HOME/miniforge3/envs/saige-build   # the build environment (Home, Installation)
export CUDA_HOME=/usr/local/cuda
export OUT_ROOT=/fast/disk/saige_tests              # everything is written under here
export GENO=/path/chr21                             # step 2: PLINK prefix of ONE chromosome (section 4)
```

Record the machine (hardware, OS, GPU, CUDA, filesystems, tools) and fill in the
few lines it cannot detect:

```bash
bash docs/collab/env_record.sh            # writes $OUT_ROOT/env/environment.txt
nano $OUT_ROOT/env/manual.txt              # site, machine type, disk type, N, chromosome, ...
```

`manual.txt` asks for: site, machine type, whether the node was exclusive, the
disk holding the genotypes (local NVMe / SSD / HDD / Lustre / GPFS / NFS), N,
the chromosome and its number of markers, number of traits and the range of case
fractions, how the sparse GRM was made.

## 2. Build

```bash
bash docs/collab/build_bins.sh 80          # compute capability: 70 V100, 80 A100, 89 L4, 90 H100
```

It runs [`examples/01_build.sh`](examples/01_build.sh) `gpu <SM>` and then a second
step-2 build with stage timers (`make PHASE_TIMING=1`), and copies to
`collab_bin/`:

| Binary | Used for |
|---|---|
| `saige-null` | step 1 (CUDA build) |
| `saige-step2` | every step-2 run except block `stage`'s CPU run. CPU-path runs use this same binary with `useGPU: false` |
| `saige-step2.phase` | block `stage`, CPU run (CPU-only build with per-stage timers) |

`collab_bin/BUILD_INFO.txt` records the commit, compiler, CUDA and the binaries'
md5s. Other build questions: [Home, Build](index.md#3-build).

## 3. Install self-check

About two minutes, simulated data only (2,000 samples x 3,000 markers from
`plink2 --dummy` with a fixed seed, 4 binary traits):

```bash
bash docs/collab/selfcheck.sh             # writes $OUT_ROOT/selfcheck/selfcheck_result.txt
```

It fits full-GRM and sparse-GRM null models (step 1 on the CPU, one thread),
runs step 2 on the CPU path and the GPU path with Firth on, and checks:

- **(a)** the GPU was used and its result files are byte-identical to the CPU
  path's. This must pass on every machine; if it fails, stop and send us
  `selfcheck_result.txt` and the logs under `$OUT_ROOT/selfcheck/s2/`.
- **(b)** md5s of the data and of the CPU results against the ones recorded when
  this page was written (GCP n1-standard-8, Intel Xeon 2.3 GHz with AVX2 and FMA, no AVX-512,
  Tesla V100, CUDA 12.9, plink2 v2.0.0-a.6.5LM). The data md5s need the same
  plink2 version. The result md5s need the same source commit and the same CPU
  instruction set: the build uses `-march=native`, so a CPU with different vector
  units (e.g. AVX-512, or an AMD CPU) may change the last digits
  of the results. A (b) mismatch with (a) passing is expected in that case and is
  not an error.

Expected end of `selfcheck_result.txt` on the reference machine:

```
PASS (a) full: CPU and GPU outputs byte-identical
PASS (a) sparse: CPU and GPU outputs byte-identical
PASS (b) 12 of 12 md5s match the recorded ones
SELF-CHECK OK
```

Step 1 runs single-threaded here on purpose: a multi-threaded CPU fit (and a GPU
fit) sums in an order that varies between runs, which changes the last digits of
the null model and therefore of the step-2 results. Step 2 itself gives the same
bytes for any `nThreads`.

### Optional: rehearse the whole protocol on simulated data

A few minutes (about 3 on the reference machine). This runs sections 5 to 7 on a simulated cohort
(4,000 samples, half in sibling pairs, 6,000 markers on two chromosomes, 30% of
them rare, 8 binary traits) so you see every output before touching real data.
Set `R_STEP2` and `RSCRIPT` as in [section 6](#6-run-the-test-matrix) first, or
add `BLOCKS="main rare stage pgen precision"` to skip R.

```bash
D=$OUT_ROOT/rehearsal; mkdir -p $D/data
python3 docs/collab/sim_cohort.py geno $D/data/cohort --n 4000 --m 6000 --chroms 2 --sib-frac 0.5
python3 docs/collab/sim_cohort.py pheno $D/data/cohort.fam $D/data/pheno.txt --nbin 8 \
    --score $D/data/cohort.pheno_score --gscale 1.2
plink2 --bfile $D/data/cohort --chr 1 --make-bed --out $D/data/chr1
plink2 --bfile $D/data/chr1 --make-pgen --out $D/data/chr1
( export OUT_ROOT=$D/run STEP1_BFILE=$D/data/cohort PHENO=$D/data/pheno.txt \
         BIN_TRAITS="b1 b2 b3 b4 b5 b6 b7 b8" COVARS="x1 x2" RELATEDNESS_CUTOFF=0.125 \
         GENO=$D/data/chr1 PGEN=$D/data/chr1
  bash docs/collab/step1_models.sh && bash docs/collab/run_matrix.sh && bash docs/collab/collect_results.sh )
```

With 8 traits the P levels are 1 and 8 (53 runs). Expected: every
`cpu_vs_gpu` and the `pgen_vs_bed` row of `compare.csv` has `identical` True; the
`cpp_vs_R` rows are identical for the full GRM and differ slightly for the sparse
GRM ([why](#things-you-may-see)). The relatedness cutoff is raised to 0.125 here
because with 3,500 GRM markers the default 0.05 keeps many chance entries between
unrelated samples (see the first item of [Things you may see](#things-you-may-see)).

## 4. Prepare the inputs

| What | Variable | Notes |
|---|---|---|
| genome-wide PLINK genotypes for the GRM | `STEP1_BFILE` | the usual step-1 input; hard calls ([Step 1](step1.md#example-genotype-input)) |
| phenotype / covariate file | `PHENO`, `IID_COL` (default `IID`) | binary traits coded 0/1 ([Step 1](step1.md#example-phenotype-file)) |
| binary trait names | `BIN_TRAITS` | space separated, or `@file` with one name per line. **The order matters**: the P = 8 runs use the first 8 |
| covariates | `COVARS`, `QCOVARS` (categorical ones among them) | space separated |
| step-2 genotypes, one chromosome | `GENO` | PLINK `.bed/.bim/.fam`; IIDs must match the step-1 samples |
| same chromosome as PGEN (optional) | `PGEN` | hard calls only; only for block `pgen` |
| existing sparse GRM (optional) | `SPARSE_GRM`, `SPARSE_GRM_IDS` | e.g. from FastSparseGRM; otherwise one is built |

Pick a small chromosome (see [run time](#9-run-time-and-what-to-drop-first)) and cut
it out with plink2:

```bash
plink2 --bfile /path/all_chrs --chr 21 --make-bed --out /path/chr21
plink2 --bfile /path/chr21 --make-pgen --out /path/chr21           # optional, for block pgen
```

If your step-2 genotypes are imputed dosages (BGEN), make hard calls first, e.g.
`plink2 --bgen chr21.bgen ref-first --sample chr21.sample --mac 1 --make-bed --out chr21`
(plink2 rounds dosages to hard calls); the GPU path reads hard calls only
([GPU](gpu.md#what-runs-on-the-gpu)).

To order traits so that the first 8 and 32 are a sensible mix, put a few common
and a few rare-case traits first. Any order is fine as long as it is recorded
(the trait names go into the logs; nothing else about them does).

## 5. Step 1: the null models

**Two step-1 runs are needed**, each fitting **all** binary traits at once
(multi-trait mode, `design.y_cols`):

1. a **full-GRM** null model (from `STEP1_BFILE`, on the GPU),
2. a **sparse-GRM** null model. It needs a sparse GRM file built beforehand.
   If you already have one (e.g. from FastSparseGRM: a MatrixMarket `.mtx` and
   its sample-ID list, format in [Step 1, Sparse GRM](step1.md#sparse-grm)), set
   `SPARSE_GRM` and `SPARSE_GRM_IDS` and it is used exactly like one we build.
   Otherwise the script builds it first, on all `.fam` samples, with
   [step 1's `make_sparse_grm_only`](step1.md#sparse-grm) (relatedness cutoff
   `RELATEDNESS_CUTOFF`, default 0.05).

The step-2 runs with P = 1, 8, 32 reuse **subsets of the same models** (the first P
traits). **Fast test on/off needs no separate null model**: in this code step 2
reads the fast-test switch from the model (`isFastTest` in `nullmodel.json`) and
step 1 only stores it, so the script makes a second view of the sparse models
(`step1/sparse_nofast/`: links to the same files, `nullmodel.json` with
`isFastTest: false`). We checked on simulated data that step 2 on this view gives
byte-identical results to step 2 on a separate step-1 fit with `fast_test: false`.

Models are fitted without LOCO (`loco: false`); step 2 tests one chromosome with
the whole-genome model.

```bash
export STEP1_BFILE=/path/all_chrs PHENO=/path/pheno.txt IID_COL=IID
export BIN_TRAITS=@/path/binary_traits.txt          # or "t1 t2 t3 ..."
export COVARS="age sex PC1 PC2 PC3 PC4 PC5" QCOVARS="sex"
# export SPARSE_GRM=/path/grm.mtx SPARSE_GRM_IDS=/path/grm.ids   # if you have one
bash docs/collab/step1_models.sh
```

Writes:

```
$OUT_ROOT/step1/
├── grm/                      sparse GRM built here (sparse_grm.mtx, .ids, step1.log), unless given
├── sparse_grm.paths          which sparse GRM files were used
├── full/                     step1.yaml step1.log models/<trait>/ vr_<trait>.varianceRatio.txt
├── sparse/                   the same for the sparse-GRM fit (fast_test: true)
└── sparse_nofast/            the sparse models with isFastTest: false (links, no refit)
```

Each fit's log ends with `/usr/bin/time -v` (wall time, peak memory). A fit is
skipped when its `step1.done` exists, so the script can be re-run. The sparse-GRM
fit is single-threaded ([Step 1](step1.md#important-inputs)) and can take hours
at biobank N. The full-GRM fit on the GPU needs about M x N / 4 bytes of GPU
memory for the GRM markers ([GPU, Step 1](gpu.md#step-1)).

## 6. Run the test matrix

```bash
export R_STEP2=/path/to/SAIGE/extdata/step2_SPAtests.R    # block vsr
export RSCRIPT=Rscript                                     # see below
export PGEN=/path/chr21                                    # optional
bash docs/collab/run_matrix.sh 2>&1 | tee -a $OUT_ROOT/matrix.log
```

Run it inside `tmux`/`screen` or as a batch job. **It resumes**: a finished run
has a `DONE` file and is skipped, so after an interruption (preemption, time
limit) start the same command again.

| Variable | Default | Meaning |
|---|---|---|
| `OUT_ROOT`, `GENO`, `BIN_TRAITS` | (required) | as above; `BIN_TRAITS` in the same order as for step 1 |
| `BLOCKS` | `main rare stage pgen vsr precision` | which blocks to run |
| `P_LEVELS` | `1 8 32 128` | P levels; see [What is tested](#what-is-tested) for the top level. `P_TOP=cap` keeps 128 as the top even with more traits |
| `NTHREADS` | all cores | `nThreads` of every C++ step-2 run (R runs use 1 thread per process) |
| `PGEN` | unset | block `pgen` is skipped without it |
| `RARE_GENO` | made from `GENO` | PLINK prefix of the rare set; by default `plink2 --bfile $GENO --mac 1 --max-mac 19` into `$OUT_ROOT/data/rare` |
| `R_STEP2`, `RSCRIPT`, `R_PAR` | -, `Rscript`, 8 | R SAIGE's step-2 script, the Rscript to use, R processes at once |
| `CACHE_MODE` | `auto` | cold page cache before every run, [below](#cold-page-cache) |
| `ONLY`, `SKIP` | unset | extended regex on run names (e.g. `SKIP='main_cpu_P128_.*_firth1'`) |
| `KEEP_OUTPUTS` | 0 | 1 keeps every result file (needs much more disk) |
| `MIN_FREE_GB` | 20 | stop when less is free on `OUT_ROOT` |
| `GPU_ID` | 0 | the GPU nvidia-smi samples (the runs use the first visible GPU) |
| `DRY_RUN` | 0 | 1 prints the list of runs and exits |

Check the plan first:

```bash
DRY_RUN=1 bash docs/collab/run_matrix.sh | head -3
# traits: 60   P levels: 1 8 32 60   cells selected: 77   threads: 64   commit: 39bc03e3
```

**R SAIGE.** `RSCRIPT` must start an Rscript that can `library(SAIGE)` (1.5.2).
If R SAIGE is in its own conda environment, a wrapper works:

```bash
cat > $OUT_ROOT/rscript.sh <<'EOF'
#!/bin/bash
source $HOME/miniforge3/etc/profile.d/conda.sh && conda activate r-saige-1.5.2
exec Rscript "$@"
EOF
chmod +x $OUT_ROOT/rscript.sh; export RSCRIPT=$OUT_ROOT/rscript.sh
```

R reads our null models after a one-time conversion to `.rda`
([`collab/arma_to_rda.R`](collab/arma_to_rda.R), into `$OUT_ROOT/rmodels/`, not
timed), so both sides test against the **same** null model and differences come
from step 2 only. R runs one process per trait with `--nThreads=1`; P = 8 runs 8
processes at once. R's own `--nThreads` is not used (in 1.5.2 it only acts together
with `--idstoIncludeFile`, splits that marker list over forked processes, and the
concatenated output can lose rows and run two rows together). The R options match ours: `minMAC`
1, `maxMissing` 0.15, `SPAcutoff` 2, `max_MAC_for_ER` 4, `pCutoffforFirth` 0.01,
`is_fastTest=FALSE`, `relatednessCutoff` 0, `LOCO=FALSE`, alt-first; sparse-GRM
runs also get `--sparseGRMFile/--sparseGRMSampleIDFile` (the files step 1 used).

### What each run does

Per run, one directory `$OUT_ROOT/cells/<run>/`:

1. writes `cfg.yaml` (step-2 config; GPU runs: `useGPU: true` with all defaults;
   the `stage` GPU run adds `gpuOverlap: false`, the `precision` runs their
   `gpuPrecision*` keys) or `rjobs.sh` (R);
2. evicts the run's input files from the page cache (`cache.json`);
3. runs `/usr/bin/time -v <binary> cfg.yaml` with `SAIGE_STEP2_ROUTE_DUMP` set
   (per-pair route records, used for the counts), and for GPU runs samples
   `nvidia-smi` memory every 200 ms (`gpu_mem.txt`);
4. summarises into `cell.json` and `md5.txt`, deletes the route records (block
   `precision`: after its comparison);
5. compares: a GPU run against its CPU twin, an R run against our CPU run, the
   PGEN run against the `.bed` GPU run (into `$OUT_ROOT/compare/`), a
   `precision` run against the all-fp64 `precision` run (into
   `$OUT_ROOT/compare_precision/`, [below](#block-precision)); result files
   are then deleted unless `KEEP_OUTPUTS=1` or they differ (kept for us to look at;
   `precision` results always differ a little and are deleted once compared).

At the end `summary.csv` and `compare.csv` are written to `$OUT_ROOT`.

### Cold page cache

**Every timed run starts with a cold page cache.** Before each run,
[`cache_evict.py`](collab/cache_evict.py) evicts exactly the files the run
reads: the genotype files (`.bed/.bim/.fam` or `.pgen/.pvar/.psam`), every file
in the null-model directories, the variance-ratio files and, for R, the `.rda`
models and the sparse GRM. Method (`CACHE_MODE=auto` takes the first that works):

| Method | When | What |
|---|---|---|
| `drop_caches` | passwordless `sudo` works | `sync; echo 3 \| sudo tee /proc/sys/vm/drop_caches` (whole machine) |
| `vmtouch` | `vmtouch` is installed | `vmtouch -e <files>` |
| `fadvise` | always (no root needed) | per file: `fsync`, then `posix_fadvise(POSIX_FADV_DONTNEED)` |

Then it measures how much of the files is still resident (`mincore`, the number
`fincore` / `vmtouch` print) and records in `cache.json` and in
`summary.csv`: `cache_method`, `cache_resident_after_MB`, `cache_cold`
(true when < 1% is still resident), the input size and the **filesystem type**.
`summary.csv` also has `fs_inputs_MB` from `/usr/bin/time` (data read from disk
during the run).

On network filesystems (Lustre, GPFS, NFS) eviction on your node may not reach
the file server's cache, so a "cold" run can still read from server memory.
The filesystem type is recorded automatically; please also say in `manual.txt`
what the storage is.

### Block precision

Step 2 on the GPU can run each of its four stages in a lower precision
([GPU, Precision](gpu.md#precision)): the marker scan (`gpuPrecisionScan`
fp64 / fp32 / int8), SPA, the exact test (ER) and Firth (`gpuPrecisionSPA`,
`gpuPrecisionER`, `gpuPrecisionFirth`, fp64 / fp32). The default is fp64
everywhere, which is what every other block runs. Whether the other modes are
faster depends on the card: on V100 / A100 / H100 fp64 runs at half the fp32
rate, on most other cards (L4, A10, T4, consumer cards) at 1/32 to 1/64, and
int8 tensor cores exist from Turing (T4) on. On V100 the fp32 / int8 modes are
not faster (SPA fp32 is slower than fp64); they are meant for weak-fp64 cards.
So we need the timings from your GPU.

The block runs the same configuration (full GRM, Firth on, P = top level,
`.bed` genotypes) eight times:

| Run | Keys on top of the defaults |
|---|---|
| `prec_fp64_...` | none (all fp64; the reference) |
| `prec_scan_fp32_...`, `prec_scan_int8_...` | `gpuPrecisionScan: fp32` / `int8` |
| `prec_spa_fp32_...`, `prec_er_fp32_...`, `prec_firth_fp32_...` | that stage in fp32 |
| `prec_all_fp32_...` | all four fp32 |
| `prec_int8_fp32_...` | scan int8, SPA / ER / Firth fp32 |

Each run is timed like every other run (cold page cache, `/usr/bin/time -v`),
and each non-fp64 run is compared with `prec_fp64_...` by
`step2_saige-step2/tools/precision_compare.py --no-ids` (rows matched on CHR,
POS, MarkerID, Allele1, Allele2; route records compared pair by pair). The
comparison writes aggregate numbers only (no marker IDs) to
`$OUT_ROOT/compare_precision/prec_vs_fp64__<mode>.{json,txt}`, and its main
numbers go into the `vs_fp64_*` columns of `summary.csv`. These results are
**not** byte-identical to the fp64 run (or to R); what we look at is the size of
the differences and whether any p-value crosses 5e-8 or 1e-5.

### Comparing two runs yourself

```bash
python3 docs/collab/compare_outputs.py dirA/out dirB/out --label test --json test.json --rows
```

prints md5 agreement and, with `--rows` or when files differ, a row-level
comparison (rows matched on CHR, POS, MarkerID, Allele1, Allele2): rows only in
one file, rows whose text differs, for p.value / BETA / SE the number of
differing rows and max absolute and relative difference, max |Δ log10 p|, how
many p-values cross 5e-8 (and 1e-5) in one file but not the other, Is.SPA
agreement. The JSON holds aggregates only.

## 7. Collect and send the results

```bash
bash docs/collab/collect_results.sh
# wrote $OUT_ROOT/saige_collab_results_<date>.tar.gz (...); contents: ....tar.gz.contents.txt
```

It copies an allow-list of files into a staging directory, scans every file for
your sample IDs (from the `.fam`, `.psam` and the phenotype file) and marker IDs
(`.bim`, `.pvar`), removes any line containing one (counts in
`REDACTION_REPORT.txt`), refuses to continue if it finds a file over 50 MB or a
result-table header, and writes the tar plus a listing of its contents.

| In the tar | Not in the tar |
|---|---|
| `env/` environment record and `manual.txt` | result files (`out/*.txt`) |
| `build/BUILD_INFO.txt` | route records |
| `selfcheck/selfcheck_result.txt` | null models, `.rda`, sparse GRM |
| `step1/*/step1.yaml`, `step1.log` | genotype, phenotype files |
| `cells/<run>/` cfg, log, `cell.json`, `cache.json`, `md5.txt`, `gpu_mem.txt`, `stage.txt`; for R runs `rjobs.sh` and each process's `/usr/bin/time` record | sample or marker ID lists |
| | R SAIGE's own logs (for sparse-GRM runs they print a table of sample IDs); only for a failed R run their last 30 lines |
| `compare/*.json` (aggregates), `compare_precision/*.json`, `*.txt` (aggregates, no marker IDs), `summary.csv`, `compare.csv` | |

Please read `<tar>.contents.txt` and `REDACTION_REPORT.txt` before sending.
Sample IDs shorter than 5 characters are not searched for (short numbers occur
everywhere in logs); if your IDs are that short, set e.g. `REDACT_MIN_LEN=2` and
check `REDACTION_REPORT.txt` for lines removed by accident.

## 8. Result files we receive

### `summary.csv`: one row per run

| Column | Meaning |
|---|---|
| `block`, `cell` | block and run name, e.g. `main_gpu_P32_sparse_nofast_firth1` |
| `path` | `cpu`, `gpu` or `R` |
| `trait_type` | `binary` |
| `P` | traits in the run |
| `grm` | `full` or `sparse` |
| `fast_test` | `true`, `false`, or `na` (full GRM) |
| `firth` | 0 / 1 (`is_Firth_beta`, p < 0.01) |
| `geno`, `binary`, `extra` | `bed` / `rare` / `pgen`; `prod` / `phase` / `R`; extra config line |
| `nthreads` | threads per process |
| `rc` | exit code (0 = success) |
| `wall_s`, `user_s`, `sys_s`, `cpu_pct` | `/usr/bin/time -v`: wall clock, CPU seconds, CPU utilisation |
| `peak_rss_kb` | peak resident memory (for R: of the largest process) |
| `r_processes`, `r_peak_rss_sum_kb`, `r_proc_wall_max_s` | R runs: processes, sum of their peak RSS, slowest process |
| `gpu_peak_mem_mib`, `gpu_mem_base_mib` | peak `nvidia-smi` memory.used during the run minus the value before it |
| `gpu_buffers_mib_log` | the step-2 log's "MiB on the device" |
| `gpu_name`, `gpu_coverage` | GPU name, `pairs on the GPU / all pairs` from the log |
| `gpu_refused` | the reason when a GPU run fell back to the CPU (empty otherwise) |
| `startup_s`, `main_loop_s` | from the step-2 `[TIMING]` marks |
| `n_markers_tested_max`, `n_result_files` | markers tested (largest over traits), result files |
| `n_pairs` | (marker, trait) pairs with an output row |
| `n_need_spa`, `n_spa` | pairs over the SPA cutoff; pairs where SPA ran and converged (not ER) |
| `n_lowmac`, `n_er` | pairs with MAC ≤ 4; pairs where the exact test ran |
| `n_need_firth`, `n_firth`, `n_firth_conv` | pairs with p < 0.01; Firth fits; converged fits |
| `n_need_fast` | sparse fast-test recomputations |
| `n_firth_log` | Firth fits according to the log |
| `n_spa_is_spa_col`, `counts_source` | when there is no route record (R; our CPU path at P = 1, which takes the single-trait code path): pairs and Is.SPA = true counted from the result files |
| `cache_method`, `cache_cold`, `cache_resident_after_MB`, `cache_input_MB`, `fstypes`, `fs_inputs_MB` | [cold page cache](#cold-page-cache) record |
| `log_errors` | lines matching error / terminate / segfault in the log |
| `cpu_vs_gpu_identical` | `main`/`rare`: CPU and GPU result files byte-identical |
| `compare_label` | the row of `compare.csv` for this run |
| `prec_scan`, `prec_spa`, `prec_er`, `prec_firth` | GPU runs: the precision of each stage, from the log line `GPU precision: scan=... SPA=... ER=... Firth=...` |
| `vs_fp64_p_max_rel`, `vs_fp64_max_abs_dlog10p` | block `precision`: against the all-fp64 run, max relative difference of p.value, max abs(Δ log10 p) |
| `vs_fp64_p_max_rel_p_lt_1e-5`, `vs_fp64_n_p_lt_1e-5` | the same over the rows with p < 1e-5 on either side, and their number |
| `vs_fp64_beta_max_rel`, `vs_fp64_se_max_rel` | max relative difference of BETA, SE |
| `vs_fp64_cross_5e-8`, `vs_fp64_cross_1e-5` | rows significant at 5e-8 (1e-5) in only one of the two runs |
| `vs_fp64_is_spa_differs` | rows whose Is.SPA differs |
| `vs_fp64_firth_route_changes`, `vs_fp64_spa_nonconv_changes`, `vs_fp64_firth_nonconv_changes` | from the route records: pairs Firth-fitted in only one run; pairs where SPA (Firth) ran in both and converged in only one |
| `vs_fp64_rows_one_side`, `vs_fp64_note` | rows present in one run only; why a comparison is missing |
| `md5_all`, `commit`, `start` | md5 over the run's `md5.txt`, source commit, start time |

### `compare.csv`: one row per comparison

`label` is `cpu_vs_gpu__<run>`, `cpp_vs_R__<run>` (A = ours, B = R) or
`pgen_vs_bed__<run>`.

| Column | Meaning |
|---|---|
| `identical`, `n_paired`, `n_md5_identical`, `files_only_a/b` | md5 agreement per trait file |
| `rows_a`, `rows_b`, `only_a`, `only_b` | rows in each, rows present in one but not the other |
| `rows_text_differ` | matched rows whose text differs in any column |
| `p_*`, `beta_*`, `se_*` | for p.value, BETA, SE: rows that differ, max absolute and max relative difference |
| `na_mismatch` | NA on one side only |
| `max_abs_dlog10p` | max abs(log10 pA - log10 pB) |
| `n_p_lt_5e-8_a/b`, `cross_5e-8_a_only/b_only` | p < 5e-8 in A and in B; significant in only one of them (same for 1e-5) |
| `is_spa_true_a/b`, `is_spa_disagree` | Is.SPA counts and disagreements |

Row-level columns are empty when all files were byte-identical (only md5 was
compared); `cpp_vs_R` always has them.

## 9. Run time and what to drop first

Step-2 time grows roughly with N x M x P (N samples, M markers on the chromosome,
P traits). Rough figures from our earlier single-variant binary-trait runs on 8
vCPUs (4 cores) and one V100, simulated data (yours will differ, especially in how
many pairs need SPA and Firth):

| Path | seconds per (sample x marker x trait) |
|---|---|
| our CPU, 8 vCPU, Firth off | 0.7–1.8 x 10⁻⁹ |
| our CPU, 8 vCPU, Firth on | about 3.5 x 10⁻⁹ |
| our GPU (V100), P ≥ 32 | about 0.02–0.05 x 10⁻⁹ (+ about 10 s start-up) |
| R SAIGE 1.5.2, one process, one trait | about 40 x 10⁻⁹ |

Example, N = 400,000 and M = 50,000 (N x M = 2 x 10¹⁰), 8 vCPUs: one CPU run at
P = 128 with Firth on is about 128 x 3.5 x 10⁻⁹ x 2 x 10¹⁰ ≈ 9,000 s ≈ 2.5 h; the 24 CPU
runs of `main` add up to about 10–15 h, of which the six P = 128 runs are about
three quarters. The 24 GPU runs take well under an hour. Block `vsr` (R) is about
2–3 h: P = 1 is ~15 min per R process, P = 8 runs 8 processes side by side. More
cores shorten the CPU runs, but not linearly; please run a pilot first:

```bash
ONLY='main_(cpu|gpu)_P8_full_firth[01]' bash docs/collab/run_matrix.sh
```

and scale its `wall_s` by P / 8 to estimate the larger runs.

If time is limited, drop runs in this order (examples are `SKIP` patterns, which
can be combined with `|`):

1. CPU, P = 128, Firth on: `SKIP='main_cpu_P128_.*_firth1'`
2. CPU, P = 128, Firth off: `main_cpu_P128_.*_firth0`
3. CPU, P = 32, sparse with fast test off (the slowest sparse variant on the CPU): `main_cpu_P32_sparse_nofast`
4. the rest of CPU P = 32: `main_cpu_P32_`
5. R at P = 8 with the sparse GRM: `vsr_.*_P8_sparse`

Keep the GPU runs, the P = 1 and 8 runs, and the `stage` block; when a CPU run is
dropped its GPU twin is still run, it is just not compared.

### Things you may see

- A GPU run with `gpu_refused` set: the run fell back to the CPU and says why. With
  the sparse GRM and fast test off this happens when a connected component of the
  sparse GRM is too large for the block inverse (log: `block inverse refused by
  the cost gate`); the results are then those of the CPU path. That is an
  outcome we want to know about, not a failure.
- Block `vsr`, sparse GRM: small differences from R are expected. Our step 2 by
  default solves the sparse-GRM variance with an exact block inverse
  (`blockSparseSigma: true`; `saige-step2` without arguments describes it), R 1.5.2
  with an iterative solver (PCG) that stops at a tolerance.
  On simulated data: p-values differed in about 37% of rows, max relative
  difference 3e-4, max abs(Δ log10 p) 1.3e-4, no p-value crossed 5e-8. Full-GRM
  runs are expected to be byte-identical to R.
- `cache_cold` false: eviction did not work for some file (e.g. a network
  filesystem with `CACHE_MODE=fadvise`); the run still counts, please mention it.
- P = 1 runs on the CPU take a different code path (single-trait) from P ≥ 2, so
  their SPA / ER / Firth counts come from the result files and are partly empty.

## What was verified

On GCP n1-standard-8 (8 vCPU, 30 GB), Tesla V100 16 GB, CUDA 12.9, commit
`39bc03e3`, simulated data only:

- `build_bins.sh 70`; `selfcheck.sh` (the md5s above, the same with `NTHREADS`
  8 and 3); `env_record.sh`.
- The [rehearsal](#optional-rehearse-the-whole-protocol-on-simulated-data) as
  written: `step1_models.sh`, `run_matrix.sh` with all five blocks (P = 1 and 8,
  45 runs, R SAIGE 1.5.2 for all 8 `vsr` runs, `CACHE_MODE=fadvise`) and
  `collect_results.sh`. Results: all 13 CPU/GPU pairs byte-identical (including
  the rare-marker set: 6,060 pairs, 31 with the exact test); PGEN = `.bed`;
  ours = R byte-for-byte for the full GRM; sparse GRM within the differences
  quoted above. The tar was unpacked and searched for sample IDs with
  `REDACT_MIN_LEN=2`: none.
- The `sparse_nofast` view against a separate step-1 fit with `fast_test: false`:
  model files identical apart from the flag, step-2 results byte-identical.
- With the default relatedness cutoff 0.05 on that cohort (one large connected
  component), the GPU runs with fast test off fell back to the CPU path
  (`gpu_refused`), and our sparse-GRM results were then byte-identical to R.
- `cache_evict.py` with `drop_caches`, `fadvise` and `none`, checked with `fincore`.
- Block `precision` alone (`BLOCKS=precision`, commits `186c15fb` and `63a94267`, the
  rehearsal cohort's step-1 models, P = 8, `CACHE_MODE=fadvise`): 8 runs,
  7 comparisons with all-fp64 in `compare_precision/` and the `vs_fp64_*`
  columns of `summary.csv`, result files and route records deleted after the
  comparisons; `collect_results.sh` packed the comparison summaries and no
  result rows (no marker or sample IDs found).

Not verified here: P = 32 and 128, biobank-scale N and M, other GPUs, network
filesystems, `CACHE_MODE=vmtouch`.
