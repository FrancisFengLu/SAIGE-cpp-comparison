---
title: Troubleshooting
nav_order: 6
---

# Troubleshooting

Messages were reproduced on the tutorial data except those marked †, which are
quoted from the source code. Step 1 errors end the
run with `terminate called after throwing ... what(): <message>`; step 2 errors
print `ERROR: <message>` and exit with status 1. Check the input before a long
run with `saige-null -c step1.yaml --dry-run`.

## Build

| Symptom | Cause / fix |
|---|---|
| `activate-gcc_linux-64.sh: line 114: SYS_SYSROOT: unbound variable` | `conda activate` inside a script with `set -u`; activate before `set -u` (as `examples/env.sh` does) |
| `USE_CUDA=1 but nvcc is not on PATH` | step 2 GPU build: add `$CUDA_HOME/bin` to `PATH` or pass `NVCC=/path/to/nvcc` |
| step 1 built, but never prints `GPU tier=` | step 1's Makefile found no `nvcc` and built the CPU version; pass `NVCC=/path/to/nvcc` |
| `undefined reference to pcre2_*@PCRE2_10.47` (step 1 link) † | the environment's `pcre2` is older than R's; install `pcre2>=10.47` |
| segfaults after changing branches or build type † | stale object files; `make clean` and rebuild |
| `nvcc warning : Support for offline compilation for architectures prior to ... _75 ...` | harmless for `SM=70` |

## Step 1

| Message | Cause / fix |
|---|---|
| `Fatal error: R home directory is not defined` | set `R_HOME` (e.g. `export R_HOME=$CONDA_PREFIX/lib/R`) |
| `IID in design not found in FAM: extra1` | the phenotype file has a sample that is not in the `.fam`. Add `design.whitelist_ids:` with the `.fam` IIDs (`cut -f2 geno.fam > ids.txt`) or remove the row |
| `ERROR: binary phenotype value must be 0 or 1, found: 2.000000 at sample per3` | recode cases/controls as 1/0 |
| `ERROR: variance of the phenotype (0.000125) is much smaller than 1. Please consider setting inv_normalize: true in config.` | quantitative trait on a small scale; set `fit.inv_normalize: true` or rescale |
| `Covariate column not found: zz` | a name in `covar_cols` is not in the header (names are matched case-insensitively) |
| `Design file: phenotype column 'nosuch' not found. ...` | set `design.y_col` / `y_cols` |
| `Refusing to overwrite existing variance-ratio file: ... (set paths.overwrite_varratio=true to allow).` | output exists; change `out_prefix_vr` or set `paths.overwrite_varratio: true` |
| `[design] dropped 10 row(s) with missing phenotype or covariates (complete.cases)` | information: rows with `NA`/empty trait or covariate are left out |
| `Converged: NO` † | the fit hit `maxiter`; do not use the model. Check the trait (case count, scale) and covariates |
| no `[parallelCrossProd] GPU tier=` line with `use_gpu: true` | step 1 ran on the CPU: CPU build, no visible GPU, or a sparse-GRM fit ([GPU](gpu.md#step-1)) |
| LOCO not applied (`LOCO: off`) † | the `.bim` has fewer than 2 autosomes, or the fit uses a sparse GRM |

Overrides with `-o` take scalars only (`-o fit.loco=false`,
`-o paths.out_prefix=dir`); edit the YAML for lists such as `covar_cols`.

## Step 2

| Message | Cause / fix |
|---|---|
| `ERROR: Cannot open output file: .../out/b1.txt` | the output directory does not exist; `mkdir -p` it first |
| `ERROR: genotype is in bgen, please set AlleleOrder=ref-first ...` | BGEN and PGEN need `AlleleOrder: ref-first` (or leave the key out) |
| `ERROR: The models were fitted on different sample sets; multi-trait testing on different sample sets currently supports genoType: plink or a hard-call pgen only (this config uses bgen). Run one config per sample set.` | BGEN / VCF / dosage PGEN with models on different samples: put the models with the same samples in one config each, or use PLINK / hard-call PGEN |
| `ERROR: model 'b4' (models[3]) lists 4511 sample IDs but 'b1' lists 5000. mtRequireSameSamples: true requires identical sample IDs in identical order.` | `mtRequireSameSamples: true` is set; remove it to allow different sample sets |
| `ERROR: outputFormat: sgs is implemented on the multi-trait single-variant paths, and this run takes the single-trait path (one model and useGPU is off). ...` | one model with `outputFormat: sgs`: add `useGPU: true` or use `outputFormat: text` |
| `ERROR: chrom needs to be specified in order to apply leave-one-chromosome-out. ...` | `LOCO: true` without `chrom:` |
| `ERROR: No markers on chrom 3 are found` | `chrom` does not match any chromosome code in the genotype file |
| `LOCO=true with genoType=bgen requires a .bgi index ...` † | create it: `bgenix -g file.bgen -index` |
| `useGPU: refused, running on the CPU (<reason>)` | not an error; see the reasons in [GPU](gpu.md#reading-the-log) |
| `gpuSparse: not active (...)` | not an error; sparse-GRM variances stay on the CPU |
| VCF results not in genomic order | VCF input with `nThreads` > 1; sort, or use `nThreads: 1` |
| BGEN results have Allele1/Allele2 and the BETA sign swapped vs. PLINK | expected: BGEN is read ref-first; `Allele2` is always the tested allele |
| `isFirth: false` but Firth still applied | Firth is controlled by `is_Firth_beta` (model default from step 1, or set in the step-2 config) |

## GPU out of memory

Step 2 reports a device failure at start-up (`device setup failed` †) and runs on
the CPU. Reduce `gpuOverlapSets` (minimum 3), `gpuBlockSize`, or the number of
traits per run; see [GPU memory](gpu.md#gpu-memory-step-2).
