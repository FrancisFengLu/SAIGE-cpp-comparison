---
title: GPU
nav_order: 4
---

# GPU

Both steps can use one NVIDIA GPU. The GPU is optional: a run that cannot use it
runs on the CPU and says so in the log (step 2) or silently (step 1, see below).

## Requirements

- The GPU build ([Home, Build](index.md#3-build)): `01_build.sh gpu <SM>` with the
  compute capability of your card (70 V100, 75 T4, 80 A100, 86 A10, 89 L4, 90 H100).
- CUDA 12.x runtime and an NVIDIA driver that supports it. Tested: CUDA 12.9,
  driver 580, Tesla V100-SXM2 16 GB.
- One GPU per process. Step 2 uses device 0 unless `gpuDevice: N` is set;
  `CUDA_VISIBLE_DEVICES` also works.

## Step 1

Switch: `fit.use_gpu: true` in the config, or `--gpu` on the command line.

| Fit | Uses the GPU |
|---|---|
| full GRM from the PLINK genotypes (`use_sparse_grm_to_fit: false`), with or without LOCO, one or several traits | yes |
| sparse GRM (`use_sparse_grm_to_fit: true`), `make_sparse_grm_only` | no (these do not use the full GRM) |

Check the log: when the GPU is used it prints

```
[parallelCrossProd] GPU tier=4 enabled (source=packed_flat_, batch K·U yes).
```

once per sample-set group. **If this line is missing, step 1 ran on the CPU** —
with a CPU build or no visible device it falls back without a message.

The GPU holds the GRM markers packed at 2 bits per genotype: about
M × N / 4 bytes (M = markers passing `min_maf_grm`, N = samples), e.g.
100,000 markers x 400,000 samples ≈ 10 GB. When that does not fit, step 1
streams the genotypes from host memory instead.

GPU and CPU fits agree to rounding; the GPU sums in a different order, so the
last digits of the estimates (and of the step-2 results computed from them)
can differ between a GPU and a CPU step 1.

## Step 2

Switch: `useGPU: true`. With it, every sub-switch below defaults to on; write
`false` to turn one off. Results are the same as on the CPU (the tutorial's text
files are byte-identical with `useGPU: true` and `false`).

| Key | Default with `useGPU: true` | What it moves to the GPU |
|---|---|---|
| `gpuBinary` | `true` | binary traits (without it only quantitative traits run on the GPU) |
| `gpuSpa` | `true` | saddle-point approximation (needs `gpuBinary`) |
| `gpuFirth` | `true` | Firth fits (needs `gpuSpa`) |
| `gpuER` | `true` | exact test for MAC ≤ `MACCutoffforER` (needs `gpuBinary`); a start-up self-check keeps it on the CPU if the GPU cannot reproduce the CPU's arithmetic exactly |
| `gpuSparse` | `true` | variance of binary traits fitted on a sparse GRM (needs `gpuBinary`) |
| `gpuOwnSampleSets` | `true` | models with different sample lists |
| `gpuPgen` | `true` | hard-call `.pgen` input |
| `gpuPrefetch` | `true` | reading the next block while the current one computes (`gpuPrefetchThreads`, default 4) |
| `gpuOverlap` | `true` | GPU and CPU work on different blocks at once |
| `gpuDecodeX2`, `gpuSpaFused`, `gpuSpaDynamic`, `gpuSpaOrder: trait`, `gpuSpaMinBlocks: 3`, `gpuSpaImpl: lib` | as shown | kernel variants; no effect on results |

Other GPU keys: `gpuDevice` (0), `gpuPrecision` (`fp64`; `fp32` disables
`gpuSparse`), `gpuBlockSize` (16384, but a batch never exceeds `marker_chunksize`,
default 10,000, so with defaults the batch is 9,984 markers), `gpuOverlapSets` (6), `gpuOverlapLag` (3),
`gpuPrefetchSets` (3), `gpuFirthMaxStep` (15), `gpuSparseMaxPairs` (5e7).

A dependency that is not met turns the dependent switch off; when the config
wrote that switch explicitly, the log says so (e.g. `gpuFirth: ignored, it needs gpuSpa: true`).

### What runs on the GPU

| Input / setting | GPU | Otherwise |
|---|---|---|
| `genoType: plink` | yes | |
| `genoType: pgen`, hard calls only | yes | |
| `genoType: pgen` with dosages | no | CPU, log: `... the file is not hard-call only: it stores dosages ...` |
| `genoType: bgen`, `vcf` | no | CPU, log: `genoType is not plink or pgen` |
| binary and quantitative traits, mixed | yes | |
| models with different sample lists | yes (PLINK, hard-call PGEN) | |
| binary trait fitted on a sparse GRM | yes | |
| quantitative trait fitted on a sparse GRM with `fast_test: false` | no | CPU, log: `... sparseGRM first pass (gpuSparse covers binary traits only)` |
| `LOCO: true` (one chromosome per run) | yes | |
| one model | yes (it is run like a one-entry `models:` list) | |
| `outputFormat: text` or `sgs` | yes | |
| conditional analysis (`condition`) | no | CPU |
| `isnoadjCov: true` | no | CPU |
| region / gene-based tests (`groupFile`) | no | CPU |

### Reading the log

Used (per trait, then totals):

```
  useGPU: Tesla V100-SXM2-16GB, 16144 MiB, sm_70; fp64, ... 9984 markers per device batch (39 x 256), 391 MiB on the device, decode x2
[b1] 5000 markers were tested (4769 on the GPU + 231 with the device SPA, 0 via the scalar CPU path; 49 Firth fits on the device).
  GPU coverage: 19999 / 20000 pairs (99.995%)
```

A few markers of a GPU run can go "via the scalar CPU path"; that is normal.

Not used — one line with the reason, then the run continues on the CPU:

```
  useGPU: refused, running on the CPU (<reason>)
```

| Reason | Meaning / fix |
|---|---|
| `built without CUDA (rebuild with: make USE_CUDA=1)` | CPU build of `saige-step2`; rebuild with `01_build.sh gpu <SM>` |
| `cudaGetDeviceCount: no CUDA-capable device is detected` | no GPU visible (check `nvidia-smi`, `CUDA_VISIBLE_DEVICES`, the driver) |
| `gpuDevice 3 out of range (1 device(s) present)` | `gpuDevice` names a GPU that does not exist |
| `genoType is not plink or pgen` | BGEN / VCF input |
| `genoType is pgen but the file is not hard-call only: ...` | PGEN with dosages or multiallelic variants |
| `genoType is pgen and gpuPgen is false` | remove `gpuPgen: false` |
| `binary traits present (gpuBinary: true runs them on the device)` | remove `gpuBinary: false` |
| `the models do not share one sample list (gpuOwnSampleSets: true ...)` | remove `gpuOwnSampleSets: false` |
| `mtBatch is false` | remove `mtBatch: false` |
| `trait '<name>': sparseGRM first pass ...` | quantitative sparse-GRM model with `fast_test: false` |
| `trait '<name>': isnoadjCov=true`, `... runs conditional analysis` | not supported on the GPU |
| `gpuSparse: ...` | the sparse-GRM variance could not be set up on the GPU (e.g. `gpuPrecision: fp32`) |
| `device setup failed` | GPU out of memory or CUDA error at start; see memory below |

Also not a refusal, but worth knowing:
`gpuSparse: not active (... refused by the cost gate ...)` means the sparse GRM
has components too large for the block inverse. With `isFastTest: true` the
affected variances are computed on the CPU and the rest of the run stays on the
GPU; with `isFastTest: false` the whole run falls back to the CPU.

### GPU memory (step 2)

Rough planning estimate for the buffers the log reports, N = samples in the
genotype file, P = traits:

```
trait-side constants   N × 7P × 8 bytes
genotype ring          6 × 9,984 × N / 4 bytes     (gpuOverlapSets × batch × packed row)
decode buffers         up to about 1 GB
```

The log line `... MiB on the device` gives the actual figure for these buffers.
The process uses more than that in total (CUDA context, cuBLAS, SPA / Firth / ER
work space). Measured on this guide's V100:

| Run | log: MiB on the device | peak in `nvidia-smi` |
|---|---|---|
| N = 5,000, P = 4 binary (tutorial) | 391 | — |
| N = 50,000, P = 10 binary, Firth on | 1,778 | 4,333 |

So leave a few GB of headroom above the estimate. To reduce memory:
`gpuOverlapSets: 3`, a smaller `gpuBlockSize`, or fewer traits per run. Host
memory is not reduced by the GPU.
