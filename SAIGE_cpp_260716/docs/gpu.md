---
title: GPU
nav_order: 4
---

# GPU

Both steps can use one NVIDIA GPU. The GPU is optional: a run that cannot use it
runs on the CPU and says so in the log (step 2) or silently (step 1, see below).

## In short

| | CPU | GPU |
|---|---|---|
| Build | `make` | `make USE_CUDA=1 SM=<compute capability>` |
| Step 1 | no flag (config `fit.use_gpu: false`) | `saige-gpu-cpp step1 --useGPU` (config `fit.use_gpu: true`) |
| Step 2 | no flag (config `useGPU: false`) | `saige-gpu-cpp step2 --useGPU` (config `useGPU: true`) |

Nothing else changes, and the step-2 output is the same (with the default
fp64 [precision](#precision)).
Recommended cards: V100, A30, A100, H100 (full-speed fp64); see the table in
[Home, Build](index.md#3-build).

## Requirements

- The GPU build ([Home, Build](index.md#3-build)): `make USE_CUDA=1 SM=<SM>` with the
  compute capability of your card (70 V100, 75 T4, 80 A100, 86 A10, 89 L4, 90 H100),
  or several at once (`SM="70 80 90"`).
- CUDA 12.x runtime and an NVIDIA driver that supports it. Tested: CUDA 12.9,
  driver 580, Tesla V100-SXM2 16 GB.
- One GPU per process. Both steps use device 0 unless `--gpuDevice N` is given
  (step 2 config `gpuDevice: N`; step 1 environment `SAIGE_GPU_DEVICE=N`);
  `CUDA_VISIBLE_DEVICES` also works.

## Step 1

Switch: `saige-gpu-cpp step1 --useGPU` (config `fit.use_gpu: true`, or `saige-null --gpu`).

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

Switch: `saige-gpu-cpp step2 --useGPU` (config `useGPU: true`). With it, every
sub-switch below defaults to on; turn one off with `--set <key>=false` (config
`<key>: false`). Results are the same as on the CPU (the tutorial's text
files are byte-identical with `useGPU: true` and `false`) as long as every
stage runs in fp64, the default ([Precision](#precision)).

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

Other GPU keys: `gpuDevice` (0), `gpuBlockSize` (16384, but a batch never exceeds `marker_chunksize`,
default 10,000, so with defaults the batch is 9,984 markers), `gpuOverlapSets` (6), `gpuOverlapLag` (3),
`gpuPrefetchSets` (3), `gpuFirthMaxStep` (15), `gpuSparseMaxPairs` (5e7).

A dependency that is not met turns the dependent switch off; when the config
wrote that switch explicitly, the log says so (e.g. `gpuFirth: ignored, it needs gpuSpa: true`).

### Precision

Step 2 on the GPU has four stages whose arithmetic precision can be chosen
separately. The default is fp64 everywhere: that is the output described above,
byte-identical to the CPU path and to R where those agree. Any combination of the
modes below may be given; a mode is never replaced by another one.

| Config key | Flag (`step2`) | Values (default first) | Stage |
|---|---|---|---|
| `gpuPrecisionScan` (old name `gpuPrecision`) | `--gpuPrecisionScan` | `fp64`, `fp32`, `int8` | genotype decode + GEMMs of the marker scan, and the sparse-GRM variance |
| `gpuPrecisionSPA` | `--gpuPrecisionSPA` | `fp64`, `fp32` | saddlepoint approximation (`gpuSpaImpl` `lib` and `own`) |
| `gpuPrecisionER` | `--gpuPrecisionER` | `fp64`, `fp32` | exact test for MAC ≤ `MACCutoffforER` |
| `gpuPrecisionFirth` | `--gpuPrecisionFirth` | `fp64`, `fp32` | Firth fits |
| `gpuInt8Slices` | `--gpuInt8Slices` | `7` (1..8) | `int8` scan: int8 slices of the trait-side matrix |

The log prints one line, e.g.
`GPU precision: scan=int8 SPA=fp32 ER=fp32 Firth=fp32 (int8 slices 7)`. The CPU
path ignores these keys (it says so once in the log). The ER start-up self-check
runs only with ER `fp64`.

What the modes do (every stage hands its results back in double):

- **Scan `fp32`**: SGEMMs over 4,096-sample chunks (16,384 for the second
  GEMM), each chunk's sum added in fp64; genotypes are shifted by the rounded
  marker mean and trait-side columns with a large mean are centred, and both
  are undone exactly in fp64. Sparse-GRM cross terms in fp32 with fp64 partial
  sums.
- **Scan `int8`**: the trait-side matrix is split into `gpuInt8Slices` int8
  digits of 7 bits each (Ozaki scheme); genotypes are exact in int8 (a
  mean-imputed missing call gets its own 0/1 column applied in fp64); one int8
  GEMM with int32 accumulation per slice, recombined in fp64. With 7 slices
  the results matched fp64 to the printed digits on every test set. The
  sparse-GRM cross terms run in fp64. Needs N < 8,388,608 samples.
- **SPA `fp32`**: the per-sample sums of the cumulant generating function and
  its derivatives, the projection passes (fp32 copies of the trait-side
  matrices, float-float for the carriers) and the block reductions in fp32
  (centred, compensated, overflow-safe forms); the Newton scalars and the whole
  tail probability in fp64, so p-values do not underflow earlier than in fp64.
  Needs 8N + 16 x traitStride bytes more device memory per trait (traitStride
  = the padded trait-side row length), on top of the fp64 tables.
- **ER `fp32`**: the 2^k enumeration in the log domain in fp32, with the final
  p-value formed in double.
- **Firth `fp32`**: the per-sample work of each Newton step in fp32 with
  compensated sums; block reductions and the 2 x 2 solves in fp64; the
  convergence rule is unchanged.

Accuracy against all-fp64 on the same build (V100; six simulated sets with
binary traits, 8 traits each, Firth on: bingpu_test bt and bm with the full and
the sparse GRM, a rare-variant set with `MACCutoffforER: 20`, and 200,000
markers x 50,000 samples; 4.4 million (marker, trait) pairs, 638 with p < 1e-5).
Largest value over all sets:

| Mode (rest fp64) | p.value max rel. diff | same, p < 1e-5 | BETA max abs diff | SE max rel. diff | p crossing 5e-8 / 1e-5 | Is.SPA, Firth route, convergence changes |
|---|---|---|---|---|---|---|
| scan `fp32` | 1.7e-3 | 1.7e-3 | 1e-5 | 3.7e-5 | 0 / 0 | 0 |
| scan `int8` | 0 | 0 | 0 | 0 | 0 / 0 | 0 |
| SPA `fp32` | 2.7e-4 | 1.2e-4 | 0 | 3.5e-5 | 0 / 0 | 0 |
| ER `fp32` | 3.9e-6 | 0 | 0 | 9.8e-6 | 0 / 0 | 0 |
| Firth `fp32` | 0 | 0 | 8e-5 | 1.1e-4 | 0 / 0 | 0 |
| all `fp32` | 1.6e-3 | 1.6e-3 | 8e-5 | 1.1e-4 | 0 / 0 | 0 |
| scan `int8`, rest `fp32` | 2.7e-4 | 1.2e-4 | 8e-5 | 1.1e-4 | 0 / 0 | 0 |

Differences are measured on the printed text output (6-7 significant digits).
The scan `fp32` maximum sits on SPA-adjusted markers whose score is at the edge
of its support (one bm marker, p = 6e-11), where the p-value moves 1e-3 for a
1e-7 relative change of the score; on the other five sets scan `fp32` stays at
or below 1.5e-4. Results in any non-fp64 mode are **not**
byte-identical to the fp64 GPU run, the CPU path or R SAIGE.

Speed: on V100, A100 and H100 fp64 runs at half the fp32 rate, so the low
precision modes do not pay off there: on V100 SPA `fp32` is 30-55% slower than
SPA `fp64` (it was built to move fp64 work off weak-fp64 cards). **The fp32 and
int8 modes are meant for weak-fp64 cards only.**

V100 has no int8 tensor cores, so scan `int8` is not faster than fp64 there. Cards with weak fp64 (L4,
A10, T4, RTX: 1/32 to 1/64 of fp32) are where they are meant to pay off, and
int8 tensor cores exist from Turing (T4) on; measure there before choosing
(the `precision` block of the [collaborator tests](collaborator_tests.md#block-precision)
does that).

When to use which:

- Datacenter cards with full-rate fp64 (V100, A100, H100, A30): keep the
  default fp64.
- Cards with weak fp64: scan `int8` (results equal to fp64 in our tests) or
  scan `fp32`, with SPA, ER and Firth `fp32`:
  `--gpuPrecisionScan int8 --gpuPrecisionSPA fp32 --gpuPrecisionER fp32 --gpuPrecisionFirth fp32`.
- One stage only, e.g. Firth in fp32 with everything else fp64:
  `--gpuPrecisionFirth fp32`.

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
| `built without CUDA (rebuild with: make USE_CUDA=1)` | CPU build of `saige-step2`; rebuild with `make USE_CUDA=1 SM=<SM>` |
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
| `gpuSparse: ...` | the sparse-GRM variance could not be set up on the GPU |
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
