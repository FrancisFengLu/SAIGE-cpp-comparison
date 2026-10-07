---
title: Home
nav_order: 1
---

# SAIGE C++ / GPU

A standalone C++ implementation of SAIGE for genome-wide association tests of
binary and quantitative traits with a generalized linear mixed model. It runs
on the CPU, and optionally on an NVIDIA GPU. Many traits can be analysed in one
run.

Pages:

| Page | Contents |
|---|---|
| [Home](index.md) | workflow, installation, tutorial |
| [Step 1: fit the null model](step1.md) | `saige-null` inputs, config keys, outputs |
| [Step 2: test genetic variants](step2.md) | `saige-step2` inputs, config keys, outputs |
| [GPU](gpu.md) | GPU build, switches, what runs on the GPU, memory |
| [HPC example](hpc_example.md) | one job per chromosome, concatenating results |
| [Troubleshooting](troubleshooting.md) | error and log messages and what they mean |

Every command and config on these pages was run on simulated data
(`plink2 --dummy`, 5,000 samples x 5,000 markers) with the scripts in
[`examples/`](examples/).

## Analysis workflow

1. **Step 1** (`saige-null`, YAML config): fit one null model per trait from the
   phenotype/covariate table and the PLINK genotypes (full GRM) or a sparse GRM.
   Output: one model directory and one variance-ratio file per trait.
2. **Step 2** (`saige-step2`, YAML config): test every marker of a genotype file
   (PLINK, PGEN, BGEN or VCF) against the step-1 models. Output: one result
   table per trait.

## Installation

No container image is published. Build from source.

### Requirements

- Linux x86-64. The build uses `-march=native`, so build on the
  machine type you will run on.
- A conda environment with the libraries below. Step 1 embeds R (for its random
  number generator), so it needs R and the Rcpp packages at build and run time.
- For the GPU build only: CUDA toolkit 12.x (`nvcc`, cuBLAS) and an NVIDIA GPU
  with compute capability 7.0 or newer. Tested: CUDA 12.9, Tesla V100 (sm_70).

### 1. Create the build environment

```bash
mamba create -y -p $HOME/miniforge3/envs/saige-build -c conda-forge -c bioconda \
    gxx_linux-64=12 make pkg-config armadillo openblas superlu yaml-cpp htslib \
    zstd zlib sqlite boost-cpp eigen tbb-devel pcre2 \
    r-base r-rcpp r-rcpparmadillo r-rcppparallel
```

### 2. Get the code and point the example scripts at your paths

```bash
git clone https://github.com/FrancisFengLu/SAIGE-cpp-comparison.git
cd SAIGE-cpp-comparison/SAIGE_cpp_260716
```

[`examples/env.sh`](examples/env.sh) is sourced by every script. It reads these
variables, each with a default you can override from the shell:

| Variable | Default | Meaning |
|---|---|---|
| `CONDA_ENV` | `$HOME/miniforge3/envs/saige-build` | the environment from step 1 |
| `CUDA_HOME` | `/usr/local/cuda` | CUDA toolkit root (GPU build only) |
| `WORK` | `./saige_example` | where the examples write data and results |
| `PLINK2` | `plink2` | only used to simulate the tutorial data |
| `CONDA_SH` | `<conda base>/etc/profile.d/conda.sh`, derived from `CONDA_ENV` | set it when the environment is not under `<conda base>/envs/` |

It also exports `R_HOME=$CONDA_PREFIX/lib/R`, which step 1 needs at run time.

### 3. Build

```bash
bash docs/examples/01_build.sh cpu        # CPU only, no CUDA needed
bash docs/examples/01_build.sh gpu 70     # CUDA build for compute capability 7.0
```

The script runs, in each program's directory:

| | CPU build | GPU build (`SM` = compute capability) |
|---|---|---|
| step 1 (`step1_saige-null/`) | `make clean && make -j8 NVCC=none` | `make clean && make -j8 NVCC=$CUDA_HOME/bin/nvcc GPU_SM=sm_$SM` |
| step 2 (`step2_saige-step2/`) | `make clean && make -j8` | `make clean && make -j8 USE_CUDA=1 SM=$SM NVCC=$CUDA_HOME/bin/nvcc` |

Results: `step1_saige-null/saige-null`, `step2_saige-step2/saige-step2` and the
converter `step2_saige-step2/tools/sgs2txt`. Always `make clean` before
switching between the CPU and GPU builds. A GPU build also runs on machines
without a GPU (it falls back to the CPU).

The two Makefiles spell the GPU architecture differently: step 1 takes
`GPU_SM=sm_XX`, step 2 takes `SM=XX`.

| GPU | step 1 | step 2 |
|---|---|---|
| V100 | `GPU_SM=sm_70` (default) | `SM=70` (default) |
| T4 | `GPU_SM=sm_75` | `SM=75` |
| A100 | `GPU_SM=sm_80` | `SM=80` |
| A10, RTX 30xx | `GPU_SM=sm_86` | `SM=86` |
| L4, RTX 40xx | `GPU_SM=sm_89` | `SM=89` |
| H100 | `GPU_SM=sm_90` | `SM=90` |

Only `sm_70` was built and run for this guide.

Step 1's Makefile finds `nvcc` under `/usr/local/cuda*/bin` on its own; if no
`nvcc` is found it silently builds the CPU version. Step 2 builds the CPU
version unless `USE_CUDA=1` is given.

## Tutorial

Four binary traits, one run of each step, GPU on. From
`SAIGE-cpp-comparison/SAIGE_cpp_260716`:

```bash
export WORK=$PWD/saige_example            # every file below lands under here
bash docs/examples/02_simulate.sh         # genotypes + phenotypes
bash docs/examples/03_step1_binary.sh     # step 1: 4 binary traits, full GRM
bash docs/examples/04_step2_binary.sh     # step 2: all 5,000 markers x 4 traits
```

What each script writes:

```
$WORK/
├── data/
│   ├── geno.bed .bim .fam          PLINK genotypes (5,000 x 5,000, chromosomes 1-2)
│   ├── geno.pgen .pvar .psam       same genotypes, PGEN
│   ├── geno.bgen .sample           same genotypes, BGEN 1.2
│   ├── geno.vcf.gz                 same genotypes, VCF with DS
│   ├── dosage.pgen .pvar .psam     1,000 markers with fractional dosages
│   └── pheno.txt                   IID b1 b2 b3 b4 q1 q2 x1 x2
├── step1_bin/
│   ├── step1.yaml  step1.log
│   ├── models/b1/ ... models/b4/   one null-model directory per trait
│   └── vr_b1.varianceRatio.txt ... one variance-ratio file per trait
└── step2_bin/
    ├── step2.yaml  step2.log
    └── out/b1.txt ... out/b4.txt   association results, one file per trait
```

The final results are `$WORK/step2_bin/out/<trait>.txt`, one row per marker
(columns: [Step 2, results](step2.md#results)). The step-2 log reports
whether the GPU was used:

```
[b1] 5000 markers were tested (4769 on the GPU + 231 with the device SPA, 0 via the scalar CPU path; 49 Firth fits on the device).
  GPU coverage: 19999 / 20000 pairs (99.995%)
```

On a machine without a usable GPU the same run prints
`useGPU: refused, running on the CPU (...)` and produces the same files.

More examples, each runnable after `02_simulate.sh`:

| Script | What it does |
|---|---|
| `05_step1_quant_loco.sh` | step 1, two quantitative traits, LOCO |
| `06_step2_quant_loco.sh` | step 2 for those models, one run per chromosome |
| `07_sparse_grm.sh` | build a sparse GRM, fit on it, test |
| `08_step2_sgs.sh` | step 2 with binary output, then convert to text (needs 04) |
| `09_step2_formats.sh` | step 2 on PGEN, dosage PGEN, BGEN, VCF (needs 03) |
| `10_hpc_chrom_job.sh`, `11_hpc_concat.sh` | per-chromosome jobs and merging (needs 05) |
