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
| [Step 1: fit the null model](step1.md) | `saige-gpu-cpp step1`: flags, inputs, outputs; YAML config |
| [Step 2: test genetic variants](step2.md) | `saige-gpu-cpp step2`: flags, inputs, outputs; YAML config |
| [GPU](gpu.md) | GPU build, switches, what runs on the GPU, memory |
| [HPC example](hpc_example.md) | one job per chromosome, concatenating results |
| [Troubleshooting](troubleshooting.md) | error and log messages and what they mean |
| [Collaborator tests](collaborator_tests.md) | performance/accuracy test protocol for running on your own data |

Every command and config on these pages was run on simulated data
(`plink2 --dummy`, 5,000 samples x 5,000 markers) with the scripts in
[`examples/`](examples/).

## Analysis workflow

1. **Step 1** (`saige-gpu-cpp step1`): fit one null model per trait from the
   phenotype/covariate table and the PLINK genotypes (full GRM) or a sparse GRM.
   Output: one directory with a model and a variance-ratio file per trait.
2. **Step 2** (`saige-gpu-cpp step2`): test every marker of a genotype file
   (PLINK, PGEN, BGEN or VCF) against the step-1 models. Output: one result
   table per trait.

The flags are R SAIGE's (`step1_fitNULLGLMM.R`, `step2_SPAtests.R`), with a few
additions (several traits per run, `--useGPU`, `--outDir`, `--step1Dir`). Each
run writes its settings as a YAML config into its output directory and runs the
step's engine (`saige-null`, `saige-step2`) on it; the engines can also be run
on a hand-written YAML config.

## Installation

No container image is published. Build from source.

### Requirements

- Linux x86-64. The build uses `-march=native` by default, so build on the
  machine type you will run on, or build with `make ARCH=x86-64-v3` for a binary
  that runs on any x86-64 CPU with AVX2 (Intel Haswell / AMD Zen and newer).
- A conda environment with the libraries below (any other source of the same
  libraries works too). Neither program needs R.
- For the GPU build only: CUDA toolkit 12.x (`nvcc`, cuBLAS) and an NVIDIA GPU
  with compute capability 7.0 or newer. Tested: CUDA 12.9, Tesla V100 (sm_70).

### 1. Create the build environment

```bash
mamba create -y -p $HOME/miniforge3/envs/saige-build -c conda-forge -c bioconda \
    gxx_linux-64=12 make pkg-config armadillo openblas superlu yaml-cpp htslib \
    zstd zlib sqlite boost-cpp eigen
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

### 3. Build

With the environment from step 1 active (`conda activate $HOME/miniforge3/envs/saige-build`):

```bash
make                          # CPU only, no CUDA needed
make USE_CUDA=1 SM=70         # CUDA build for compute capability 7.0
```

`bash docs/examples/01_build.sh cpu` / `bash docs/examples/01_build.sh gpu 70` do
the same after activating the environment named in `env.sh`. The result is one
directory, `bin/`:

| Program | |
|---|---|
| `bin/saige-gpu-cpp` | the command line: `saige-gpu-cpp step1 ...`, `step2 ...`, `sgs2txt ...` |
| `bin/saige-null` | step-1 engine (`saige-null -c step1.yaml`) |
| `bin/saige-step2` | step-2 engine (`saige-step2 step2.yaml`) |
| `bin/sgs2txt` | converter for the binary step-2 output |

`saige-gpu-cpp` runs the engines in its own directory, so keep the four
together (copy or link `bin/` anywhere; a symlink to `saige-gpu-cpp` works too).
The programs find their libraries without `conda activate`. A GPU build also
runs on machines without a GPU (it falls back to the CPU).

Make variables:

| Variable | Default | Meaning |
|---|---|---|
| `USE_CUDA` | `0` | `1`: GPU build (needs `nvcc`, found as `$CUDA_HOME/bin/nvcc`, `CUDA_HOME` default `/usr/local/cuda`) |
| `SM` | `70` | GPU compute capability; several as `SM="70 80 90"` |
| `ARCH` | `native` | CPU instruction set (`-march`): `native` or e.g. `x86-64-v3` |
| `PROGRAM` | `saige-gpu-cpp` | name of the command-line program |
| `JOBS` | `8` | parallel compile jobs |

Changing `ARCH`, `USE_CUDA` or `SM` rebuilds both engines from scratch.

`SM` is the card's compute capability, a fixed number per GPU model (not a
setting). Look it up on the machine you will run on:

```bash
nvidia-smi --query-gpu=name,compute_cap --format=csv    # e.g. "Tesla V100-SXM2-16GB, 7.0" -> SM=70
```

Compute capability 7.0 or newer is required.

| GPU | `SM` | fp64 speed (step 2 runs in fp64) |
|---|---|---|
| V100 | `70` (default) | full — **recommended, tested** |
| T4 | `75` | 1/32 of fp32 — slow |
| A100, A30 | `80` | full — recommended |
| A10, A10G, RTX 30xx | `86` | 1/64 of fp32 — slow |
| L4, L40S, RTX 40xx | `89` | 1/64 of fp32 — slow |
| H100 | `90` | full — recommended |

Step 2 computes in double precision (fp64) only, so it runs on every card above
but is much slower on the cards marked slow. An fp32 version is not available yet.

One build for several GPU types: `make USE_CUDA=1 SM="70 80 90"` puts device
code for each listed compute capability into the programs (plus PTX for the
highest, which newer cards compile at start-up). Measured on this guide's
machine: build 5.0 min instead of 3.3 min, `saige-step2` 25.3 MB instead of
13.1 MB, `saige-null` 5.2 MB instead of 4.7 MB; results on the V100 identical.

`make ARCH=x86-64-v3` gave byte-identical results and the same run times as the
`native` build on this guide's machine (Intel Haswell-class CPU) for the
tutorial and a 50,000-sample test (`tests/cli/arch_compare.sh`).

Only `SM=70` was run for this guide.

## Tutorial

Four binary traits, one run of each step, GPU on. From
`SAIGE-cpp-comparison/SAIGE_cpp_260716`, after the GPU build:

```bash
export WORK=$PWD/saige_example            # every file below lands under here
bash docs/examples/02_simulate.sh         # genotypes + phenotypes
bash docs/examples/03_step1_binary.sh     # step 1: 4 binary traits, full GRM
bash docs/examples/04_step2_binary.sh     # step 2: all 5,000 markers x 4 traits
```

The two steps, as the scripts run them in `$WORK`:

```bash
saige-gpu-cpp step1 \
  --plinkFile data/geno \
  --phenoFile data/pheno.txt \
  --phenoCol b1,b2,b3,b4 \
  --covarColList x1,x2 \
  --traitType binary \
  --LOCO=FALSE \
  --nThreads 8 \
  --useGPU \
  --IsOverwriteVarianceRatioFile=TRUE \
  --outDir step1_bin > step1_bin.log 2>&1

saige-gpu-cpp step2 \
  --step1Dir step1_bin \
  --plinkFile data/geno \
  --minMAF 0 \
  --minMAC 1 \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth 0.01 \
  --nThreads 8 \
  --useGPU \
  --outDir step2_bin > step2_bin.log 2>&1
```

(`saige-gpu-cpp` is `bin/saige-gpu-cpp`; the scripts call it as `$SAIGE`.)
Flags: [Step 1](step1.md#flags), [Step 2](step2.md#flags);
`saige-gpu-cpp step1 --help` and `step2 --help` list them with their defaults.

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
├── step1_bin.log                   step-1 log
├── step1_bin/
│   ├── step1.yaml                  the settings of this run (YAML), written by step1
│   ├── models/b1/ ... models/b4/   one null-model directory per trait
│   └── vr_b1.varianceRatio.txt ... one variance-ratio file per trait
├── step2_bin.log                   step-2 log
└── step2_bin/
    ├── step2.yaml                  the settings of this run (YAML), written by step2
    └── b1.txt ... b4.txt           association results, one file per trait
```

The final results are `$WORK/step2_bin/<trait>.txt`, one row per marker
(columns: [Step 2, results](step2.md#results)). The step-2 log reports
whether the GPU was used:

```
[b1] 5000 markers were tested (4769 on the GPU + 231 with the device SPA, 0 via the scalar CPU path; 49 Firth fits on the device).
  GPU coverage: 19999 / 20000 pairs (99.995%)
```

On a machine without a usable GPU the same run prints
`useGPU: refused, running on the CPU (...)` and produces the same files.

Each `stepN.yaml` repeats the run without the command line:
`bin/saige-null -c step1_bin/step1.yaml` and `bin/saige-step2 step2_bin/step2.yaml`
(the exact commands are in the file's header).

More examples, each runnable after `02_simulate.sh`:

| Script | What it does |
|---|---|
| `05_step1_quant_loco.sh` | step 1, two quantitative traits, LOCO |
| `06_step2_quant_loco.sh` | step 2 for those models, one run per chromosome |
| `07_sparse_grm.sh` | build a sparse GRM, fit on it, test |
| `08_step2_sgs.sh` | step 2 with binary output, then convert to text (needs 04) |
| `09_step2_formats.sh` | step 2 on PGEN, dosage PGEN, BGEN, VCF (needs 03) |
| `10_hpc_chrom_job.sh`, `11_hpc_concat.sh` | per-chromosome jobs and merging (needs 05) |
