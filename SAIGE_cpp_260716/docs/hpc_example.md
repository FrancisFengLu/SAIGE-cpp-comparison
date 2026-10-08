---
title: HPC example
nav_order: 5
---

# HPC example

Step 1 runs once per set of traits. Step 2 is split into one job per chromosome
(required with LOCO models, useful anyway); every job writes its own files, which
are concatenated at the end. All traits of the step-1 run go into one step-2 job
(`--step1Dir`); the genotypes are then read once for all of them.

## 1. Step 1 (one job)

```bash
bash docs/examples/05_step1_quant_loco.sh      # writes $WORK/step1_qt/
```

## 2. Step 2, one job per chromosome

[`examples/10_hpc_chrom_job.sh`](examples/10_hpc_chrom_job.sh) takes the
chromosome from its argument or from `SLURM_ARRAY_TASK_ID` and runs step 2 on
that chromosome:

```bash
#!/bin/bash
#SBATCH --job-name=saige-step2
#SBATCH --array=1-2                  # one task per chromosome (1-22 for a real genome)
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --gres=gpu:1                 # drop this line (and --useGPU) for CPU nodes
#SBATCH --output=saige-step2_%a.log
set -euo pipefail
source "${EXAMPLES:-$(dirname "$0")}/env.sh"
CHR=${1:-${SLURM_ARRAY_TASK_ID:?give a chromosome}}
cd "$WORK"
mkdir -p hpc
$SAIGE step2 \
  --step1Dir step1_qt \
  --plinkFile data/geno \
  --minMAF 0.01 \
  --LOCO=TRUE \
  --chrom "$CHR" \
  --nThreads "${SLURM_CPUS_PER_TASK:-8}" \
  --useGPU \
  --outDir hpc/chr$CHR > hpc/chr$CHR.log 2>&1
echo "chr$CHR done: $(grep -h 'GPU coverage\|useGPU: refused' hpc/chr$CHR.log)"
```

Without a scheduler, loop over the chromosomes:

```bash
for c in 1 2; do bash docs/examples/10_hpc_chrom_job.sh $c; done
```

```
chr1 done:   GPU coverage: 4910 / 4910 pairs (100%)
chr2 done:   GPU coverage: 4890 / 4890 pairs (100%)
```

On SLURM, submit the same file as an array. SLURM runs a copy of the script, so
tell it where `env.sh` is and pass `WORK`:

```bash
sbatch --export=ALL,EXAMPLES=$PWD/docs/examples,WORK=$WORK docs/examples/10_hpc_chrom_job.sh
```

(The SLURM path was checked by running the script with `SLURM_ARRAY_TASK_ID=2`
and `SLURM_CPUS_PER_TASK=4` set by hand; no SLURM cluster was used.)

With a genotype file per chromosome, change `--plinkFile` to that chromosome's
prefix (e.g. `data/geno_chr$CHR`); everything else stays the same.

## 3. Concatenate

After all jobs have finished, [`examples/11_hpc_concat.sh`](examples/11_hpc_concat.sh)
writes one genome-wide file per trait, keeping the header once:

```bash
bash docs/examples/11_hpc_concat.sh
```

```bash
for t in q1 q2; do
  first=1
  for CHR in 1 2; do                 # 1..22 for a real genome
    if [ $first = 1 ]; then cat "$H/chr$CHR/$t.txt"; first=0
    else tail -n +2 "$H/chr$CHR/$t.txt"; fi
  done > "$H/all/$t.txt"
done
```

Results: `$WORK/hpc/all/q1.txt`, `$WORK/hpc/all/q2.txt`.

## Saving storage and GPU time with binary output

With many traits, `--outputFormat sgs` in the step-2 jobs writes compact binary
files; convert them to text later on a CPU-only machine with `saige-gpu-cpp sgs2txt`
([Step 2, Binary output](step2.md#binary-output-sgs)). Each chromosome job writes
its own `.sgs` files, so convert per chromosome, then concatenate as above.

## Resource notes

- Threads: `--nThreads` in both steps. Match it to `--cpus-per-task`.
- One GPU per job; several jobs on one node each need their own GPU
  (SLURM sets `CUDA_VISIBLE_DEVICES`).
- GPU memory: see [GPU memory](gpu.md#gpu-memory-step-2).
- Step 1 with a sparse GRM runs single-threaded (it sets `nthreads: 1` itself).
- Each job writes its resolved config next to its results (`hpc/chr<N>/step2.yaml`).
