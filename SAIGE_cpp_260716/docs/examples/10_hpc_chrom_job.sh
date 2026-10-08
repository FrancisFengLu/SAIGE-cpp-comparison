#!/bin/bash
#SBATCH --job-name=saige-step2
#SBATCH --array=1-2                  # one task per chromosome (1-22 for a real genome)
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --gres=gpu:1                 # drop this line (and --useGPU) for CPU nodes
#SBATCH --output=saige-step2_%a.log
#
# 10_hpc_chrom_job.sh [CHR] -- step 2 for one chromosome, all traits of one step-1
# run (the LOCO models of 05_step1_quant_loco.sh). The chromosome comes from the
# argument or from SLURM_ARRAY_TASK_ID. Writes $WORK/hpc/chr<CHR>/<trait>.txt.
#   sbatch 10_hpc_chrom_job.sh              (SLURM array, one task per chromosome)
#   bash   10_hpc_chrom_job.sh 1            (one chromosome, any machine)
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
