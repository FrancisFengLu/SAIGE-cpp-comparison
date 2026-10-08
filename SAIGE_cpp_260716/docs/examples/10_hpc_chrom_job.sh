#!/bin/bash
#SBATCH --job-name=saige-step2
#SBATCH --array=1-2                  # one task per chromosome (1-22 for a real genome)
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --gres=gpu:1                 # drop this line (and set useGPU: false) for CPU nodes
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
D=$WORK/data
M=$WORK/step1_qt
O=$WORK/hpc/chr$CHR
mkdir -p "$O"
{
cat <<YAML
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first
minMAF: 0.01
LOCO: true
chrom: "$CHR"
nThreads: ${SLURM_CPUS_PER_TASK:-8}
useGPU: true
outputFormat: text
models:
YAML
for t in q1 q2; do
cat <<YAML
  - traitName: $t
    modelFile: $M/models/$t
    varianceRatioFile: $M/vr_$t.varianceRatio.txt
    outputFile: $O/$t.txt
YAML
done
} > "$O/step2.yaml"
"$S2" "$O/step2.yaml" > "$O/step2.log" 2>&1
echo "chr$CHR done: $(grep -h 'GPU coverage\|useGPU: refused' "$O/step2.log")"
