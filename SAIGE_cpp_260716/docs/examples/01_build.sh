#!/bin/bash
# 01_build.sh [cpu|gpu] [SM]
#   cpu       CPU-only build of step 1 and step 2 (no CUDA toolkit needed)
#   gpu SM    CUDA build for compute capability SM (70 = V100, 80 = A100,
#             86 = A10/RTX 30xx, 89 = L4/RTX 40xx, 90 = H100); default 70
# Builds in place: step1_saige-null/saige-null, step2_saige-step2/saige-step2,
# step2_saige-step2/tools/sgs2txt.
set -euo pipefail
source "$(dirname "$0")/env.sh"
MODE=${1:-cpu}
SM=${2:-70}
JOBS=${JOBS:-8}

cd "$SAIGE_HOME/step1_saige-null"
make clean > /dev/null
if [ "$MODE" = gpu ]; then
  make -j"$JOBS" NVCC="$CUDA_HOME/bin/nvcc" GPU_SM=sm_"$SM"
else
  make -j"$JOBS" NVCC=none                 # NVCC=none: link the CPU stub, no CUDA
fi

cd "$SAIGE_HOME/step2_saige-step2"
make clean > /dev/null
if [ "$MODE" = gpu ]; then
  make -j"$JOBS" USE_CUDA=1 SM="$SM" NVCC="$CUDA_HOME/bin/nvcc"
else
  make -j"$JOBS"
fi
ls -l "$S1" "$S2" "$SGS2TXT"
