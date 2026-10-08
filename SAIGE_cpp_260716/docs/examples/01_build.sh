#!/bin/bash
# 01_build.sh [cpu|gpu] [SM]
#   cpu       CPU-only build (no CUDA toolkit needed)
#   gpu SM    CUDA build for compute capability SM (70 = V100, 80 = A100,
#             86 = A10/RTX 30xx, 89 = L4/RTX 40xx, 90 = H100); default 70.
#             Several at once: gpu "70 80 90"
# Builds $SAIGE_HOME/bin/: saige-gpu-cpp (the command line), saige-null and
# saige-step2 (the two engines it runs), sgs2txt.
set -euo pipefail
source "$(dirname "$0")/env.sh"
MODE=${1:-cpu}
SM=${2:-70}
cd "$SAIGE_HOME"
if [ "$MODE" = gpu ]; then
  make -j"${JOBS:-8}" USE_CUDA=1 SM="$SM" NVCC="$CUDA_HOME/bin/nvcc"
else
  make -j"${JOBS:-8}"
fi
ls -l bin/
