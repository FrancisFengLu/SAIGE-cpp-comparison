#!/bin/bash
# build_bins.sh SM  -- the three binaries the test protocol needs, copied to $BIN_DIR:
#   saige-null          step 1, CUDA build          (= docs/examples/01_build.sh gpu SM)
#   saige-step2         step 2, CUDA build          (runs both the CPU path, useGPU: false,
#                                                    and the GPU path, useGPU: true)
#   saige-step2.phase   step 2, CPU build with PHASE_TIMING=1 (per-stage timers; block 3b)
#   sgs2txt             converter (not used by the protocol, copied for completeness)
# SM = compute capability without the dot (70 V100, 80 A100, 86 A10, 89 L4, 90 H100).
# Builds in place in the source tree (make clean between configurations), like 01_build.sh.
set -euo pipefail
source "$(dirname "$0")/common.sh"
SM=${1:?usage: build_bins.sh SM   (70 V100, 80 A100, 89 L4, 90 H100)}
JOBS=${JOBS:-$NPROC}
mkdir -p "$BIN_DIR"
LOG=$BIN_DIR/build.log; : > "$LOG"

JOBS=$JOBS bash "$COLLAB/../examples/01_build.sh" gpu "$SM" >> "$LOG" 2>&1
cp "$S1" "$S2" "$SGS2TXT" "$BIN_DIR/"

cd "$SAIGE_HOME/step2_saige-step2"
make clean > /dev/null
make -j"$JOBS" PHASE_TIMING=1 TARGET=saige-step2.phase saige-step2.phase >> "$LOG" 2>&1
cp saige-step2.phase "$BIN_DIR/"
# leave the tree as 01_build.sh left it (GPU build)
make clean > /dev/null
make -j"$JOBS" USE_CUDA=1 SM="$SM" NVCC="$CUDA_HOME/bin/nvcc" >> "$LOG" 2>&1

{ echo "commit $(git -C "$SAIGE_HOME" rev-parse HEAD)"
  echo "dirty  $(git -C "$SAIGE_HOME" status --short --untracked-files=no | wc -l) tracked file(s) modified"
  echo "SM     $SM"
  echo "cxx    $(${CXX:-x86_64-conda-linux-gnu-c++} --version 2>/dev/null | head -1)"
  echo "nvcc   $("$CUDA_HOME/bin/nvcc" --version | tail -1)"
  ( cd "$BIN_DIR" && md5sum saige-null saige-step2 saige-step2.phase sgs2txt )
} > "$BIN_DIR/BUILD_INFO.txt"
cat "$BIN_DIR/BUILD_INFO.txt"
echo BUILD_OK
