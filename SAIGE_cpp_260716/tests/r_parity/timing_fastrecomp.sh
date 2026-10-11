#!/bin/bash
# g200k P = 128 (the 32 balanced binary models cycled, one sample list), R defaults + isFastTest: true,
# useGPU, sgs output, 8 threads: before (gpuprep binary) vs after (fastrecomp), one run each, hot cache,
# NOT timing grade -- the host tail and the fast-test lines only.
set -u
F=/opt/saige/logs/fastrecomp/g200k
BASE=/opt/saige/worktrees/gpuprep/SAIGE_cpp_260716/step2_saige-step2/saige-step2
NEW=/opt/saige/worktrees/fastrecomp/SAIGE_cpp_260716/step2_saige-step2/saige-step2
BAL=/opt/saige/logs/binsplit/models/spa
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
mkdir -p $F
log(){ echo "$(date '+%F %T') $*"; }
gen(){  # gen <outdir>
  echo "genoType: plink"; echo "plinkFile: /opt/saige/logs/binsplit/data/g200k"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isMoreOutput: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "useGPU: true"; echo "outputFormat: sgs"; echo "isFastTest: true"
  echo "models:"
  for k in $(seq 1 128); do m=$(( (k-1) % 32 + 1 ))
    echo "  - traitName: y$k"; echo "    modelFile: $BAL/m/y$m"; echo "    varianceRatioFile: $BAL/mvr_y$m.varianceRatio.txt"
    echo "    outputFile: $1/y$k.txt"; done
}
for which in before after; do
  [ $which = before ] && BIN=$BASE || BIN=$NEW
  while [ -e /opt/saige/logs/TIMING_LOCK ] || [ -e /opt/saige/logs/TIMING_PAUSE ]; do log "waiting for TIMING_LOCK/PAUSE"; sleep 60; done
  D=$F/$which; rm -rf $D; mkdir -p $D/out
  gen $D/out > $D/cfg.yaml
  log "START $which"
  ( cd $D && /usr/bin/time -v -o time.txt $BIN cfg.yaml > log.txt 2>&1 ); rc=$?
  log "DONE $which rc=$rc wall=$(grep -oP 'Elapsed \(wall clock\).*: \K.*' $D/time.txt)"
  grep -h "\[gpu breakdown\]\|\[gpu pipeline\]\|\[gpu busy\]\|gate:\|fast-test recompute\|gpuDeviceStats: [0-9]\|GPU coverage\|device SPA:" $D/log.txt | cut -c1-300
  rm -rf $D/out
done
