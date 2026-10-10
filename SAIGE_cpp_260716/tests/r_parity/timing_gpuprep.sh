#!/bin/bash
# timing_gpuprep.sh: g200k (N = 50,000, 200,000 markers, no missing calls) at P = 128, sgs output, the R
# defaults, 8 threads, useGPU: true -- before (origin/ownstats binary) vs after (gpuprep):
#   same    the 32 balanced binary models (binsplit/models/spa) cycled to 128: one sample list
#   own4    the bingpu_test bm_full models (8 binary, 4 distinct missing-phenotype patterns) cycled to 128
#   own128  128 models with 128 distinct sample lists (subset_model.py of the bm_full models: 1-5% of the
#           samples dropped per model, every list different -- a timing-only construction)
# plus qsp: the 8 qm_sparse_fast quantitative models (own lists) cycled to 32, isFastTest false (every
# marker on the sparse variance), before = refused by the old binary (runs on the CPU).
# Timing grade: TIMING_LOCK per cell, drop_caches before each cell, n = 2.
set -u
S=/opt/saige/logs/gpuprep
LOCK=/opt/saige/logs/TIMING_LOCK
TAG=gpuprep
BASE=${BASE:-/opt/saige/worktrees/ownstats/SAIGE_cpp_260716/step2_saige-step2/saige-step2}
NEW=${NEW:-/opt/saige/worktrees/gpuprep/SAIGE_cpp_260716/step2_saige-step2/saige-step2}
D=/opt/saige/data/bingpu_test
BAL=/opt/saige/logs/binsplit/models/spa
OWN=$S/models/own128
Q=$S/timing; mkdir -p $Q
LOG=$Q/driver.log
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
log(){ echo "$(date '+%F %T') $*" >> $LOG; }
idle(){
  [ -e /opt/saige/logs/TIMING_PAUSE ] && return 1
  [ -s $LOCK ] && ! grep -q "$TAG" $LOCK 2>/dev/null && return 1
  if ! grep -q "$TAG" $LOCK 2>/dev/null; then
    ps -eo comm= | grep -qE '^(saige-step2|saige-step2.cuda|saige-null|spadrv|spa_gpu_check|sgs2txt|nvcc|cc1plus|make|Rscript|R)$' && return 1
  fi
  return 0
}
acquire_lock(){
  local n=0
  while [ $n -lt 6 ]; do if idle; then n=$((n+1)); else n=0; fi; sleep 5; done
  printf '%s: %s (%s)\n' "$TAG" "$1" "$(date '+%F %T')" > $LOCK
}
release_lock(){ grep -q "$TAG" $LOCK 2>/dev/null && rm -f $LOCK; true; }
freegb(){ df --output=avail -BG / | tail -1 | tr -dc 0-9; }
models(){  # models <kind> <P> <outdir>
  local kind=$1 P=$2 OD=$3 k m
  for k in $(seq 1 $P); do
    case $kind in
      same)   m=$(( (k-1) % 32 + 1 )); echo "  - traitName: y$k"; echo "    modelFile: $BAL/m/y$m"; echo "    varianceRatioFile: $BAL/mvr_y$m.varianceRatio.txt";;
      own4)   m=$(( (k-1) % 8 + 1 ));  echo "  - traitName: y$k"; echo "    modelFile: $D/models/bm_full/out/m/bm$m"; echo "    varianceRatioFile: $D/models/bm_full/out/mvr_bm$m.varianceRatio.txt";;
      own128) m=$(( (k-1) % 8 + 1 ));  echo "  - traitName: y$k"; echo "    modelFile: $OWN/y$k"; echo "    varianceRatioFile: $D/models/bm_full/out/mvr_bm$m.varianceRatio.txt";;
      qsp)    m=$(( (k-1) % 8 + 1 ));  echo "  - traitName: y$k"; echo "    modelFile: $D/qm/models/qm_sparse_fast/out/m/q$m"; echo "    varianceRatioFile: $D/qm/models/qm_sparse_fast/out/mvr_q$m.varianceRatio.txt";;
    esac
    echo "    outputFile: $OD/y$k.txt"
  done
}
gen(){  # gen <kind> <P> <outdir> <extra...>
  local kind=$1 P=$2 OD=$3; shift 3
  echo "genoType: plink"; echo "plinkFile: /opt/saige/logs/binsplit/data/g200k"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isMoreOutput: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "useGPU: true"; echo "outputFormat: sgs"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  models $kind $P $OD
}
drop(){ sync; sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'; }
cell(){  # cell <name> <bin> <kind> <P> <extra...>
  local NAME=$1 BIN=$2 KIND=$3 P=$4; shift 4
  local DD=$Q/$NAME
  [ -f $DD/wall ] && return
  rm -rf "$DD"; mkdir -p "$DD/out"
  gen $KIND $P $DD/out "$@" > $DD/cfg.yaml
  acquire_lock "timed cell $NAME running; do not start anything"
  while [ $(freegb) -lt 15 ]; do log "only $(freegb)G free"; sleep 60; done
  drop
  log "START $NAME load=$(cut -d' ' -f1-3 /proc/loadavg)"
  ( cd $DD && /usr/bin/time -v -o time.txt $BIN cfg.yaml > log.txt 2>&1 ); local rc=$?
  local w=$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K.*' $DD/time.txt | awk -F: '{if(NF==3) printf "%.2f",$1*3600+$2*60+$3; else printf "%.2f",$1*60+$2}')
  echo "$w" > $DD/wall
  log "DONE  $NAME rc=$rc wall=$w $(grep -h 'useGPU: refused\|gpuDeviceStats: [0-9]' $DD/log.txt | head -2 | cut -c1-160 | tr '\n' ';')"
  grep -h "\[gpu breakdown\]\|\[gpu pipeline\]\|\[gpu busy\]\|\[gpu device time\]\|\[gpu overlap\]\|device SPA:\|stats kernel" $DD/log.txt | cut -c1-300 >> $LOG
  rm -rf $DD/out
  release_lock
}
for r in 1 2; do
  for kind in same own4 own128; do
    cell ${kind}_before_r$r $BASE $kind 128
    cell ${kind}_after_r$r  $NEW  $kind 128
  done
  cell qsp_before_r$r $BASE qsp 32 "blockSparseSigma: true"
  cell qsp_after_r$r  $NEW  qsp 32 "blockSparseSigma: true"
done
log "ALL DONE"
