#!/bin/bash
# timing.sh: g200k (N = 50,000, 200,000 markers, no missing calls) x the bingpu_test bm_full models (8 binary
# traits, 4 distinct sample lists) cycled to P = 128, sgs output, the R defaults (+ Firth cutoff 0.05 in the
# F=1 cells), 8 threads. before = origin/rdefaults binary (own-sample-set runs: host tail), after = ownstats.
# Timing grade: TIMING_LOCK per cell, drop_caches before each cell, n = 2.
# Then one non-timing pair on bt (26,060 markers, 3,640 with missing calls) under isnoadjCov: false +
# impute_method: mean, where every pair carries a non-zero d and the missing-cell sums run on the device.
set -u
S=/opt/saige/logs/ownstats
LOCK=/opt/saige/logs/TIMING_LOCK
TAG=ownstats
BASE=$S/bin_base/saige-step2
NEW=/opt/saige/worktrees/ownstats/SAIGE_cpp_260716/bin/saige-step2
D=/opt/saige/data/bingpu_test
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
gen(){  # gen <plink> <F> <outdir> <extra...>
  local PF=$1 F=$2 OD=$3; shift 3
  echo "genoType: plink"; echo "plinkFile: $PF"
  echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isMoreOutput: false"
  if [ "$F" = 1 ]; then echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.05"; fi
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  echo "useGPU: true"; echo "outputFormat: sgs"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  local k
  for k in $(seq 1 128); do
    local m=$(( (k-1) % 8 + 1 ))
    echo "  - traitName: y$k"; echo "    modelFile: $D/models/bm_full/out/m/bm$m"
    echo "    varianceRatioFile: $D/models/bm_full/out/mvr_bm$m.varianceRatio.txt"; echo "    outputFile: $OD/y$k.txt"
  done
}
drop(){ sync; sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'; }
cell(){  # cell <name> <bin> <plink> <F> <cold 1|0> <timed 1|0> <extra...>
  local NAME=$1 BIN=$2 PF=$3 F=$4 COLD=$5 TIMED=$6; shift 6
  local DD=$Q/$NAME
  [ -f $DD/wall ] && return
  rm -rf "$DD"; mkdir -p "$DD/out"
  gen $PF $F $DD/out "$@" > $DD/cfg.yaml
  if [ "$TIMED" = 1 ]; then acquire_lock "timed cell $NAME running; do not start anything"; fi
  while [ $(freegb) -lt 15 ]; do log "only $(freegb)G free"; sleep 60; done
  if [ "$COLD" = 1 ]; then drop; else cat $PF.bed > /dev/null; fi
  log "START $NAME load=$(cut -d' ' -f1-3 /proc/loadavg)"
  ( cd $DD && /usr/bin/time -v -o time.txt $BIN cfg.yaml > log.txt 2>&1 ); local rc=$?
  local w=$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K.*' $DD/time.txt | awk -F: '{if(NF==3) printf "%.2f",$1*3600+$2*60+$3; else printf "%.2f",$1*60+$2}')
  echo "$w" > $DD/wall
  log "DONE  $NAME rc=$rc wall=$w $(grep -h 'gpuDeviceStats: [0-9]\|gpuDeviceStats: off' $DD/log.txt | cut -c1-160)"
  grep -h "\[gpu breakdown\]\|\[gpu pipeline\]\|\[gpu busy\]\|device SPA:\|device Firth:\|stats kernel" $DD/log.txt | cut -c1-240 >> $LOG
  rm -rf $DD/out
  if [ "$TIMED" = 1 ]; then release_lock; fi
}
G=/opt/saige/logs/binsplit/data/g200k
for r in 1 2; do
  for F in 0 1; do
    cell g200k_before_f${F}_r$r $BASE $G $F 1 1
    cell g200k_after_f${F}_r$r  $NEW  $G $F 1 1
  done
done
# bt, mean imputation: the missing-cell sums run for every slot with missing calls
cell bt_adj_before $BASE $D/geno/bt 0 0 0 "isnoadjCov: false" "impute_method: mean"
cell bt_adj_after  $NEW  $D/geno/bt 0 0 0 "isnoadjCov: false" "impute_method: mean"
cell bt_def_before $BASE $D/geno/bt 0 0 0
cell bt_def_after  $NEW  $D/geno/bt 0 0 0
log "ALL DONE"
