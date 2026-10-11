#!/bin/bash
# Quick (not timing-grade) device SPA / Firth per-pair times by covariate count: covt
# (gen_covp.py at N = M = 20,000, b1 + b2, one sample list; R step 1 at p = 4 9 13 24 40,
# r_runs_covp.sh), configuration defF (R defaults + Firth), useGPU, nThreads 4.
#   timing_covp.sh <new binary> <base binary> <OUT>
# base: p = 4 only (it refuses p > 8); new: every p, lib fp64 (default), own fp64
# (gpuSpaImpl: own) and lib fp32 (gpuPrecisionSPA / gpuPrecisionFirth: fp32). Each cell
# runs twice; the us/pair come from the binary's "device SPA:" / "device Firth:" lines
# (CUDA event time of the kernels / pairs; CUDA_MODULE_LOADING=EAGER so the module load is
# not inside it). Waits for TIMING_LOCK; not a lock holder.
set -u
NEW=$1; BASE=$2; OUT=$3
B=${COVLIM:-/opt/saige/logs/covlimit}
PS=${PS:-"4 9 13 24 40"}
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
export OPENBLAS_NUM_THREADS=1
# The kernels' first launch otherwise loads the module inside the timed region (the SPA
# library's 80 instantiations take ~0.13 s on this V100, more than its kernels at this N).
export CUDA_MODULE_LOADING=EAGER
RDA2ARMA=$(cd "$(dirname "$0")/../../step2_saige-step2/tools" && pwd)/rda_to_arma.R
log(){ echo "$(date '+%F %T') $*"; }
( source /opt/saige/SAIGE-work/optimization/comp8/env_r.sh
  for p in $PS; do for t in b1 b2; do
    M=$B/models/covt/p$p/$t
    [ -f $M/nullmodel.json ] && continue
    mkdir -p $M; Rscript $RDA2ARMA $B/R/covt/s1/p$p/$t.rda $M > $M/convert.log 2>&1 || { log "convert covt p=$p $t FAILED"; cat $M/convert.log; }
  done; done )
run(){  # run <label> <binary> <p> <rep> [extra keys...]
  local lab=$1 bin=$2 p=$3 rep=$4; shift 4
  local O=$OUT/$lab/p$p/r$rep; mkdir -p $O
  { echo "genoType: plink"; echo "plinkFile: $B/data/covt/g"; echo "AlleleOrder: alt-first"
    echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"; echo "LOCO: false"; echo "nThreads: 4"
    echo "is_Firth_beta: true"; echo "pCutoffforFirth: 0.01"; echo "useGPU: true"
    for k in "$@"; do echo "$k"; done
    echo "models:"
    for t in b1 b2; do
      echo "  - traitName: $t"; echo "    modelFile: $B/models/covt/p$p/$t"
      echo "    varianceRatioFile: $B/R/covt/s1/p$p/$t.varianceRatio.txt"; echo "    outputFile: $O/$t.txt"
    done; } > $O/cfg.yaml
  while [ -e /opt/saige/logs/TIMING_LOCK ]; do sleep 10; done
  local t0=$(date +%s.%N)
  $bin $O/cfg.yaml > $O/log.txt 2>&1; local rc=$?
  log "$lab p=$p r$rep rc=$rc $(printf %.1f $(echo "$(date +%s.%N) - $t0" | bc))s | $(grep -h 'device SPA: \|device Firth: ' $O/log.txt | sed 's/ of kernel time//; s/ kept the device result.*//; s/ of the .* Firth fits.*//' | tr '\n' ';' | cut -c1-200)"
}
for rep in 1 2; do
  run base $BASE 4 $rep
  for p in $PS; do
    run lib64 $NEW $p $rep
    run own64 $NEW $p $rep "gpuSpaImpl: own"
    run lib32 $NEW $p $rep "gpuPrecisionSPA: fp32" "gpuPrecisionFirth: fp32"
  done
done
log "timing runs done"
