#!/bin/bash
# sgsPrecision f64 vs f32: the same P sweep as bench3.sh (10^6 markers x 50,000
# samples, quantitative, GPU path, page cache dropped before every timed run,
# one job at a time).  For each P: size, write time, end-to-end wall.
# Fidelity is measured separately by sgsfid on the fp64 files.
set -u
S2=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2
BIN=$S2/saige-step2.cuda
CVT=$S2/tools/sgs2txt
FID=$S2/tools/sgs_fidelity
D=/opt/saige/logs/gpu_step2/writer/f32/bench4
M=/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out
BED=/opt/saige/logs/tg2_step2/data/g1m
mkdir -p $D

gen() { local P=$1 OD=$2; shift 2; mkdir -p $OD
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  for k in $(seq 1 $P); do
    echo "  - traitName: y$k"; echo "    modelFile: $M/m/y$k"
    echo "    varianceRatioFile: $M/mvr_y$k.varianceRatio.txt"
    echo "    outputFile: $OD/y$k.txt"
  done
}

one() { # tag cfg reps
  local tag=$1 cfg=$2 reps=$3 best=999999 r rc sec
  for r in $(seq 1 $reps); do
    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches >/dev/null; sleep 3
    /usr/bin/time -v $BIN $cfg > $D/$tag.log 2> $D/$tag.r$r.time
    rc=$?
    sec=$(awk -F': ' '/Elapsed \(wall/{print $NF}' $D/$tag.r$r.time |
          awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else if(NF==2)print $1*60+$2; else print $1}')
    echo "RUN $tag rep$r rc=$rc wall=${sec}s"
    if (( $(echo "$sec < $best" | bc -l) )); then best=$sec; fi
  done
  echo "RESULT $tag min_wall=${best}s (n=$reps)"
  grep -E "^  outputFormat|^  \[gpu breakdown\]" $D/$tag.log
}

for P in 1 8 32 128; do
  echo "######## P=$P  $(date +%T)"
  R=2; [ $P -ge 128 ] && R=1
  gen $P $D/s64_P$P "useGPU: true" "outputFormat: sgs"                      > $D/s64_P$P.yaml
  gen $P $D/s32_P$P "useGPU: true" "outputFormat: sgs" "sgsPrecision: fp32" > $D/s32_P$P.yaml
  one s64_P$P $D/s64_P$P.yaml $R
  one s32_P$P $D/s32_P$P.yaml $R
  b64=$(du -sb $D/s64_P$P | cut -f1); b32=$(du -sb $D/s32_P$P | cut -f1)
  echo "  bytes fp64 = $b64   fp32 = $b32   ratio = $(echo "scale=3; $b64/$b32" | bc)"
  for w in 64 32; do
    t0=$(date +%s.%N)
    $CVT -j 8 $D/s${w}_P$P/*.sgs >/dev/null 2>$D/cvt${w}_P$P.log
    echo "  sgs2txt fp$w rc=$? wall=$(echo "$(date +%s.%N) - $t0" | bc)s"
  done
  echo "---- fidelity from the fp64 file, P=$P ----"
  $FID -j 8 -m $D/s64_P$P/y1.txt.markers.sgs $D/s64_P$P/y*.txt.sgs
  rm -rf $D/s64_P$P $D/s32_P$P
  echo
done
echo "BENCH4 DONE $(date +%T)"
