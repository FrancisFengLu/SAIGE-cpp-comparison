#!/bin/bash
# Writer re-measure: the same P sweep as bench2.sh (10^6 markers x 50,000
# samples, quantitative, GPU path) with the new writer.
#   gtxt_P*  GPU, text output, fast formatter + per-trait parallel write
# (copy of logs/gpu_step2/writer/bench3.sh; the run it produced is recorded in
# out_fast.hpp. bench_writer_walltimes.sh pulls the real wall times out of the
# /usr/bin/time -v files, because this script prints min_wall=0s -- its awk
# splits the "Elapsed (wall clock) time (h:mm:ss or m:ss): M:SS" line on the
# wrong separator. bench2.sh had the same bug.)
#   gsgs_P*  GPU, outputFormat: sgs, then sgs2txt and cmp against gtxt
# One timed job at a time; page cache dropped before every run.
set -u
BIN=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/saige-step2.cuda
CVT=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/tools/sgs2txt
D=/opt/saige/logs/gpu_step2/writer/bench3
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

one() { # one tag cfg reps
  local tag=$1 cfg=$2 reps=$3
  local best=999999 r rc w sec
  for r in $(seq 1 $reps); do
    sync; echo 3 | sudo tee /proc/sys/vm/drop_caches >/dev/null; sleep 3
    /usr/bin/time -v $BIN $cfg > $D/$tag.log 2> $D/$tag.r$r.time
    rc=$?
    w=$(awk -F'clock) ' '/Elapsed .wall/{print $2}' $D/$tag.r$r.time)
    sec=$(echo "$w" | awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else print $1*60+$2}')
    echo "RUN $tag rep$r rc=$rc wall=${sec}s cpu=$(awk -F': ' '/Percent of CPU/{print $2}' $D/$tag.r$r.time) rssKB=$(awk -F': ' '/Maximum resident/{print $2}' $D/$tag.r$r.time)"
    if (( $(echo "$sec < $best" | bc -l) )); then best=$sec; fi
  done
  echo "RESULT $tag min_wall=${best}s (n=$reps)"
  grep -E "^  \[gpu breakdown\]|^  \[gpu device time\]|^  outputFormat|^  GPU coverage" $D/$tag.log
}

for P in 1 8 32 128; do
  echo "######## P=$P  $(date +%T)"
  R=2; [ $P -ge 128 ] && R=1
  gen $P $D/gtxt_P$P "useGPU: true"                        > $D/gtxt_P$P.yaml
  gen $P $D/gsgs_P$P "useGPU: true" "outputFormat: sgs"    > $D/gsgs_P$P.yaml
  one gtxt_P$P $D/gtxt_P$P.yaml $R
  one gsgs_P$P $D/gsgs_P$P.yaml $R
  echo "  bytes text = $(du -sb $D/gtxt_P$P | cut -f1)"
  echo "  bytes sgs  = $(du -sb $D/gsgs_P$P | cut -f1)"
  echo "---- round trip P=$P ----"
  t0=$(date +%s)
  $CVT -j 8 $D/gsgs_P$P/*.sgs 2> $D/cvt_P$P.log
  echo "  sgs2txt rc=$? wall=$(( $(date +%s) - t0 ))s"
  bad=0; rows=0
  for k in $(seq 1 $P); do
    cmp -s $D/gtxt_P$P/y$k.txt $D/gsgs_P$P/y$k.txt || { echo "  DIFF y$k"; bad=1; }
  done
  rows=$(( $(wc -l < $D/gtxt_P$P/y1.txt) - 1 ))
  if [ $bad -eq 0 ]; then
    echo "  ROUND TRIP OK: $P traits x $rows markers = $(( P * rows )) rows byte-identical"
  else
    echo "  ROUND TRIP FAILED"
  fi
  rm -rf $D/gtxt_P$P $D/gsgs_P$P
  echo
done
echo "BENCH3 DONE $(date +%T)"
