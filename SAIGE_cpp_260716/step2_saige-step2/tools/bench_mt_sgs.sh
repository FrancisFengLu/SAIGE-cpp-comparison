#!/bin/bash
# bench_mt_sgs.sh -- outputFormat: sgs vs text on the CPU multi-trait path
# (mainMarkerMT). The GPU-path equivalent is tools/bench_writer.sh.
#
# One run per invocation. g1m = 50,000 samples x 1,000,000 markers,
# quantitative, sparse-GRM null models. The page cache is dropped first (needs
# passwordless sudo) so the .bed read costs the same in both arms; the arms are
# then compared on wall clock, on the `[mt breakdown] output write` line, and
# on output bytes.
#
# usage: bench_mt_sgs.sh <P> <fmt: text|sgs> [tag]
# results: optimization/torchgwas2/MT_SGS.md section 3
set -u
P=$1; FMT=$2; TAG=${3:-}
BIN=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/saige-step2
DATA=/opt/saige/logs/tg2_step2
M=$DATA/runs/s1_cpp_P128/out
W=/opt/saige/logs/mtsgs/bench/P${P}_${FMT}${TAG}
rm -rf "$W"; mkdir -p "$W/out"
{
  echo "genoType: plink"; echo "plinkFile: $DATA/data/g1m"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  [ "$FMT" = sgs ] && echo "outputFormat: sgs"
  echo "models:"
  for k in $(seq 1 "$P"); do
    echo "  - traitName: y$k"
    echo "    modelFile: $M/m/y$k"
    echo "    varianceRatioFile: $M/mvr_y$k.varianceRatio.txt"
    echo "    outputFile: $W/out/y$k.txt"
  done
} > "$W/cfg.yaml"
pgrep -x 'saige-step2|saige-null|R|Rscript' > "$W/busy.txt"
sync; sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches'
/usr/bin/time -v "$BIN" "$W/cfg.yaml" > "$W/log.txt" 2> "$W/time.txt"
echo "rc=$?" > "$W/rc"
du -sb "$W/out" | cut -f1 > "$W/bytes"
grep -E "Elapsed|Maximum resident" "$W/time.txt"
grep -E "\[mt breakdown\]|outputFormat: sgs wrote" "$W/log.txt"
echo "bytes: $(cat "$W/bytes")"
