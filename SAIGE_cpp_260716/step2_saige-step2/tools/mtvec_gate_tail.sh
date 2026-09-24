#!/bin/bash
# mtvec_gate_tail.sh -- the gate the null data cannot give: p-values small
# enough to exercise mtVecQuantStats' boost fallback AND the "%.1fE%d"
# underflow branch of format_score_result.
#
# tests/make_extreme_model.py scales a quantitative null model's res and S_a by
# k, which multiplies the chi-square statistic by k^2 (see its docstring). k=6
# and k=15 drag the 5,000 g5k markers across p = 1 .. 1e-716, so both the
# fallback (p < 1e-5) and the p == 0 log form are hit thousands of times.
#
# Passes when every text file is byte identical between switch off and on, and
# between the pre-change binary and both.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
S2="$(dirname "$HERE")"
W=${W:-/opt/saige/logs/mtvec/tail}
BIN=${BIN:-/opt/saige/logs/mtvec/bin/saige-step2}
BASEBIN=${BASEBIN:-/opt/saige/logs/mtvec/bin/saige-step2.base}
M=/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out
BED=${BED:-/opt/saige/logs/tg2_step2/data/g5k}
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh

rm -rf "$W"; mkdir -p "$W/m" "$W/off" "$W/on" "$W/base"
cd "$S2/tests"
python3 make_extreme_model.py $M/m/y1 "$W/m/e6"  6   || exit 1
python3 make_extreme_model.py $M/m/y2 "$W/m/e15" 15  || exit 1
python3 make_extreme_model.py $M/m/y3 "$W/m/e3"  3   || exit 1

gen() {  # gen OUTDIR VEC
  local OD=$1 VEC=$2
  { echo "genoType: plink"; echo "plinkFile: $BED"
    echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
    echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
    echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
    echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
    echo "mtFoldQuantProj: true"
    echo "mtVecQuantStats: $VEC"
    echo "models:"
    local i=1
    for t in e3 e6 e15; do
      echo "  - traitName: $t"; echo "    modelFile: $W/m/$t"
      echo "    varianceRatioFile: $M/mvr_y$i.varianceRatio.txt"
      echo "    outputFile: $OD/$t.txt"
      i=$((i+1))
    done
    echo "  - traitName: plain"; echo "    modelFile: $M/m/y4"
    echo "    varianceRatioFile: $M/mvr_y4.varianceRatio.txt"
    echo "    outputFile: $OD/plain.txt"; }
}
gen "$W/off"  false > "$W/off.yaml"
gen "$W/on"   true  > "$W/on.yaml"
gen "$W/base" false > "$W/base.yaml"

"$BIN"     "$W/off.yaml"  > "$W/off.log"  2>&1 || { echo FAIL off;  tail -5 "$W/off.log";  exit 1; }
"$BIN"     "$W/on.yaml"   > "$W/on.log"   2>&1 || { echo FAIL on;   tail -5 "$W/on.log";   exit 1; }
"$BASEBIN" "$W/base.yaml" > "$W/base.log" 2>&1 || { echo FAIL base; tail -5 "$W/base.log"; exit 1; }

echo "===== p-value coverage (switch off run) ====="
for t in e3 e6 e15 plain; do
  f="$W/off/$t.txt"
  tot=$(( $(wc -l < "$f") - 1 ))
  lt5=$(awk -F'\t' 'NR>1 && $13 !~ /E[0-9]/ && $13+0 < 1e-5 {c++} END{print c+0}' "$f")
  logf=$(awk -F'\t' 'NR>1 && $13 ~ /^[0-9]\.[0-9]E-[0-9]+$/ {c++} END{print c+0}' "$f")
  mn=$(awk -F'\t' 'NR>1{print $13}' "$f" | sort -g | head -1)
  mx=$(awk -F'\t' 'NR>1{print $13}' "$f" | sort -g | tail -1)
  printf "  %-6s rows %5d   p<1e-5 (boost fallback) %5d   \"%%.1fE%%d\" log form %5d   p range %s .. %s\n" \
         "$t" "$tot" "$lt5" "$logf" "$mx" "$mn"
done

echo "===== byte identity ====="
ok=1
for t in e3 e6 e15 plain; do
  a="$W/off/$t.txt"; b="$W/on/$t.txt"; c="$W/base/$t.txt"
  s1=$(cmp -s "$a" "$b" && echo same || echo DIFF)
  s2=$(cmp -s "$c" "$a" && echo same || echo DIFF)
  s3=$(cmp -s "$c" "$b" && echo same || echo DIFF)
  printf "  %-6s off vs on: %-4s   base vs off: %-4s   base vs on: %-4s\n" "$t" "$s1" "$s2" "$s3"
  [ "$s1$s2$s3" = "samesamesame" ] || ok=0
done
[ "$ok" = 1 ] && echo "TAIL GATE PASS" || echo "TAIL GATE FAIL"
