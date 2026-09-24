#!/bin/bash
# mtvec_gates.sh -- acceptance for mtVecQuantStats.
#
#   gate 1  noise: switch OFF, same binary, twice
#   gate 2  switch OFF vs ON, quantitative, fp64 (.sgs) and text
#   gate 3  a mixed binary + quantitative run: every binary trait byte for byte
#           what the pre-change binary wrote, switch either way
#   gate 4  the fp64 tail accuracy on THIS run's own stat values
#
# env: DATA (plink prefix), NQ, NB, BIN (new), BASEBIN (pre-change), RUNROOT
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DATA=${DATA:-/opt/saige/data/mid}
NQ=${NQ:-8}
NB=${NB:-8}
BIN=${BIN:-/opt/saige/logs/mtvec/bin/saige-step2}
BASEBIN=${BASEBIN:-/opt/saige/logs/mtvec/bin/saige-step2.base}
export RUNROOT=${RUNROOT:-/opt/saige/logs/mtvec/runs}
R=$RUNROOT

run() {  # run TAG VEC SPEC...   (env SGS, BIN already exported)
  local tag=$1 vec=$2; shift 2
  BIN=$BIN bash "$HERE/mtvec_run_one.sh" "$tag" "$DATA" "$vec" "$@" > "$R/$tag.run.log" 2>&1 \
    || { echo "RUN FAILED: $tag"; tail -5 "$R/$tag.run.log"; exit 1; }
}

cmpdir() {  # cmpdir A B PREFIX  -> count of byte-identical / total
  local a=$1 b=$2 p=$3 same=0 tot=0
  for f in "$a"/$p*.txt; do
    [ -e "$f" ] || continue
    tot=$((tot+1))
    cmp -s "$f" "$b/$(basename "$f")" && same=$((same+1))
  done
  echo "$same/$tot"
}

echo "===== data $DATA   NQ=$NQ NB=$NB ====="
echo "new  $(md5sum "$BIN" | cut -d' ' -f1)"
echo "base $(md5sum "$BASEBIN" | cut -d' ' -f1)"

# ---------------- gate 1: noise ----------------
echo
echo "----- gate 1: switch OFF twice (noise floor) -----"
SGS=1 run g1_off_a 0 q:$NQ
SGS=1 run g1_off_b 0 q:$NQ
python3 "$HERE/mtvec_cmp_sgs.py" "$R/g1_off_a/out" "$R/g1_off_b/out" q_
run g1t_off_a 0 q:$NQ
run g1t_off_b 0 q:$NQ
echo "  text byte-identical: $(cmpdir "$R/g1t_off_a/out" "$R/g1t_off_b/out" q_)"

# ---------------- gate 2: OFF vs ON ----------------
echo
echo "----- gate 2: switch OFF vs ON, quantitative -----"
SGS=1 run g2_on 1 q:$NQ
python3 "$HERE/mtvec_cmp_sgs.py" "$R/g1_off_a/out" "$R/g2_on/out" q_ \
        --dump-stat "$R/g2_stat.f64"
run g2t_on 1 q:$NQ
echo "  text byte-identical: $(cmpdir "$R/g1t_off_a/out" "$R/g2t_on/out" q_)"
if [ -x "$HERE/out_fieldcmp" ]; then
  for f in "$R/g1t_off_a/out"/q_*.txt; do
    "$HERE/out_fieldcmp" "$f" "$R/g2t_on/out/$(basename "$f")"
  done | awk '{print "  "$0}'
fi

# ---------------- gate 3: binary untouched ----------------
echo
echo "----- gate 3: mixed run, binary traits -----"
BIN=$BASEBIN run g3_base 0 q:$NQ b:$NB
run g3_off 0 q:$NQ b:$NB
run g3_on  1 q:$NQ b:$NB
echo "  base vs new(off)   binary $(cmpdir "$R/g3_base/out" "$R/g3_off/out" b_)   quant $(cmpdir "$R/g3_base/out" "$R/g3_off/out" q_)"
echo "  base vs new(on)    binary $(cmpdir "$R/g3_base/out" "$R/g3_on/out"  b_)   quant $(cmpdir "$R/g3_base/out" "$R/g3_on/out"  q_)"
echo "  new(off) vs new(on) binary $(cmpdir "$R/g3_off/out" "$R/g3_on/out" b_)   quant $(cmpdir "$R/g3_off/out" "$R/g3_on/out" q_)"
echo "  quantitative results independent of the batch's binary company:"
echo "    mixed(on) vs quant-only(on) $(cmpdir "$R/g3_on/out" "$R/g2t_on/out" q_)"

# ---------------- gate 4: fp64 tail accuracy on this run ----------------
echo
echo "----- gate 4: fp64 tail accuracy on this run's stat values -----"
if [ -x "$HERE/mtvec_tail_accuracy" ] && [ -f "$R/g2_stat.f64" ]; then
  "$HERE/mtvec_tail_accuracy" "$R/g2_stat.f64" | sed -n '/run stat values/,/^$/p'
else
  echo "  (build tools/mtvec_tail_accuracy first)"
fi
