#!/bin/bash
# run_cli_identity.sh WORKDIR [REF_COMMIT]
#
# The command-line examples (docs/examples/03..11, saige-gpu-cpp flags) against
# the YAML examples they replaced (docs/examples of REF_COMMIT, default
# ef674945): both sets run on the same simulated data with the same engines
# (bin/saige-null, bin/saige-step2, bin/sgs2txt), then every result file is
# compared with cmp. Prints one line per compared file group and exits non-zero
# on any difference.
#
# Needs: `make` (bin/) done, plink2 (PLINK2=...), the build environment of
# docs/examples/env.sh (CONDA_ENV=...).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
export SAIGE_HOME=$(cd "$HERE/../.." && pwd)
W=$(mkdir -p "${1:?usage: run_cli_identity.sh WORKDIR [REF_COMMIT]}" && cd "$1" && pwd)
REF=${2:-ef674945}
export S1=$SAIGE_HOME/bin/saige-null S2=$SAIGE_HOME/bin/saige-step2 SGS2TXT=$SAIGE_HOME/bin/sgs2txt
export SAIGE=$SAIGE_HOME/bin/saige-gpu-cpp
for b in "$S1" "$S2" "$SGS2TXT" "$SAIGE"; do [ -x "$b" ] || { echo "missing $b (run make)"; exit 2; }; done

# the YAML examples as they were at REF
rm -rf "$W/yaml_examples"; mkdir -p "$W/yaml_examples" "$W/yaml" "$W/cli"
TOP=$(git -C "$SAIGE_HOME" rev-parse --show-toplevel)
PREFIX=$(git -C "$SAIGE_HOME" rev-parse --show-prefix)
git -C "$TOP" archive "$REF:${PREFIX}docs/examples" | tar -x -C "$W/yaml_examples"

run() {  # run <dir with the scripts> <WORK> <script> [args]
  local D=$1 WK=$2 S=$3; shift 3
  echo "  $(basename "$D")/$S $*"
  ( cd "$WK" && WORK=$WK bash "$D/$S" "$@" ) > "$WK/$S.out" 2>&1 || { echo "FAILED: $D/$S (see $WK/$S.out)"; exit 1; }
}

echo "== YAML examples ($REF)"
Y=$W/yaml_examples
run "$Y" "$W/yaml" 02_simulate.sh
for s in 03_step1_binary.sh 04_step2_binary.sh 05_step1_quant_loco.sh 06_step2_quant_loco.sh \
         07_sparse_grm.sh 08_step2_sgs.sh 09_step2_formats.sh; do run "$Y" "$W/yaml" $s; done
run "$Y" "$W/yaml" 10_hpc_chrom_job.sh 1
run "$Y" "$W/yaml" 10_hpc_chrom_job.sh 2
run "$Y" "$W/yaml" 11_hpc_concat.sh

echo "== command-line examples (this tree)"
C=$SAIGE_HOME/docs/examples
ln -sfn "$W/yaml/data" "$W/cli/data"       # the same simulated data
for s in 03_step1_binary.sh 04_step2_binary.sh 05_step1_quant_loco.sh 06_step2_quant_loco.sh \
         07_sparse_grm.sh 08_step2_sgs.sh 09_step2_formats.sh; do run "$C" "$W/cli" $s; done
run "$C" "$W/cli" 10_hpc_chrom_job.sh 1
run "$C" "$W/cli" 10_hpc_chrom_job.sh 2
run "$C" "$W/cli" 11_hpc_concat.sh

echo "== cmp (YAML run vs command-line run)"
nfail=0
cmpgroup() {  # cmpgroup <label> <yaml dir> <cli dir> [find args]: every file of the yaml dir
  local L=$1 A=$W/yaml/$2 B=$W/cli/$3; shift 3
  local n=0 bad=0 f
  while IFS= read -r f; do
    n=$((n+1))
    if [ ! -e "$B/$f" ]; then echo "    missing: $B/$f"; bad=$((bad+1));
    elif ! cmp -s "$A/$f" "$B/$f"; then echo "    differs: $f"; bad=$((bad+1)); fi
  done < <(cd "$A" && find . -type f ! -name '*.log' ! -name '*.yaml' ! -name '*.out' "$@" | sort)
  if [ "$n" = 0 ]; then echo "  $L: no files"; nfail=$((nfail+1));
  elif [ "$bad" = 0 ]; then echo "  $L: $n files identical";
  else echo "  $L: $bad of $n files differ"; nfail=$((nfail+1)); fi
}
cmpgroup "03 step 1, binary (models + VR)"       step1_bin           step1_bin
cmpgroup "04 step 2, binary"                     step2_bin/out       step2_bin
cmpgroup "05 step 1, quantitative LOCO"          step1_qt            step1_qt
cmpgroup "06 step 2, LOCO chr1"                  step2_qt/chr1       step2_qt/chr1
cmpgroup "06 step 2, LOCO chr2"                  step2_qt/chr2       step2_qt/chr2
cmpgroup "07 sparse GRM (grm.mtx, grm.ids)"      sparse              sparse          -maxdepth 1 -name 'grm.*'
cmpgroup "07 step 1 on the sparse GRM (models)"  sparse/models       sparse/step1/models
cmpgroup "07 step 1 on the sparse GRM (VR)"      sparse              sparse/step1    -maxdepth 1 -name 'vr_*'
cmpgroup "07 step 2 on the sparse GRM"           sparse/out          sparse/step2
cmpgroup "08 sgs -> text"                        step2_sgs/text      step2_sgs/text
for f in pgen pgen_dosage bgen; do
  cmpgroup "09 step 2, $f"                       formats/$f          formats/$f
done
# VCF with nThreads > 1 writes the rows in completion order (docs/step2.md), so the
# 8-thread files are compared as sorted rows, and a 1-thread run of each must be
# byte-identical.
n=0; bad=0
for t in b1 b2 b3; do
  n=$((n+1))
  cmp -s <(sort "$W/yaml/formats/vcf/$t.txt") <(sort "$W/cli/formats/vcf/$t.txt") || bad=$((bad+1))
done
[ $bad = 0 ] && echo "  09 step 2, vcf (8 threads, rows sorted): $n files identical" \
             || { echo "  09 step 2, vcf (8 threads, rows sorted): $bad of $n differ"; nfail=$((nfail+1)); }
mkdir -p "$W/yaml/formats/vcf1" "$W/cli/formats/vcf1"
sed -e 's/^nThreads: .*/nThreads: 1/' -e "s#/formats/vcf/#/formats/vcf1/#" \
    "$W/yaml/formats/vcf/step2.yaml" > "$W/yaml/formats/vcf1/step2.yaml"
"$S2" "$W/yaml/formats/vcf1/step2.yaml" > "$W/yaml/formats/vcf1/step2.log" 2>&1
( cd "$W/cli" && "$SAIGE" step2 --vcfFile data/geno.vcf.gz --vcfField DS --step1Dir step1_bin \
    --phenoCol b1,b2,b3 --minMAC 1 --nThreads 1 --useGPU --outDir formats/vcf1 > formats_vcf1.log 2>&1 )
cmpgroup "09 step 2, vcf (1 thread)"            formats/vcf1        formats/vcf1
cmpgroup "10/11 per-chromosome jobs + concat"    hpc                 hpc             -name '*.txt'
# the sgs text also equals the text writer's output
for t in b1 b2 b3 b4; do
  cmp -s "$W/cli/step2_sgs/text/$t.txt" "$W/cli/step2_bin/$t.txt" || { echo "  08: $t sgs text != text output"; nfail=$((nfail+1)); }
done
if [ "$nfail" = 0 ]; then echo "ALL IDENTICAL"; else echo "$nfail group(s) differ"; exit 1; fi
