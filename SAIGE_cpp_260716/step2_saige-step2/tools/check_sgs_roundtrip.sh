#!/bin/bash
# The fp32 round-trip check, done the only way it can be done: `cmp` says
# "differs" and that is all it can say, so this counts fields instead.
#   - text run vs sgs/fp64 -> sgs2txt : cmp, must be byte-identical (unchanged)
#   - text run vs sgs/fp32 -> sgs2txt : fieldcmp, per-column mismatch counts
#   - and the same counts predicted by sgsfid from the fp64 file alone, so the
#     P=128 numbers (which would need 24 GB of text to check this way) can be
#     trusted.
# 10^6 markers x 50,000 samples, P=8.
set -u
W=/opt/saige/logs/gpu_step2/writer/f32/cross
S2=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2
BIN=$S2/saige-step2.cuda
CVT=$S2/tools/sgs2txt
FID=$S2/tools/sgs_fidelity
FC=$S2/tools/out_fieldcmp
M=/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out
BED=/opt/saige/logs/tg2_step2/data/g1m
P=8
rm -rf $W; mkdir -p $W/txt $W/f64 $W/f32
gen(){ local OD=$1; shift
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  for k in $(seq 1 $P); do echo "  - traitName: y$k"; echo "    modelFile: $M/m/y$k"
    echo "    varianceRatioFile: $M/mvr_y$k.varianceRatio.txt"; echo "    outputFile: $OD/y$k.txt"; done; }
gen $W/txt "useGPU: true"                                          > $W/txt.yaml
gen $W/f64 "useGPU: true" "outputFormat: sgs"                      > $W/f64.yaml
gen $W/f32 "useGPU: true" "outputFormat: sgs" "sgsPrecision: fp32" > $W/f32.yaml
source /opt/saige/SAIGE-work/optimization/torchgwas2/scripts/env_cpp.sh
for t in txt f64 f32; do
  echo "### run $t $(date +%T)"
  $BIN $W/$t.yaml > $W/$t.log 2>&1 || { echo "$t RUN FAILED"; tail -5 $W/$t.log; exit 1; }
done
$CVT -j 8 $W/f64/*.sgs >/dev/null 2>$W/cvt64.log
$CVT -j 8 $W/f32/*.sgs >/dev/null 2>$W/cvt32.log
bad=0; for k in $(seq 1 $P); do cmp -s $W/txt/y$k.txt $W/f64/y$k.txt || { echo "fp64 DIFF y$k"; bad=1; }; done
rows=$(( $(wc -l < $W/txt/y1.txt) - 1 ))
[ $bad -eq 0 ] && echo "fp64 ROUND TRIP OK: $P x $rows = $((P*rows)) rows byte-identical" \
               || echo "fp64 ROUND TRIP FAILED"
echo
echo "=== fp32 measured: text run vs sgs/fp32 -> sgs2txt, all $P traits ==="
for k in $(seq 1 $P); do
  cmp -s $W/txt/y$k.txt $W/f32/y$k.txt && echo "y$k: identical (unexpected)" || true
  echo "--- y$k ---"; $FC $W/txt/y$k.txt $W/f32/y$k.txt
done > $W/fieldcmp.out 2>&1
grep -E "^rows|^total differing" $W/fieldcmp.out | paste - - | head -$P
echo
echo "per-column totals over the $P traits (measured):"
awk '/^(AC_Allele2|AF_Allele2|MissingRate|BETA|SE|Tstat|var|p\.value|N|CHR|POS|MarkerID|Allele1|Allele2) /{c[$1]+=$2} END{for(k in c) printf "  %-12s %12d\n",k,c[k]}' $W/fieldcmp.out | sort
echo
echo "=== the same thing predicted by sgsfid from the fp64 file alone ==="
$FID -j 8 -m $W/f64/y1.txt.markers.sgs $W/f64/y*.txt.sgs
echo
echo "=== worst p-value rows, all traits ==="
grep -h "^p.value " $W/fieldcmp.out
