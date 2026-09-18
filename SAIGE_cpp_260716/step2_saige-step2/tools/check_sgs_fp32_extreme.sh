#!/bin/bash
# The case the benchmark data cannot show: p-values small enough that a float
# cannot hold them.
#
# The 10^6-marker set is a null simulation, so its smallest p is around 1e-9 and
# every p sits comfortably inside float's normal range. A real GWAS does not:
# p = 1e-50, 1e-120, 1e-300 are ordinary for a strong locus, the score test
# prints all of them with "%.6E" (it only leaves that form below ~1e-308), and
# float's smallest normal is 1.2e-38 and its smallest subnormal 1.4e-45. So
# under sgsPrecision: fp32 every p below ~1.2e-38 loses digits and every p below
# ~1.4e-45 becomes exactly 0.000000E+00.
#
# tests/make_extreme_model.py scales a null model's residuals by k, which
# multiplies the chi-square statistic by k^2 and drags every p-value down. k=6
# and k=15 put the 5,000 markers across 1 .. 1e-300.
set -u
W=/opt/saige/logs/gpu_step2/writer/f32/extreme
S2=/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2
BIN=$S2/saige-step2.cuda
CVT=$S2/tools/sgs2txt
FC=$S2/tools/out_fieldcmp
M=/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out
BED=/opt/saige/logs/tg2_step2/data/g5k
rm -rf $W; mkdir -p $W/m $W/txt $W/f32
source /opt/saige/SAIGE-work/optimization/torchgwas2/scripts/env_cpp.sh
cd $S2/tests
python3 make_extreme_model.py $M/m/y1 $W/m/e6  6  || exit 1
python3 make_extreme_model.py $M/m/y2 $W/m/e15 15 || exit 1
cd $W
gen(){ local OD=$1; shift
  echo "genoType: plink"; echo "plinkFile: $BED"
  echo "minMAF: 0"; echo "minMAC: 1"; echo "maxMissRate: 0.15"
  echo "AlleleOrder: alt-first"; echo "LOCO: false"; echo "isnoadjCov: false"
  echo "isMoreOutput: false"; echo "isFirth: false"; echo "is_Firth_beta: false"
  echo "MACCutoffforER: 4"; echo "relatednessCutoff: 0"; echo "nThreads: 8"
  for L in "$@"; do echo "$L"; done
  echo "models:"
  echo "  - traitName: e6";  echo "    modelFile: $W/m/e6"
  echo "    varianceRatioFile: $M/mvr_y1.varianceRatio.txt"; echo "    outputFile: $OD/e6.txt"
  echo "  - traitName: e15"; echo "    modelFile: $W/m/e15"
  echo "    varianceRatioFile: $M/mvr_y2.varianceRatio.txt"; echo "    outputFile: $OD/e15.txt"; }
gen $W/txt "useGPU: true"                                          > $W/txt.yaml
gen $W/f32 "useGPU: true" "outputFormat: sgs" "sgsPrecision: fp32" > $W/f32.yaml
$BIN $W/txt.yaml > $W/txt.log 2>&1 || { echo "text run FAILED"; tail -5 $W/txt.log; exit 1; }
$BIN $W/f32.yaml > $W/f32.log 2>&1 || { echo "fp32 run FAILED"; tail -5 $W/f32.log; exit 1; }
$CVT -j 2 $W/f32/*.sgs >/dev/null 2>$W/cvt.log
for t in e6 e15; do
  echo "######## $t"
  echo "p-value range in the text run: $(awk -F'\t' 'NR>1{print $13}' $W/txt/$t.txt | sort -g | tail -1) .. $(awk -F'\t' 'NR>1{print $13}' $W/txt/$t.txt | sort -g | head -1)"
  echo "rows with p < 1.2e-38 (float's smallest normal): $(awk -F'\t' 'NR>1 && $13+0 < 1.2e-38 && $13 !~ /E[0-9]/ {c++} END{print c+0}' $W/txt/$t.txt)"
  echo "rows the fp32 file turned into 0.000000E+00: $(awk -F'\t' 'NR>1 && $13=="0.000000E+00"{c++} END{print c+0}' $W/f32/$t.txt)"
  echo "  (the text run wrote 0.000000E+00 on $(awk -F'\t' 'NR>1 && $13=="0.000000E+00"{c++} END{print c+0}' $W/txt/$t.txt) rows)"
  echo "rows kept verbatim as the %.1fE%d underflow form: $(awk -F'\t' 'NR>1 && $13 ~ /^[0-9]\.[0-9]E-/{c++} END{print c+0}' $W/txt/$t.txt)"
  $FC $W/txt/$t.txt $W/f32/$t.txt
  echo "--- five worst p-value rows (text -> fp32) ---"
  paste <(cut -f3,13 $W/txt/$t.txt) <(cut -f13 $W/f32/$t.txt) | tail -n +2 |
    awk '$2!=$3' | sort -t$'\t' -k2,2g | head -5
  echo
done
