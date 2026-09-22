#!/bin/bash
# Gate for the "sparse GRM in use but the variance-ratio file has no sparse row"
# refusal added in main.cpp (model construction).
#   A control : the VR file as step 1 wrote it -> step 2 runs, rc=0
#   B break   : same run with the 'sparse' row stripped -> must exit non-zero
#               and name the missing row. Before the patch this run succeeded
#               and handed -1.0 out as the variance ratio for that MAC bin.
set -u
G=/opt/saige/logs/vrgate
S2=${S2:-/opt/saige/SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/saige-step2}
SRC=/opt/saige/logs/comp8qs/runs/base_single_q1_r1
MODEL=$SRC/out/m/q1
VR=$SRC/out/mvr_q1.varianceRatio.txt
SGRM=/opt/saige/logs/missing_mt/step1/data/mid.sgrm2.mtx
SGRM_IDS=/opt/saige/logs/missing_mt/step1/data/mid.sgrm2.ids
rm -rf $G/A $G/B; mkdir -p $G/A/out $G/B/out
cp "$VR" $G/A/vr.txt
grep -v $'\tsparse\t' "$VR" > $G/B/vr.txt
echo "--- A has $(wc -l < $G/A/vr.txt) rows, B has $(wc -l < $G/B/vr.txt) (sparse row stripped)"
for d in $G/A $G/B; do cat > $d/cfg.yaml <<YAML
genoType: plink
plinkFile: /opt/saige/data/mid
minMAF: 0
minMAC: 1
maxMissRate: 0.15
AlleleOrder: alt-first
LOCO: false
isnoadjCov: false
isMoreOutput: false
isFirth: false
is_Firth_beta: false
MACCutoffforER: 4
relatednessCutoff: 0
nThreads: 4
sparseGRMFile: $SGRM
sparseGRMSampleIDFile: $SGRM_IDS
models:
  - traitName: q1
    modelFile: $MODEL
    varianceRatioFile: $d/vr.txt
    outputFile: $d/out/q1.txt
YAML
done
for d in A B; do
  ( cd $G/$d && timeout 900 $S2 cfg.yaml > log.txt 2>&1; echo $? > rc )
  echo "=== $d rc=$(cat $G/$d/rc)"; tail -2 $G/$d/log.txt | sed 's/^/    /'
done
a=$(cat $G/A/rc); b=$(cat $G/B/rc)
if [ "$a" = 0 ] && [ "$b" != 0 ] && grep -q "no 'sparse' row" $G/B/log.txt; then
  echo "GATE PASS: control rc=0, stripped file refused rc=$b with the new message"
else
  echo "GATE FAIL: control rc=$a, stripped rc=$b"; grep -i "sparse" $G/B/log.txt | tail -3; exit 1
fi
