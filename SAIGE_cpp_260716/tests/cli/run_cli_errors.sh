#!/bin/bash
# run_cli_errors.sh WORKDIR -- the front end's error messages and flag parsing.
# Every case runs bin/saige-gpu-cpp (with --dryRun where a run would start), checks
# the exit status and greps the message. Needs no data: the step-1 directory it
# uses is a stub with empty files (only their existence is checked).
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
SAIGE=${SAIGE:-$(cd "$HERE/../.." && pwd)/bin/saige-gpu-cpp}
W=$(mkdir -p "${1:?usage: run_cli_errors.sh WORKDIR}" && cd "$1" && pwd)
cd "$W"
rm -rf s1 s1empty out o o1 o2 o3 o4 rg rg2
mkdir -p s1/models/b1 s1/models/b2 s1empty
: > s1/models/b1/nullmodel.json; : > s1/models/b2/nullmodel.json
: > s1/vr_b1.varianceRatio.txt; : > s1/vr_b2.varianceRatio.txt

npass=0; nfail=0
expect() {  # expect <exit status> <regex> -- <args...>
  local want=$1 re=$2; shift 3
  local out rc
  out=$("$SAIGE" "$@" 2>&1); rc=$?
  if [ "$rc" = "$want" ] && grep -qE -- "$re" <<<"$out"; then
    npass=$((npass+1)); printf '  ok    %s\n        -> %s\n' "$*" "$(grep -E -m1 -- "$re" <<<"$out")"
  else
    nfail=$((nfail+1)); printf '  FAIL  %s (exit %s, wanted %s /%s/)\n%s\n' "$*" "$rc" "$want" "$re" "$out"
  fi
}
G="--plinkFile g --phenoFile p.txt"
echo "== unknown and refused flags"
expect 2 "unknown flag --phenoColl .*did you mean --phenoCol" -- step1 $G --phenoColl b1 --outDir o
expect 2 "unknown flag --useGpu .*did you mean --useGPU" -- step2 --useGpu
expect 2 "unknown flag --frobnicate for saige-gpu-cpp step2" -- step2 --frobnicate=1
expect 2 "--memoryChunk is an R SAIGE flag that .* step1 does not support" -- step1 $G --phenoCol b1 --outDir o --memoryChunk 2
expect 2 "--skipModelFitting is an R SAIGE flag" -- step1 $G --phenoCol b1 --outDir o --skipModelFitting=TRUE
expect 2 "--SPAcutoff is an R SAIGE flag .*give it to .saige-gpu-cpp step1. --SPAcutoff" -- step2 --step1Dir s1 --plinkFile g --outDir o --SPAcutoff 2
expect 2 "--is_fastTest is an R SAIGE flag .*step1. --is_fastTest" -- step2 --is_fastTest=TRUE
expect 2 "--idstoIncludeFile is an R SAIGE flag .*plink2 --extract" -- step2 --idstoIncludeFile ids.txt
expect 2 "--savFile is an R SAIGE flag" -- step2 --savFile x.sav
expect 2 "unexpected argument 'b1'" -- step1 --phenoCol b1 b1
expect 2 "unknown command 'step3'" -- step3
echo "== R defaults of refused flags are accepted"
expect 0 "config written" -- step1 $G --phenoCol b1 --outDir o --tauInit 0,0 --skipModelFitting=FALSE --dryRun
expect 0 "config written" -- step2 --step1Dir s1 --plinkFile g --outDir out --maxMAC_in_groupTest 0 --is_no_weight_in_groupTest=FALSE --dryRun
echo "== values"
expect 2 "--LOCO expects TRUE or FALSE, got 'maybe'" -- step1 $G --phenoCol b1 --outDir o --LOCO=maybe
expect 2 "--tol expects a number, got 'abc'" -- step1 $G --phenoCol b1 --outDir o --tol abc
expect 2 "--maxiter expects an integer, got '2.5'" -- step1 $G --phenoCol b1 --outDir o --maxiter 2.5
expect 2 "--outDir needs a value" -- step1 $G --phenoCol b1 --outDir
expect 2 "--traitType must be binary, quantitative or survival" -- step1 $G --phenoCol b1 --outDir o --traitType qt
expect 2 "--phenoCol: empty column name" -- step1 $G --phenoCol b1,,b2 --outDir o
echo "== missing required flags"
expect 2 "missing --outDir" -- step1 $G --phenoCol b1
expect 2 "missing --phenoCol" -- step1 $G --outDir o
expect 2 "missing --plinkFile" -- step1 --phenoFile p.txt --phenoCol b1 --outDir o
expect 2 "missing --phenoFile" -- step1 --plinkFile g --phenoCol b1 --outDir o
expect 2 "--makeSparseGRMOnly needs --sparseGRMFile" -- step1 $G --phenoCol b1 --outDir o --makeSparseGRMOnly
expect 2 "missing genotypes" -- step2 --step1Dir s1 --outDir out
expect 2 "missing models: --step1Dir DIR" -- step2 --plinkFile g --outDir out
expect 2 "--step1Dir needs --outDir" -- step2 --step1Dir s1 --plinkFile g
expect 2 "--GMMATmodelFile needs --varianceRatioFile" -- step2 --plinkFile g --GMMATmodelFile s1/models/b1 --outDir out
expect 2 "--bgenFile needs --sampleFile" -- step2 --step1Dir s1 --outDir out --bgenFile g.bgen
echo "== --step1Dir"
expect 2 "--step1Dir .*s1empty: no step-1 models found" -- step2 --step1Dir s1empty --plinkFile g --outDir out
expect 2 "--step1Dir: not a directory" -- step2 --step1Dir nosuchdir --plinkFile g --outDir out
expect 2 "--phenoCol b9: no such trait in .*s1 \(it has: b1, b2\)" -- step2 --step1Dir s1 --phenoCol b9 --plinkFile g --outDir out
expect 2 "--phenoCol in step 2 picks traits out of --step1Dir" -- step2 --phenoCol b1 --plinkFile g --GMMATmodelFile s1/models/b1 --varianceRatioFile s1/vr_b1.varianceRatio.txt --outDir out
expect 0 "config written" -- step2 --step1Dir s1 --phenoCol b2 --plinkFile g --outDir out --dryRun
echo "== genotype inputs"
expect 2 "give the genotypes in one format" -- step2 --step1Dir s1 --outDir out --plinkFile g --vcfFile g.vcf.gz
expect 2 "must share one prefix" -- step2 --step1Dir s1 --outDir out --bedFile a.bed --bimFile b.bim --famFile a.fam
expect 2 "reads the index from <bgenFile>.bgi" -- step2 --step1Dir s1 --outDir out --bgenFile g.bgen --sampleFile g.sample --bgenFileIndex other.bgi
expect 0 "--vcfFileIndex is not used" -- step2 --step1Dir s1 --outDir out --vcfFile g.vcf.gz --vcfFileIndex g.vcf.gz.csi --dryRun
echo "== output"
mkdir -p out && : > out/b1.txt
expect 2 "--is_overwrite_output=FALSE and .*out/b1.txt exists" -- step2 --step1Dir s1 --plinkFile g --outDir out --is_overwrite_output=FALSE
expect 2 "give either --outDir or --outputPrefix" -- step1 $G --phenoCol b1 --outDir o --outputPrefix p
echo "== logical flags: bare, =TRUE, separate TRUE/FALSE, R's T/F"
expect 0 "config written" -- step1 $G --phenoCol b1 --outDir o1 --LOCO FALSE --useGPU --invNormalize=T --isCovariateOffset F --dryRun
grep -q "loco: false" o1/step1.yaml && grep -q "use_gpu: true" o1/step1.yaml && grep -q "inv_normalize: true" o1/step1.yaml \
  && grep -q "covariate_offset: false" o1/step1.yaml && { npass=$((npass+1)); echo "  ok    o1/step1.yaml has loco false, use_gpu true, inv_normalize true, covariate_offset false"; } \
  || { nfail=$((nfail+1)); echo "  FAIL  o1/step1.yaml:"; cat o1/step1.yaml; }
echo "== --config + flags (flags win)"
cat > base.yaml <<YAML
paths:
  plinkFile: g
design:
  csv: p.txt
  y_col: b1
fit:
  trait: quantitative
  nthreads: 2
YAML
expect 0 "config written" -- step1 --config base.yaml --nThreads 4 --outDir o2 --dryRun
grep -q "nthreads: 4" o2/step1.yaml && grep -q "trait: quantitative" o2/step1.yaml && grep -q "plinkFile: $W/g" o2/step1.yaml \
  && { npass=$((npass+1)); echo "  ok    o2/step1.yaml: nthreads 4 (flag), trait quantitative (config), paths absolute"; } \
  || { nfail=$((nfail+1)); echo "  FAIL  o2/step1.yaml:"; cat o2/step1.yaml; }
echo "== the written config holds only what was set (no defaults filled in)"
"$SAIGE" step1 --plinkFile g --phenoFile p.txt --phenoCol b1 --outDir o3 --dryRun > /dev/null
k=$(awk '!/^#/' o3/step1.yaml | grep -E '^ *[A-Za-z_]+:' | sed -E 's/:.*//; s/^ +/  /' | tr '\n' ' ')
want="paths   out_prefix   out_prefix_vr   plinkFile design   csv   y_cols "
[ "$k" = "$want" ] && { npass=$((npass+1)); echo "  ok    step1 minimal: only $k"; } \
  || { nfail=$((nfail+1)); echo "  FAIL  step1 minimal wrote: $k"; cat o3/step1.yaml; }
"$SAIGE" step2 --step1Dir s1 --phenoCol b1 --plinkFile g --outDir o4 --groupFile grp.txt --dryRun > /dev/null
k=$(awk '!/^#/' o4/step2.yaml | grep -E '^[A-Za-z_]+:' | sed -E 's/:.*//' | tr '\n' ' ')
want="groupFile genoType plinkFile models "
[ "$k" = "$want" ] && { npass=$((npass+1)); echo "  ok    step2 minimal (+groupFile): only $k"; } \
  || { nfail=$((nfail+1)); echo "  FAIL  step2 minimal wrote: $k"; cat o4/step2.yaml; }
echo
echo "$npass passed, $nfail failed"
[ "$nfail" = 0 ]
