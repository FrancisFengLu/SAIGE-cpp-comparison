#!/bin/bash
# arch_compare.sh WORKDIR BIN_A BIN_B -- two builds of the engines (e.g. make ARCH=native
# vs make ARCH=x86-64-v3 BINDIR=bin-v3) on the same runs: wall time of each step and
# whether the outputs are byte-identical. Single runs, not timing-grade.
#
# Data: the tutorial data (docs/examples/02_simulate.sh, 5,000 x 5,000) and a
# simulated N = 50,000 x 20,000 cohort (plink2 --dummy, two chromosomes, the
# tutorial's phenotype recipe). Runs, each with BIN_A then BIN_B:
#   tut/gpu   step 1 GPU 8 threads, step 2 GPU 8 threads
#   tut/cpu   step 1 CPU 1 thread (bit-reproducible), step 2 CPU 8 threads
#   n50k/gpu  step 1 GPU 8 threads, step 2 GPU 8 threads
#   n50k/cpu  step 2 CPU 8 threads on the n50k/gpu models of the same build
# The front end is BIN_A/saige-gpu-cpp when present, else bin/saige-gpu-cpp; it
# runs the engines next to it, so each build's engines are linked into a
# private directory with a copy of the front end.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
export SAIGE_HOME=$(cd "$HERE/../.." && pwd)
W=$(mkdir -p "${1:?usage: arch_compare.sh WORKDIR BIN_A BIN_B}" && cd "$1" && pwd)
A=$(cd "${2:?}" && pwd); B=$(cd "${3:?}" && pwd)
FRONT=$SAIGE_HOME/bin/saige-gpu-cpp
PLINK2=${PLINK2:-plink2}

# data
if [ ! -f "$W/tut/data/geno.bed" ]; then
  mkdir -p "$W/tut"; ( cd "$W/tut" && WORK=$W/tut bash "$SAIGE_HOME/docs/examples/02_simulate.sh" > /dev/null )
fi
if [ ! -f "$W/n50k/data/geno.bed" ]; then
  mkdir -p "$W/n50k/data"; cd "$W/n50k/data"
  $PLINK2 --dummy 50000 20000 0.01 acgt --seed 3 --make-bed --out geno0 > /dev/null
  awk 'BEGIN{OFS="\t"} {if (NR > 10000) $1 = 2; $2 = "snp" NR; $4 = 1000 * NR; print}' geno0.bim > geno.bim
  mv geno0.bed geno.bed; mv geno0.fam geno.fam; rm -f geno0.*
  awk 'BEGIN{srand(11); OFS="\t"; print "IID","b1","b2","b3","b4","x1","x2"}
       function gauss(){ return sqrt(-2*log(1-rand()))*cos(6.283185307*rand()) }
       { x1 = gauss(); x2 = (rand() < 0.5) ? 1 : 0
         printf "%s\t%d\t%d\t%d\t%d\t%.4f\t%d\n", $2, (rand() < 0.10 + 0.03*x1), (rand() < 0.30),
                (rand() < 0.05), (rand() < 0.20), x1, x2 }' geno.fam > pheno.txt
fi

declare -A T
for L in A B; do
  D=$W/bin_$L; BIN=$([ $L = A ] && echo "$A" || echo "$B")
  rm -rf "$D"; mkdir -p "$D"
  for e in saige-null saige-step2 sgs2txt; do ln -s "$BIN/$e" "$D/$e"; done
  cp "$FRONT" "$D/saige-gpu-cpp"
  S=$D/saige-gpu-cpp
  timed() {  # timed <key> <cmd...>
    local k=$1; shift; local t0=$(date +%s.%N)
    "$@" > "$W/run_$L.$(echo "$k" | tr / _).log" 2>&1 || { echo "FAILED $k ($L)"; exit 1; }
    T[$L,$k]=$(echo "$(date +%s.%N) - $t0" | bc)
  }
  for c in tut n50k; do
    cd "$W/$c"
    timed "$c/gpu/step1" $S step1 --plinkFile data/geno --phenoFile data/pheno.txt --phenoCol b1,b2,b3,b4 \
        --covarColList x1,x2 --LOCO=FALSE --nThreads 8 --useGPU --IsOverwriteVarianceRatioFile --outDir $L/gpu/step1
    timed "$c/gpu/step2" $S step2 --step1Dir $L/gpu/step1 --plinkFile data/geno --minMAC 1 --nThreads 8 \
        --useGPU --outDir $L/gpu/step2
    if [ $c = tut ]; then
      timed "$c/cpu/step1" $S step1 --plinkFile data/geno --phenoFile data/pheno.txt --phenoCol b1,b2,b3,b4 \
          --covarColList x1,x2 --LOCO=FALSE --nThreads 1 --IsOverwriteVarianceRatioFile --outDir $L/cpu/step1
      M=$L/cpu/step1
    else
      M=$L/gpu/step1
    fi
    timed "$c/cpu/step2" $S step2 --step1Dir $M --plinkFile data/geno --minMAC 1 --nThreads 8 --outDir $L/cpu/step2
  done
done

printf '\n%-16s %10s %10s   %s\n' "run" "A (s)" "B (s)" "outputs A vs B"
for k in tut/gpu/step1 tut/gpu/step2 tut/cpu/step1 tut/cpu/step2 n50k/gpu/step1 n50k/gpu/step2 n50k/cpu/step2; do
  c=${k%%/*}; r=${k#*/}; d=$W/$c/A/$r; e=$W/$c/B/$r
  n=0; bad=0
  while IFS= read -r f; do n=$((n+1)); cmp -s "$d/$f" "$e/$f" || bad=$((bad+1)); done \
    < <(cd "$d" && find . -type f ! -name '*.yaml' ! -name '*.log' | sort)
  res=$([ $bad = 0 ] && echo "$n files identical" || echo "$bad of $n files differ")
  printf '%-16s %10.1f %10.1f   %s\n' "$k" "${T[A,$k]}" "${T[B,$k]}" "$res"
done
