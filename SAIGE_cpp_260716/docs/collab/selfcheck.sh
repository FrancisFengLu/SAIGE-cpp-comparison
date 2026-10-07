#!/bin/bash
# selfcheck.sh -- install self-check on simulated data (about 2 minutes, no real data).
#   1. simulate 2,000 samples x 3,000 markers (plink2 --dummy, fixed seed) and 4 binary traits
#   2. step 1: full-GRM null model (CPU), sparse GRM + sparse-GRM null model
#   3. step 2: CPU path and GPU path, Firth on, both models
#   4. checks: (a) the GPU path was used and its output is byte-identical to the CPU path
#              (b) md5s against the ones recorded when this protocol was written
# (a) must pass on any machine. (b) is expected to pass only with the same plink2 version (for the
# data) and the same build flags and CPU instruction set (the build uses -march=native; a CPU
# with different vector units, e.g. AVX-512 vs AVX2, can change last digits of the results).
# A (b) mismatch with (a) passing is not an error: send the result file, we compare.
# Environment: SELFCHECK_DIR (default $OUT_ROOT/selfcheck, or ./saige_selfcheck), NTHREADS.
set -uo pipefail
source "$(dirname "$0")/common.sh"
W=${SELFCHECK_DIR:-${OUT_ROOT:+$OUT_ROOT/selfcheck}}; W=${W:-$PWD/saige_selfcheck}
NTHREADS=${NTHREADS:-$NPROC}
mkdir -p "$W/data" "$W/s1" "$W/s2"; cd "$W"
RES=$W/selfcheck_result.txt; : > "$RES"
say() { echo "$*" | tee -a "$RES"; }

# Recorded on GCP n1-standard-8 (Intel Xeon 2.30 GHz, AVX2 + FMA, no AVX-512), Tesla V100, CUDA 12.9,
# plink2 v2.0.0-a.6.5LM (22 Dec 2024), source commit given in BUILD_INFO.txt of that run.
declare -A EXP=(
  [data/geno.bed]=e4df542264c9da921ba567eb5a1f222a
  [data/geno.bim]=25b7ffdd21722371648c178c0c4cc3a4
  [data/geno.fam]=da01ac9cced30a5662c271617625a604
  [data/pheno.txt]=94e389390318a76ba78549d27f6e49a4
  [full/b1.txt]=6653d123713351a76b3fc41456950483
  [full/b2.txt]=af3645de0220b40721494f8ac099188d
  [full/b3.txt]=587c813fb838e725aadc74df714ddfed
  [full/b4.txt]=58fd0310ce45647ae5f6be0c5b599909
  [sparse/b1.txt]=3bac2d845a7526f2523f467b16ada31d
  [sparse/b2.txt]=d6f0f8cda6b627308786fe27a4c3e659
  [sparse/b3.txt]=b2f4a18abfc80fe39ab08ebb99dc0bb5
  [sparse/b4.txt]=865dccbc46446ddee259611335adca78
)

# 1. data
"$PLINK2" --dummy 2000 3000 0.01 acgt --seed 11 --make-bed --out data/geno > data/plink2.log 2>&1 \
  || { say "FAIL plink2 --dummy (see $W/data/plink2.log)"; exit 1; }
awk 'BEGIN{OFS="\t"} {$2 = "m" NR; $4 = 1000 * NR; print}' data/geno.bim > data/g.bim && mv data/g.bim data/geno.bim
$PYTHON "$COLLAB/sim_cohort.py" pheno data/geno.fam data/pheno.txt --nbin 4 > /dev/null
say "plink2: $("$PLINK2" --version | head -1)"

# 2. step 1 on the CPU with one thread: the full-GRM fit with several CPU threads (or on the GPU)
#    sums in a varying order, which changes last digits of the model; step 2 does not depend on nThreads
cat > s1/full.yaml <<YAML
paths: {plinkFile: $W/data/geno, out_prefix: $W/s1/full/models, out_prefix_vr: $W/s1/full/vr, overwrite_varratio: true}
design: {csv: $W/data/pheno.txt, iid_col: IID, y_cols: [b1, b2, b3, b4], covar_cols: [x1, x2]}
fit: {trait: binary, loco: false, nthreads: 1, use_gpu: false, firth_beta: true, fast_test: true}
YAML
cat > s1/grm.yaml <<YAML
paths: {plinkFile: $W/data/geno, out_prefix: $W/s1/grm/run, sparse_grm: $W/s1/grm.mtx, sparse_grm_ids: $W/s1/grm.ids}
design: {csv: $W/data/pheno.txt, iid_col: IID, y_col: b1}
fit: {trait: binary, use_sparse_grm_to_fit: true, make_sparse_grm_only: true, relatedness_cutoff: 0.05, min_maf_grm: 0.01}
YAML
cat > s1/sparse.yaml <<YAML
paths: {plinkFile: $W/data/geno, out_prefix: $W/s1/sparse/models, out_prefix_vr: $W/s1/sparse/vr, sparse_grm: $W/s1/grm.mtx, sparse_grm_ids: $W/s1/grm.ids, overwrite_varratio: true}
design: {csv: $W/data/pheno.txt, iid_col: IID, y_cols: [b1, b2, b3, b4], covar_cols: [x1, x2]}
fit: {trait: binary, use_sparse_grm_to_fit: true, use_sparse_grm_for_vr: true, firth_beta: true, fast_test: true}
YAML
rm -f s1/grm.mtx s1/grm.ids
for c in full grm sparse; do
  "$NULLBIN" -c s1/$c.yaml > s1/$c.log 2>&1 || { say "FAIL step 1 ($c), see $W/s1/$c.log"; exit 1; }
done
say "step 1: $(grep -c '^Converged: yes' s1/full.log)/4 full-GRM and $(grep -c '^Converged: yes' s1/sparse.log)/4 sparse-GRM models converged"

# 3. step 2
for g in full sparse; do for p in cpu gpu; do
  d=s2/${g}_$p; mkdir -p "$d/out"
  { printf 'genoType: plink\nplinkFile: %s\nAlleleOrder: alt-first\nminMAF: 0\nminMAC: 1\n' "$W/data/geno"
    printf 'isFirth: true\nis_Firth_beta: true\npCutoffforFirth: 0.01\nnThreads: %s\nuseGPU: %s\nmodels:\n' "$NTHREADS" "$([ $p = gpu ] && echo true || echo false)"
    for t in b1 b2 b3 b4; do
      printf '  - {traitName: %s, modelFile: %s, varianceRatioFile: %s, outputFile: %s}\n' \
        $t "$W/s1/$g/models/$t" "$W/s1/$g/vr_$t.varianceRatio.txt" "$W/$d/out/$t.txt"
    done; } > "$d/cfg.yaml"
  "$S2BIN" "$d/cfg.yaml" > "$d/log.txt" 2>&1 || { say "FAIL step 2 $g $p, see $W/$d/log.txt"; exit 1; }
done; done

# 4. checks
ok=1
for g in full sparse; do
  if grep -q "useGPU: refused" s2/${g}_gpu/log.txt; then
    say "FAIL (a) $g: GPU not used: $(grep -o 'useGPU: refused.*' s2/${g}_gpu/log.txt)"; ok=0
  else say "     (a) $g: $(grep -m1 -o 'GPU coverage: .*' s2/${g}_gpu/log.txt)"; fi
  if diff -rq s2/${g}_cpu/out s2/${g}_gpu/out > /dev/null; then say "PASS (a) $g: CPU and GPU outputs byte-identical"
  else say "FAIL (a) $g: CPU and GPU outputs differ"; ok=0; fi
done
nm=0; nb=0
for k in $(printf '%s\n' "${!EXP[@]}" | sort); do
  case $k in data/*) f=$k;; *) f=s2/${k%%/*}_cpu/out/${k#*/};; esac
  got=$(md5sum "$f" | cut -d' ' -f1)
  if [ "$got" = "${EXP[$k]}" ]; then nm=$((nm + 1)); else nb=$((nb + 1)); say "     (b) $k: md5 $got, recorded ${EXP[$k]}"; fi
done
say "$([ $nb = 0 ] && echo PASS || echo NOTE) (b) $nm of $((nm + nb)) md5s match the recorded ones"
[ $ok = 1 ] && say "SELF-CHECK OK" || say "SELF-CHECK FAILED"
echo "result: $RES"
