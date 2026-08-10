#!/usr/bin/env bash
# =============================================================================
# Head-to-head wall-clock benchmark: C++ SAIGE port vs R SAIGE 1.5.1.
#
#   bench_cpp_vs_r.sh <data_dir> <out_dir> [phase]
#
#     data_dir   output of gen_bench_data.py (grm.*, wes.*, group.txt, pheno.txt,
#                design.csv)
#     out_dir    where the CSV + raw per-run artefacts are written
#     phase      A   OPENBLAS_NUM_THREADS=1  (isolates our own OMP scaling)
#                B   BLAS threading left at the machine default (what a user gets)
#                all (default) = A then B
#
# Env knobs
#   REPS=3            measured repetitions per cell (a warm-up run is always
#                     done first and recorded with warmup=1; drop it)
#   STAGES="step1 region single"
#   THREADS_A="4 2 1" / THREADS_B="4 1"   (cheapest first: the warm-up uses the first)
#   REPS_B=2
#   SKIP_PREP=1       reuse an existing r_null.rda / r_null_arma
#
# Everything measured with /usr/bin/time -v; wall, user, sys and max RSS are
# recorded for every single run.  Raw time/stdout files are kept under
# <out_dir>/raw so the numbers can be audited later.
#
# NOTE on fairness
#   * Both engines run the SAME null model in step 2: the R .rda is converted to
#     the C++ .arma directory by rda_to_arma.R, so step-2 differences are step-2
#     only.
#   * LOCO is off on both sides so the two do identical work.
#   * The R SAIGE used here is the DE-INSTRUMENTED build (R_LIBS below); the
#     copy in code_copy/SAIGE_isolated writes debug .rds/.csv dumps on every
#     step-1 iteration and would make R look artificially slow.
# =============================================================================
set -uo pipefail

DATA="${1:?usage: bench_cpp_vs_r.sh <data_dir> <out_dir> [phase]}"
OUT="${2:?}"
PHASE="${3:-all}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"                 # SAIGE_cpp_260716
REPO="$(cd "$ROOT/.." && pwd)"                    # SAIGE-cpp-comparison
S1BIN="${S1BIN:-$ROOT/step1_saige-null/saige-null}"
S2BIN="${S2BIN:-$ROOT/step2_saige-step2/saige-step2}"
RLIB="${RLIB:-$REPO/bench_data/Rlib_clean}"

REPS="${REPS:-3}"
REPS_B="${REPS_B:-2}"
STAGES="${STAGES:-step1 region single}"
THREADS_A="${THREADS_A:-4 2 1}"
THREADS_B="${THREADS_B:-4 1}"

[[ -x "$S1BIN" ]] || { echo "no saige-null at $S1BIN" >&2; exit 2; }
[[ -x "$S2BIN" ]] || { echo "no saige-step2 at $S2BIN" >&2; exit 2; }
DATA="$(cd "$DATA" && pwd)"
mkdir -p "$OUT/raw" "$OUT/results"
OUT="$(cd "$OUT" && pwd)"
CSV="$OUT/timings.csv"
[[ -f "$CSV" ]] || echo "phase,blas,stage,engine,threads,rep,warmup,exit,wall_s,user_s,sys_s,max_rss_kb,n_units" > "$CSV"

# shellcheck disable=SC1091
source "${CONDA_ROOT:-$HOME/miniforge3}/etc/profile.d/conda.sh"
conda activate "${ENV_NAME:-saige-build}"
export R_HOME="$(Rscript -e 'cat(R.home())')"
export R_LIBS="$RLIB"

log() { echo "[bench $(date +%H:%M:%S)] $*" >&2; }

run_cell() {  # run_cell <phase> <blas_env> <stage> <engine> <threads> <rep> <warmup>
  local ph=$1 blas=$2 st=$3 en=$4 th=$5 rp=$6 wu=$7
  local tag="${ph}_${st}_${en}_T${th}_rep${rp}"
  [[ "$wu" == 1 ]] && tag="${ph}_${st}_${en}_T${th}_warmup"
  local tf="$OUT/raw/${tag}.time" lf="$OUT/raw/${tag}.log"
  local -a cmd
  local outfile="$OUT/results/${tag}.out"

  case "$st:$en" in
    step1:cpp)
      local cfg="$OUT/raw/${tag}.yaml"
      cat > "$cfg" <<EOF
paths:
  plinkFile:     $DATA/grm
  out_prefix:    $OUT/results/${tag}_null
  out_prefix_vr: $OUT/results/${tag}_null
  overwrite_varratio: true
design:
  csv: $DATA/design.csv
  iid_col: IID
  y_col: y_quantitative
  covar_cols:
    - x1
    - x2
fit:
  trait: quantitative
  loco: false
  nthreads: $th
  maxiter: 20
  tol: 0.02
  tolPCG: 1e-5
  maxiterPCG: 500
  nrun: 30
  num_markers_for_vr: 30
  min_maf_grm: 0.01
  overwrite_vr: true
EOF
      cmd=("$S1BIN" -c "$cfg") ;;
    step1:r)
      cmd=(Rscript "$HERE/bench_r_step1.R" "$DATA/grm" "$DATA/pheno.txt" \
           y_quantitative quantitative "$OUT/results/${tag}_null" "$th") ;;
    region:cpp)
      local cfg="$OUT/raw/${tag}.yaml"
      cat > "$cfg" <<EOF
modelFile:         $DATA/r_null_arma
varianceRatioFile: $DATA/r_null.varianceRatio.txt
genoType:  plink
plinkFile: $DATA/wes
outputFile: $outfile
minMAF: 0
minMAC: 0.5
maxMissRate: 0.15
AlleleOrder: alt-first
isMoreOutput: false
groupFile: $DATA/group.txt
annotationList:
  - "lof"
  - "missense;lof"
  - "missense;lof;synonymous"
maxMAFList:
  - 0.0001
  - 0.001
  - 0.01
r_corr: 0
MACCutoff_to_CollapseUltraRare: 10
markers_per_chunk_in_groupTest: 500
max_markers_region: 100000
groups_per_chunk: 100
weights_beta:
  - 1
  - 25
isSingleInGroupTest: false
isOutputMarkerList: false
isFirth: false
nThreads: $th
EOF
      cmd=("$S2BIN" "$cfg") ;;
    region:r)
      cmd=(Rscript "$HERE/bench_r_step2.R" region "$DATA/r_null" "$DATA/wes" \
           "$outfile" "$DATA/group.txt") ;;
    single:cpp)
      local cfg="$OUT/raw/${tag}.yaml"
      cat > "$cfg" <<EOF
modelFile:         $DATA/r_null_arma
varianceRatioFile: $DATA/r_null.varianceRatio.txt
genoType:  plink
plinkFile: $DATA/wes
outputFile: $outfile
minMAF: 0
minMAC: 0.5
maxMissRate: 0.15
AlleleOrder: alt-first
isMoreOutput: true
isFirth: false
nThreads: $th
EOF
      cmd=("$S2BIN" "$cfg") ;;
    single:r)
      cmd=(Rscript "$HERE/bench_r_step2.R" single "$DATA/r_null" "$DATA/wes" "$outfile") ;;
    *) echo "bad cell $st:$en" >&2; return 1 ;;
  esac

  log "$tag  (blas=$blas)"
  # OPENBLAS_NUM_THREADS only touches the BLAS; our own OMP parallelism is
  # driven by the nThreads/nthreads key in the config, so the two are separable.
  env OPENBLAS_NUM_THREADS="$blas" MKL_NUM_THREADS="$blas" \
      /usr/bin/time -v "${cmd[@]}" > "$lf" 2> "$tf"
  local ec=$?

  local wall user sys rss
  wall=$(grep -oP 'Elapsed \(wall clock\).*: \K[0-9:.]+' "$tf" | tail -1)
  user=$(grep -oP 'User time \(seconds\): \K[0-9.]+' "$tf" | tail -1)
  sys=$( grep -oP 'System time \(seconds\): \K[0-9.]+' "$tf" | tail -1)
  rss=$( grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$tf" | tail -1)
  case "$wall" in
    *:*:*) wall=$(awk -F: '{print $1*3600+$2*60+$3}' <<<"$wall") ;;
    *:*)   wall=$(awk -F: '{print $1*60+$2}' <<<"$wall") ;;
  esac
  local units=NA
  [[ -f "$outfile" ]] && units=$(( $(wc -l < "$outfile") - 1 ))
  echo "$ph,$blas,$st,$en,$th,$rp,$wu,$ec,${wall:-NA},${user:-NA},${sys:-NA},${rss:-NA},$units" >> "$CSV"
  log "  -> exit=$ec wall=${wall}s user=${user}s sys=${sys}s rss=${rss}kb units=$units"
}

# ---- prep: R step-1 null + .arma conversion shared by all step-2 cells ------
prep() {
  if [[ "${SKIP_PREP:-0}" == 1 && -f "$DATA/r_null.rda" && -d "$DATA/r_null_arma" ]]; then
    log "prep: reusing $DATA/r_null.rda"
    return
  fi
  log "prep: fitting the shared step-2 null model (R, 4 threads)"
  Rscript "$HERE/bench_r_step1.R" "$DATA/grm" "$DATA/pheno.txt" y_quantitative \
          quantitative "$DATA/r_null" 4 > "$OUT/raw/prep_r_step1.log" 2>&1
  Rscript "$HERE/rda_to_arma.R" "$DATA/r_null.rda" "$DATA/r_null_arma" \
          > "$OUT/raw/prep_rda_to_arma.log" 2>&1
}

phase_run() {  # phase_run <name> <blas> <reps> <threads list>
  local ph=$1 blas=$2 reps=$3 thl=$4
  # The warm-up run is done ONCE per (stage, engine) -- on the FIRST thread
  # value in the list, which is why THREADS_* is ordered cheapest-first.  Page
  # cache is per-file, not per-thread-count, so warming it once is sufficient;
  # doing a full warm-up at every thread count would triple the wall budget.
  for st in $STAGES; do
    local first=1
    for th in $thl; do
      if [[ $first == 1 ]]; then run_cell "$ph" "$blas" "$st" cpp "$th" 0 1; first=0; fi
      for r in $(seq 1 "$reps"); do run_cell "$ph" "$blas" "$st" cpp "$th" "$r" 0; done
    done
    # --- R: step 2 has no thread knob; step 1 does (nThreads) ---
    if [[ "$st" == step1 ]]; then
      first=1
      for th in $thl; do
        if [[ $first == 1 ]]; then run_cell "$ph" "$blas" "$st" r "$th" 0 1; first=0; fi
        for r in $(seq 1 "$reps"); do run_cell "$ph" "$blas" "$st" r "$th" "$r" 0; done
      done
    else
      run_cell "$ph" "$blas" "$st" r 1 0 1
      for r in $(seq 1 "$reps"); do run_cell "$ph" "$blas" "$st" r 1 "$r" 0; done
    fi
  done
}

prep
NCPU=$(nproc)
case "$PHASE" in
  A)   phase_run A 1 "$REPS" "$THREADS_A" ;;
  B)   phase_run B "$NCPU" "$REPS_B" "$THREADS_B" ;;
  all) phase_run A 1 "$REPS" "$THREADS_A"
       phase_run B "$NCPU" "$REPS_B" "$THREADS_B" ;;
  *) echo "unknown phase $PHASE" >&2; exit 2 ;;
esac

log "done -> $CSV"
column -s, -t "$CSV"
