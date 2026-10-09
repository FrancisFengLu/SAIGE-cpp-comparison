#!/bin/bash
# run_matrix.sh -- the step-2 test matrix of docs/collaborator_tests.md (binary traits).
#
# Blocks (BLOCKS, space separated; default "main rare stage pgen vsr precision"):
#   main   path {cpu,gpu} x P x GRM {full, sparse_nofast, sparse_fast} x Firth {0,1}
#   rare   rare-variant marker set (MAC < 20), cpu + gpu, full GRM, Firth on, P = top level
#   stage  per-stage timing: CPU PHASE_TIMING build + GPU with gpuOverlap false, full, Firth on, P = top
#   pgen   GPU on a hard-call PGEN of the same chromosome, full, Firth on, P = top (needs PGEN)
#   vsr    our CPU vs R SAIGE 1.5.2: P {1,8} x GRM {full, sparse_nofast} x Firth {0,1}
#   precision  GPU per-stage precision, full GRM, Firth on, P = top: all fp64, scan fp32, scan int8,
#          SPA fp32, ER fp32, Firth fp32, all fp32, scan int8 + SPA/ER/Firth fp32; each one
#          compared with the all-fp64 run by precision_compare.py (aggregate numbers only)
#
# Required environment:
#   OUT_ROOT     output root (the step-1 models from step1_models.sh are under $OUT_ROOT/step1)
#   GENO         PLINK prefix (.bed/.bim/.fam) of ONE chromosome, for step 2
#   BIN_TRAITS   the binary trait names, same order as for step1_models.sh (or @file)
# Optional:
#   P_LEVELS     "1 8 32 128". Levels >= the number of traits are dropped and the number of traits
#                becomes the top level (P_TOP=all, default; with more traits than the largest
#                level, the largest level is replaced by all traits). P_TOP=cap keeps the levels.
#   NTHREADS     step-2 threads (default: all cores)
#   PGEN         PLINK 2 prefix (.pgen/.pvar/.psam), hard calls, same chromosome (block pgen)
#   RARE_GENO    PLINK prefix of the rare-marker set; default: made from GENO with plink2
#                (--mac 1 --max-mac 19) into $OUT_ROOT/data/rare
#   PLINK2       plink2 binary (only to make RARE_GENO)
#   R_STEP2      path to R SAIGE 1.5.2's extdata/step2_SPAtests.R (block vsr)
#   RSCRIPT      Rscript to use (default: Rscript); R_PAR (8): R processes run at once
#   CACHE_MODE   auto (default) | drop_caches | vmtouch | fadvise | none   (see cache_evict.py)
#   GPU_ID       GPU index for nvidia-smi sampling (0); the run itself uses CUDA_VISIBLE_DEVICES
#   MIN_MAC      step-2 minMAC (1)
#   KEEP_OUTPUTS 0 (default): delete result files once compared; 1: keep them
#   ONLY / SKIP  extended regex on cell names: run only matching / skip matching cells
#   MIN_FREE_GB  stop when the output filesystem has less free space (20)
#   DRY_RUN=1    print the cell list and exit
# Each cell: $OUT_ROOT/cells/<cell>/ {cfg.yaml, log.txt, cell.json, cache.json, md5.txt, DONE}.
# A cell with DONE is skipped, so the script can simply be started again after an interruption.
set -uo pipefail
source "$(dirname "$0")/common.sh"
: "${OUT_ROOT:?}" "${GENO:?}" "${BIN_TRAITS:?}"
BLOCKS=${BLOCKS:-main rare stage pgen vsr precision}
P_LEVELS=${P_LEVELS:-1 8 32 128}; P_TOP=${P_TOP:-all}
NTHREADS=${NTHREADS:-$NPROC}; CACHE_MODE=${CACHE_MODE:-auto}; GPU_ID=${GPU_ID:-0}
MIN_MAC=${MIN_MAC:-1}; KEEP_OUTPUTS=${KEEP_OUTPUTS:-0}; R_PAR=${R_PAR:-8}; RSCRIPT=${RSCRIPT:-Rscript}
ONLY=${ONLY:-}; SKIP=${SKIP:-}; MIN_FREE_GB=${MIN_FREE_GB:-20}
[[ $BIN_TRAITS == @* ]] && BIN_TRAITS=$(grep -v '^\s*$' "${BIN_TRAITS#@}" | tr '\n' ' ')
read -r -a TR <<< "$BIN_TRAITS"; NT=${#TR[@]}
M=$OUT_ROOT/step1; C=$OUT_ROOT/cells; CMP=$OUT_ROOT/compare; CMPP=$OUT_ROOT/compare_precision
mkdir -p "$C" "$CMP" "$CMPP" "$OUT_ROOT/data"
COMMIT=$(git -C "$SAIGE_HOME" rev-parse --short HEAD 2>/dev/null || echo unknown)

# ---------- P levels ----------
LEVELS=()
for p in $(tr ' ' '\n' <<< "$P_LEVELS" | sort -n); do [ "$p" -lt "$NT" ] && LEVELS+=("$p"); done
NLEV=$(wc -w <<< "$P_LEVELS")
if [ "$P_TOP" = all ] && [ "${#LEVELS[@]}" -eq "$NLEV" ]; then unset 'LEVELS[-1]'; fi
if [ "$P_TOP" = all ] || [ "$(tr ' ' '\n' <<< "$P_LEVELS" | sort -n | tail -1)" -ge "$NT" ]; then LEVELS+=("$NT"); fi
PTOP=${LEVELS[-1]}

# ---------- cell list: name|block|path|P|grm|firth|geno|bin|extra ----------
CELLS=()
add() { CELLS+=("$1"); }
for b in $BLOCKS; do case $b in
  main)
    for P in "${LEVELS[@]}"; do for g in full sparse_nofast sparse_fast; do for f in 0 1; do
      for pth in cpu gpu; do add "main_${pth}_P${P}_${g}_firth${f}|main|$pth|$P|$g|$f|bed|prod|"; done
    done; done; done ;;
  rare)
    for pth in cpu gpu; do add "rare_${pth}_P${PTOP}_full_firth1|rare|$pth|$PTOP|full|1|rare|prod|"; done ;;
  stage)
    add "stage_cpu_P${PTOP}_full_firth1|stage|cpu|$PTOP|full|1|bed|phase|"
    add "stage_gpu_P${PTOP}_full_firth1|stage|gpu|$PTOP|full|1|bed|prod|gpuOverlap: false" ;;
  pgen)
    if [ -n "${PGEN:-}" ]; then add "pgen_gpu_P${PTOP}_full_firth1|pgen|gpu|$PTOP|full|1|pgen|prod|"
    else echo "block pgen: PGEN not set, skipped"; fi ;;
  vsr)
    for P in 1 8; do [ "$P" -gt "$NT" ] && continue; for g in full sparse_nofast; do for f in 0 1; do
      add "vsr_cpp_P${P}_${g}_firth${f}|vsr|cpu|$P|$g|$f|bed|prod|"
      add "vsr_R_P${P}_${g}_firth${f}|vsr|R|$P|$g|$f|bed|R|"
    done; done; done ;;
  precision)
    # name|...|extra: extra YAML lines separated by ';'. prec_fp64 first: the others compare with it.
    for pc in "fp64|" \
              "scan_fp32|gpuPrecisionScan: fp32" \
              "scan_int8|gpuPrecisionScan: int8" \
              "spa_fp32|gpuPrecisionSPA: fp32" \
              "er_fp32|gpuPrecisionER: fp32" \
              "firth_fp32|gpuPrecisionFirth: fp32" \
              "all_fp32|gpuPrecisionScan: fp32;gpuPrecisionSPA: fp32;gpuPrecisionER: fp32;gpuPrecisionFirth: fp32" \
              "int8_fp32|gpuPrecisionScan: int8;gpuPrecisionSPA: fp32;gpuPrecisionER: fp32;gpuPrecisionFirth: fp32"; do
      add "prec_${pc%%|*}_gpu_P${PTOP}_full_firth1|precision|gpu|$PTOP|full|1|bed|prod|${pc#*|}"
    done ;;
  *) echo "unknown block $b"; exit 1 ;;
esac; done
SEL=()
for c in "${CELLS[@]}"; do n=${c%%|*}
  [ -n "$ONLY" ] && ! grep -Eq "$ONLY" <<< "$n" && continue
  [ -n "$SKIP" ] && grep -Eq "$SKIP" <<< "$n" && continue
  SEL+=("$c")
done
echo "traits: $NT   P levels: ${LEVELS[*]}   cells selected: ${#SEL[@]}   threads: $NTHREADS   commit: $COMMIT"
if [ "${DRY_RUN:-0}" = 1 ]; then printf '%s\n' "${SEL[@]%%|*}"; exit 0; fi

# ---------- inputs ----------
need() { [ -e "$1" ] || { echo "missing $1"; exit 1; }; }
need "$GENO.bed"; need "$S2BIN"
for g in full sparse sparse_nofast; do need "$M/$g/models/${TR[0]}/nullmodel.json"; done
geno_files() { case $1 in
  bed)  echo "$GENO.bed $GENO.bim $GENO.fam" ;;
  rare) echo "$RARE_GENO.bed $RARE_GENO.bim $RARE_GENO.fam" ;;
  pgen) echo "$PGEN.pgen $PGEN.pvar $PGEN.psam" ;; esac; }
if grep -q "|rare|" <<< "${SEL[*]}"; then
  if [ -z "${RARE_GENO:-}" ]; then
    RARE_GENO=$OUT_ROOT/data/rare
    if [ ! -s "$RARE_GENO.bed" ]; then
      "$PLINK2" --bfile "$GENO" --mac 1 --max-mac 19 --make-bed --out "$RARE_GENO" > "$OUT_ROOT/data/rare.plink2.log" 2>&1 \
        || { echo "plink2 failed, see $OUT_ROOT/data/rare.plink2.log"; exit 1; }
    fi
  fi
  echo "rare-variant set: $(wc -l < "$RARE_GENO.bim") markers ($RARE_GENO)"
fi
if grep -q "|R|" <<< "${SEL[*]}"; then
  : "${R_STEP2:?R_STEP2 (path to step2_SPAtests.R of R SAIGE 1.5.2) is needed for block vsr}"
  mapfile -t SG < "$M/sparse_grm.paths"
  # R reads .rda null models: convert ours once (not timed)
  for g in full sparse; do for t in "${TR[@]:0:8}"; do
    rda=$OUT_ROOT/rmodels/$g/$t.rda
    [ -s "$rda" ] && continue
    mkdir -p "$(dirname "$rda")"
    "$RSCRIPT" "$COLLAB/arma_to_rda.R" "$M/$g/models/$t" "$rda" > "$rda.log" 2>&1 || { echo "arma_to_rda failed: $rda.log"; exit 1; }
  done; done
fi

# ---------- helpers ----------
mdir() { case $1 in sparse_fast) echo sparse;; *) echo "$1";; esac; }   # model directory of a GRM label
free_gb() { df -Pk "$OUT_ROOT" | awk 'NR==2{print int($4/1048576)}'; }

write_cfg() {   # dir P grm firth geno path extra
  local d=$1 P=$2 g=$3 f=$4 geno=$5 pth=$6 extra=$7 F=false GPUV=false
  [ "$f" = 1 ] && F=true; [ "$pth" = gpu ] && GPUV=true
  {
    case $geno in
      bed)  printf 'genoType: plink\nplinkFile: %s\nAlleleOrder: alt-first\n' "$GENO" ;;
      rare) printf 'genoType: plink\nplinkFile: %s\nAlleleOrder: alt-first\n' "$RARE_GENO" ;;
      pgen) printf 'genoType: pgen\npgenFile: %s.pgen\npvarFile: %s.pvar\npsamFile: %s.psam\nAlleleOrder: ref-first\n' "$PGEN" "$PGEN" "$PGEN" ;;
    esac
    cat <<YAML
minMAF: 0
minMAC: $MIN_MAC
maxMissRate: 0.15
LOCO: false
isFirth: $F
is_Firth_beta: $F
pCutoffforFirth: 0.01
MACCutoffforER: 4
relatednessCutoff: 0
nThreads: $NTHREADS
useGPU: $GPUV
outputFormat: text
YAML
    [ -n "$extra" ] && tr ';' '\n' <<< "$extra"
    echo "models:"
    for t in "${TR[@]:0:$P}"; do
      printf '  - traitName: %s\n    modelFile: %s\n    varianceRatioFile: %s\n    outputFile: %s\n' \
        "$t" "$M/$(mdir "$g")/models/$t" "$M/$(mdir "$g")/vr_$t.varianceRatio.txt" "$d/out/$t.txt"
    done
  } > "$d/cfg.yaml"
}

write_rjobs() {   # dir P grm firth
  local d=$1 P=$2 g=$3 f=$4 F=FALSE sp="" rg=$g
  [ "$f" = 1 ] && F=TRUE
  [ "$g" != full ] && { sp="--sparseGRMFile=${SG[0]} --sparseGRMSampleIDFile=${SG[1]}"; rg=sparse; }
  {
    echo "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1"
    local i=0
    for t in "${TR[@]:0:$P}"; do
      echo "/usr/bin/time -v -o $d/rlogs/$t.time $RSCRIPT $R_STEP2 --bedFile=$GENO.bed --bimFile=$GENO.bim --famFile=$GENO.fam --AlleleOrder=alt-first --GMMATmodelFile=$OUT_ROOT/rmodels/$rg/$t.rda --varianceRatioFile=$M/$g/vr_$t.varianceRatio.txt $sp --SAIGEOutputFile=$d/out/$t.txt --minMAF=0 --minMAC=$MIN_MAC --maxMissing=0.15 --LOCO=FALSE --is_noadjCov=FALSE --is_output_moreDetails=FALSE --is_Firth_beta=$F --pCutoffforFirth=0.01 --impute_method=mean --is_fastTest=FALSE --SPAcutoff=2 --max_MAC_for_ER=4 --relatednessCutoff=0 --nThreads=1 > $d/rlogs/$t.log 2>&1 &"
      i=$((i + 1)); [ $((i % R_PAR)) -eq 0 ] && echo "wait"
    done
    echo "wait"
    echo "! grep -l '^Error' $d/rlogs/*.log"
  } > "$d/rjobs.sh"
}

cache_files() {   # P grm geno path -> files the run reads
  local P=$1 g=$2 geno=$3 pth=$4 t
  geno_files "$geno"
  for t in "${TR[@]:0:$P}"; do
    if [ "$pth" = R ]; then local rg=$g; [ "$g" != full ] && rg=sparse
      echo "$OUT_ROOT/rmodels/$rg/$t.rda"; else echo "$M/$(mdir "$g")/models/$t"; fi
    echo "$M/$(mdir "$g")/vr_$t.varianceRatio.txt"
  done
  [ "$pth" = R ] && [ "$g" != full ] && echo "${SG[0]} ${SG[1]}"
  return 0
}

run_cell() {   # spec
  IFS='|' read -r name blk pth P g f geno bin extra <<< "$1"
  local d=$C/$name
  [ -e "$d/DONE" ] && { echo "done  $name"; return 0; }
  [ "$(free_gb)" -lt "$MIN_FREE_GB" ] && { echo "less than $MIN_FREE_GB GB free on $OUT_ROOT; stopping"; exit 2; }
  rm -rf "$d"; mkdir -p "$d/out"
  local fast=na; case $g in sparse_fast) fast=true;; sparse_nofast) fast=false;; esac
  local grm=full; [ "$g" != full ] && grm=sparse
  local B=$S2BIN; [ "$bin" = phase ] && B=$S2PHASE
  printf '{"cell":"%s","block":"%s","path":"%s","trait_type":"binary","P":%s,"grm":"%s","fast_test":"%s","firth":%s,"geno":"%s","binary":"%s","extra":"%s","nthreads":%s,"commit":"%s","start":"%s"}\n' \
    "$name" "$blk" "$pth" "$P" "$grm" "$fast" "$f" "$geno" "$bin" "$extra" "$([ "$pth" = R ] && echo 1 || echo "$NTHREADS")" "$COMMIT" "$(date -Is)" > "$d/meta.json"
  echo "run   $name"
  $PYTHON "$COLLAB/cache_evict.py" "$CACHE_MODE" "$d/cache.json" $(cache_files "$P" "$g" "$geno" "$pth") | sed 's/^/      /'
  local smp=""
  if [ "$pth" = gpu ] && command -v nvidia-smi > /dev/null; then
    nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$GPU_ID" > "$d/gpu_mem_base.txt" 2>/dev/null
    nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "$GPU_ID" -lms 200 > "$d/gpu_mem.txt" 2>/dev/null &
    smp=$!
  fi
  if [ "$pth" = R ]; then
    mkdir -p "$d/rlogs"; write_rjobs "$d" "$P" "$g" "$f"
    ( cd "$d" && /usr/bin/time -v bash rjobs.sh > log.txt 2>&1 ); echo $? > "$d/rc"
  else
    write_cfg "$d" "$P" "$g" "$f" "$geno" "$pth" "$extra"
    mkdir -p "$d/routes"
    ( cd "$d" && SAIGE_STEP2_ROUTE_DUMP=$d/routes /usr/bin/time -v "$B" cfg.yaml > log.txt 2>&1 ); echo $? > "$d/rc"
  fi
  [ -n "$smp" ] && { kill "$smp" 2>/dev/null; wait "$smp" 2>/dev/null; }
  $PYTHON "$COLLAB/cell_post.py" "$d" | sed 's/^/      /'
  if [ "$blk" = stage ]; then
    $PYTHON "$COLLAB/stage_table.py" "$d" --csv "$d/stage.csv" > "$d/stage.txt" 2>&1 || true
  fi
  [ "$blk" = precision ] || rm -rf "$d/routes"     # precision: kept until compared with all-fp64
  [ "$(cat "$d/rc")" = 0 ] || echo "      WARNING: $name exited with rc=$(cat "$d/rc"); see $d/log.txt"
  touch "$d/DONE"
}

compare_pair() {   # label dirA dirB [--rows]
  local lab=$1 a=$2 b=$3; shift 3
  [ -e "$CMP/$lab.json" ] && return 0
  [ -e "$a/DONE" ] && [ -e "$b/DONE" ] || return 0
  if [ "$(cat "$a/rc")" != 0 ] || [ "$(cat "$b/rc")" != 0 ]; then
    echo "{\"label\":\"$lab\",\"identical\":null,\"note\":\"a run failed (rc $(cat "$a/rc") / $(cat "$b/rc"))\"}" > "$CMP/$lab.json"
    echo "      $lab: not compared, a run failed"; return 0
  fi
  if [ -z "$(ls "$a/out" 2>/dev/null)" ] || [ -z "$(ls "$b/out" 2>/dev/null)" ]; then
    # outputs already removed: md5 lists only
    if cmp -s "$a/md5.txt" "$b/md5.txt"; then echo "{\"label\":\"$lab\",\"identical\":true,\"from\":\"md5.txt\"}" > "$CMP/$lab.json"
    else echo "{\"label\":\"$lab\",\"identical\":false,\"from\":\"md5.txt\"}" > "$CMP/$lab.json"; fi
    echo "      $lab: md5 lists $(cmp -s "$a/md5.txt" "$b/md5.txt" && echo identical || echo DIFFER)"
    return 0
  fi
  $PYTHON "$COLLAB/compare_outputs.py" "$a/out" "$b/out" --label "$lab" --json "$CMP/$lab.json" "$@" | sed 's/^/      /'
}

compare_prec() {   # label dirA(all fp64) dirB -> $CMPP/label.{json,txt}; aggregate numbers only
  local lab=$1 a=$2 b=$3
  [ -e "$CMPP/$lab.json" ] && return 0
  [ -e "$a/DONE" ] && [ -e "$b/DONE" ] || return 0
  if [ "$(cat "$a/rc")" != 0 ] || [ "$(cat "$b/rc")" != 0 ]; then
    echo "{\"label\":\"$lab\",\"note\":\"a run failed (rc $(cat "$a/rc") / $(cat "$b/rc"))\"}" > "$CMPP/$lab.json"
    echo "      $lab: not compared, a run failed"; return 0
  fi
  if [ -z "$(ls "$a/out" 2>/dev/null)" ] || [ -z "$(ls "$b/out" 2>/dev/null)" ]; then
    echo "{\"label\":\"$lab\",\"note\":\"result files already removed\"}" > "$CMPP/$lab.json"
    echo "      $lab: not compared, result files already removed"; return 0
  fi
  $PYTHON "$PREC_CMP" "$a" "$b" --quiet --no-ids --json "$CMPP/$lab.json" > "$CMPP/$lab.txt" 2>&1
  sed -n 's/^   //; /^p.value \|^BETA \|^SE \|crossings\|p < 1e-5\|Is.SPA\|routes:/p' "$CMPP/$lab.txt" | sed "s/^/      $lab: /"
}

cleanup_outputs() {   # dirs... : remove result files when every comparison involving them is done
  [ "$KEEP_OUTPUTS" = 1 ] && return 0
  local d; for d in "$@"; do rm -f "$d"/out/*; done     # results and R's <output>.index files
}

# ---------- run ----------
for spec in "${SEL[@]}"; do
  run_cell "$spec"
  name=${spec%%|*}
  case $name in
    main_gpu_*|rare_gpu_*)
      cpu=$C/${name/_gpu_/_cpu_}; gpu=$C/$name
      if [ -e "$cpu/DONE" ]; then
        compare_pair "cpu_vs_gpu__${name/_gpu_/_}" "$cpu" "$gpu"
        grep -q '"identical": *true' "$CMP/cpu_vs_gpu__${name/_gpu_/_}.json" && cleanup_outputs "$cpu" "$gpu"
      fi ;;
    vsr_R_*)
      cpp=$C/${name/_R_/_cpp_}
      compare_pair "cpp_vs_R__${name#vsr_R_}" "$cpp" "$C/$name" --rows
      cleanup_outputs "$cpp" "$C/$name" ;;
    pgen_gpu_*)
      ref=$C/main_gpu_P${PTOP}_full_firth1
      [ -e "$ref/DONE" ] && compare_pair "pgen_vs_bed__P${PTOP}_full_firth1" "$ref" "$C/$name"
      cleanup_outputs "$C/$name" ;;
    stage_*) cleanup_outputs "$C/$name" ;;
    prec_fp64_*) ;;
    prec_*)
      ref=$C/prec_fp64_gpu_${name##*_gpu_}; pl=${name%%_gpu_*}; pl=prec_vs_fp64__${pl#prec_}
      compare_prec "$pl" "$ref" "$C/$name"
      [ -e "$CMPP/$pl.json" ] && { cleanup_outputs "$C/$name"; [ "$KEEP_OUTPUTS" = 1 ] || rm -rf "$C/$name/routes"; } ;;
  esac
done
# precision: the all-fp64 reference once every selected comparison with it is written
for spec in "${SEL[@]}"; do
  name=${spec%%|*}
  case $name in prec_fp64_*)
    left=0
    for s2 in "${SEL[@]}"; do n2=${s2%%|*}
      pl=${n2%%_gpu_*}; pl=prec_vs_fp64__${pl#prec_}
      case $n2 in prec_fp64_*) ;; prec_*) [ -e "$CMPP/$pl.json" ] || left=1 ;; esac
    done
    [ "$left" = 0 ] && { cleanup_outputs "$C/$name"; [ "$KEEP_OUTPUTS" = 1 ] || rm -rf "$C/$name/routes"; } ;;
  esac
done
# CPU cells whose GPU partner was not selected: nothing left to compare against
for spec in "${SEL[@]}"; do
  name=${spec%%|*}
  case $name in main_cpu_*|rare_cpu_*)
    grep -q "|${name/_cpu_/_gpu_}|" <<< "|$(printf '%s|' "${SEL[@]%%|*}")" || cleanup_outputs "$C/$name" ;;
  esac
done
$PYTHON "$COLLAB/make_summary.py" "$OUT_ROOT"
echo "summary: $OUT_ROOT/summary.csv  $OUT_ROOT/compare.csv"
