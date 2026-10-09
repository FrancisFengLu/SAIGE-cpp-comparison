#!/bin/bash
# collect_results.sh -- pack what we need back into ONE tar, and nothing else.
#
# Included (allow-list): environment record, build info, self-check result, step-1 configs and
# logs, per-run configs / logs / timing / cache records / GPU-memory samples / route-count
# summaries / md5s of result files, comparison summaries (aggregate numbers only),
# precision comparisons with all-fp64 (compare_precision/*.json and .txt, aggregate numbers,
# written with --no-ids), summary.csv, compare.csv.
# Never included: result files (out/*.txt), route dumps, null models, genotype or phenotype
# files, sparse GRM, sample ID lists. As a last check every included text file is scanned for
# sample IDs and marker IDs (redact.py); lines that contain one are removed and counted.
#
# Environment: OUT_ROOT, GENO, PHENO (+ IID_COL), STEP1_BFILE; optional PGEN, RARE_GENO,
#   SELFCHECK_DIR (default $OUT_ROOT/selfcheck), TAR (default $OUT_ROOT/saige_collab_results_<date>.tar.gz),
#   REDACT_MIN_LEN (5: sample IDs shorter than this are not searched for)
set -euo pipefail
source "$(dirname "$0")/common.sh"
: "${OUT_ROOT:?}" "${GENO:?}" "${PHENO:?}" "${STEP1_BFILE:?}"
IID_COL=${IID_COL:-IID}
SELFCHECK_DIR=${SELFCHECK_DIR:-$OUT_ROOT/selfcheck}
TAR=${TAR:-$OUT_ROOT/saige_collab_results_$(date +%Y%m%d_%H%M).tar.gz}
ST=$OUT_ROOT/.collect/saige_collab_results
rm -rf "$OUT_ROOT/.collect"; mkdir -p "$ST"

cp_if() { [ -f "$1" ] && { mkdir -p "$(dirname "$2")"; cp "$1" "$2"; }; return 0; }
cp_if "$OUT_ROOT/env/environment.txt" "$ST/env/environment.txt"
cp_if "$OUT_ROOT/env/manual.txt"      "$ST/env/manual.txt"
cp_if "$BIN_DIR/BUILD_INFO.txt"       "$ST/build/BUILD_INFO.txt"
cp_if "$SELFCHECK_DIR/selfcheck_result.txt" "$ST/selfcheck/selfcheck_result.txt"
for d in grm full sparse; do
  cp_if "$OUT_ROOT/step1/$d/step1.yaml" "$ST/step1/$d/step1.yaml"
  cp_if "$OUT_ROOT/step1/$d/step1.log"  "$ST/step1/$d/step1.log"
done
for d in "$OUT_ROOT"/cells/*/; do
  n=$(basename "$d")
  for f in meta.json cell.json cache.json md5.txt cfg.yaml rjobs.sh log.txt rc gpu_mem.txt gpu_mem_base.txt stage.txt stage.csv; do
    cp_if "$d/$f" "$ST/cells/$n/$f"
  done
  # R: the /usr/bin/time records only. R SAIGE's own logs print a table of sample IDs for
  # sparse-GRM runs, so they are left out; for a failed R run the last 30 lines are kept.
  for f in "$d"/rlogs/*.time; do cp_if "$f" "$ST/cells/$n/rlogs/$(basename "$f")"; done
  if [ -d "$d/rlogs" ] && [ "$(cat "$d/rc" 2>/dev/null)" != 0 ]; then
    for f in "$d"/rlogs/*.log; do tail -30 "$f" > "$ST/cells/$n/rlogs/$(basename "$f").tail"; done
  fi
done
for f in "$OUT_ROOT"/compare/*.json; do cp_if "$f" "$ST/compare/$(basename "$f")"; done
for f in "$OUT_ROOT"/compare_precision/*.json "$OUT_ROOT"/compare_precision/*.txt; do
  cp_if "$f" "$ST/compare_precision/$(basename "$f")"
done
cp_if "$OUT_ROOT/summary.csv" "$ST/summary.csv"
cp_if "$OUT_ROOT/compare.csv" "$ST/compare.csv"

# last check: sample / marker IDs
args=(--fam "$GENO.fam" --fam "$STEP1_BFILE.fam" --bim "$GENO.bim" --pheno "$PHENO:$IID_COL")
[ -n "${PGEN:-}" ] && args+=(--psam "$PGEN.psam" --pvar "$PGEN.pvar")
RARE_GENO=${RARE_GENO:-$OUT_ROOT/data/rare}
[ -s "$RARE_GENO.bim" ] && args+=(--bim "$RARE_GENO.bim")
$PYTHON "$COLLAB/redact.py" "$ST" "${args[@]}" --min-len "${REDACT_MIN_LEN:-5}"

# nothing big, nothing that looks like a result table
big=$(find "$ST" -type f -size +50M)
[ -z "$big" ] || { echo "refusing: unexpectedly large file(s): $big"; exit 1; }
if grep -rl -E '^CHR	POS	MarkerID' "$ST" > /dev/null; then echo "refusing: a result table header was found"; exit 1; fi

tar -C "$OUT_ROOT/.collect" -czf "$TAR" saige_collab_results
tar -tzvf "$TAR" > "$TAR.contents.txt"
echo "wrote $TAR ($(du -h "$TAR" | cut -f1), $(wc -l < "$TAR.contents.txt") entries); contents: $TAR.contents.txt"
echo "file types in the tar:"
tar -tzf "$TAR" | grep -v '/$' | sed 's/.*\///' | sed -E 's/^.*(\.[a-z]+)$/\1/; s/^(rc|[a-z_]+)$/\1/' | sort | uniq -c
