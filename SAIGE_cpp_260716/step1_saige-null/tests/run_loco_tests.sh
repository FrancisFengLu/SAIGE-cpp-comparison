#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# LOCO regression tests for saige-null (step 1).
#
#   make test                 # from step1_saige-null/
#   tests/run_loco_tests.sh   # equivalent
#
# Environment:
#   SAIGE_NULL_BIN   path to the binary under test (default ../saige-null)
#   SAIGE_NULL_REF   path to a PRE-CHANGE binary. When set, test 1 asserts that
#                    a non-LOCO run is byte-identical between the two builds.
#                    Without it, test 1 degrades to a self-reproducibility check
#                    (same binary, two runs) and says so.
#   SAIGE_TEST_DATA  input directory (default: repo test_data/SAIGE/extdata/input)
#   SAIGE_TEST_TMP   scratch directory (default: mktemp -d)
#
# NOTE on determinism (measured 2026-08-10, predates the LOCO work):
#   nthreads > 1  -> NOT reproducible. Two back-to-back runs of the identical
#                    binary+config on plinkforGRM_1000samples_10kMarkers with
#                    nthreads: 4 gave tau[1] = 0.233564 and 0.279550 (a 20%
#                    swing), and mu.arma / V.arma differed. Cause: the TBB
#                    parallelReduce in parallelCrossProd accumulates marker
#                    contributions in a nondeterministic order in float32, and
#                    AI-REML amplifies that into tau.
#   nthreads == 1 -> fully reproducible, byte-for-byte, including nullmodel.json.
# Every config in this suite therefore pins nthreads: 1 and trace_seed: 10.
# Test 1 asserts byte identity outright; it does NOT downgrade to a tolerance.
# ---------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP1_DIR="$(cd "$HERE/.." && pwd)"
REPO_ROOT="$(cd "$STEP1_DIR/../.." && pwd)"

BIN="${SAIGE_NULL_BIN:-$STEP1_DIR/saige-null}"
DATA="${SAIGE_TEST_DATA:-$REPO_ROOT/test_data/SAIGE/extdata/input}"
TMP="${SAIGE_TEST_TMP:-$(mktemp -d /tmp/saige_loco_test.XXXXXX)}"
mkdir -p "$TMP"
TPL="$HERE/loco/config_template.yaml"

PLINK_2CHR="$DATA/plinkforGRM_1000samples_10kMarkers"
PLINK_22CHR="$DATA/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr_random1000"
PHENO="$DATA/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"

PASS=0; FAIL=0
ok()   { echo "  PASS: $*"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL: $*"; FAIL=$((FAIL+1)); }
head1() { echo; echo "=== $* ==="; }

if [[ ! -x "$BIN" ]]; then echo "binary not found/executable: $BIN" >&2; exit 2; fi
if [[ -z "${R_HOME:-}" ]] && command -v Rscript >/dev/null 2>&1; then
  R_HOME="$(Rscript -e 'cat(R.home())')"; export R_HOME
fi

mkcfg() { # mkcfg <out.yaml> <plink> <outdir> <pheno> <loco> <sparse_fit> [ycol]
  local out="$1" plink="$2" dir="$3" pheno="$4" loco="$5" sparse="$6" ycol="${7:-y_binary}"
  sed -e "s|@PLINK@|$plink|g" -e "s|@OUT@|$dir|g" -e "s|@PHENO@|$pheno|g" \
      -e "s|@LOCO@|$loco|g" -e "s|@SPARSE_FIT@|$sparse|g" \
      -e "s|@SPARSE_GRM@||g" -e "s|@SPARSE_GRM_IDS@||g" \
      -e "s|y_col: y_binary|y_col: $ycol|" "$TPL" > "$out"
}

run() { # run <cfg> <logfile>
  ( cd "$TMP" && "$BIN" -c "$1" ) > "$2" 2>&1
}

# manifest of the artifacts we compare for byte identity (excludes logs, which
# contain timings)
STABLE_FILES="nullmodel.json mu.arma res.arma y.arma V.arma S_a.arma X.arma \
XV.arma XVX.arma XVX_inv.arma XXVX_inv.arma XVX_inv_XV.arma offset.arma"

cmp_dirs() { # cmp_dirs <a> <b> -> prints differing files
  local a="$1" b="$2" f rc=0
  for f in $STABLE_FILES; do
    if ! cmp -s "$a/$f" "$b/$f"; then echo "    differs: $f"; rc=1; fi
  done
  return $rc
}

PY=python3

# ---------------------------------------------------------------------------
head1 "Test 1: non-LOCO run is byte-identical (regression guard)"
mkcfg "$TMP/t1a.yaml" "$PLINK_2CHR" "$TMP/t1a" "$PHENO" false false
run "$TMP/t1a.yaml" "$TMP/t1a.log" || { bad "run t1a failed (see $TMP/t1a.log)"; }

if [[ -n "${SAIGE_NULL_REF:-}" && -x "${SAIGE_NULL_REF}" ]]; then
  mkcfg "$TMP/t1b.yaml" "$PLINK_2CHR" "$TMP/t1b" "$PHENO" false false
  ( cd "$TMP" && "$SAIGE_NULL_REF" -c "$TMP/t1b.yaml" ) > "$TMP/t1b.log" 2>&1
  echo "  comparing new build vs reference build $SAIGE_NULL_REF"
else
  echo "  SAIGE_NULL_REF not set -> comparing two runs of the SAME binary"
  echo "  (this checks reproducibility, NOT pre/post-change equivalence)"
  mkcfg "$TMP/t1b.yaml" "$PLINK_2CHR" "$TMP/t1b" "$PHENO" false false
  run "$TMP/t1b.yaml" "$TMP/t1b.log"
fi

if cmp_dirs "$TMP/t1a" "$TMP/t1b" > "$TMP/t1.diff"; then
  ok "non-LOCO artifacts identical"
else
  # nullmodel.json legitimately gains one additive line ("loco_chroms": []) vs a
  # pre-LOCO reference build; that is part of the LOCO_FORMAT.md contract. Every
  # .arma payload must still be byte-identical. Anything else is a regression.
  onlyjson=1
  while read -r _ f; do [[ "$f" == "nullmodel.json" ]] || onlyjson=0; done < "$TMP/t1.diff"
  if [[ $onlyjson -eq 1 ]] \
     && diff <(grep -v '"loco_chroms"' "$TMP/t1a/nullmodel.json") \
             <(grep -v '"loco_chroms"' "$TMP/t1b/nullmodel.json") > /dev/null; then
    ok "all .arma payloads byte-identical; nullmodel.json differs only by the additive \"loco_chroms\" field"
  else
    bad "non-LOCO artifacts DIFFER"
    cat "$TMP/t1.diff"
  fi
fi
if [[ -d "$TMP/t1a/chr1" ]]; then bad "non-LOCO run created chr1/"; else ok "no chr<j>/ dirs in a non-LOCO run"; fi
if grep -q '"loco": false' "$TMP/t1a/nullmodel.json"; then ok 'nullmodel.json has "loco": false'; else bad 'nullmodel.json loco flag wrong'; fi
if grep -q '"loco_chroms": \[\]' "$TMP/t1a/nullmodel.json"; then ok 'loco_chroms is []'; else bad 'loco_chroms not empty for non-LOCO run'; fi

# ---------------------------------------------------------------------------
head1 "Test 2: LOCO on 22-chromosome input produces the full chr<j>/ set"
mkcfg "$TMP/t2.yaml" "$PLINK_22CHR" "$TMP/t2" "$PHENO" true false
run "$TMP/t2.yaml" "$TMP/t2.log" || bad "run t2 failed (see $TMP/t2.log)"

missing=0
for j in $(seq 1 22); do
  d="$TMP/t2/chr$j"
  [[ -d "$d" ]] || { echo "    missing dir chr$j"; missing=1; continue; }
  for f in mu res V offset XV XVX XVX_inv XVX_inv_XV XXVX_inv S_a; do
    [[ -s "$d/$f.arma" ]] || { echo "    missing $j/$f.arma"; missing=1; }
  done
done
[[ $missing -eq 0 ]] && ok "all 22 chr<j>/ dirs contain the 10-file per-chromosome set" \
                     || bad "chr<j>/ set incomplete"

if grep -q '"loco": true' "$TMP/t2/nullmodel.json"; then ok '"loco": true'; else bad '"loco" not true'; fi
if grep -q '"loco_chroms": \[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22\]' "$TMP/t2/nullmodel.json"; then
  ok "loco_chroms == 1..22"
else
  bad "loco_chroms wrong: $(grep loco_chroms "$TMP/t2/nullmodel.json")"
fi

# ---------------------------------------------------------------------------
head1 "Test 3: chromosome-invariant quantities are NOT duplicated into chr<j>/"
dup=0
for j in $(seq 1 22); do
  [[ -e "$TMP/t2/chr$j/X.arma" ]] && { echo "    chr$j/X.arma exists"; dup=1; }
  [[ -e "$TMP/t2/chr$j/y.arma" ]] && { echo "    chr$j/y.arma exists"; dup=1; }
done
[[ $dup -eq 0 ]] && ok "no X.arma / y.arma under chr<j>/" || bad "X/y duplicated into chr<j>/"
[[ -s "$TMP/t2/X.arma" && -s "$TMP/t2/y.arma" ]] && ok "top-level X.arma and y.arma present" \
                                                 || bad "top-level X.arma / y.arma missing"

# ---------------------------------------------------------------------------
head1 "Test 4: numeric sanity — per-chromosome mu really differs"
# Uses the 2-chromosome dataset, where tau[1] > 0 so the GRM actually matters.
# On the 22-chr random-1000 dataset AI-REML drives tau[1] to 0, which makes
# Sigma diagonal and LOCO a mathematical no-op; that dataset can therefore
# only be used for the structural checks above.
mkcfg "$TMP/t4.yaml" "$PLINK_2CHR" "$TMP/t4" "$PHENO" true false
run "$TMP/t4.yaml" "$TMP/t4.log" || bad "run t4 failed (see $TMP/t4.log)"

"$PY" - "$TMP/t4" <<'EOF'
import sys, os, numpy as np
d = sys.argv[1]
def rd(p):
    with open(p,'rb') as f:
        f.readline(); dim=f.readline().split(); r,c=int(dim[0]),int(dim[1])
        a=np.fromfile(f,dtype=np.float64,count=r*c)
        return a.reshape((c,r)).T if c>1 else a
tau1 = None
import re, json
txt = open(os.path.join(d,'nullmodel.json')).read()
m = re.search(r'"theta":\s*\[([^\]]*)\]', txt)
theta = [float(x) for x in m.group(1).split(',')]
print(f"    theta = {theta}")
if theta[1] <= 0:
    print("    FAIL: tau[1] == 0, LOCO is a no-op on this dataset — test is vacuous")
    sys.exit(1)
full = rd(os.path.join(d,'mu.arma'))
chrs = sorted(int(x[3:]) for x in os.listdir(d) if x.startswith('chr') and os.path.isdir(os.path.join(d,x)))
mus = {j: rd(os.path.join(d,f'chr{j}','mu.arma')) for j in chrs}
rc = 0
for j in chrs:
    delta = np.abs(mus[j]-full).max()
    rel   = delta/np.abs(full).max()
    print(f"    chr{j}: max|mu_chr - mu_full| = {delta:.3e}  (rel {rel:.3e})")
    if rel < 1e-6:
        print(f"    FAIL: chr{j} mu is indistinguishable from the full-genome mu "
              f"-> the chromosome exclusion did nothing")
        rc = 1
for a in range(len(chrs)):
    for b in range(a+1, len(chrs)):
        ja, jb = chrs[a], chrs[b]
        delta = np.abs(mus[ja]-mus[jb]).max()
        if delta/np.abs(full).max() < 1e-6:
            print(f"    FAIL: chr{ja} and chr{jb} mu are identical")
            rc = 1
# res must equal y - mu, and V must equal mu(1-mu), per chromosome
y = rd(os.path.join(d,'y.arma'))
for j in chrs:
    res = rd(os.path.join(d,f'chr{j}','res.arma'))
    V   = rd(os.path.join(d,f'chr{j}','V.arma'))
    if np.abs(res-(y-mus[j])).max() > 1e-6:
        print(f"    FAIL: chr{j} res != y - mu"); rc = 1
    if np.abs(V-mus[j]*(1-mus[j])).max() > 1e-6:
        print(f"    FAIL: chr{j} V != mu(1-mu)"); rc = 1
    off_top = rd(os.path.join(d,'offset.arma'))
    off_chr = rd(os.path.join(d,f'chr{j}','offset.arma'))
    if np.abs(off_top-off_chr).max() != 0.0:
        print(f"    FAIL: chr{j} offset differs from the top-level offset"); rc = 1
sys.exit(rc)
EOF
if [[ $? -eq 0 ]]; then ok "per-chromosome mu differs from full-genome mu and from each other"
else bad "LOCO numerics look like a no-op (or res/V/offset inconsistent)"; fi

# ---------------------------------------------------------------------------
head1 "Test 5: auto-disable — single-autosome input"
# Synthesized from the 2-chromosome plink set by relabelling every marker to
# chr1 (same BED/FAM, so the fit itself is cheap and the only thing that changes
# is the number of autosomes LOCO can see).
mkdir -p "$TMP/onechr"
ln -sf "$PLINK_2CHR.bed" "$TMP/onechr/g.bed"
ln -sf "$PLINK_2CHR.fam" "$TMP/onechr/g.fam"
awk 'BEGIN{OFS="\t"}{$1=1; print}' "$PLINK_2CHR.bim" > "$TMP/onechr/g.bim"
mkcfg "$TMP/t5.yaml" "$TMP/onechr/g" "$TMP/t5" "$PHENO" true false
run "$TMP/t5.yaml" "$TMP/t5.log"
if grep -q "number of autosomal chromosomes is 1" "$TMP/t5.log"; then
  ok "logged the <2-autosome auto-disable"
else
  bad "no auto-disable message (see $TMP/t5.log)"
fi
if grep -q '"loco": false' "$TMP/t5/nullmodel.json" 2>/dev/null; then ok 'single-autosome run reports "loco": false'
else bad 'single-autosome run did not report loco:false'; fi
if compgen -G "$TMP/t5/chr*" > /dev/null; then bad "chr<j>/ written despite auto-disable"; else ok "no chr<j>/ written"; fi

# ---------------------------------------------------------------------------
head1 "Test 6: auto-disable — sparse GRM used to fit the null model"
SPARSE_MTX="$DATA/../output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
SPARSE_IDS="$SPARSE_MTX.sampleIDs.txt"
if [[ -f "$SPARSE_MTX" && -f "$SPARSE_IDS" ]]; then
  mkcfg "$TMP/t6.yaml" "$PLINK_2CHR" "$TMP/t6" "$PHENO" true true
  sed -i "s|sparse_grm:     \"\"|sparse_grm:     \"$SPARSE_MTX\"|" "$TMP/t6.yaml"
  sed -i "s|sparse_grm_ids: \"\"|sparse_grm_ids: \"$SPARSE_IDS\"|" "$TMP/t6.yaml"
  run "$TMP/t6.yaml" "$TMP/t6.log"
  if [[ -f "$TMP/t6/nullmodel.json" ]] && grep -q '"loco": false' "$TMP/t6/nullmodel.json"; then
    ok 'sparse-GRM fit reports "loco": false'
  else
    bad "sparse-GRM fit did not report loco:false (see $TMP/t6.log)"
  fi
  if compgen -G "$TMP/t6/chr*" > /dev/null; then bad "chr<j>/ written under sparse-GRM fit"; else ok "no chr<j>/ written"; fi
else
  echo "  SKIP: sparse GRM fixture not found ($SPARSE_MTX)"
fi

# ---------------------------------------------------------------------------
echo
echo "==================================="
echo "  passed: $PASS   failed: $FAIL"
echo "  artifacts: $TMP"
echo "==================================="
[[ $FAIL -eq 0 ]] || exit 1
