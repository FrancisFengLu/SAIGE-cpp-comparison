#!/usr/bin/env bash
# LOCO (leave-one-chromosome-out) tests for saige-step2.
#
# Usage:  tests/run_loco_tests.sh [path/to/saige-step2] [path/to/baseline/saige-step2]
#
# The baseline binary is optional; when given, test 1 asserts that a LOCO=false
# run of the new binary is BYTE-IDENTICAL to the pre-change build.
#
# Test data comes from the sibling test_data/ tree; the tests skip (exit 0 with
# a SKIP message) if it is not present.

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP2_DIR="$(dirname "$HERE")"
BIN="${1:-$STEP2_DIR/saige-step2}"
BASELINE_BIN="${2:-}"

ROOT="$(cd "$STEP2_DIR/../.." && pwd)"
DATA="$ROOT/test_data/Step_2_Feb_11/test/data"
EXTDATA="$ROOT/test_data/SAIGE/extdata/input"
SRC_MODEL="$DATA/nullmodel_from_rda_binary"
VR="$DATA/varianceRatio.txt"
PLINK="$EXTDATA/genotype_100markers_2chr"

WORK="${LOCO_TEST_WORK:-/tmp/saige_loco_tests}"

PASS=0; FAIL=0
ok()   { echo "  PASS: $1"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

for f in "$SRC_MODEL/nullmodel.json" "$VR" "$PLINK.bed" "$PLINK.bim" "$PLINK.fam"; do
    if [ ! -f "$f" ]; then
        echo "SKIP: test data not found: $f"
        exit 0
    fi
done
if [ ! -x "$BIN" ]; then echo "FAIL: binary not found: $BIN"; exit 1; fi

rm -rf "$WORK"; mkdir -p "$WORK"
echo "=== Building synthetic LOCO model ==="
python3 "$HERE/make_loco_model.py" "$SRC_MODEL" "$WORK/model_loco" || exit 1
# A copy with loco:false but chr<j>/ directories still present, to prove that
# LOCO:false ignores them and that LOCO:true on a non-LOCO model errors.
cp -r "$WORK/model_loco" "$WORK/model_noloco"
python3 - "$WORK/model_noloco/nullmodel.json" <<'PY'
import sys
p = sys.argv[1]
s = open(p).read().replace('"loco": true', '"loco": false')
open(p, 'w').write(s)
PY

mkcfg() {  # mkcfg <out.yaml> <modeldir> <outputfile> <LOCO> <chrom> [extra...]
    local cfg=$1 model=$2 out=$3 loco=$4 chrom=$5; shift 5
    {
        echo "modelFile:         $model"
        echo "varianceRatioFile: $VR"
        echo "genoType:   plink"
        echo "plinkFile:  $PLINK"
        echo "outputFile: $out"
        echo "minMAF: 0"
        echo "minMAC: 0.5"
        echo "maxMissRate: 0.15"
        echo "AlleleOrder: alt-first"
        echo "isMoreOutput: true"
        echo "isFirth: false"
        echo "nThreads: 1"
        echo "LOCO: $loco"
        echo "chrom: \"$chrom\""
        for e in "$@"; do echo "$e"; done
    } > "$cfg"
}

run() {  # run <cfg> <logfile> -> echoes exit code
    "$BIN" "$1" > "$2" 2>&1
    echo $?
}

##############################################################################
echo
echo "=== Test 1: LOCO off is byte-identical to the pre-change build ==="
mkcfg "$WORK/t1.yaml" "$WORK/model_noloco" "$WORK/t1_new.txt" false ""
rc=$(run "$WORK/t1.yaml" "$WORK/t1_new.log")
if [ "$rc" != 0 ]; then bad "LOCO=false run exited $rc"; sed -n '$p' "$WORK/t1_new.log"; fi

if [ -n "$BASELINE_BIN" ] && [ -x "$BASELINE_BIN" ]; then
    # Baseline binary predates the LOCO keys; it ignores them, but feed it a
    # config without them anyway so the comparison is unambiguous.
    grep -v -e '^LOCO:' -e '^chrom:' "$WORK/t1.yaml" | \
        sed "s#$WORK/t1_new.txt#$WORK/t1_base.txt#" > "$WORK/t1_base.yaml"
    rcb=$("$BASELINE_BIN" "$WORK/t1_base.yaml" > "$WORK/t1_base.log" 2>&1; echo $?)
    if [ "$rcb" != 0 ]; then
        bad "baseline run exited $rcb"
    elif cmp -s "$WORK/t1_new.txt" "$WORK/t1_base.txt"; then
        ok "LOCO=false output byte-identical to baseline build"
    else
        bad "LOCO=false output DIFFERS from baseline build"
        diff "$WORK/t1_base.txt" "$WORK/t1_new.txt" | head -5
    fi
else
    echo "  (no baseline binary given -- skipping byte-identity comparison)"
fi

# LOCO=false must ignore chr<j>/ even when they exist: run against the model
# that HAS chr1/chr2 and confirm the result equals the LOCO=false run above.
mkcfg "$WORK/t1b.yaml" "$WORK/model_loco" "$WORK/t1b.txt" false ""
rc=$(run "$WORK/t1b.yaml" "$WORK/t1b.log")
if [ "$rc" = 0 ] && cmp -s "$WORK/t1_new.txt" "$WORK/t1b.txt"; then
    ok "LOCO=false ignores chr<j>/ subdirectories (guard 4)"
else
    bad "LOCO=false did not ignore chr<j>/ subdirectories (rc=$rc)"
fi

##############################################################################
echo
echo "=== Test 2: LOCO on differs from LOCO off, and differs per chromosome ==="
mkcfg "$WORK/t2c1.yaml" "$WORK/model_loco" "$WORK/t2_chr1.txt" true 1
rc1=$(run "$WORK/t2c1.yaml" "$WORK/t2_chr1.log")
mkcfg "$WORK/t2c2.yaml" "$WORK/model_loco" "$WORK/t2_chr2.txt" true 2
rc2=$(run "$WORK/t2c2.yaml" "$WORK/t2_chr2.log")

if [ "$rc1" != 0 ]; then bad "LOCO chrom=1 exited $rc1"; tail -3 "$WORK/t2_chr1.log"; fi
if [ "$rc2" != 0 ]; then bad "LOCO chrom=2 exited $rc2"; tail -3 "$WORK/t2_chr2.log"; fi

if [ "$rc1" = 0 ]; then
    # Compare p-values on the chr1 markers common to the LOCO and non-LOCO runs.
    if python3 "$HERE/cmp_pval.py" differ "$WORK/t1_new.txt" "$WORK/t2_chr1.txt"; then
        ok "LOCO chrom=1 results differ from LOCO=off (chr1/ was actually read)"
    else
        bad "LOCO chrom=1 results are IDENTICAL to LOCO=off -- chr1/ was not read"
    fi
fi
if [ "$rc2" = 0 ]; then
    if python3 "$HERE/cmp_pval.py" differ "$WORK/t1_new.txt" "$WORK/t2_chr2.txt"; then
        ok "LOCO chrom=2 results differ from LOCO=off (chr2/ was actually read)"
    else
        bad "LOCO chrom=2 results are IDENTICAL to LOCO=off -- chr2/ was not read"
    fi
fi
if [ "$rc1" = 0 ] && [ "$rc2" = 0 ]; then
    # chr1 and chr2 outputs cover disjoint marker sets, so compare the models
    # instead: rerun chrom=2's model against chr1's marker set is not possible
    # without recompiling, so we assert the two runs used different mu vectors
    # via the logged per-chromosome directory.
    d1=$(grep -c "chr1$" "$WORK/t2_chr1.log")
    d2=$(grep -c "chr2$" "$WORK/t2_chr2.log")
    if [ "$d1" -ge 1 ] && [ "$d2" -ge 1 ]; then
        ok "per-chromosome model directories chr1/ and chr2/ were selected"
    else
        bad "per-chromosome directory selection not logged (d1=$d1 d2=$d2)"
    fi
fi

##############################################################################
echo
echo "=== Test 3: marker filtering restricts to chrom ==="
if [ "$rc1" = 0 ]; then
    if python3 "$HERE/cmp_pval.py" onlychr 1 "$WORK/t2_chr1.txt"; then
        ok "chrom=1 output contains no chromosome-2 markers"
    else
        bad "chrom=1 output contains markers from another chromosome"
    fi
fi
if [ "$rc2" = 0 ]; then
    if python3 "$HERE/cmp_pval.py" onlychr 2 "$WORK/t2_chr2.txt"; then
        ok "chrom=2 output contains only chromosome-2 markers"
    else
        bad "chrom=2 output contains markers from another chromosome"
    fi
fi
# "chr1"-style label must normalize to the bare "1" used in the .bim.
mkcfg "$WORK/t3.yaml" "$WORK/model_loco" "$WORK/t3.txt" true chr1
rc=$(run "$WORK/t3.yaml" "$WORK/t3.log")
if [ "$rc" = 0 ] && cmp -s "$WORK/t2_chr1.txt" "$WORK/t3.txt"; then
    ok "chrom='chr1' normalizes to the .bim's '1'"
else
    bad "chrom='chr1' did not normalize (rc=$rc)"
fi

##############################################################################
echo
echo "=== Test 4: guards ==="
# Guard 1: LOCO true, model loco false -> error, non-zero exit.
mkcfg "$WORK/t4a.yaml" "$WORK/model_noloco" "$WORK/t4a.txt" true 1
rc=$(run "$WORK/t4a.yaml" "$WORK/t4a.log")
if [ "$rc" != 0 ] && grep -qi "does not contain LOCO results" "$WORK/t4a.log"; then
    ok "guard 1: LOCO=true on a non-LOCO model errors (exit $rc)"
else
    bad "guard 1: expected non-zero exit with an explicit message (rc=$rc)"
fi

# Guard 2: LOCO true, chrom empty -> error, non-zero exit.
mkcfg "$WORK/t4b.yaml" "$WORK/model_loco" "$WORK/t4b.txt" true ""
rc=$(run "$WORK/t4b.yaml" "$WORK/t4b.log")
if [ "$rc" != 0 ] && grep -qi "chrom needs to be specified" "$WORK/t4b.log"; then
    ok "guard 2: LOCO=true with empty chrom errors (exit $rc)"
else
    bad "guard 2: expected non-zero exit with an explicit message (rc=$rc)"
fi

# Guard 3: chrom not in loco_chroms -> SILENT fallback, exit 0, full-genome fit.
# Chromosome 2 exists in the .bim, so use a model that only has chr1/ so the
# run still has markers to test.
cp -r "$WORK/model_loco" "$WORK/model_chr1only"
rm -rf "$WORK/model_chr1only/chr2"
python3 - "$WORK/model_chr1only/nullmodel.json" <<'PY'
import sys
p = sys.argv[1]
s = open(p).read().replace('"loco_chroms": [1, 2]', '"loco_chroms": [1]')
open(p, 'w').write(s)
PY
mkcfg "$WORK/t4c.yaml" "$WORK/model_chr1only" "$WORK/t4c.txt" true 2
rc=$(run "$WORK/t4c.yaml" "$WORK/t4c.log")
if [ "$rc" = 0 ] && ! grep -qiE "^(error|terminate|what\(\))" "$WORK/t4c.log"; then
    # Must have fallen back to the full-genome fit: results on the chr2 markers
    # must equal the LOCO=off run, and must NOT equal the chr2-LOCO run.
    if python3 "$HERE/cmp_pval.py" same "$WORK/t1_new.txt" "$WORK/t4c.txt" && \
       python3 "$HERE/cmp_pval.py" differ "$WORK/t2_chr2.txt" "$WORK/t4c.txt"; then
        ok "guard 3: chrom absent from loco_chroms silently falls back (exit 0)"
    else
        bad "guard 3: exited 0 but did not use the full-genome fit"
    fi
else
    bad "guard 3: expected a silent exit 0, got rc=$rc"
fi

# Guard 3b: a genuine non-autosome label (X).  No X markers exist in the .bim,
# so the run legitimately has nothing to test; assert only that the loader did
# not raise the LOCO guards.
mkcfg "$WORK/t4d.yaml" "$WORK/model_loco" "$WORK/t4d.txt" true X
rc=$(run "$WORK/t4d.yaml" "$WORK/t4d.log")
if grep -q "using the full-genome fit" "$WORK/t4d.log" && \
   ! grep -qi "does not contain LOCO results" "$WORK/t4d.log"; then
    ok "guard 3b: non-autosome 'X' falls back to the full-genome fit"
else
    bad "guard 3b: 'X' did not take the silent fallback path"
fi

##############################################################################
echo
echo "=== Test 5: region tests work under LOCO ==="
# group.txt holds GENE1 (all chr1 markers) and GENE2 (all chr2 markers).
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group.txt"
region_cfg() {  # region_cfg <cfg> <model> <out> <LOCO> <chrom> [groupfile]
    local grp="${6:-$WORK/group.txt}"
    mkcfg "$1" "$2" "$3" "$4" "$5" \
        "groupFile:   $grp" \
        "annotationList:" \
        "  - \"lof\"" \
        "maxMAFList:" \
        "  - 0.5" \
        "r_corr: 0" \
        "MACCutoff_to_CollapseUltraRare: 10" \
        "is_output_moreDetails: true"
}
region_cfg "$WORK/t5off.yaml" "$WORK/model_noloco" "$WORK/t5off.txt" false ""
rco=$(run "$WORK/t5off.yaml" "$WORK/t5off.log")
region_cfg "$WORK/t5on.yaml" "$WORK/model_loco" "$WORK/t5on.txt" true 1
rcn=$(run "$WORK/t5on.yaml" "$WORK/t5on.log")
if [ "$rco" = 0 ] && [ "$rcn" = 0 ]; then
    if [ -s "$WORK/t5on.txt" ] && ! cmp -s "$WORK/t5off.txt" "$WORK/t5on.txt"; then
        ok "region test runs under LOCO and picks up the chr1/ model"
    elif [ -s "$WORK/t5on.txt" ]; then
        bad "region test under LOCO produced identical output to LOCO=off"
    else
        bad "region test under LOCO produced no output"
    fi
else
    bad "region test exited non-zero (off=$rco on=$rcn)"
    tail -3 "$WORK/t5on.log"
fi

##############################################################################
echo
echo "=== Test 6: LOCO restricts REGIONS to the LOCO chromosome ==="
# The load-bearing test. group.txt spans chr1 and chr2; under LOCO chrom=1 the
# run must (a) emit GENE1 and not GENE2, and (b) be bit-for-bit the same as a
# run over a chr1-ONLY group file. (b) is what distinguishes "filtered" from
# "filtered but corrupted the region assembly / ultra-rare collapsing chunk".
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_chr1.txt" chr1
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_chr2.txt" chr2
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_mixed.txt" mixed

# 6a: only chr1 genes in the output of the both-chromosomes group file.
if [ "$rcn" = 0 ]; then
    if grep -q "^GENE1[[:space:]]" "$WORK/t5on.txt" && \
       ! grep -q "^GENE2[[:space:]]" "$WORK/t5on.txt"; then
        ok "6a: LOCO chrom=1 emits GENE1 and drops the chr2 gene GENE2"
    else
        bad "6a: expected GENE1 only in $WORK/t5on.txt"
        cut -f1 "$WORK/t5on.txt" | head -5
    fi
    if grep -q "not on LOCO chromosome 1" "$WORK/t5on.log" && \
       grep -q "regions skipped as not on chromosome 1: 1 of 2" "$WORK/t5on.log"; then
        ok "6b: the skipped region and the count are logged to stdout"
    else
        bad "6b: off-chromosome region skip was not clearly logged"
    fi
fi

# 6c: identical to a chr1-only group file -- both the region table and the
# per-variant table (the latter also covers the ultra-rare pseudo-marker chunk).
region_cfg "$WORK/t6c.yaml" "$WORK/model_loco" "$WORK/t6c.txt" true 1 "$WORK/group_chr1.txt"
rc=$(run "$WORK/t6c.yaml" "$WORK/t6c.log")
if [ "$rc" = 0 ] && [ "$rcn" = 0 ] && cmp -s "$WORK/t5on.txt" "$WORK/t6c.txt"; then
    ok "6c: both-chromosome group file == chr1-only group file (region p-values)"
else
    bad "6c: both-chromosome group file differs from chr1-only group file (rc=$rc)"
    diff "$WORK/t6c.txt" "$WORK/t5on.txt" | head -5
fi
if [ "$rc" = 0 ] && [ "$rcn" = 0 ] && \
   cmp -s "$WORK/t5on.txt.singleAssoc.txt" "$WORK/t6c.txt.singleAssoc.txt"; then
    ok "6d: single-variant-in-group output also identical (region assembly intact)"
else
    bad "6d: singleAssoc output differs between the two group files"
    diff "$WORK/t6c.txt.singleAssoc.txt" "$WORK/t5on.txt.singleAssoc.txt" | head -5
fi

# 6e: a region whose variants span two chromosomes must be a HARD ERROR under
# LOCO. Silently dropping half a gene would change the burden statistic with no
# signal in the output (R does exactly that; we refuse instead).
region_cfg "$WORK/t6e.yaml" "$WORK/model_loco" "$WORK/t6e.txt" true 1 "$WORK/group_mixed.txt"
rc=$(run "$WORK/t6e.yaml" "$WORK/t6e.log")
if [ "$rc" != 0 ] && grep -qi "spans more than one chromosome" "$WORK/t6e.log"; then
    ok "6e: chromosome-spanning region errors under LOCO (exit $rc)"
else
    bad "6e: expected non-zero exit for a chromosome-spanning region (rc=$rc)"
fi

# 6f: the same mixed group file with LOCO off must still run -- the check is
# LOCO-only and must not break existing non-LOCO region analyses.
region_cfg "$WORK/t6f.yaml" "$WORK/model_noloco" "$WORK/t6f.txt" false "" "$WORK/group_mixed.txt"
rc=$(run "$WORK/t6f.yaml" "$WORK/t6f.log")
if [ "$rc" = 0 ] && grep -q "^GENEMIX[[:space:]]" "$WORK/t6f.txt"; then
    ok "6f: LOCO=off still analyses a chromosome-spanning region"
else
    bad "6f: LOCO=off run over the mixed group file failed (rc=$rc)"
fi

# 6g: every region off-chromosome -> non-zero exit rather than an empty table.
region_cfg "$WORK/t6g.yaml" "$WORK/model_loco" "$WORK/t6g.txt" true 1 "$WORK/group_chr2.txt"
rc=$(run "$WORK/t6g.yaml" "$WORK/t6g.log")
if [ "$rc" != 0 ] && grep -qi "No regions on chrom 1 are found" "$WORK/t6g.log"; then
    ok "6g: a group file with nothing on the LOCO chromosome errors (exit $rc)"
else
    bad "6g: expected non-zero exit when no region is on the LOCO chromosome (rc=$rc)"
fi

# 6h: region test with LOCO off is byte-identical to the pre-LOCO build.
if [ -n "$BASELINE_BIN" ] && [ -x "$BASELINE_BIN" ]; then
    grep -v -e '^LOCO:' -e '^chrom:' "$WORK/t5off.yaml" | \
        sed "s#$WORK/t5off.txt#$WORK/t6h_base.txt#" > "$WORK/t6h_base.yaml"
    rcb=$("$BASELINE_BIN" "$WORK/t6h_base.yaml" > "$WORK/t6h_base.log" 2>&1; echo $?)
    if [ "$rcb" = 0 ] && [ "$rco" = 0 ] && \
       cmp -s "$WORK/t5off.txt" "$WORK/t6h_base.txt" && \
       cmp -s "$WORK/t5off.txt.singleAssoc.txt" "$WORK/t6h_base.txt.singleAssoc.txt"; then
        ok "6h: LOCO=off region output byte-identical to the pre-LOCO build"
    else
        bad "6h: LOCO=off region output DIFFERS from the pre-LOCO build (rc=$rcb)"
        diff "$WORK/t6h_base.txt" "$WORK/t5off.txt" | head -5
    fi
else
    echo "  (no baseline binary given -- skipping region byte-identity comparison)"
fi

##############################################################################
echo
echo "======================================"
echo "  LOCO tests: $PASS passed, $FAIL failed"
echo "  artifacts in $WORK"
echo "======================================"
[ "$FAIL" = 0 ]
