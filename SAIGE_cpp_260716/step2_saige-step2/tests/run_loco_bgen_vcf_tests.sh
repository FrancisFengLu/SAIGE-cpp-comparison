#!/usr/bin/env bash
# LOCO chromosome-restriction tests for the BGEN and VCF genotype readers.
#
# Usage:  tests/run_loco_bgen_vcf_tests.sh [path/to/saige-step2] [path/to/baseline/saige-step2]
#
# Companion to run_loco_tests.sh (which covers PLINK). Test data is
# genotype_100markers_2chr.{bgen,vcf.gz}: 90 markers on chr1, 10 on chr2.
#
# What each test proves:
#   1  non-LOCO bgen/vcf output is byte-identical to the baseline build
#   2  LOCO chrom=N emits ONLY chromosome-N markers
#   3  the retained markers are the RIGHT ones (identical to the PLINK run's
#      marker set / p-values), i.e. we are not just reading the first K
#      variants sequentially
#   4  a bgen with no .bgi errors out under LOCO instead of warning

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
BGEN="$EXTDATA/genotype_100markers_2chr.bgen"
SAMPLE="$EXTDATA/genotype_100markers_2chr.sample"
VCF="$EXTDATA/genotype_100markers_2chr.vcf.gz"

WORK="${LOCO_TEST_WORK:-/tmp/saige_loco_bgen_vcf}"

PASS=0; FAIL=0
ok()  { echo "  PASS: $1"; PASS=$((PASS+1)); }
bad() { echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

for f in "$SRC_MODEL/nullmodel.json" "$VR" "$BGEN" "$BGEN.bgi" "$SAMPLE" "$VCF"; do
    if [ ! -f "$f" ]; then echo "SKIP: test data not found: $f"; exit 0; fi
done
if [ ! -x "$BIN" ]; then echo "FAIL: binary not found: $BIN"; exit 1; fi

rm -rf "$WORK"; mkdir -p "$WORK"
echo "=== Building synthetic LOCO model ==="
python3 "$HERE/make_loco_model.py" "$SRC_MODEL" "$WORK/model_loco" || exit 1
cp -r "$WORK/model_loco" "$WORK/model_noloco"
python3 - "$WORK/model_noloco/nullmodel.json" <<'PY'
import sys
p = sys.argv[1]
s = open(p).read().replace('"loco": true', '"loco": false')
open(p, 'w').write(s)
PY

# mkcfg <out.yaml> <fmt> <modeldir> <outputfile> <LOCO|omit> <chrom> [bgen_override]
mkcfg() {
    local cfg=$1 fmt=$2 model=$3 out=$4 loco=$5 chrom=$6 bgenov="${7:-$BGEN}"
    {
        echo "modelFile:         $model"
        echo "varianceRatioFile: $VR"
        echo "genoType:   $fmt"
        case "$fmt" in
            bgen) echo "bgenFile: $bgenov"; echo "bgenSampleFile: $SAMPLE";;
            vcf)  echo "vcfFile: $VCF"; echo "vcfField: GT";;
        esac
        echo "outputFile: $out"
        echo "minMAF: 0"
        echo "minMAC: 0.5"
        echo "maxMissRate: 0.15"
        echo "AlleleOrder: alt-first"
        echo "isMoreOutput: true"
        echo "isFirth: false"
        echo "nThreads: 1"
        if [ "$loco" != "omit" ]; then
            echo "LOCO: $loco"
            echo "chrom: \"$chrom\""
        fi
    } > "$cfg"
}

run() { "$BIN" "$1" > "$2" 2>&1; echo $?; }

##############################################################################
echo
echo "=== Test 1: non-LOCO output is byte-identical to the baseline build ==="
for fmt in bgen vcf; do
    mkcfg "$WORK/t1_$fmt.yaml" "$fmt" "$WORK/model_noloco" "$WORK/t1_$fmt.txt" false ""
    rc=$(run "$WORK/t1_$fmt.yaml" "$WORK/t1_$fmt.log")
    if [ "$rc" != 0 ]; then bad "$fmt LOCO=false run exited $rc"; tail -3 "$WORK/t1_$fmt.log"; continue; fi
    if [ -n "$BASELINE_BIN" ] && [ -x "$BASELINE_BIN" ]; then
        mkcfg "$WORK/t1_${fmt}_base.yaml" "$fmt" "$WORK/model_noloco" "$WORK/t1_${fmt}_base.txt" omit ""
        rcb=$("$BASELINE_BIN" "$WORK/t1_${fmt}_base.yaml" > "$WORK/t1_${fmt}_base.log" 2>&1; echo $?)
        if [ "$rcb" != 0 ]; then
            bad "$fmt baseline run exited $rcb"
        elif cmp -s "$WORK/t1_$fmt.txt" "$WORK/t1_${fmt}_base.txt"; then
            ok "$fmt LOCO=false output byte-identical to baseline build"
        else
            bad "$fmt LOCO=false output DIFFERS from baseline build"
            diff "$WORK/t1_${fmt}_base.txt" "$WORK/t1_$fmt.txt" | head -5
        fi
    else
        echo "  (no baseline binary given -- skipping $fmt byte-identity comparison)"
    fi
done

##############################################################################
echo
echo "=== Test 2: LOCO chrom=N emits only chromosome-N markers ==="
for fmt in bgen vcf; do
    for j in 1 2; do
        mkcfg "$WORK/t2_${fmt}_c$j.yaml" "$fmt" "$WORK/model_loco" "$WORK/t2_${fmt}_c$j.txt" true $j
        rc=$(run "$WORK/t2_${fmt}_c$j.yaml" "$WORK/t2_${fmt}_c$j.log")
        if [ "$rc" != 0 ]; then
            bad "$fmt LOCO chrom=$j exited $rc"; tail -3 "$WORK/t2_${fmt}_c$j.log"; continue
        fi
        if python3 "$HERE/cmp_pval.py" onlychr $j "$WORK/t2_${fmt}_c$j.txt"; then
            ok "$fmt chrom=$j output contains only chromosome-$j markers"
        else
            bad "$fmt chrom=$j output contains markers from another chromosome"
        fi
    done
done

##############################################################################
echo
echo "=== Test 3: the retained markers match the PLINK LOCO run exactly ==="
# The PLINK path is the already-verified reference implementation. bgen/vcf must
# select the SAME marker set -- this is what catches "read the first K variants
# sequentially" (which would keep chr1 markers under chrom=2).
PLINK="$EXTDATA/genotype_100markers_2chr"
for j in 1 2; do
    {
        echo "modelFile:         $WORK/model_loco"
        echo "varianceRatioFile: $VR"
        echo "genoType:   plink"
        echo "plinkFile:  $PLINK"
        echo "outputFile: $WORK/t3_plink_c$j.txt"
        echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"
        echo "AlleleOrder: alt-first"; echo "isMoreOutput: true"; echo "isFirth: false"
        echo "nThreads: 1"; echo "LOCO: true"; echo "chrom: \"$j\""
    } > "$WORK/t3_plink_c$j.yaml"
    run "$WORK/t3_plink_c$j.yaml" "$WORK/t3_plink_c$j.log" > /dev/null
done
for fmt in bgen vcf; do
    for j in 1 2; do
        python3 - "$WORK/t3_plink_c$j.txt" "$WORK/t2_${fmt}_c$j.txt" "$fmt" "$j" <<'PY'
import csv, sys
def rows(f):
    return {(r['CHR'], r['POS']): r for r in csv.DictReader(open(f), delimiter='\t')}
a, b, fmt, j = rows(sys.argv[1]), rows(sys.argv[2]), sys.argv[3], sys.argv[4]
if set(a) != set(b):
    print("    %s chr%s: marker set differs (plink=%d %s=%d, sym-diff=%d)"
          % (fmt, j, len(a), fmt, len(b), len(set(a) ^ set(b))))
    sys.exit(1)
worst = max(abs(float(a[k]['p.value']) - float(b[k]['p.value'])) for k in a)
print("    %s chr%s: %d markers, same set as plink, max |dp| = %.3g" % (fmt, j, len(a), worst))
sys.exit(0 if worst < 1e-8 else 2)
PY
        rc=$?
        if [ "$rc" = 0 ]; then
            ok "$fmt chrom=$j marker set and p-values match the PLINK LOCO run"
        else
            bad "$fmt chrom=$j does not match the PLINK LOCO run (rc=$rc)"
        fi
    done
done

##############################################################################
echo
echo "=== Test 4: bgen without a .bgi errors under LOCO (never warn-and-run) ==="
cp "$BGEN" "$WORK/nobgi.bgen"   # deliberately no .bgi alongside
mkcfg "$WORK/t4.yaml" bgen "$WORK/model_loco" "$WORK/t4.txt" true 1 "$WORK/nobgi.bgen"
rc=$(run "$WORK/t4.yaml" "$WORK/t4.log")
if [ "$rc" != 0 ] && grep -qi "requires a .bgi index" "$WORK/t4.log"; then
    ok "guard: LOCO bgen without .bgi errors out (exit $rc)"
else
    bad "guard: LOCO bgen without .bgi did not error (rc=$rc)"
fi
# ...and the same file WITHOUT LOCO still runs fine.
mkcfg "$WORK/t4b.yaml" bgen "$WORK/model_noloco" "$WORK/t4b.txt" false "" "$WORK/nobgi.bgen"
rc=$(run "$WORK/t4b.yaml" "$WORK/t4b.log")
if [ "$rc" = 0 ] && cmp -s "$WORK/t4b.txt" "$WORK/t1_bgen.txt"; then
    ok "bgen without .bgi still runs unchanged when LOCO is off"
else
    bad "bgen without .bgi regressed on the non-LOCO path (rc=$rc)"
fi

##############################################################################
echo
echo "=== Test 5: REGION tests are restricted to the LOCO chromosome ==="
# Group files are built from the PLINK .bim, whose chr:pos:ref:alt IDs are the
# same markers as in the .bgen/.vcf.gz.
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group.txt"
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_chr1.txt" chr1
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_mixed.txt" mixed
python3 "$HERE/make_group_file.py" "$PLINK.bim" "$WORK/group_chr1_nonstd.txt" nonstdchr1

# region_cfg <cfg> <fmt> <model> <out> <LOCO|omit> <chrom> <groupfile>
region_cfg() {
    mkcfg "$1" "$2" "$3" "$4" "$5" "$6"
    {
        echo "groupFile:   $7"
        echo "annotationList:"
        echo "  - \"lof\""
        echo "maxMAFList:"
        echo "  - 0.5"
        echo "r_corr: 0"
        echo "MACCutoff_to_CollapseUltraRare: 10"
        echo "is_output_moreDetails: true"
    } >> "$1"
}

# plink_region_cfg <cfg> <model> <out> <chrom> <groupfile>
plink_region_cfg() {
    {
        echo "modelFile:         $2"
        echo "varianceRatioFile: $VR"
        echo "genoType:   plink"
        echo "plinkFile:  $PLINK"
        echo "outputFile: $3"
        echo "minMAF: 0"; echo "minMAC: 0.5"; echo "maxMissRate: 0.15"
        echo "AlleleOrder: alt-first"; echo "isMoreOutput: true"; echo "isFirth: false"
        echo "nThreads: 1"; echo "LOCO: true"; echo "chrom: \"$4\""
        echo "groupFile:   $5"
        echo "annotationList:"; echo "  - \"lof\""
        echo "maxMAFList:"; echo "  - 0.5"
        echo "r_corr: 0"
        echo "MACCutoff_to_CollapseUltraRare: 10"
        echo "is_output_moreDetails: true"
    } > "$1"
}

# PLINK region reference: both-chromosome group file under LOCO chrom=1.
plink_region_cfg "$WORK/t5_plink_both.yaml" "$WORK/model_loco" \
    "$WORK/t5_plink_both.txt" 1 "$WORK/group.txt"
rcp=$(run "$WORK/t5_plink_both.yaml" "$WORK/t5_plink_both.log")
if [ "$rcp" != 0 ]; then
    bad "plink region LOCO reference run exited $rcp"; tail -3 "$WORK/t5_plink_both.log"
fi

for fmt in bgen vcf; do
    # both-chromosome group file under LOCO chrom=1
    region_cfg "$WORK/t5_${fmt}_both.yaml" "$fmt" "$WORK/model_loco" \
        "$WORK/t5_${fmt}_both.txt" true 1 "$WORK/group.txt"
    rcb=$(run "$WORK/t5_${fmt}_both.yaml" "$WORK/t5_${fmt}_both.log")
    # chr1-only group file under LOCO chrom=1
    region_cfg "$WORK/t5_${fmt}_c1.yaml" "$fmt" "$WORK/model_loco" \
        "$WORK/t5_${fmt}_c1.txt" true 1 "$WORK/group_chr1.txt"
    rc1=$(run "$WORK/t5_${fmt}_c1.yaml" "$WORK/t5_${fmt}_c1.log")

    # The marker-ID map must be non-empty for every format. (It used to be
    # empty for VCF -- VcfClass::getMarkerIDToIndex() read m_chr/m_pd/... which
    # only the streaming reader filled -- so every VCF region run silently
    # found 0 markers. prescanMarkerCount() now populates that metadata.)
    nmap=$(sed -n 's/^  Built map with \([0-9]*\) entries.*/\1/p' "$WORK/t5_${fmt}_c1.log" | head -1)
    if [ -n "$nmap" ] && [ "$nmap" -gt 0 ]; then
        ok "$fmt region marker-ID map is non-empty ($nmap entries)"
    else
        bad "$fmt region marker-ID map is empty (region tests would find 0 markers)"
    fi

    if [ "$rcb" = 0 ] && [ "$rc1" = 0 ]; then
        if grep -q "^GENE1[[:space:]]" "$WORK/t5_${fmt}_both.txt" && \
           ! grep -q "^GENE2[[:space:]]" "$WORK/t5_${fmt}_both.txt" && \
           cmp -s "$WORK/t5_${fmt}_both.txt" "$WORK/t5_${fmt}_c1.txt"; then
            ok "$fmt region LOCO chrom=1 drops the chr2 gene and matches a chr1-only group file"
        else
            bad "$fmt region LOCO chrom=1 did not match the chr1-only group file"
            diff "$WORK/t5_${fmt}_c1.txt" "$WORK/t5_${fmt}_both.txt" | head -5
        fi
    else
        bad "$fmt region LOCO run exited non-zero (both=$rcb chr1=$rc1)"
        tail -3 "$WORK/t5_${fmt}_both.log"
    fi

    # ...and the surviving chr1 gene must be the SAME result the PLINK path
    # produces from the same group file. This is the decisive check that the
    # marker-ID map points at the right variants (an off-by-one index space
    # would still produce plausible p-values).
    if [ "$rcb" = 0 ] && [ "$rcp" = 0 ]; then
        if cmp -s "$WORK/t5_${fmt}_both.txt" "$WORK/t5_plink_both.txt"; then
            ok "$fmt region LOCO chrom=1 output is cmp-identical to the PLINK run"
        elif [ "$fmt" = bgen ] && python3 - "$WORK/t5_plink_both.txt" "$WORK/t5_${fmt}_both.txt" <<'PY'
# BGEN stores 8-bit probabilities, so its dosages are not bit-identical to
# PLINK hard calls; require agreement to ~1e-12 relative instead of cmp.
import csv, sys
a = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
b = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
if [r['Region'] for r in a] != [r['Region'] for r in b]:
    print("    region sets differ"); sys.exit(1)
worst = 0.0
for x, y in zip(a, b):
    for k in ('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT', 'BETA_Burden', 'SE_Burden'):
        u, v = float(x[k]), float(y[k])
        worst = max(worst, abs(u - v) / max(abs(u), 1e-300))
print("    %d regions, worst relative deviation vs plink = %.3g" % (len(a), worst))
sys.exit(0 if worst < 1e-12 else 1)
PY
        then
            ok "$fmt region LOCO chrom=1 matches the PLINK run to 1e-12"
        else
            bad "$fmt region LOCO chrom=1 output differs from the PLINK run"
            diff "$WORK/t5_plink_both.txt" "$WORK/t5_${fmt}_both.txt" | head -5
        fi
    else
        bad "$fmt/plink region comparison skipped (both=$rcb plink=$rcp)"
    fi

    # Defect 2: variant IDs that do NOT parse as chr:pos:ref:alt must not be a
    # hard error. The chromosome now comes from the genotype file, so unknown
    # IDs are simply dropped (as R does) and the gene is scored on the rest.
    region_cfg "$WORK/t5_${fmt}_nonstd.yaml" "$fmt" "$WORK/model_loco" \
        "$WORK/t5_${fmt}_nonstd.txt" true 1 "$WORK/group_chr1_nonstd.txt"
    rcn=$(run "$WORK/t5_${fmt}_nonstd.yaml" "$WORK/t5_${fmt}_nonstd.log")
    if [ "$rcn" = 0 ] && \
       ! grep -q "chr:pos:ref:alt' form" "$WORK/t5_${fmt}_nonstd.log" && \
       grep -q "are not in 'GenoFile'" "$WORK/t5_${fmt}_nonstd.log" && \
       cmp -s "$WORK/t5_${fmt}_nonstd.txt" "$WORK/t5_${fmt}_c1.txt"; then
        ok "$fmt LOCO region with non-parseable variant IDs runs and matches the clean run"
    else
        bad "$fmt LOCO region with non-parseable variant IDs failed (rc=$rcn)"
        tail -3 "$WORK/t5_${fmt}_nonstd.log"
    fi

    # chromosome-spanning region -> hard error
    region_cfg "$WORK/t5_${fmt}_mix.yaml" "$fmt" "$WORK/model_loco" \
        "$WORK/t5_${fmt}_mix.txt" true 1 "$WORK/group_mixed.txt"
    rc=$(run "$WORK/t5_${fmt}_mix.yaml" "$WORK/t5_${fmt}_mix.log")
    if [ "$rc" != 0 ] && grep -qi "spans more than one chromosome" "$WORK/t5_${fmt}_mix.log"; then
        ok "$fmt chromosome-spanning region errors under LOCO (exit $rc)"
    else
        bad "$fmt chromosome-spanning region did not error (rc=$rc)"
    fi
done

##############################################################################
echo
echo "======================================"
echo "  LOCO bgen/vcf tests: $PASS passed, $FAIL failed"
echo "  artifacts in $WORK"
echo "======================================"
[ "$FAIL" -eq 0 ]
