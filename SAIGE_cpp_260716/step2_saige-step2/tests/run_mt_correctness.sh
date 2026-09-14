#!/bin/bash
# Multi-trait correctness gate (MULTITRAIT_DESIGN.md section 11).
#
# The judgement is always the same one: for every trait in a P > 1 run, the
# output file must be BYTE-IDENTICAL to the file the SAME binary produces when
# that model is run alone through the single-trait path. The single-trait path
# is in turn pinned to the pre-change binary by tests/run_p1_regression.sh, so
# this transitively pins multi-trait output to the original implementation.
#
# usage: run_mt_correctness.sh <saige-step2 binary> [workdir]
#
# env overrides:
#   MT_PLINK   default /opt/saige/data/mid2k    (50k samples x 2k markers)
#   MT_ARMA    default /opt/saige/logs/step2mt/arma
#   MT_VR_B    default /opt/saige/logs/step2mt/models64      (w0_<name>.varianceRatio.txt)
#   MT_VR_Q    default /opt/saige/logs/step2mt/models_q8     (<name>.varianceRatio.txt)
#   MT_THREADS default 1
#
# Cases:
#   C1  4 binary traits
#   C2  4 quantitative traits
#   C3  mixed 4 binary + 4 quantitative (also checks the binary-first internal
#       reordering does not leak into any output)
#   C4  mixed, isMoreOutput + Firth on
#   C5  mixed, per-model override (one trait isnoadjCov) -> that trait is gated
#       out of any batch path and must still match its own golden
#   C6  mixed, marker_chunksize smaller than the marker count -> exercises the
#       chunked writer; output must not change
#   C7  C3 again at nThreads=8 -> must be byte-identical to the 1-thread run

set -u

BIN="${1:?usage: run_mt_correctness.sh <saige-step2 binary> [workdir]}"
WORK="${2:-$(mktemp -d)}"
PLINK="${MT_PLINK:-/opt/saige/data/mid2k}"
ARMA="${MT_ARMA:-/opt/saige/logs/step2mt/arma}"
VR_B="${MT_VR_B:-/opt/saige/logs/step2mt/models64}"
VR_Q="${MT_VR_Q:-/opt/saige/logs/step2mt/models_q8}"
THREADS="${MT_THREADS:-1}"

BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
[ -x "$BIN" ] || { echo "FAIL: binary not found: $BIN"; exit 1; }

BTRAITS="b1 b2 b3 b4"
QTRAITS="q1 q2 q3 q4"
# Reduced-covariate variants (tests/make_reduced_p_model.py). Their whole job is
# to make p differ between traits inside one run, which is what the tl_X1 /
# tl_A1 grow test in scoreTestFast has to survive (design section 9.1). Skipped
# silently when they have not been generated.
PTRAITS=""
for cand in b3p2 q3p1; do
    [ -e "$ARMA/$cand/nullmodel.json" ] && PTRAITS="$PTRAITS $cand"
done

vrfile () {  # vrfile <trait>
    case "$1" in
        b3p2) echo "$VR_B/w0_y3.varianceRatio.txt" ;;
        q3p1) echo "$VR_Q/q3.varianceRatio.txt" ;;
        b*)   echo "$VR_B/w0_y${1#b}.varianceRatio.txt" ;;
        q*)   echo "$VR_Q/$1.varianceRatio.txt" ;;
    esac
}

for t in $BTRAITS $QTRAITS; do
    [ -e "$ARMA/$t/nullmodel.json" ] || { echo "SKIP: no model $ARMA/$t"; exit 0; }
    [ -e "$(vrfile "$t")" ]          || { echo "SKIP: no VR $(vrfile "$t")"; exit 0; }
done
[ -e "$PLINK.bed" ] || { echo "SKIP: no genotypes $PLINK.bed"; exit 0; }

mkdir -p "$WORK/cfg" "$WORK/golden" "$WORK/out"
pass=0; fail=0

common () {  # common <nThreads> <moreOutput> <firth> <chunksize>
    cat <<EOF
genoType:   plink
plinkFile:  $PLINK
minMAF: 0
minMAC: 1
maxMissRate: 0.15
AlleleOrder: alt-first
MACCutoffforER: 4
nThreads: $1
isMoreOutput: $2
isFirth: $3
marker_chunksize: $4
EOF
}

# ---- golden: every trait on its own, through the single-trait path ---------
# Two variants, because isMoreOutput / isFirth change the columns written.
golden_one () {  # golden_one <trait> <tag> <moreOutput> <firth> <extra yaml>
    local t="$1" tag="$2" mo="$3" fi="$4" extra="${5:-}"
    local cfg="$WORK/cfg/golden_${tag}_$t.yaml"
    {
        echo "modelFile:         $ARMA/$t"
        echo "varianceRatioFile: $(vrfile "$t")"
        echo "outputFile:        $WORK/golden/${tag}_$t.txt"
        common 1 "$mo" "$fi" 10000
        [ -n "$extra" ] && echo "$extra"
    } > "$cfg"
    "$BIN" "$cfg" > "$WORK/golden/${tag}_$t.log" 2>&1 ||
        { echo "  FAIL golden $tag/$t exited $?"; fail=$((fail+1)); }
}

echo "=== building goldens (single-trait path) ==="
for t in $BTRAITS $QTRAITS $PTRAITS; do golden_one "$t" plain  false false; done
for t in $BTRAITS $QTRAITS;          do golden_one "$t" more   true  true;  done
# C5's odd trait out: q1 with isnoadjCov
golden_one q1 noadj false false "isnoadjCov: true"

# ---- multi-trait configs --------------------------------------------------
mt_cfg () {  # mt_cfg <name> <nThreads> <moreOutput> <firth> <chunksize> <tag> <traits...>
    local name="$1" thr="$2" mo="$3" fi="$4" chunk="$5" tag="$6"; shift 6
    local cfg="$WORK/cfg/$name.yaml"
    { common "$thr" "$mo" "$fi" "$chunk"; echo "models:"; } > "$cfg"
    for t in "$@"; do
        cat >> "$cfg" <<EOF
  - traitName: $t
    modelFile:         $ARMA/$t
    varianceRatioFile: $(vrfile "$t")
    outputFile: $WORK/out/$name.$t.txt
EOF
    done
    "$BIN" "$cfg" > "$WORK/out/$name.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  FAIL $name exited $rc"; sed -n '$p' "$WORK/out/$name.log" | sed 's/^/        /'
        fail=$((fail+1)); return
    fi
    for t in "$@"; do
        if cmp -s "$WORK/out/$name.$t.txt" "$WORK/golden/$tag.$t.txt" 2>/dev/null ||
           cmp -s "$WORK/out/$name.$t.txt" "$WORK/golden/${tag}_$t.txt"; then
            pass=$((pass+1))
        else
            echo "  FAIL $name/$t differs from golden/${tag}_$t.txt"
            diff <(head -5 "$WORK/golden/${tag}_$t.txt") <(head -5 "$WORK/out/$name.$t.txt") |
                sed 's/^/        /' | head -12
            fail=$((fail+1))
        fi
    done
    echo "  ran $name ($# traits)"
}

echo "=== multi-trait runs ==="
mt_cfg C1_bin4   1 false false 10000 plain $BTRAITS
mt_cfg C2_qnt4   1 false false 10000 plain $QTRAITS
mt_cfg C3_mix8   1 false false 10000 plain $BTRAITS $QTRAITS
mt_cfg C4_more8  1 true  true  10000 more  $BTRAITS $QTRAITS
mt_cfg C6_chunk  1 false false   256 plain $BTRAITS $QTRAITS
mt_cfg C7_t8     8 false false 10000 plain $BTRAITS $QTRAITS
if [ -n "$PTRAITS" ]; then
    # Different p in one run, single-threaded and then multi-threaded: the
    # multi-threaded one is the case that actually corrupts memory if the
    # tl_X1 / tl_A1 grow test regresses to n_rows.
    mt_cfg C8_pmix   1 false false 10000 plain $BTRAITS $QTRAITS $PTRAITS
    mt_cfg C9_pmix8  8 false false 10000 plain $BTRAITS $QTRAITS $PTRAITS
fi

# C5: per-model override -- q1 gets isnoadjCov, everyone else does not.
{
    common 1 false false 10000
    echo "models:"
    for t in $BTRAITS $QTRAITS; do
        cat <<EOF
  - traitName: $t
    modelFile:         $ARMA/$t
    varianceRatioFile: $(vrfile "$t")
    outputFile: $WORK/out/C5_ov.$t.txt
EOF
        [ "$t" = q1 ] && echo "    isnoadjCov: true"
    done
} > "$WORK/cfg/C5_ov.yaml"
"$BIN" "$WORK/cfg/C5_ov.yaml" > "$WORK/out/C5_ov.log" 2>&1 || {
    echo "  FAIL C5_ov exited $?"; sed -n '$p' "$WORK/out/C5_ov.log"; fail=$((fail+1)); }
for t in $BTRAITS $QTRAITS; do
    want="$WORK/golden/plain_$t.txt"
    [ "$t" = q1 ] && want="$WORK/golden/noadj_q1.txt"
    if cmp -s "$WORK/out/C5_ov.$t.txt" "$want"; then pass=$((pass+1))
    else echo "  FAIL C5_ov/$t differs from $(basename "$want")"; fail=$((fail+1)); fi
done
echo "  ran C5_ov (per-model isnoadjCov on q1)"

# C7 must also equal C3 exactly (thread count may not move a byte).
for t in $BTRAITS $QTRAITS; do
    if cmp -s "$WORK/out/C7_t8.$t.txt" "$WORK/out/C3_mix8.$t.txt"; then pass=$((pass+1))
    else echo "  FAIL nThreads=8 != nThreads=1 for $t"; fail=$((fail+1)); fi
done

echo
echo "=== $pass passed, $fail failed  (workdir $WORK) ==="
[ $fail -eq 0 ]
