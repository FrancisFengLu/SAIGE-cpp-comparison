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
#       out of the batch path and must still match its own golden
#   C6  mixed, marker_chunksize smaller than the marker count -> exercises the
#       chunked writer
#   C7  C3 again at nThreads=8
#   C8  mixed p (1 / 2 / 3) in one run -- the tl_X1 / tl_A1 grow test
#   C9  C8 at nThreads=8, where a wrong grow test corrupts the heap
#   CA  mtBatch: false -- every pair down the scalar path (the A/B reference)
#   CB  mtBlockSize 32      CC  mtBlockSize 512      CD  mtBlockSize 7 + 8 threads
#   CE  a quantitative trait whose p-values all underflow, so the batch kernel's
#       fixed-point "%.1fE%d" branch is exercised on non-fallback pairs
#   CF  LOCO with the models disagreeing about whether it applied
#
# C7 / CA / CB / CC must additionally be byte-identical to C3: neither the
# thread count, nor the batch kernel, nor the block width may move a byte.

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
# Quantitative trait with scaled residuals (tests/make_extreme_model.py): every
# marker underflows to p == 0, so the batch kernel's "%.1fE%d" branch is
# exercised on pairs that are NOT routed to the fallback (design section 9.4).
XTRAITS=""
[ -e "$ARMA/q2x/nullmodel.json" ] && XTRAITS="q2x"

vrfile () {  # vrfile <trait>
    case "$1" in
        b3p2) echo "$VR_B/w0_y3.varianceRatio.txt" ;;
        q3p1) echo "$VR_Q/q3.varianceRatio.txt" ;;
        q2x)  echo "$VR_Q/q2.varianceRatio.txt" ;;
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
ok ()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad () { echo "  FAIL  $1"; fail=$((fail+1)); }

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
for t in $BTRAITS $QTRAITS $PTRAITS $XTRAITS; do golden_one "$t" plain  false false; done
for t in $BTRAITS $QTRAITS;          do golden_one "$t" more   true  true;  done
# C5's odd trait out: q1 with isnoadjCov
golden_one q1 noadj false false "isnoadjCov: true"

# ---- multi-trait configs --------------------------------------------------
MT_EXTRA=""   # extra top-level yaml lines for the next mt_cfg call
mt_cfg () {  # mt_cfg <name> <nThreads> <moreOutput> <firth> <chunksize> <tag> <traits...>
    local name="$1" thr="$2" mo="$3" fi="$4" chunk="$5" tag="$6"; shift 6
    local cfg="$WORK/cfg/$name.yaml"
    { common "$thr" "$mo" "$fi" "$chunk"
      [ -n "$MT_EXTRA" ] && printf '%s\n' "$MT_EXTRA"
      echo "models:"; } > "$cfg"
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
            # cmp is the gate; the field-by-field report is only here to say
            # WHAT moved, so a last-bit rounding artefact of the printed
            # 6-significant-figure format can be told apart from a real change.
            echo "  FAIL $name/$t differs from golden/${tag}_$t.txt"
            python3 "$(dirname "${BASH_SOURCE[0]}")/cmp_assoc.py" \
                "$WORK/golden/${tag}_$t.txt" "$WORK/out/$name.$t.txt" 2>&1 |
                sed 's/^/      /' | head -20
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

if [ -n "$XTRAITS" ]; then
    # p == 0 on a quantitative trait: the fixed-point "%.1fE%d" branch, reached
    # through the batch kernel (quantitative pairs never fall back).
    mt_cfg CE_extreme 1 false false 10000 plain $BTRAITS $QTRAITS $XTRAITS
    n0=$(awk -F'\t' 'NR>1 && $13 ~ /^[0-9][.][0-9]E-/ {c++} END{print c+0}' \
             "$WORK/out/CE_extreme.q2x.txt" 2>/dev/null || echo 0)
    if [ "${n0:-0}" -gt 0 ]; then
        ok "CE_extreme: $n0 q2x rows took the log-scale p-value branch"
    else
        bad "CE_extreme: no q2x row reached the log-scale p-value branch -- the case is not being tested"
    fi
fi

# Batch kernel off: every pair takes the scalar path. This is the Phase 1
# behaviour and the A/B reference for the kernel -- it must match the goldens
# too, otherwise a kernel failure and a plumbing failure look the same.
MT_EXTRA="mtBatch: false"
mt_cfg CA_nobatch 1 false false 10000 plain $BTRAITS $QTRAITS
MT_EXTRA=""
# Block width may not move a byte: it only changes where the block boundaries
# fall, never any pair's arithmetic.
MT_EXTRA="mtBlockSize: 32"
mt_cfg CB_blk32  1 false false 10000 plain $BTRAITS $QTRAITS
MT_EXTRA="mtBlockSize: 512"
mt_cfg CC_blk512 1 false false 10000 plain $BTRAITS $QTRAITS
MT_EXTRA="mtBlockSize: 7"
mt_cfg CD_blk7   8 true  true  1000  more  $BTRAITS $QTRAITS
MT_EXTRA=""

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

# Cross-checks: none of these knobs may move a byte.
xcheck () {  # xcheck <label> <run A> <run B>
    for t in $BTRAITS $QTRAITS; do
        if cmp -s "$WORK/out/$2.$t.txt" "$WORK/out/$3.$t.txt"; then pass=$((pass+1))
        else echo "  FAIL $1 for $t ($2 != $3)"; fail=$((fail+1)); fi
    done
}
xcheck "nThreads=8 != nThreads=1" C7_t8     C3_mix8
xcheck "mtBatch=false != batched" CA_nobatch C3_mix8
xcheck "mtBlockSize=32 moved"     CB_blk32  C3_mix8
xcheck "mtBlockSize=512 moved"    CC_blk512 C3_mix8

# ---------------------------------------------------------------------------
# LOCO, including the mixed state (design section 4.6)
# ---------------------------------------------------------------------------
# loco_ab / loco_aq carry chr1/ and chr2/ and list both, so chrom=1 really
# swaps their per-chromosome files in. loco_gw lists only chr2, so chrom=1 hits
# the loader's silent fallback to the genome-wide fit -- which is legal but has
# to be announced when the P models disagree about it.
HELPER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/make_loco_model.py"
# Same construction as tests/run_p1_regression.sh: copy the model, write a
# perturbed chr<N>/ set, and declare which chromosomes it holds. The
# perturbation is what makes the test meaningful -- without it a loader that
# silently read the top-level files would still look correct.
mk_loco () {  # mk_loco <src model> <dst dir> <chrom> <factor> <loco_chroms json>
    python3 - "$1" "$2" "$3" "$4" "$5" "$HELPER" <<'PY2'
import importlib.util, json, os, shutil, sys
src, dst, chrom, factor, chroms, helper = sys.argv[1:7]
spec = importlib.util.spec_from_file_location("mlm", helper)
mlm = importlib.util.module_from_spec(spec); spec.loader.exec_module(mlm)
if os.path.exists(dst): shutil.rmtree(dst)
shutil.copytree(src, dst)
cdir = os.path.join(dst, "chr" + chrom); os.makedirs(cdir, exist_ok=True)
for name in mlm.PER_CHROM:
    header, rows, cols, vals = mlm.read_arma(os.path.join(src, name + ".arma"))
    mlm.write_arma(os.path.join(cdir, name + ".arma"), header, rows, cols,
                   mlm.perturb(name, vals, float(factor)))
j = json.load(open(os.path.join(dst, "nullmodel.json")))
j["loco"] = True; j["loco_chroms"] = json.loads(chroms)
json.dump(j, open(os.path.join(dst, "nullmodel.json"), "w"))
PY2
}
if [ -e "$HELPER" ]; then
    # loco_ab / loco_aq hold chr1 and declare it, so chrom=1 really swaps files
    # in. loco_gw holds only chr2, so chrom=1 hits the loader's silent fallback
    # to the genome-wide fit -- legal, but it has to be announced when the P
    # models disagree about it.
    mk_loco "$ARMA/b1" "$WORK/loco_ab" 1 0.90 '[1]' &&
    mk_loco "$ARMA/q1" "$WORK/loco_aq" 1 0.85 '[1]' &&
    mk_loco "$ARMA/b2" "$WORK/loco_gw" 2 1.15 '[2]' || {
        echo "  FAIL could not build LOCO models"; fail=$((fail+1)); }

    loco_vr () { case "$1" in loco_ab) vrfile b1;; loco_aq) vrfile q1;; loco_gw) vrfile b2;; esac; }
    for t in loco_ab loco_aq loco_gw; do
        cfg="$WORK/cfg/golden_loco_$t.yaml"
        { echo "modelFile:         $WORK/$t"
          echo "varianceRatioFile: $(loco_vr $t)"
          echo "outputFile:        $WORK/golden/loco_$t.txt"
          echo "LOCO: true"
          echo "chrom: \"1\""
          common 1 false false 10000; } > "$cfg"
        "$BIN" "$cfg" > "$WORK/golden/loco_$t.log" 2>&1 ||
            { echo "  FAIL golden loco/$t exited $?"; fail=$((fail+1)); }
    done
    cfg="$WORK/cfg/CF_loco.yaml"
    { common 1 false false 10000
      echo "LOCO: true"
      echo "chrom: \"1\""
      echo "models:"
      for t in loco_ab loco_aq loco_gw; do
        echo "  - traitName: $t"
        echo "    modelFile:         $WORK/$t"
        echo "    varianceRatioFile: $(loco_vr $t)"
        echo "    outputFile: $WORK/out/CF_loco.$t.txt"
      done; } > "$cfg"
    if ! "$BIN" "$cfg" > "$WORK/out/CF_loco.log" 2>&1; then
        echo "  FAIL CF_loco exited non-zero"; sed -n '$p' "$WORK/out/CF_loco.log"
        fail=$((fail+1))
    else
        for t in loco_ab loco_aq loco_gw; do
            if cmp -s "$WORK/out/CF_loco.$t.txt" "$WORK/golden/loco_$t.txt"; then pass=$((pass+1))
            else echo "  FAIL CF_loco/$t differs from its single-trait LOCO golden"; fail=$((fail+1)); fi
        done
        # b1's chr1 fit must differ from b1 without LOCO, or the chr<N>/ files
        # were never read and the whole case proves nothing.
        if cmp -s "$WORK/golden/loco_loco_ab.txt" "$WORK/golden/plain_b1.txt" 2>/dev/null ||
           cmp -s "$WORK/golden/loco_ab.txt" "$WORK/golden/plain_b1.txt"; then
            bad "CF_loco: chr1 fit == non-LOCO fit; the chr1/ files were not read"
        else
            ok "CF_loco: chr1 fit differs from the non-LOCO fit"
        fi
        if grep -q "fell back to the genome-wide fit" "$WORK/out/CF_loco.log"; then
            ok "CF_loco: mixed-LOCO warning printed"
        else
            bad "CF_loco: mixed LOCO state was NOT announced"
        fi
    fi
    echo "  ran CF_loco (3 traits, mixed LOCO)"
fi

echo
echo "=== $pass passed, $fail failed  (workdir $WORK) ==="
[ $fail -eq 0 ]
