#!/bin/bash
# Multi-trait with DIFFERENT sample sets per model (MULTITRAIT_DESIGN.md 4.7).
#
# Judgement, as in run_mt_correctness.sh: every trait's output file from a P > 1
# run must be BYTE-IDENTICAL to the file the same binary writes when that model
# runs alone (P = 1, single-trait path, nThreads 1). Any difference is printed
# field by field (cmp_assoc.py) and counts as a failure.
#
# usage: run_mt_subset_tests.sh <saige-step2 binary> [workdir]
#
# env:
#   MS_NULL   single-trait step 1 models, one fit per trait
#             (default /opt/saige/logs/missing_mt/step2/null):
#               block16/m/y1..y16  y1-y8 on everyone, y9-y16 on one 60% subset
#               indep16/m/y1..y16  each missing an independent 5% (all n = 47500)
#               qmiss/m/qm1..qm4   quantitative; qm1-3 independent 8%, qm4 a 70% block
#   MS_MID    common-variant PLINK set (default /opt/saige/data/mid2k)
#   MS_RARE   rare-variant PLINK set, same samples (default /opt/saige/data/rare)
#   MS_JOBS   concurrent golden runs (default 3)
#   MS_CASES  subset of cases to run (default: all)
#
# Cases (P, what it is there to catch):
#   D1  block16, 16  two sample sets; y1-y8 are exactly the union
#   D2  indep16, 16  16 sample sets, every n equal -- "same n => same samples" bugs
#   D3  7 interleaved binary/quantitative traits from three families, first model
#       not the full set (no trait equals the union, union not in .fam order)
#   D4  5 with sample-order permutations: block16 y1 and a permuted copy of it
#       (same set, different order), indep16 y4 and its permuted copy, a permuted
#       quantitative model
#   D5  6 with SPA_Cutoff 0.5 on two models: most binary pairs fall back to SPA
#   D6  6 on the rare set: ER (MAC <= 4), dosage-zeroing gate (MAC <= 10), missing
#   D7  D6 traits, minMAC 5: markers passing QC for one trait and not another
#   D8  D6 traits, MACCutoffforER 20
#   D9  D3 traits, isMoreOutput + Firth (pCutoffforFirth 0.05)
#   D10 D2 at nThreads 8        == D2 byte for byte
#   D11 D2 with mtBatch: false  == D2 byte for byte (all pairs scalar)
#   D12 D3 with mtBlockSize 7 and nThreads 4 == D3
#   D13 D6 under SAIGE_STEP2_SCALAR_DECODE=1 (goldens under the same env)
#   D14 D3 with a per-model isnoadjCov on one own-sample trait (not batchable)
#   D15 LOCO chrom 1, P=4: two models with a perturbed chr1/ fit (one full set,
#       one subset), one subset model whose loco_chroms lacks chr1 (silent
#       genome-wide fallback), one quantitative subset model with a chr1/ fit
#   E1-E3 refusals: mtRequireSameSamples: true, genoType bgen, conditional analysis

set -u

BIN="${1:?usage: run_mt_subset_tests.sh <saige-step2 binary> [workdir]}"
WORK="${2:-$(mktemp -d)}"
NULLD="${MS_NULL:-/opt/saige/logs/missing_mt/step2/null}"
MID="${MS_MID:-/opt/saige/data/mid2k}"
RARE="${MS_RARE:-/opt/saige/data/rare}"
JOBS="${MS_JOBS:-3}"
CASES="${MS_CASES:-D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 D13 D14 D15 E1 E2 E3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"
[ -x "$BIN" ] || { echo "FAIL: binary not found: $BIN"; exit 1; }
for f in "$NULLD/block16/m/y16/nullmodel.json" "$NULLD/indep16/m/y16/nullmodel.json" \
         "$NULLD/qmiss/m/qm4/nullmodel.json" "$MID.bed" "$RARE.bed"; do
    [ -e "$f" ] || { echo "SKIP: missing $f"; exit 0; }
done

mkdir -p "$WORK/cfg" "$WORK/golden" "$WORK/out" "$WORK/models"
pass=0; fail=0
ok ()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad () { echo "  FAIL  $1"; fail=$((fail+1)); }
has_case () { case " $CASES " in *" $1 "*) return 0;; esac; return 1; }

# ---- model keys ------------------------------------------------------------
# B<k> block16 y<k>, I<k> indep16 y<k>, Q<k> qmiss qm<k>; anything else is a
# derived model under $WORK/models whose VR is its base's.
model_dir () {
    case "$1" in
        B[0-9]*) echo "$NULLD/block16/m/y${1#B}" ;;
        I[0-9]*) echo "$NULLD/indep16/m/y${1#I}" ;;
        Q[0-9]*) echo "$NULLD/qmiss/m/qm${1#Q}" ;;
        *)       echo "$WORK/models/$1" ;;
    esac
}
vr_file () {
    case "$1" in
        Bp1)  vr_file B1 ;;  Ip4) vr_file I4 ;;  Qp2) vr_file Q2 ;;
        Bs10) vr_file B10 ;; Is6) vr_file I6 ;;
        Bl1)  vr_file B1 ;;  Bl9) vr_file B9 ;;  Il2) vr_file I2 ;;  Ql3) vr_file Q3 ;;
        B[0-9]*) echo "$NULLD/block16/mvr_y${1#B}.varianceRatio.txt" ;;
        I[0-9]*) echo "$NULLD/indep16/mvr_y${1#I}.varianceRatio.txt" ;;
        Q[0-9]*) echo "$NULLD/qmiss/mvr_qm${1#Q}.varianceRatio.txt" ;;
    esac
}

echo "=== derived models ==="
MK="$HERE/make_model_variant.py"
python3 "$MK" "$(model_dir B1)"  "$WORK/models/Bp1"  --permute 11 &&
python3 "$MK" "$(model_dir I4)"  "$WORK/models/Ip4"  --permute 12 &&
python3 "$MK" "$(model_dir Q2)"  "$WORK/models/Qp2"  --permute 13 &&
python3 "$MK" "$(model_dir B10)" "$WORK/models/Bs10" --json SPA_Cutoff=0.5 &&
python3 "$MK" "$(model_dir I6)"  "$WORK/models/Is6"  --json SPA_Cutoff=0.5 ||
    { echo "FAIL: could not build derived models"; exit 1; }
# LOCO variants: same construction as run_mt_correctness.sh (copy the model,
# write a perturbed chr<N>/ set, declare loco_chroms). The perturbation is what
# makes a silently ignored chr1/ visible.
mk_loco () {  # mk_loco <src model> <dst dir> <chrom> <factor> <loco_chroms json>
    python3 - "$1" "$2" "$3" "$4" "$5" "$HERE/make_loco_model.py" <<'PY2'
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
if has_case D15; then
    mk_loco "$(model_dir B1)" "$WORK/models/Bl1" 1 0.90 '[1]' &&
    mk_loco "$(model_dir B9)" "$WORK/models/Bl9" 1 1.10 '[1]' &&
    mk_loco "$(model_dir I2)" "$WORK/models/Il2" 2 0.85 '[2]' &&
    mk_loco "$(model_dir Q3)" "$WORK/models/Ql3" 1 1.20 '[1]' ||
        { echo "FAIL: could not build LOCO models"; exit 1; }
fi

# ---- config pieces ---------------------------------------------------------
# tag -> genotype file + the keys that change the columns or the numbers
tag_yaml () {  # tag_yaml <tag> <nThreads>
    local plink="$MID" minmac=1 er=4 more=false firth=false extra=""
    case "$1" in
        mid_plain) ;;
        mid_more)  more=true; firth=true; extra=$'is_Firth_beta: true\npCutoffforFirth: 0.05' ;;
        mid_noadj) extra="isnoadjCov: true" ;;
        rare_plain|rare_scalar) plink="$RARE" ;;
        rare_mac5) plink="$RARE"; minmac=5 ;;
        rare_er20) plink="$RARE"; er=20 ;;
        mid_loco)  extra=$'LOCO: true\nchrom: "1"' ;;
        *) echo "unknown tag $1" >&2; return 1 ;;
    esac
    cat <<EOF
genoType: plink
plinkFile: $plink
minMAF: 0
minMAC: $minmac
maxMissRate: 0.15
AlleleOrder: alt-first
MACCutoffforER: $er
nThreads: $2
isMoreOutput: $more
isFirth: $firth
EOF
    [ -n "$extra" ] && printf '%s\n' "$extra"
    return 0
}
tag_env () { [ "$1" = rare_scalar ] && echo "SAIGE_STEP2_SCALAR_DECODE=1" || echo "SAIGE_STEP2_SCALAR_DECODE=0"; }

# ---- goldens: every (model, tag) a case needs, run once, in parallel -------
GLIST="$WORK/golden/needed.txt"; : > "$GLIST"
need () { local tag="$1"; shift; for m in "$@"; do echo "$m $tag" >> "$GLIST"; done; }

D1T="B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16"
D2T="I1 I2 I3 I4 I5 I6 I7 I8 I9 I10 I11 I12 I13 I14 I15 I16"
D3T="I3 B10 I5 B2 Q2 I7 B12"
D4T="B1 Bp1 Ip4 I4 Qp2"
D5T="Bs10 Is6 B1 I9 Q1 B11"
D6T="B1 B9 I1 I2 Q1 Q4"

has_case D1  && need mid_plain $D1T
has_case D2  && need mid_plain $D2T
{ has_case D3 || has_case D12 || has_case D14; } && need mid_plain $D3T
has_case D14 && need mid_noadj I3
has_case D4  && need mid_plain $D4T
has_case D5  && need mid_plain $D5T
has_case D6  && need rare_plain $D6T
has_case D7  && need rare_mac5 $D6T
has_case D8  && need rare_er20 $D6T
has_case D9  && need mid_more $D3T
has_case D13 && need rare_scalar $D6T
D15T="Bl1 Bl9 Il2 Ql3"
has_case D15 && need mid_loco $D15T
has_case D15 && need mid_plain B9
sort -u "$GLIST" -o "$GLIST"

golden_one () {  # golden_one <model key> <tag>
    local m="$1" tag="$2" cfg="$WORK/cfg/golden_$2_$1.yaml"
    [ -s "$WORK/golden/${tag}_$m.txt" ] && return 0
    { echo "modelFile:         $(model_dir "$m")"
      echo "varianceRatioFile: $(vr_file "$m")"
      echo "outputFile:        $WORK/golden/${tag}_$m.txt"
      tag_yaml "$tag" 1; } > "$cfg"
    env "$(tag_env "$tag")" "$BIN" "$cfg" > "$WORK/golden/${tag}_$m.log" 2>&1 ||
        echo "golden $tag/$m exited $?" >> "$WORK/golden/errors.txt"
}
export -f golden_one tag_yaml tag_env model_dir vr_file
export WORK BIN NULLD MID RARE
echo "=== goldens: $(wc -l < "$GLIST") single-trait runs ($JOBS at a time) ==="
: > "$WORK/golden/errors.txt"
xargs -P "$JOBS" -L 1 bash -c 'golden_one "$0" "$1"' < "$GLIST"
if [ -s "$WORK/golden/errors.txt" ]; then
    sed 's/^/  FAIL /' "$WORK/golden/errors.txt"; fail=$((fail+$(wc -l < "$WORK/golden/errors.txt")))
fi

# ---- multi-trait runs ------------------------------------------------------
MT_EXTRA=""; MT_OVERRIDE=""; MT_ENVV="SAIGE_STEP2_SCALAR_DECODE=0"
mt_case () {  # mt_case <name> <tag> <nThreads> <model keys...>
    local name="$1" tag="$2" thr="$3"; shift 3
    local cfg="$WORK/cfg/$name.yaml"
    { tag_yaml "$tag" "$thr"
      [ -n "$MT_EXTRA" ] && printf '%s\n' "$MT_EXTRA"
      echo "models:"
      for m in "$@"; do
          echo "  - traitName: $m"
          echo "    modelFile:         $(model_dir "$m")"
          echo "    varianceRatioFile: $(vr_file "$m")"
          echo "    outputFile: $WORK/out/$name.$m.txt"
          [ "$m" = "$MT_OVERRIDE" ] && echo "    isnoadjCov: true"
      done; } > "$cfg"
    env "$MT_ENVV" "$BIN" "$cfg" > "$WORK/out/$name.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        bad "$name exited $rc: $(tail -1 "$WORK/out/$name.log")"; return
    fi
    for m in "$@"; do
        local gtag="$tag"
        [ "$m" = "$MT_OVERRIDE" ] && gtag=mid_noadj
        local g="$WORK/golden/${gtag}_$m.txt"
        if [ -s "$g" ] && cmp -s "$WORK/out/$name.$m.txt" "$g"; then
            pass=$((pass+1))
        else
            bad "$name/$m differs from golden ${gtag}_$m"
            python3 "$HERE/cmp_assoc.py" "$g" "$WORK/out/$name.$m.txt" 2>&1 | sed 's/^/      /' | head -20
        fi
    done
    echo "  ran $name (P=$#)"
}
xsame () {  # xsame <label> <run A> <run B> <keys...>
    local label="$1" a="$2" b="$3"; shift 3
    for m in "$@"; do
        if cmp -s "$WORK/out/$a.$m.txt" "$WORK/out/$b.$m.txt"; then pass=$((pass+1))
        else bad "$label for $m ($a != $b)"; fi
    done
}

echo "=== multi-trait runs ==="
has_case D1 && mt_case D1_block16 mid_plain 1 $D1T
has_case D2 && mt_case D2_indep16 mid_plain 1 $D2T
has_case D3 && mt_case D3_interleaved mid_plain 1 $D3T
has_case D4 && mt_case D4_permuted mid_plain 1 $D4T
has_case D5 && mt_case D5_spa mid_plain 1 $D5T
has_case D6 && mt_case D6_rare rare_plain 1 $D6T
has_case D7 && mt_case D7_minmac5 rare_mac5 1 $D6T
has_case D8 && mt_case D8_er20 rare_er20 1 $D6T
has_case D9 && mt_case D9_more_firth mid_more 1 $D3T
if has_case D10 && has_case D2; then
    mt_case D10_indep16_t8 mid_plain 8 $D2T
    xsame "nThreads=8 moved a byte" D10_indep16_t8 D2_indep16 $D2T
fi
if has_case D11 && has_case D2; then
    MT_EXTRA="mtBatch: false"
    mt_case D11_indep16_nobatch mid_plain 1 $D2T
    MT_EXTRA=""
    xsame "mtBatch=false moved a byte" D11_indep16_nobatch D2_indep16 $D2T
fi
if has_case D12 && has_case D3; then
    MT_EXTRA="mtBlockSize: 7"
    mt_case D12_interleaved_blk7 mid_plain 4 $D3T
    MT_EXTRA=""
    xsame "mtBlockSize=7/nThreads=4 moved a byte" D12_interleaved_blk7 D3_interleaved $D3T
fi
if has_case D13; then
    MT_ENVV="SAIGE_STEP2_SCALAR_DECODE=1"
    mt_case D13_rare_scalar_decode rare_scalar 1 $D6T
    MT_ENVV="SAIGE_STEP2_SCALAR_DECODE=0"
fi
if has_case D14; then
    MT_OVERRIDE="I3"
    mt_case D14_noadj_override mid_plain 1 $D3T
    MT_OVERRIDE=""
    if grep -q "isnoadjCov=true" "$WORK/out/D14_noadj_override.log"; then
        ok "D14: I3 gated out of the batch path (isnoadjCov=true)"
    else
        bad "D14: the override did not reach the gate table"
    fi
fi
if has_case D15; then
    mt_case D15_loco mid_loco 1 $D15T
    if grep -q "fell back to the genome-wide fit" "$WORK/out/D15_loco.log"; then
        ok "D15: mixed LOCO state announced"
    else
        bad "D15: mixed LOCO state NOT announced"
    fi
    if cmp -s "$WORK/golden/mid_loco_Bl9.txt" "$WORK/golden/mid_plain_B9.txt"; then
        bad "D15: Bl9 chr1 fit == B9 genome-wide fit; chr1/ was not read"
    else
        ok "D15: Bl9 chr1 fit differs from the genome-wide fit"
    fi
fi

# ---- coverage: the cases have to actually reach what they claim to test ----
echo "=== coverage ==="
cov_sum () {  # cov_sum <case> <extended regex ending in a number> <label>: sums that number
    local f="$WORK/out/$1.log"
    [ -e "$f" ] || return 0
    local v
    v=$(grep -oE "$2" "$f" | sed -E 's/^.*[^0-9]([0-9]+)$/\1/' | awk '{s+=$1} END{print s+0}')
    if [ "${v:-0}" -gt 0 ]; then ok "$1: $3 = $v"; else bad "$1: $3 = 0 -- not exercised"; fi
}
cov_sum D1_block16 "flip opposite to the union column: [0-9]+" "pairs whose flip is opposite to the union column"
cov_sum D1_block16 "column: [0-9]+ pairs \([0-9]+" "  ...of which stayed batched"
cov_sum D2_indep16 "flip opposite to the union column: [0-9]+" "pairs whose flip is opposite to the union column"
cov_sum D2_indep16 "column: [0-9]+ pairs \([0-9]+" "  ...of which stayed batched"
cov_sum D3_interleaved "flip opposite to the union column: [0-9]+" "pairs whose flip is opposite to the union column"
cov_sum D7_minmac5 "only some traits: [0-9]+" "markers passing QC for only some traits"
cov_sum D6_rare "only some traits: [0-9]+" "markers passing QC for only some traits"
cov_sum D5_spa "batched, [0-9]+" "fallback pairs"
if [ -e "$WORK/out/D2_indep16.log" ]; then
    ns=$(grep -cE "\] n=47500, flip opposite" "$WORK/out/D2_indep16.log")
    [ "$ns" -eq 16 ] && ok "D2: 16 traits, each n=47500 with its own sample list" ||
        bad "D2: expected 16 own-sample traits with n=47500, saw $ns"
fi
if [ -e "$WORK/golden/rare_plain_B9.txt" ]; then
    # ER pairs: binary rows with MAC <= 4 (n = 30000 for B9)
    ner=$(awk -F'\t' 'NR>1 {ac=$6; m=(ac<60000-ac)?ac:60000-ac; if (m<=4) c++} END{print c+0}' "$WORK/golden/rare_plain_B9.txt")
    [ "$ner" -gt 0 ] && ok "D6: B9 has $ner ER rows (MAC <= 4)" || bad "D6: no ER rows"
fi
if [ -e "$WORK/golden/mid_plain_Bs10.txt" ]; then
    nspa=$(awk -F'\t' 'NR>1 && $15=="true" {c++} END{print c+0}' "$WORK/golden/mid_plain_Bs10.txt")
    [ "$nspa" -gt 100 ] && ok "D5: Bs10 has $nspa SPA-converged rows" || bad "D5: only $nspa SPA rows"
fi

# ---- refusals ---------------------------------------------------------------
echo "=== refusals ==="
expect_error () {  # expect_error <name> <substring>   (config on stdin)
    local name="$1" want="$2" cfg="$WORK/cfg/$1.yaml" log="$WORK/out/$1.log"
    cat > "$cfg"
    "$BIN" "$cfg" > "$log" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then bad "$name: exited 0, expected a hard error"
    elif grep -qF "$want" "$log"; then ok "$name: rejected -- $want"
    else bad "$name: exited $rc without '$want': $(tail -1 "$log")"; fi
}
two_models () {
    echo "models:"
    for m in B1 B9; do
        echo "  - traitName: $m"
        echo "    modelFile: $(model_dir $m)"
        echo "    varianceRatioFile: $(vr_file $m)"
        echo "    outputFile: $WORK/out/err.$m.txt"
    done
}
if has_case E1; then
    { tag_yaml mid_plain 1; echo "mtRequireSameSamples: true"; two_models; } |
        expect_error E1_require_same "mtRequireSameSamples: true requires identical sample IDs"
fi
if has_case E2; then
    { printf 'genoType: bgen\nbgenFile: /opt/saige/data/readers/mid2kr_reffirst.bgen\n'
      printf 'bgenSampleFile: /opt/saige/data/readers/mid2kr_reffirst.sample\n'
      printf 'AlleleOrder: ref-first\nnThreads: 1\n'; two_models; } |
        expect_error E2_bgen "currently supports genoType: plink only"
fi
if has_case E3; then
    cond=$(awk 'NR==10 {print $2}' "$MID.bim")
    { tag_yaml mid_plain 1; echo "condition: \"$cond\""; two_models; } |
        expect_error E3_condition "Conditional analysis is not supported in a multi-trait run"
fi

echo
echo "=== $pass passed, $fail failed  (workdir $WORK) ==="
[ $fail -eq 0 ]
