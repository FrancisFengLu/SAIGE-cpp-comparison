#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# End-to-end LOCO integration test: step 1 -> step 2, plus the R golden
# reference at test_data/loco_reference/.
#
#   tests/run_loco_e2e.sh
#
# Environment:
#   SAIGE_NULL_BIN    step-1 binary        (default ../step1_saige-null/saige-null)
#   SAIGE_STEP2_BIN   step-2 binary        (default ../step2_saige-step2/saige-step2)
#   SAIGE_PRE_NULL    PRE-LOCO step-1 binary  (no longer used; step-1 byte
#                     identity was deliberately given up -- see PRIORITY 4b)
#   SAIGE_PRE_STEP2   PRE-LOCO step-2 binary  (optional; enables check 4c)
#   SAIGE_E2E_TMP     scratch dir (default mktemp -d)
#
# To build the pre-LOCO binaries for check 4:
#   git archive <pre-loco-rev> SAIGE_cpp_260716/step1_saige-null \
#       SAIGE_cpp_260716/step2_saige-step2 | tar -x -C /tmp/preloco
#   (cd /tmp/preloco/SAIGE_cpp_260716/step1_saige-null && make -j8)
#   (cd /tmp/preloco/SAIGE_cpp_260716/step2_saige-step2 && make -j8)
#
# NOTES that this test depends on and that you must not "optimise away":
#  * nthreads MUST be 1 for step 1. At nthreads>1 the TBB float32 reduction in
#    parallelCrossProd is order-nondeterministic and tau swings ~20%.
#  * The R reference's tau is NOT reproducible unless skipVarianceRatioEstimation
#    is on. Only the *_novr exports are usable as a tight oracle, and those have
#    no step-2 output. So:
#       - checks against export/grm2chr_binary_novr/  are TIGHT (reproducible)
#       - checks against export/grm2chr_binary/step2/ are LOOSE (tau noise)
#  * The C++ full-genome fit DOES now re-solve Get_Coef at the final tau, as R
#    does (SAIGE_fitGLMM_fast.R:967 binary / :963 quantitative). What remains is
#    the trace-estimator RNG difference: the C++ tau lands ~10% below R's, and
#    that tau gap is what drives the residual full-genome mu error. Thresholds
#    below are set accordingly. See PRIORITY 4b for the recorded numbers.
# ---------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CPP_ROOT="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$CPP_ROOT/.." && pwd)"

S1="${SAIGE_NULL_BIN:-$CPP_ROOT/step1_saige-null/saige-null}"
S2="${SAIGE_STEP2_BIN:-$CPP_ROOT/step2_saige-step2/saige-step2}"
DATA="$REPO/test_data/SAIGE/extdata/input"
REF="$REPO/test_data/loco_reference/export/grm2chr_binary"
REFN="${REF}_novr"
TMP="${SAIGE_E2E_TMP:-$(mktemp -d /tmp/saige_loco_e2e.XXXXXX)}"
mkdir -p "$TMP"

PLINK_GRM="$DATA/plinkforGRM_1000samples_10kMarkers"
PLINK_G2="$DATA/genotype_100markers_2chr"
PHENO="$DATA/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"

PASS=0; FAIL=0; SKIP=0
ok()   { echo "  PASS  $*"; PASS=$((PASS+1)); }
bad()  { echo "  FAIL  $*"; FAIL=$((FAIL+1)); }
skip() { echo "  SKIP  $*"; SKIP=$((SKIP+1)); }
hdr()  { echo; echo "=== $* ==="; }

for b in "$S1" "$S2"; do
  [[ -x "$b" ]] || { echo "binary not found/executable: $b" >&2; exit 2; }
done
[[ -d "$REFN" ]] || { echo "R reference not found: $REFN" >&2; exit 2; }
if [[ -z "${R_HOME:-}" ]] && command -v Rscript >/dev/null 2>&1; then
  R_HOME="$(Rscript -e 'cat(R.home())')"; export R_HOME
fi

echo "scratch: $TMP"

# ---------------------------------------------------------------------------
# step-1 configs
# ---------------------------------------------------------------------------
mk_s1cfg() { # mk_s1cfg <out.yaml> <outdir> <loco>
  cat > "$1" <<EOF
paths:
  plinkFile:     $PLINK_GRM
  out_prefix:    $2
  out_prefix_vr: $2/vr
  sparse_grm:     ""
  sparse_grm_ids: ""
  overwrite_varratio: true
design:
  csv: $PHENO
  iid_col: IID
  y_col: y_binary
  covar_cols:
    - x1
    - x2
fit:
  trait: binary
  loco: $3
  nthreads: 1
  maxiter: 20
  tol: 0.02
  tolPCG: 1e-5
  maxiterPCG: 500
  nrun: 30
  trace_seed: 10
  num_markers_for_vr: 0
  use_sparse_grm_to_fit: false
  use_pcg_with_sparse_grm: false
  relatedness_cutoff: 0.125
  overwrite_vr: true
EOF
}

mk_s2cfg() { # mk_s2cfg <out.yaml> <modeldir> <outfile> <LOCO> <chrom> [omit_loco_keys]
  cat > "$1" <<EOF
modelFile:         $2
varianceRatioFile: $TMP/vr.txt
genoType:   plink
plinkFile:  $PLINK_G2
outputFile: $3
minMAF: 0
minMAC: 0.5
maxMissRate: 0.15
AlleleOrder: alt-first
isMoreOutput: true
isFirth: true
MACCutoffforER: 4
nThreads: 1
EOF
  if [[ "${6:-}" != "omit" ]]; then
    printf 'LOCO: %s\nchrom: "%s"\n' "$4" "$5" >> "$1"
  fi
}

# Variance ratio: reuse R's single-VR file so step-2 scaling matches the
# reference. (The C++ step 1 here is run with num_markers_for_vr: 0.)
cp "$REF/varianceRatio.txt" "$TMP/vr.txt"

hdr "step 1: LOCO and non-LOCO runs (nthreads=1)"
mk_s1cfg "$TMP/s1_loco.yaml"   "$TMP/m_loco"   true
mk_s1cfg "$TMP/s1_noloco.yaml" "$TMP/m_noloco" false
"$S1" -c "$TMP/s1_loco.yaml"   > "$TMP/s1_loco.log"   2>&1 || { bad "step1 LOCO run"; }
"$S1" -c "$TMP/s1_noloco.yaml" > "$TMP/s1_noloco.log" 2>&1 || { bad "step1 non-LOCO run"; }
grep -q "LOCO: on" "$TMP/s1_loco.log" && ok "step1 reports LOCO on" || bad "step1 did not report LOCO on"

# ---------------------------------------------------------------------------
hdr "PRIORITY 1 - format contract (LOCO_FORMAT.md)"

PERCHR="mu res V offset XV XVX XVX_inv XVX_inv_XV XXVX_inv S_a"
for j in 1 2; do
  d="$TMP/m_loco/chr$j"
  if [[ -d "$d" ]]; then ok "chr$j/ directory exists"; else bad "chr$j/ directory missing"; continue; fi
  miss=""
  for f in $PERCHR; do [[ -f "$d/$f.arma" ]] || miss="$miss $f"; done
  [[ -z "$miss" ]] && ok "chr$j/ has all 10 per-chromosome .arma files" \
                   || bad "chr$j/ missing:$miss"
  extra=""
  for f in X y; do [[ -f "$d/$f.arma" ]] && extra="$extra $f"; done
  [[ -z "$extra" ]] && ok "chr$j/ does not duplicate X/y" || bad "chr$j/ wrongly contains:$extra"
done
[[ -d "$TMP/m_noloco/chr1" ]] && bad "non-LOCO run wrote chr1/" || ok "non-LOCO run writes no chr<j>/"

python3 - "$TMP/m_loco" "$TMP/m_noloco" <<'PY'
import json,sys,os
lo=json.load(open(os.path.join(sys.argv[1],'nullmodel.json')))
nl=json.load(open(os.path.join(sys.argv[2],'nullmodel.json')))
bad=[]
if lo.get('loco') is not True: bad.append('loco != true in LOCO model')
if lo.get('loco_chroms') != [1,2]: bad.append('loco_chroms=%r, expected [1,2]'%lo.get('loco_chroms'))
if nl.get('loco') is not False: bad.append('loco != false in non-LOCO model')
if nl.get('loco_chroms') not in (None,[]): bad.append('non-LOCO loco_chroms=%r'%nl.get('loco_chroms'))
print('\n'.join('  FAIL  json: '+b for b in bad) if bad else '  PASS  nullmodel.json loco / loco_chroms fields')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

# dimensions of every per-chromosome .arma vs the top-level copy
python3 - "$TMP/m_loco" <<'PY'
import sys,os,numpy as np
def dims(p):
    with open(p,'rb') as f:
        f.readline(); d=f.readline().split()
        return (int(d[0]), int(d[1]) if len(d)>1 else 1)
m=sys.argv[1]; bad=[]
for j in (1,2):
    for f in "mu res V offset XV XVX XVX_inv XVX_inv_XV XXVX_inv S_a".split():
        a=dims(os.path.join(m,f+'.arma')); b=dims(os.path.join(m,'chr%d'%j,f+'.arma'))
        if a!=b: bad.append('chr%d/%s dims %s != top-level %s'%(j,f,b,a))
print('\n'.join('  FAIL  '+b for b in bad) if bad else '  PASS  chr<j>/*.arma dimensions match the top-level copies')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

# step 2 consumes it
mk_s2cfg "$TMP/s2_l1.yaml" "$TMP/m_loco" "$TMP/s2_loco_chr1.txt" true  1
mk_s2cfg "$TMP/s2_l2.yaml" "$TMP/m_loco" "$TMP/s2_loco_chr2.txt" true  2
mk_s2cfg "$TMP/s2_nl.yaml" "$TMP/m_loco" "$TMP/s2_noloco.txt"    false ""
"$S2" "$TMP/s2_l1.yaml" > "$TMP/s2_l1.log" 2>&1 && ok "step2 LOCO chrom=1 ran"  || bad "step2 LOCO chrom=1 failed"
"$S2" "$TMP/s2_l2.yaml" > "$TMP/s2_l2.log" 2>&1 && ok "step2 LOCO chrom=2 ran"  || bad "step2 LOCO chrom=2 failed"
"$S2" "$TMP/s2_nl.yaml" > "$TMP/s2_nl.log" 2>&1 && ok "step2 LOCO=false ran"    || bad "step2 LOCO=false failed"
grep -q "per-chromosome fit will be read from .*chr1" "$TMP/s2_l1.log" \
  && ok "step2 reports reading chr1/" || bad "step2 never reported reading chr1/"
grep -q "restricting to chromosome 1: 90 of 100" "$TMP/s2_l1.log" \
  && ok "step2 LOCO marker filter: 90/100 on chr1" || bad "step2 chr1 marker filter wrong"
grep -q "restricting to chromosome 2: 10 of 100" "$TMP/s2_l2.log" \
  && ok "step2 LOCO marker filter: 10/100 on chr2" || bad "step2 chr2 marker filter wrong"

hdr "PRIORITY 1b - guards (LOCO_FORMAT.md 'Step 2 behaviour')"
mk_s2cfg "$TMP/g1.yaml" "$TMP/m_noloco" "$TMP/g1.txt" true 1
"$S2" "$TMP/g1.yaml" > "$TMP/g1.log" 2>&1
[[ $? -ne 0 ]] && grep -qi "does not contain LOCO results" "$TMP/g1.log" \
  && ok "guard: LOCO=true against loco:false model -> non-zero exit" \
  || bad "guard: LOCO=true against loco:false model did not error"

mk_s2cfg "$TMP/g2.yaml" "$TMP/m_loco" "$TMP/g2.txt" true ""
"$S2" "$TMP/g2.yaml" > "$TMP/g2.log" 2>&1
[[ $? -ne 0 ]] && grep -qi "chrom needs to be specified" "$TMP/g2.log" \
  && ok "guard: LOCO=true with empty chrom -> non-zero exit" \
  || bad "guard: LOCO=true with empty chrom did not error"

mk_s2cfg "$TMP/g3.yaml" "$TMP/m_loco" "$TMP/g3.txt" true X
"$S2" "$TMP/g3.yaml" > "$TMP/g3.log" 2>&1
grep -q "has no LOCO result in the null model; using the full-genome fit" "$TMP/g3.log" \
  && ok "guard: chrom=X falls back to the full-genome fit" \
  || bad "guard: chrom=X did not fall back"
# (the run then dies with 'No markers on chrom X' because the test bed has none;
#  that is a property of the data, not of the guard.)

# ---------------------------------------------------------------------------
hdr "PRIORITY 2 - LOCO actually changes the answer"

python3 - "$TMP/m_loco" <<'PY'
import sys,os,numpy as np
def rd(p):
    with open(p,'rb') as f:
        f.readline(); d=f.readline().split(); r=int(d[0]); c=int(d[1]) if len(d)>1 else 1
        a=np.frombuffer(f.read(r*c*8),dtype='<f8')
        return a.reshape((c,r)).T if c>1 else a
m=sys.argv[1]; bad=[]
full=rd(m+'/mu.arma'); c1=rd(m+'/chr1/mu.arma'); c2=rd(m+'/chr2/mu.arma')
d1=np.abs(c1-full).max(); d2=np.abs(c2-full).max(); d12=np.abs(c1-c2).max()
print('        mu: max|chr1-full|=%.4g  max|chr2-full|=%.4g  max|chr1-chr2|=%.4g'%(d1,d2,d12))
if d1 < 1e-4: bad.append('chr1 mu is identical to the full-genome mu')
if d2 < 1e-6: bad.append('chr2 mu is identical to the full-genome mu')
if d12< 1e-4: bad.append('chr1 mu is identical to chr2 mu')
print('\n'.join('  FAIL  '+b for b in bad) if bad else '  PASS  chr1/chr2/full mu are all genuinely different')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

python3 - "$TMP" <<'PY'
import sys,csv,os,numpy as np
T=sys.argv[1]
def rd(f): return {x['MarkerID']:x for x in csv.DictReader(open(f),delimiter='\t')}
nl=rd(T+'/s2_noloco.txt')
bad=[]
for j,thr in ((1,1e-3),(2,1e-4)):
    a=rd(T+'/s2_loco_chr%d.txt'%j); ks=sorted(set(a)&set(nl))
    if not ks: bad.append('no common markers for chr%d'%j); continue
    rel=max(abs(float(a[k]['p.value'])-float(nl[k]['p.value']))/float(nl[k]['p.value']) for k in ks)
    print('        chr%d: %d markers, max relative p-value change vs non-LOCO = %.4g'%(j,len(ks),rel))
    if rel < thr: bad.append('chr%d LOCO p-values are effectively identical to non-LOCO (%.3g)'%(j,rel))
print('\n'.join('  FAIL  '+b for b in bad) if bad else '  PASS  step-2 p-values move when LOCO is on')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

# The chr<j>/ directory really is the one being read: swap chr2/ for chr1/ and
# the chr2 p-values must change.
rm -rf "$TMP/m_swap"; cp -r "$TMP/m_loco" "$TMP/m_swap"
rm -rf "$TMP/m_swap/chr2"; cp -r "$TMP/m_loco/chr1" "$TMP/m_swap/chr2"
mk_s2cfg "$TMP/s2_swap.yaml" "$TMP/m_swap" "$TMP/s2_swap_chr2.txt" true 2
"$S2" "$TMP/s2_swap.yaml" > "$TMP/s2_swap.log" 2>&1
cmp -s "$TMP/s2_swap_chr2.txt" "$TMP/s2_loco_chr2.txt" \
  && bad "chr2/ content is ignored (swapping chr1/ in changed nothing)" \
  || ok  "chr<j>/ content is genuinely read (chr1/-for-chr2/ swap changes p-values)"

# ---------------------------------------------------------------------------
hdr "PRIORITY 3 - agreement with the R reference"

echo "  -- TIGHT (vs ${REFN##*/}, reproducible: skipVarianceRatioEstimation=TRUE) --"
grep -q "leave chromosome 1 out  \[0,9637\]" "$TMP/s1_loco.log" \
  && ok "chr1 QC'd marker range [0,9637] matches R chrom_index.csv" \
  || bad "chr1 QC'd marker range does not match R"
grep -q "leave chromosome 2 out  \[9638,9649\]" "$TMP/s1_loco.log" \
  && ok "chr2 QC'd marker range [9638,9649] matches R chrom_index.csv" \
  || bad "chr2 QC'd marker range does not match R"
grep -qE "chr1: Get_Coef_LOCO iterations = 1" "$TMP/s1_loco.log" && \
grep -qE "chr2: Get_Coef_LOCO iterations = 1" "$TMP/s1_loco.log" \
  && ok "1 Newton iteration per chromosome, as in R loco_newton_iterations.csv" \
  || bad "Newton iteration count differs from R (warm-start chain suspect)"

# mu / res / V are invariant to the covariate basis, so they are comparable even
# though the C++ X.arma is the raw design and R's X.csv is QR space.
python3 - "$TMP/m_loco" "$REFN" <<'PY'
import sys,os,csv,numpy as np
m,R=sys.argv[1],sys.argv[2]
def rd(p):
    with open(p,'rb') as f:
        f.readline(); d=f.readline().split(); r=int(d[0]); c=int(d[1]) if len(d)>1 else 1
        a=np.frombuffer(f.read(r*c*8),dtype='<f8'); return a.reshape((c,r)).T if c>1 else a
def col(f,c): return np.array([float(x[c]) for x in csv.DictReader(open(f))])
rfull=col(R+'/null_vectors.csv','fitted.values'); cfull=rd(m+'/mu.arma')
base=np.abs(rfull-cfull).max()
print('        full-genome mu: max|R-C++| = %.4g   (tau-gap-driven; see PRIORITY 4b)'%base)
bad=[]
for j,tag in ((1,'chr01'),(2,'chr02')):
    rj=col(R+'/loco/%s_vectors.csv'%tag,'fitted.values'); cj=rd(m+'/chr%d/mu.arma'%j)
    tot=np.abs(rj-cj).max()
    dd=np.abs((rj-rfull)-(cj-cfull)).max()
    print('        chr%d mu: max|R-C++| = %.4g ;  max|delta_R - delta_C++| = %.4g'%(j,tot,dd))
    # the LOCO solve must not add materially on top of the known full-fit gap
    if tot > base*3 + 5e-3: bad.append('chr%d mu error %.3g is far above the full-fit gap %.3g'%(j,tot,base))
print('\n'.join('  FAIL  '+b for b in bad) if bad
      else '  PASS  per-chromosome mu tracks R to within the pre-existing full-fit gap')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

echo "  -- LOOSE (vs ${REF##*/}/step2, R tau is NOT reproducible; see README 2a) --"
python3 - "$TMP" "$REF" <<'PY'
import sys,csv,numpy as np
T,R=sys.argv[1],sys.argv[2]
def rd(f): return {x['MarkerID']:x for x in csv.DictReader(open(f),delimiter='\t')}
bad=[]
for j in (1,2):
    r=rd(R+'/step2/step2_chr%d.txt'%j); c=rd(T+'/s2_loco_chr%d.txt'%j)
    if set(r)!=set(c): bad.append('chr%d marker set differs: R=%d C++=%d'%(j,len(r),len(c)))
    ks=sorted(set(r)&set(c))
    dl=[abs(np.log10(float(r[k]['p.value']))-np.log10(float(c[k]['p.value']))) for k in ks]
    print('        chr%d: %d markers, max|dlog10 p| = %.4g, median = %.4g'%(j,len(ks),max(dl),np.median(dl)))
    if max(dl) > 0.3: bad.append('chr%d max|dlog10 p| = %.3g exceeds the loose 0.3 bar'%(j,max(dl)))
print('\n'.join('  FAIL  '+b for b in bad) if bad
      else '  PASS  step-2 marker sets identical to R; p-values agree loosely')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

# ---------------------------------------------------------------------------
hdr "PRIORITY 4 - no regression on the non-LOCO path"

# 4a: within this build, a LOCO run must leave the top-level (full-genome)
#     artifacts byte-identical to a non-LOCO run.
diffs=""
for f in mu res y V S_a X XV XVX XVX_inv XXVX_inv XVX_inv_XV offset; do
  cmp -s "$TMP/m_loco/$f.arma" "$TMP/m_noloco/$f.arma" || diffs="$diffs $f"
done
[[ -z "$diffs" ]] && ok "LOCO run leaves the top-level .arma set byte-identical to non-LOCO" \
                  || bad "LOCO run perturbed the full-genome fit:$diffs"

# 4b: the non-LOCO full-genome fit vs the R golden reference.
#
# This USED to assert byte-identity against a genuinely pre-LOCO build. That
# guarantee was deliberately given up when glmm.cpp gained R's unconditional
# final Get_Coef re-solve at the final tau (SAIGE_fitGLMM_fast.R:967 / :963);
# the whole point of that change is that the non-LOCO .arma set moves. Byte
# identity to an older binary is not the ground truth -- R is. So we assert
# agreement with R instead, using the same TIGHT *_novr oracle as PRIORITY 3.
#
# Recorded numbers so a future regression is visible:
#   BEFORE the final-Get_Coef fix : full-genome mu max|R-C++| = 7.556e-3
#                                   tau = 0.242919  (R 0.271582, 10.6% low)
#   AFTER  the final-Get_Coef fix : full-genome mu max|R-C++| = 8.751e-3
#                                   tau = 0.242919  (UNCHANGED -- the re-solve
#                                   happens after the AI-REML loop, so it cannot
#                                   move tau; the residual mu gap is driven by
#                                   the tau gap, i.e. by the trace-estimator RNG,
#                                   not by the missing re-solve)
# Thresholds are set from the measured AFTER numbers with headroom.
python3 - "$TMP/m_noloco" "$REFN" <<'PY'
import sys,os,csv,numpy as np
m,R=sys.argv[1],sys.argv[2]
def rd(p):
    with open(p,'rb') as f:
        f.readline(); d=f.readline().split(); r=int(d[0]); c=int(d[1]) if len(d)>1 else 1
        a=np.frombuffer(f.read(r*c*8),dtype='<f8'); return a.reshape((c,r)).T if c>1 else a
def col(f,c): return np.array([float(x[c]) for x in csv.DictReader(open(f))])

MU_TOL   = 1.2e-2   # measured 8.751e-3 (was 7.556e-3 pre-fix)
TAU_RTOL = 0.15     # measured 0.106     (unchanged by the fix)

bad=[]
rmu=col(R+'/null_vectors.csv','fitted.values'); cmu=rd(m+'/mu.arma')
dmu=np.abs(rmu-cmu).max()
print('        non-LOCO mu:  max|R-C++| = %.4g   (bar %.3g; pre-fix 7.556e-3)'%(dmu,MU_TOL))
if dmu > MU_TOL: bad.append('non-LOCO mu error %.4g exceeds %.3g'%(dmu,MU_TOL))

rres=col(R+'/null_vectors.csv','residuals'); cres=rd(m+'/res.arma')
dres=np.abs(rres-cres).max()
print('        non-LOCO res: max|R-C++| = %.4g   (bar %.3g)'%(dres,MU_TOL))
if dres > MU_TOL: bad.append('non-LOCO res error %.4g exceeds %.3g'%(dres,MU_TOL))

rtau=[float(x) for x in open(R+'/theta.csv').read().split()[1:]]
import json
ctau=json.load(open(os.path.join(m,'nullmodel.json')))['theta']
rel=abs(ctau[-1]-rtau[-1])/abs(rtau[-1])
print('        non-LOCO tau: C++ %.6f vs R %.6f  (rel %.4g; bar %.3g)'%(ctau[-1],rtau[-1],rel,TAU_RTOL))
if rel > TAU_RTOL: bad.append('tau relative error %.4g exceeds %.3g'%(rel,TAU_RTOL))

print('\n'.join('  FAIL  '+b for b in bad) if bad
      else '  PASS  step1 non-LOCO full-genome fit agrees with the R reference')
sys.exit(1 if bad else 0)
PY
[[ $? -eq 0 ]] && PASS=$((PASS+1)) || FAIL=$((FAIL+1))

if [[ -n "${SAIGE_PRE_STEP2:-}" && -x "${SAIGE_PRE_STEP2}" ]]; then
  mk_s2cfg "$TMP/s2_preold.yaml" "$TMP/m_noloco" "$TMP/s2_pre.txt" false "" omit
  "${SAIGE_PRE_STEP2}" "$TMP/s2_preold.yaml" > "$TMP/s2_pre.log" 2>&1
  mk_s2cfg "$TMP/s2_newnl.yaml" "$TMP/m_noloco" "$TMP/s2_newnl.txt" false ""
  "$S2" "$TMP/s2_newnl.yaml" > "$TMP/s2_newnl.log" 2>&1
  cmp -s "$TMP/s2_pre.txt" "$TMP/s2_newnl.txt" \
    && ok "step2 non-LOCO output byte-identical to the pre-LOCO build" \
    || bad "step2 non-LOCO output changed vs pre-LOCO build"
else
  skip "step2 pre-LOCO byte-identity (set SAIGE_PRE_STEP2)"
fi

# ---------------------------------------------------------------------------
echo
echo "==================================================="
printf "  PASS %d   FAIL %d   SKIP %d\n" "$PASS" "$FAIL" "$SKIP"
echo "  artifacts: $TMP"
echo "==================================================="
[[ $FAIL -eq 0 ]] && exit 0 || exit 1
