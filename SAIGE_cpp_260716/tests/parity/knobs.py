"""Machine-readable inventory of every Step-2 switch that can change a number.

One row per knob.  Fields:

  name        canonical name used by parity.py (--set name=value)
  r           R step2_SPAtests.R CLI flag (None = R has no CLI flag for it)
  r_def       R CLI default, as the *CLI* defines it (NOT the SPAGMMATtest()
              function default -- they disagree for is_noadjCov, see notes)
  cpp         where the C++ side reads it:
                ("yaml", key)  top-level key in the step-2 config YAML
                ("json", key)  field in <modelFile>/nullmodel.json
                ("yaml+json", key) YAML key that overrides the JSON field
                None           C++ has no equivalent -- unreachable path
  cpp_def     C++ default when the key is absent
  kind        bool / num / str / list
  scope       single / region / both / io / perf
  effect      the code the knob steers
  note        default mismatch or other trap.  "" = defaults agree.

`MISMATCH` below is the derived list of knobs whose two defaults disagree:
those are the ones that manufacture fake cpp-vs-R differences when a run
leaves them unset on either side, so parity.py always writes both sides
explicitly rather than trusting a default.

Data-triggered branches (no flag at all) are in DATA_BRANCHES.
"""

# fmt: off
KNOBS = [
 # --- genotype input / QC -------------------------------------------------
 dict(name="genoType",        r=None,                       r_def="(implied by --bedFile/--bgenFile/--vcfFile/--pgenPrefix)",
      cpp=("yaml","genoType"),          cpp_def="plink",    kind="str",  scope="io",
      effect="PLINK/BGEN/VCF/PGEN reader selection", note=""),
 dict(name="AlleleOrder",     r="--AlleleOrder",            r_def="alt-first",
      cpp=("yaml","AlleleOrder"),       cpp_def="alt-first",kind="str",  scope="io",
      effect="which bim/bgen allele is Allele2; flips BETA sign", note=""),
 dict(name="vcfField",        r="--vcfField",               r_def="DS",
      cpp=("yaml","vcfField"),          cpp_def="GT",       kind="str",  scope="io",
      effect="dosage vs hard-call read from VCF",
      note="DEFAULT MISMATCH: R=DS, cpp=GT. DS vs GT changes every dosage."),
 dict(name="minMAF",          r="--minMAF",                 r_def=0.0,
      cpp=("yaml","minMAF"),            cpp_def=0.0,        kind="num",  scope="both",
      effect="marker drop filter (max of minMAF/minMAC applies)", note=""),
 dict(name="minMAC",          r="--minMAC",                 r_def=0.5,
      cpp=("yaml","minMAC"),            cpp_def=0.5,        kind="num",  scope="both",
      effect="marker drop filter", note=""),
 dict(name="maxMissing",      r="--maxMissing",             r_def=0.15,
      cpp=("yaml","maxMissRate"),       cpp_def=0.15,       kind="num",  scope="both",
      effect="marker drop filter", note="renamed key"),
 dict(name="minInfo",         r="--minInfo",                r_def=0.0,
      cpp=("yaml","minINFO"),           cpp_def=0.0,        kind="num",  scope="both",
      effect="imputation-info filter (only when is_imputed_data)", note="renamed key"),
 dict(name="is_imputed_data", r="--is_imputed_data",        r_def=False,
      cpp=("yaml","isImputation"),      cpp_def=False,      kind="bool", scope="both",
      effect="gates info-score output AND the dosage_zerod_* zeroing", note="renamed key"),
 dict(name="dosage_zerod_cutoff", r="--dosage_zerod_cutoff", r_def=0.2,
      cpp=("yaml","dosage_zerod_cutoff"), cpp_def=0.2,      kind="num",  scope="both",
      effect="dosages <= cutoff set to 0 on low-MAC markers", note=""),
 dict(name="dosage_zerod_MAC_cutoff", r="--dosage_zerod_MAC_cutoff", r_def=10.0,
      cpp=("yaml","dosage_zerod_MAC_cutoff"), cpp_def=10.0, kind="num",  scope="both",
      effect="MAC below which dosage zeroing applies", note=""),
 dict(name="impute_method",   r="--impute_method",          r_def="best_guess",
      cpp=("json","impute_method"),     cpp_def="mean",     kind="str",  scope="both",
      effect="missing-dosage fill: best_guess / mean / minor",
      note="DEFAULT MISMATCH and NO cpp YAML KEY. cpp reads it only from "
           "nullmodel.json; a step-2 run cannot change it from the config."),
 dict(name="idstoIncludeFile", r="--idstoIncludeFile",      r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="io",
      effect="marker subsetting before testing", note="NOT PORTED"),
 dict(name="rangestoIncludeFile", r="--rangestoIncludeFile", r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="io",
      effect="region subsetting before testing", note="NOT PORTED"),
 dict(name="vcfFilters",      r="--vcfFilters",             r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="io",
      effect="VCF FILTER-field predicates", note="NOT PORTED"),
 dict(name="subSampleFile",   r="--subSampleFile",          r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="io",
      effect="subset step-1 samples for step 2 (re-derives mu/res subset)",
      note="NOT PORTED"),

 # --- LOCO ----------------------------------------------------------------
 dict(name="LOCO",            r="--LOCO",                   r_def=True,
      cpp=("yaml","LOCO"),              cpp_def=False,      kind="bool", scope="both",
      effect="swap in per-chromosome mu/res/obj.noK; restrict markers to chrom",
      note="DEFAULT MISMATCH: R=TRUE, cpp=FALSE (documented in main.cpp:3867-3870)."),
 dict(name="chrom",           r="--chrom",                  r_def="",
      cpp=("yaml","chrom"),             cpp_def="",         kind="str",  scope="both",
      effect="which chr<N>/ block is loaded; also the marker filter", note=""),

 # --- score test / SPA / ER ----------------------------------------------
 dict(name="SPAcutoff",       r="--SPAcutoff",              r_def=2.0,
      cpp=("json","SPA_Cutoff"),        cpp_def=2.0,        kind="num",  scope="both",
      effect="|Tstat|/sqrt(var) above this -> SPA instead of normal approx",
      note="NO cpp YAML KEY; cpp reads it only from nullmodel.json."),
 dict(name="max_MAC_for_ER",  r="--max_MAC_for_ER",         r_def=4.0,
      cpp=("yaml","MACCutoffforER"),    cpp_def=4.0,        kind="num",  scope="both",
      effect="MAC <= this -> efficient-resampling exact test (binary only)",
      note="renamed key; both sides force 0 for non-binary traits"),
 dict(name="is_noadjCov",     r="--is_noadjCov",            r_def=True,
      cpp=("yaml+json","isnoadjCov"),   cpp_def=False,      kind="bool", scope="both",
      effect="scoreTestFast_noadjCov (no covariate adjustment of G)",
      note="DEFAULT MISMATCH: R CLI=TRUE, R function=FALSE, cpp=FALSE. "
           "MEASURED 2026-09-13 (tests/parity binary_single_noadjcov): with "
           "TRUE on both sides the two outputs are BIT-IDENTICAL over 2000 "
           "markers -- cpp reproduces R's arithmetic exactly, AF>0.5 defect "
           "included. So this is a DEFAULT difference, not a maths difference. "
           "Keep cpp's FALSE. Second, separate R bug: with is_noadjCov=TRUE "
           "and a VR file that has no null_noXadj row (every plain dense-fit "
           "step-1 output), R 1.5.2 aborts inside mainMarkerInCPP with "
           "'Mat::init(): requested size is not compatible with row vector "
           "layout' -- readInGLMM.R:418-427 lets a zero-length ratio vector "
           "through unguarded. Since TRUE is the CLI default, the plain "
           "documented R 1.5.2 command line does not run at all on such a "
           "model."),
 dict(name="is_fastTest",     r="--is_fastTest",            r_def=False,
      cpp=("json","isFastTest"),        cpp_def=True,       kind="bool", scope="both",
      effect="first pass with the no-GRM variance ratio, recompute markers "
             "with p < pval_cutoff_for_fastTest using the sparse-GRM VR",
      note="DEFAULT MISMATCH: R=FALSE, cpp=TRUE, and NO cpp YAML KEY. "
           "Harmless only because with no sparse GRM and isnoadjCov=FALSE "
           "R force-disables it (SAIGE_Test_main.R:349-351) and cpp detects "
           "the identical-context recompute and skips it (main.cpp:1618-1653). "
           "Verify, do not assume."),
 dict(name="pval_cutoff_for_fastTest", r=None,              r_def=0.05,
      cpp=("json","pval_cutoff_for_fastTest"), cpp_def=0.05, kind="num", scope="both",
      effect="p below which the fast test re-runs with the GRM variance ratio",
      note="no R CLI flag (function-only, always 0.05 from the CLI)"),
 dict(name="cateVarRatioMinMACVecExclude", r="--cateVarRatioMinMACVecExclude", r_def="10,20.5",
      cpp=("yaml","cateVarRatioMinMACVecExclude"), cpp_def="(from the VR file)", kind="list", scope="both",
      effect="MAC bin lower bounds selecting which categorical VR a marker gets",
      note="DEFAULT MISMATCH in source: R always passes its own CLI default; "
           "cpp derives the bins from the varianceRatio file unless overridden."),
 dict(name="cateVarRatioMaxMACVecInclude", r="--cateVarRatioMaxMACVecInclude", r_def="20.5",
      cpp=("yaml","cateVarRatioMaxMACVecInclude"), cpp_def="(from the VR file)", kind="list", scope="both",
      effect="MAC bin upper bounds",
      note="R appends nsample as a final open bin (SAIGE_Test_main.R:373); "
           "cpp does not -- it uses 1e10 from the loader instead. Same "
           "behaviour as long as MAC <= 2N, but it is not the same vector."),

 # --- effect sizes --------------------------------------------------------
 dict(name="is_Firth_beta",   r="--is_Firth_beta",          r_def=False,
      cpp=("yaml+json","is_Firth_beta"), cpp_def=False,     kind="bool", scope="both",
      effect="approximate-Firth BETA/SE for binary markers with p <= pCutoffforFirth",
      note="TRAP: the YAML key `isFirth` is NOT this switch. `isFirth` only "
           "controls the 'Firth approx was applied to N markers' log line "
           "(main.cpp:3874 -> mainMarkerInCPP t_isFirth). examples/"
           "step2_single.yaml sets `isFirth: true` and does NOT enable Firth."),
 dict(name="pCutoffforFirth", r="--pCutoffforFirth",        r_def=0.01,
      cpp=("yaml+json","pCutoffforFirth"), cpp_def=0.01,    kind="num",  scope="both",
      effect="p threshold below which Firth runs", note=""),

 # --- output shape (changes columns, not values) --------------------------
 dict(name="is_output_moreDetails", r="--is_output_moreDetails", r_def=False,
      cpp=("yaml","isMoreOutput"),      cpp_def=False,      kind="bool", scope="both",
      effect="adds hom/het counts by case/control to the output",
      note="renamed key. examples/step2_region.yaml writes "
           "`is_output_moreDetails` which the cpp parser ignores."),
 dict(name="markers_per_chunk", r="--markers_per_chunk",    r_def=10000,
      cpp=("yaml","marker_chunksize"),  cpp_def=10000,      kind="num",  scope="perf",
      effect="I/O batch size; must not change values", note="renamed key"),
 dict(name="nThreads",        r="--nThreads",               r_def=1,
      cpp=("yaml","nThreads"),          cpp_def=1,          kind="num",  scope="perf",
      effect="R>1 splits by idstoIncludeFile and concatenates; cpp uses OpenMP",
      note="ROW ORDER is nondeterministic in cpp region mode with nThreads>1. "
           "Always compare at nThreads=1."),

 # --- region / group tests ------------------------------------------------
 dict(name="groupFile",       r="--groupFile",              r_def="",
      cpp=("yaml","groupFile"),         cpp_def="",         kind="str",  scope="region",
      effect="non-empty switches the whole run from single-variant to region", note=""),
 dict(name="annotation_in_groupTest", r="--annotation_in_groupTest",
      r_def="lof,missense;lof,missense;lof;synonymous",
      cpp=("yaml","annotationList"),    cpp_def="(none; required)", kind="list", scope="region",
      effect="masks tested per gene",
      note="format differs: R = one comma-separated string with ';' unions, "
           "cpp = a YAML sequence of ';'-joined strings"),
 dict(name="maxMAF_in_groupTest", r="--maxMAF_in_groupTest", r_def="0.0001,0.001,0.01",
      cpp=("yaml","maxMAFList"),        cpp_def="(none; required)", kind="list", scope="region",
      effect="MAF cutoffs, one result row per cutoff", note="renamed key"),
 dict(name="maxMAC_in_groupTest", r="--maxMAC_in_groupTest", r_def="0",
      cpp=None,                          cpp_def=None,      kind="list", scope="region",
      effect="extra cutoffs expressed as MAC, folded into maxMAF list as MAC/(2N)",
      note="NOT PORTED. R default 0 = inactive, so this only bites when set."),
 dict(name="r.corr",          r="--r.corr",                 r_def=0.0,
      cpp=("yaml","r_corr"),            cpp_def=0.0,        kind="num",  scope="region",
      effect="0 = SKAT-O over rho grid {0,.01,.04,.09,.25,.5,1}; 1 = burden only",
      note="renamed key; cpp rejects anything but 0 or 1"),
 dict(name="MACCutoff_to_CollapseUltraRare", r="--MACCutoff_to_CollapseUltraRare", r_def=10.0,
      cpp=("yaml","MACCutoff_to_CollapseUltraRare"), cpp_def=10.0, kind="num", scope="region",
      effect="markers with MAC <= this are collapsed into one pseudo-marker", note=""),
 dict(name="markers_per_chunk_in_groupTest", r="--markers_per_chunk_in_groupTest", r_def=100,
      cpp=("yaml","markers_per_chunk_in_groupTest"), cpp_def=500, kind="num", scope="region",
      effect="chunk size for the P1/P2 covariance build",
      note="DEFAULT MISMATCH: R=100, cpp=500. Should be value-neutral (it is a "
           "blocking factor) but it changes float summation order, so set it "
           "equal on both sides before calling a residual difference real."),
 dict(name="groups_per_chunk", r="--groups_per_chunk",      r_def=100,
      cpp=("yaml","groups_per_chunk"),  cpp_def=100,        kind="num",  scope="perf",
      effect="genes read per I/O batch", note=""),
 dict(name="is_single_in_groupTest", r="--is_single_in_groupTest", r_def=False,
      cpp=("yaml","isSingleInGroupTest"), cpp_def=True,     kind="bool", scope="region",
      effect="also emit per-marker results from a region run",
      note="DEFAULT MISMATCH: R=FALSE, cpp=TRUE. Both force TRUE when r.corr=0 "
           "(SKAT-O), so it only bites with r.corr=1."),
 dict(name="is_output_markerList_in_groupTest", r="--is_output_markerList_in_groupTest", r_def=False,
      cpp=("yaml","isOutputMarkerList"), cpp_def=False,     kind="bool", scope="region",
      effect="writes the marker list per mask", note="renamed key"),
 dict(name="is_no_weight_in_groupTest", r="--is_no_weight_in_groupTest", r_def=False,
      cpp=None,                          cpp_def=None,      kind="bool", scope="region",
      effect="drop Beta(MAF,a,b) weights entirely", note="NOT PORTED"),
 dict(name="weights.beta",    r="--weights.beta",           r_def="1,25",
      cpp=("yaml","weights_beta"),      cpp_def="1,25",     kind="list", scope="region",
      effect="Beta(MAF, a, b) marker weights", note="renamed key"),
 dict(name="minGroupMAC_in_BurdenTest", r="--minGroupMAC_in_BurdenTest", r_def=5.0,
      cpp=("yaml","min_gourpmac_for_burdenonly"), cpp_def=5.0, kind="num", scope="region",
      effect="drop a burden pseudo-marker below this MAC (r.corr=1 only)",
      note="renamed key, and the cpp spelling has a typo ('gourpmac')"),

 # --- conditional analysis ------------------------------------------------
 dict(name="condition",       r="--condition",              r_def="",
      cpp=("yaml","condition"),         cpp_def="",         kind="str",  scope="both",
      effect="condition the score test on these markers; also force-disables "
             "is_fastTest on both sides",
      note="cpp accepts a comma string or a YAML sequence; cpp refuses "
           "conditional analysis for bgen/vcf input (main.cpp:4406)"),
 dict(name="weights_for_condition", r="--weights_for_condition", r_def=None,
      cpp=("yaml","weights_for_condition"), cpp_def=None,   kind="list", scope="region",
      effect="weights of the conditioning markers in region tests", note=""),

 # --- sparse GRM (OUT OF SCOPE this round, listed so it is not forgotten) --
 dict(name="sparseGRMFile",   r="--sparseGRMFile",          r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="both",
      effect="isSparseGRM=TRUE: enables the sparse-GRM variance path AND is a "
             "precondition for is_fastTest doing anything",
      note="OUT OF SCOPE this round. cpp loads a sparse GRM only from "
           "nullmodel.json-adjacent files, not from a step-2 config key."),
 dict(name="relatednessCutoff", r="--relatednessCutoff",    r_def=0.0,
      cpp=None,                          cpp_def=None,      kind="num",  scope="both",
      effect="sparse GRM thresholding", note="OUT OF SCOPE this round"),

 # --- chrX ----------------------------------------------------------------
 dict(name="sampleFile_male", r="--sampleFile_male",        r_def="",
      cpp=None,                          cpp_def=None,      kind="str",  scope="both",
      effect="male IDs for chrX non-PAR handling", note="NOT PORTED"),
 dict(name="X_PARregion",     r="--X_PARregion",            r_def="",
      cpp=("yaml","X_PARregion"),       cpp_def="",         kind="str",  scope="both",
      effect="PAR boundaries on chrX",
      note="NOT REACHABLE: main.cpp:149,152 declare g_is_rewrite_XnonPAR_forMales "
           "and g_X_PARregion_mat and nothing ever writes or reads them."),
 dict(name="is_rewrite_XnonPAR_forMales", r="--is_rewrite_XnonPAR_forMales", r_def=False,
      cpp=None,                          cpp_def=None,      kind="bool", scope="both",
      effect="double male dosages outside PAR", note="NOT PORTED"),
]
# fmt: on

# Branches with no flag at all -- taken purely on what the data looks like.
# These are why "same config" is not the same as "same code path": a parity
# run has to be shown to actually *enter* each of them.
DATA_BRANCHES = [
    dict(name="ER vs score/SPA", trigger="binary trait AND MAC <= max_MAC_for_ER",
         code="R Main.cpp/SAIGE_test.cpp getMarkerPval -> ER path; "
              "cpp main.cpp:493,1085,1279,1686 -> er_binary.cpp",
         how_to_hit="a marker with MAC <= 4; mid has none at minMAC=0.5, so "
                    "region/ultra-rare data is needed"),
    dict(name="SPA vs normal approx", trigger="|Tstat|/sqrt(var) > SPAcutoff AND trait != quantitative",
         code="SAIGE_test.cpp:545,556 / cpp saige_test.cpp:957",
         how_to_hit="any genome-wide scan on a binary trait; Is.SPA column says which"),
    dict(name="SPA root-finding fallback", trigger="SPA fails to converge",
         code="SPA_binary.cpp; reported via Is.SPA / p.value.NA",
         how_to_hit="extreme case-control imbalance at low MAC"),
    dict(name="sparse vs dense score variance", trigger="flagSparseGRM_cur, itself a "
         "function of is_fastTest and the marker's MAC bin",
         code="SAIGE_test.cpp:165-180 / cpp saige_test.cpp:1347-1510",
         how_to_hit="needs a sparse GRM -- out of scope this round"),
    dict(name="G-flip on AF>0.5", trigger="altFreq > 0.5 -> g := 2-g and alleles swap",
         code="UTIL.cpp:75-122 both sides",
         how_to_hit="always present in a real scan; this is the is_noadjCov bug site"),
    dict(name="categorical VR bin", trigger="which [minMACexcl, maxMACincl) bin the MAC falls in",
         code="saige_test.cpp:1347/1385/1463/1497 assignVarianceRatio*",
         how_to_hit="needs a cate-VR variance-ratio file (region runs)"),
    dict(name="ultra-rare collapsing", trigger="MAC <= MACCutoff_to_CollapseUltraRare inside a mask",
         code="cpp main.cpp:3061-3124; R SAIGE_SPATest_Region.R",
         how_to_hit="region run on rare variants"),
    dict(name="burden pseudo-marker drop", trigger="r.corr=1 AND group MAC < minGroupMAC_in_BurdenTest",
         code="cpp main.cpp:3330-3343",
         how_to_hit="burden-only region run on a tiny gene"),
    dict(name="Firth", trigger="binary AND is_Firth_beta AND p <= pCutoffforFirth",
         code="saige_test.cpp:1081,1088,1152 -> Firth.R / cpp Firth solver",
         how_to_hit="needs a genuinely associated marker, or raise pCutoffforFirth"),
    dict(name="missing-genotype imputation", trigger="any missing call in a marker",
         code="impute_method dispatch in the readers",
         how_to_hit="mid has ~1% missing, so every marker hits it"),
    dict(name="LOCO non-autosome fallback", trigger="chrom not in loco_chroms",
         code="R readInGLMM.R:107-113 / cpp null_model_loader.cpp guard 3",
         how_to_hit="--chrom=X on a LOCO model; both silently use the full fit"),
]

BY_NAME = {k["name"]: k for k in KNOBS}
MISMATCH = [k["name"] for k in KNOBS if "DEFAULT MISMATCH" in k["note"]]
UNPORTED = [k["name"] for k in KNOBS if k["cpp"] is None]


def _fmt(v):
    if v is True:
        return "TRUE"
    if v is False:
        return "FALSE"
    if v is None:
        return "-"
    return str(v)


def print_table():
    hdr = ("knob", "R flag", "R default", "C++ location", "C++ default", "note")
    rows = [hdr]
    for k in KNOBS:
        cpp = "-" if k["cpp"] is None else "%s:%s" % k["cpp"]
        rows.append((k["name"], k["r"] or "-", _fmt(k["r_def"]), cpp,
                     _fmt(k["cpp_def"]), k["note"].split(".")[0] if k["note"] else ""))
    w = [max(len(r[i]) for r in rows) for i in range(len(hdr))]
    for n, r in enumerate(rows):
        print("  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip())
        if n == 0:
            print("  ".join("-" * w[i] for i in range(len(hdr))))
    print("\ndefault mismatches (%d): %s" % (len(MISMATCH), ", ".join(MISMATCH)))
    print("not ported to C++ (%d): %s" % (len(UNPORTED), ", ".join(UNPORTED)))
    print("\ndata-triggered branches (no flag):")
    for b in DATA_BRANCHES:
        print("  %-32s %s" % (b["name"], b["trigger"]))


if __name__ == "__main__":
    print_table()
