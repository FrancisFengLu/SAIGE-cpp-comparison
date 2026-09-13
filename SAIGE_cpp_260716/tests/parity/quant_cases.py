#!/usr/bin/env python3
"""Single-variant QUANTITATIVE parity cases (R SAIGE 1.5.2 vs C++ step 2).

Registered as an overlay on parity.py rather than inside it so that several
agents can add case groups without colliding on one file.  Usage is the same:

    ./quant_cases.py --list
    ./quant_cases.py base_quant_single
    ./quant_cases.py --all
    ./quant_cases.py base_quant_single --desync cpp:impute_method=minor

MODEL
-----
/opt/saige/logs/step2mt/models_q8/q1.rda -- SAIGE step-1 fit on
/opt/saige/data/mid (N=50000, M=40000) with /opt/saige/data/mid.q8.pheno.txt
column q1, traitType=quantitative, covariates x1,x2, invNormalize=FALSE,
LOCO=FALSE, single (non-categorical) variance ratio, theta=[1.00222361087799,
0.0378887318074703].  V = 1/tau[1] = 0.9977813 for every sample (checked
against modglmm$obj.noK$V: identical).  The .arma directory
/opt/saige/logs/step2/null_q1_arma was produced from that same .rda by
tools/rda_to_arma.R, so both sides read one fit.

WHY THE QUANTITATIVE PATH IS NOT JUST "BINARY WITH A DIFFERENT V"
----------------------------------------------------------------
Three branches are supposed to be unreachable for a quantitative trait:
SPA (SAIGE_test.cpp / saige_test.cpp gate on traitType), efficient
resampling (binary-only), and approximate Firth (binary-only).  The cases
quant_single_spa1 / _er / _firth set those knobs to values that would visibly
move a binary run and assert that the quantitative output does not move --
that is a test of the GATE, which a "both sides agree" run alone would not
give (both could be wrong in the same direction, but they are independent
implementations of the gate).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import parity as P                                     # noqa: E402

# Registered defensively with setdefault: parity.py may already carry these
# (another case group adds the same entries).  Re-registering identical values
# is a no-op; this file stays runnable on its own either way.
P.DATASETS.setdefault("rare", dict(kind="plink", prefix="/opt/saige/data/rare"))
P.MODELS["q1_quant"] = dict(
    trait="quantitative",
    rda="/opt/saige/logs/step2mt/models_q8/q1.rda",
    arma="/opt/saige/logs/step2/null_q1_arma",
    vr="/opt/saige/logs/step2mt/models_q8/q1.varianceRatio.txt")

Q = dict(data="mid2k", model="q1_quant")

QCASES = {
    # ---- the configuration that must agree ------------------------------
    "base_quant_single": dict(
        desc="quantitative, single-variant, non-LOCO, is_noadjCov=FALSE, "
             "impute_method=mean -- the reference configuration",
        set={}, **Q),

    # ---- QC filters ------------------------------------------------------
    "quant_single_minmaf": dict(
        desc="minMAF=0.2 on both sides -- marker drop filter",
        set=dict(minMAF=0.2), **Q),
    "quant_single_minmac": dict(
        desc="minMAC=8000 on both sides -- the other drop filter",
        set=dict(minMAC=8000.0), **Q),
    "quant_single_maxmissing": dict(
        desc="maxMissing=0.0105 -- mid was generated with 1%% missing and the "
             "per-marker rate spans 0.0084..0.0117, so this cutoff splits the "
             "2000 markers instead of dropping all of them",
        set=dict(maxMissing=0.0105), **Q),

    # ---- imputation ------------------------------------------------------
    "quant_single_bestguess": dict(
        desc="impute_method=best_guess (R's own CLI default)",
        set=dict(impute_method="best_guess"), **Q),
    "quant_single_minorimpute": dict(
        desc="impute_method=minor",
        set=dict(impute_method="minor"), **Q),

    # ---- reader / orientation -------------------------------------------
    "quant_single_alleleorder": dict(
        desc="AlleleOrder=ref-first -- swaps which bim allele is Allele2 and "
             "hence the sign of BETA",
        set=dict(AlleleOrder="ref-first"), **Q),

    # ---- output shape ----------------------------------------------------
    "quant_single_moredetails": dict(
        desc="is_output_moreDetails -- the extra columns must match too",
        set=dict(is_output_moreDetails=True), **Q),
    "quant_single_imputed": dict(
        desc="is_imputed_data=TRUE on hard-call PLINK input -- adds the "
             "imputationInfo column and arms the dosage-zeroing branch",
        set=dict(is_imputed_data=True), **Q),

    # ---- perf knobs that must be value-neutral --------------------------
    "quant_single_chunk": dict(
        desc="markers_per_chunk=1000 -- an I/O blocking factor, must not move "
             "a single number. 1000 is the smallest value R accepts "
             "(checkArgs.R checkArgNumeric: 'markers_per_chunk should be a "
             "numeric value greater than or equal to 1000'); C++ has no such "
             "floor, so anything below 1000 is not comparable at all",
        set=dict(markers_per_chunk=1000), **Q),
    "quant_single_dosagezero": dict(
        desc="is_imputed_data=TRUE at R's legal MAXIMUM for the zeroing knobs "
             "(dosage_zerod_cutoff=0.5, dosage_zerod_MAC_cutoff=100; "
             "checkArgs.R:62-63 caps them at 0.5 and 100, and C++ enforces no "
             "such cap). On hard calls the only values that can be <= 0.5 and "
             "non-zero are the MEAN-imputed missing cells (2*AF ~ 2e-3 on the "
             "rarest markers), so this is the only way to enter the branch "
             "with a value change from PLINK input",
        set=dict(is_imputed_data=True, dosage_zerod_cutoff=0.5,
                 dosage_zerod_MAC_cutoff=100.0), **Q),

    # ---- branches that are supposed to be OFF for a quantitative trait ---
    "quant_single_spa1": dict(
        desc="SPAcutoff=1. On the binary model this moves 434 markers onto "
             "the SPA branch; for quantitative the gate should make it a "
             "no-op on BOTH sides",
        set=dict(SPAcutoff=1.0), **Q),
    "quant_single_er": dict(
        desc="max_MAC_for_ER=100000 -- efficient resampling is binary-only, "
             "so this must be a no-op on both sides",
        set=dict(max_MAC_for_ER=100000.0), **Q),
    "quant_single_firth": dict(
        desc="is_Firth_beta=TRUE with pCutoffforFirth=0.5 -- Firth is "
             "binary-only, so this must be a no-op on both sides",
        set=dict(is_Firth_beta=True, pCutoffforFirth=0.5), **Q),
    "quant_single_fasttest": dict(
        desc="is_fastTest=TRUE with no sparse GRM. R force-disables it "
             "(SAIGE_Test_main.R); cpp detects the identical-context "
             "recompute and skips it (main.cpp:1618-1653)",
        set=dict(is_fastTest=True), **Q),

    # ---- rare spectrum: the only place the low-MAC branches are reachable
    # `rare` is 3000 markers on the SAME 50000 individuals as `mid`
    # (tests/parity/make_rare.py), 600 of them at MAC<=4.  `mid2k`'s smallest
    # MAC is 105, so on mid2k every low-MAC branch is dead code.
    "quant_rare_single": dict(
        desc="baseline config on the rare-spectrum dataset -- 600 markers at "
             "MAC<=4, which on a BINARY trait would all take the exact "
             "resampling path; for quantitative none of them may",
        data="rare", model="q1_quant", set={}),
    "quant_rare_er": dict(
        desc="max_MAC_for_ER=1000 on the rare dataset: 2/3 of the markers are "
             "under that MAC. Efficient resampling is binary-only, so both "
             "sides must produce exactly the quant_rare_single numbers",
        data="rare", model="q1_quant", set=dict(max_MAC_for_ER=1000.0)),
    "quant_rare_dosagezero": dict(
        desc="is_imputed_data=TRUE, dosage_zerod_cutoff=0.5, "
             "dosage_zerod_MAC_cutoff=100 (R's legal maximum) on the rare "
             "dataset -- here markers with MAC<=100 DO exist, so the "
             "mean-imputed missing cells actually get zeroed",
        data="rare", model="q1_quant",
        set=dict(is_imputed_data=True, dosage_zerod_cutoff=0.5,
                 dosage_zerod_MAC_cutoff=100.0)),
    "quant_rare_bestguess": dict(
        desc="impute_method=best_guess on the rare dataset: rounding 2*AF to "
             "the nearest hard call sends every imputed cell to 0, which is a "
             "much bigger move than on mid2k's common markers",
        data="rare", model="q1_quant", set=dict(impute_method="best_guess")),

    # ---- the known R defect, reached deliberately ------------------------
    "quant_single_noadjcov": dict(
        desc="is_noadjCov=TRUE on BOTH sides. Same experiment as the binary "
             "case: does the C++ port reproduce R's noadjCov arithmetic "
             "(AF>0.5 centring defect included) for a quantitative trait too?",
        set=dict(is_noadjCov=True), vr_add_noXadj=True, **Q),

    # ---- conditional analysis -------------------------------------------
    "quant_single_condition": dict(
        desc="conditional on snp0,snp1 -- the first two markers in the .bed. "
             "MEASURED: R mis-reads the marker at file index 0 (see "
             "quant_single_condition_mid for the isolated mechanism)",
        set=dict(condition="snp0,snp1"), expect="differ", **Q),
    "quant_single_condition_mid": dict(
        desc="conditional on snp500,snp1200 -- neither is the marker at file "
             "index 0, so this separates 'conditioning marker' from 'marker 0'. "
             "R's PlinkClass::getOneMarker (src/PLINK.cpp:195-207) only seeks "
             "when t_gIndex > 0, so after the conditioning markers have moved "
             "the .bed file pointer the main loop's FIRST marker (gIndex==0) is "
             "read from wherever the pointer landed, while its CHR/POS/ID/"
             "alleles still come from index 0. The C++ port forces an absolute "
             "SEEK_SET whenever t_gIndex_prev==0 (genotype_reader.cpp:832-844) "
             "and is unaffected.",
        set=dict(condition="snp500,snp1200"), expect="differ", **Q),
}

P.CASES.update(QCASES)


def main():
    # Reuse parity.py's whole CLI, but list only this group when --list is
    # given without a case.
    argv = sys.argv[1:]
    if "--list" in argv or (not argv):
        print("quantitative single-variant cases:")
        for n, c in QCASES.items():
            print("  %-28s %s" % (n, c["desc"].split("--")[0].strip()))
        print("\nmodel: q1_quant  (%s)" % P.MODELS["q1_quant"]["rda"])
        return 0
    if "--all" in argv:
        argv = [a for a in argv if a != "--all"]
        rc = 0
        for n in QCASES:
            sys.argv = [sys.argv[0], n] + argv
            rc = max(rc, P.main() or 0)
        return rc
    sys.argv = [sys.argv[0]] + argv
    return P.main()


if __name__ == "__main__":
    sys.exit(main())
