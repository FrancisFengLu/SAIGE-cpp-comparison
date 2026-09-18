#!/usr/bin/env python3
"""fit.fused_variance_ratio A/B gate (run_ab_gate.sh is the entry point).

Runs ONE binary twice per case -- fused_variance_ratio false, then true, every
other setting identical.

What is claimed, and therefore what is checked:

  fu_*   the flag is USED and changes ONLY the variance-ratio products.
         The null model itself must be byte-identical: the variance ratio is
         computed after the fit and never feeds back into tau, so a difference
         in nullmodel.json / the model files is a bug, not a tolerance question.
         Files allowed to differ or to be missing on one side:
           *.varianceRatio.txt              (the point of the feature)
           *markers.SAIGE.results.txt       (its name carries the marker count)

  rf_*   the flag must be REFUSED, with the reason printed, and then the two
         runs are byte-identical everywhere including the variance ratio.

  ab_gate.py <saige-null> [--workdir DIR] [--cases a,b] [--scale small|mid]
"""
import argparse
import filecmp
import fnmatch
import os
import shutil
import subprocess
import sys
import time

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
MT = os.path.join(os.path.dirname(HERE), "mt_tolerance")
sys.path.insert(0, MT)
import gate  # noqa: E402

SPDATA = "/opt/saige/logs/missing_mt/step1/data"
QT_PHENO = f"{SPDATA}/small.sp.full.pheno.txt"

SPARSE = dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=False)
SPARSE_VR = dict(SPARSE, use_sparse_grm_for_vr=True)
CATE = dict(isCateVarianceRatio=True, cateVarRatioMinMACVecExclude=[10.5, 20.5],
            cateVarRatioMaxMACVecInclude=[20.5])

VR_FILES = ["*.varianceRatio.txt", "*markers.SAIGE.results.txt"]

CASES = {
    # the main path: quantitative + sparse fit, 10 markers per bin for delta
    "fu_qt":       dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=SPARSE, fused=dict(fused_vr_markers=10), expect="on"),
    # zero marker budget: pure closed form, no PCG solve and no genotype read
    "fu_qt_zero":  dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=SPARSE, fused=dict(fused_vr_markers=0), expect="on",
                        needs=["marker budget 0"]),
    # with the sparse VR row on, the per-marker sparse ratio still has to be sampled
    "fu_qt_vr":    dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=SPARSE_VR, fused=dict(fused_vr_markers=10), expect="on"),
    # MAC bins: one anchor, one delta test per bin
    "fu_qt_cate":  dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=dict(SPARSE, **CATE), fused=dict(fused_vr_markers=10),
                        expect="on"),
    # several sample-set groups: the anchor is rebuilt per group after
    # reset_step1_state_for_new_sample_set()
    "fu_qt_mt":    dict(pheno=QT_PHENO, ys=["sq1", "sq2"], trait="quantitative",
                        fit=SPARSE, fused=dict(fused_vr_markers=10), expect="on"),
    # forces the "keep delta" branch (z = 0 makes any |delta-1| > 0 significant),
    # which the data itself never triggers: the written ratio must then be the
    # sampled bin mean, not the anchor.
    "fu_qt_keep":  dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=SPARSE, fused=dict(fused_vr_markers=10, fused_vr_delta_z=0.0),
                        expect="on", needs=["keep delta (bin differs from the anchor)"]),
    # combined with the other new flag
    "fu_qt_sel":   dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=dict(SPARSE, selective_geno_load=True),
                        fused=dict(fused_vr_markers=10), expect="on"),

    "rf_bin":      dict(data="sparse", y="sA", trait="binary", fit=SPARSE,
                        fused=dict(fused_vr_markers=10), expect="off",
                        needs=["the trait is binary"]),
    "rf_dense":    dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=dict(use_sparse_grm_to_fit=False, use_sparse_grm_for_vr=True),
                        fused=dict(fused_vr_markers=10), expect="off",
                        needs=["use_sparse_grm_to_fit is off"]),
    "rf_vr0":      dict(pheno=QT_PHENO, y="sq1", trait="quantitative",
                        fit=SPARSE_VR, fused=dict(fused_vr_markers=0), expect="off",
                        needs=["use_sparse_grm_for_vr"]),
}


def config_for(spec, workdir, rundir, fused, scale):
    fit = dict(gate.FIT_DEFAULTS)
    fit.update(trait=spec["trait"], nthreads=1, multi_lockstep=False)
    fit.update(spec.get("fit", {}))
    fit["fused_variance_ratio"] = bool(fused)
    if fused:
        fit.update(spec.get("fused", {}))
    pheno = spec.get("pheno") or gate.load_manifest(workdir, scale, spec["data"])["pheno"]
    paths = {"plinkFile": f"/opt/saige/data/{scale}",
             "out_prefix": f"{rundir}/m", "out_prefix_vr": f"{rundir}/mvr",
             "sparse_grm": f"{SPDATA}/{scale}.sgrm2.mtx",
             "sparse_grm_ids": f"{SPDATA}/{scale}.sgrm2.ids",
             "overwrite_varratio": True}
    design = {"csv": pheno, "iid_col": "IID", "covar_cols": ["x1", "x2"]}
    if "ys" in spec:
        design["y_cols"] = spec["ys"]
    else:
        design["y_col"] = spec["y"]
    return {"paths": paths, "design": design, "fit": fit}


def tree(root):
    return {os.path.relpath(os.path.join(d, f), root)
            for d, _, fs in os.walk(root) for f in fs}


def is_vr(f):
    return any(fnmatch.fnmatch(os.path.basename(f), pat) for pat in VR_FILES)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("binary")
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--outdir", default="/opt/saige/logs/fusedvr_gate")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--scale", default="small")
    a = ap.parse_args(argv)
    binary = os.path.abspath(a.binary)
    cases = [c for c in a.cases.split(",") if c]
    for c in cases:
        if c not in CASES:
            raise SystemExit(f"unknown case {c}; known: {list(CASES)}")
    os.makedirs(a.outdir, exist_ok=True)
    gate.log(f"fused_variance_ratio A/B: {binary} ({gate.md5_file(binary)}); {len(cases)} cases")

    results = []
    for c in cases:
        spec = CASES[c]
        secs, logs, dirs = {}, {}, {}
        for mode in ("off", "on"):
            rd = os.path.join(a.outdir, f"{c}_{mode}")
            shutil.rmtree(rd, ignore_errors=True)
            os.makedirs(rd)
            cfg = config_for(spec, a.workdir, rd, mode == "on", a.scale)
            cfg_path = os.path.join(a.outdir, f"{c}_{mode}.yaml")
            with open(cfg_path, "w") as f:
                yaml.safe_dump(cfg, f, sort_keys=False)
            t0 = time.time()
            with open(rd + ".log", "w") as lf:
                rc = subprocess.run([binary, "-c", cfg_path], cwd=rd, stdout=lf,
                                    stderr=subprocess.STDOUT, env=gate.clean_env()).returncode
            secs[mode] = time.time() - t0
            dirs[mode] = rd
            logs[mode] = open(rd + ".log").read()
            if rc != 0:
                gate.log(f"{c}/{mode}: EXIT {rc} — see {rd}.log")
                results.append((c, False, f"exit {rc} in the {mode} run"))
                break
        else:
            used = "[fusedVR] ON" in logs["on"]
            refused = "requested but NOT used" in logs["on"]
            missing = [s for s in spec.get("needs", []) if s not in logs["on"]]
            A, B = tree(dirs["off"]), tree(dirs["on"])
            only_off, only_on = sorted(A - B), sorted(B - A)
            differ = sorted(f for f in A & B
                            if not filecmp.cmp(os.path.join(dirs["off"], f),
                                               os.path.join(dirs["on"], f), shallow=False))
            if spec["expect"] == "on":
                bad_diff = [f for f in differ if not is_vr(f)]
                bad_only = [f for f in only_off + only_on if not is_vr(f)]
                vr_moved = [f for f in differ + only_off + only_on if is_vr(f)]
                ok = (used and not refused and not missing and not bad_diff
                      and not bad_only and vr_moved)
                why = ("flag not used" if not used else
                       f"missing log lines {missing} " if missing else
                       f"non-VR files moved: differ={bad_diff} only={bad_only}"
                       if (bad_diff or bad_only) else "the variance ratio did not change")
            else:
                ok = (refused and not used and not missing and not only_on
                      and not only_off and not differ)
                why = ("flag was used but must be refused" if used else
                       f"missing log lines {missing}" if missing else
                       f"only_off={only_off} only_on={only_on} differ={differ}")
            gate.log(f"{c}: {'PASS' if ok else 'FAIL'}  "
                     f"[{len(A)} files, flag {'used' if used else 'refused'}]  "
                     f"off {secs['off']:.1f}s -> on {secs['on']:.1f}s"
                     + ("" if ok else f"  ({why})"))
            results.append((c, ok, "" if ok else why))

    bad = [c for c, ok, _ in results if not ok]
    gate.log(f"{len(results) - len(bad)}/{len(results)} cases pass"
             + (f"; FAILED: {bad}" if bad else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
