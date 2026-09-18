#!/usr/bin/env python3
"""fit.selective_geno_load A/B gate (run_ab_gate.sh is the entry point).

Runs ONE binary twice per case — selective_geno_load false, then true, every
other setting identical — and requires every output file to be byte-identical
except <out_prefix>.grm_diag.txt, which the selective run deliberately does not
write (it needs every marker; step 2 never reads it).

That is the whole correctness claim of the feature: the variance-ratio marker
draw happens before any BED byte is read and a marker's VR eligibility depends
on that marker alone, so decoding only the drawn markers must reproduce the VR
pool, the mt19937(200) shuffle over it, the variance ratio, tau and every model
artifact bit for bit. A difference here is a bug, not a tolerance question.

The rf_* cases are the other half: configurations where the flag must be
REFUSED (the genotype matrix really is needed). There the two runs must be
identical INCLUDING .grm_diag.txt, and the log must say why.

  ab_gate.py <saige-null> [--workdir DIR] [--cases a,b] [--scale small|mid]

Inputs come from tests/mt_tolerance/make_cases.py (small scale) plus the sparse
GRMs under /opt/saige/logs/missing_mt/step1/data.
"""
import argparse
import filecmp
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

SPARSE = dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=False)
SPARSE_VR = dict(SPARSE, use_sparse_grm_for_vr=True)
CATE = dict(isCateVarianceRatio=True, cateVarRatioMinMACVecExclude=[10.5, 20.5],
            cateVarRatioMaxMACVecInclude=[20.5])

# name -> spec. "expect": "on" = the flag must be used, "off" = must be refused.
CASES = {
    "sp_bin":      dict(data="sparse", y="sA", trait="binary", fit=SPARSE, expect="on"),
    "sp_qt":       dict(pheno=f"{SPDATA}/small.sp.full.pheno.txt", y="sq1",
                        trait="quantitative", fit=SPARSE, expect="on"),
    "sp_bin_vr":   dict(data="sparse", y="sA", trait="binary", fit=SPARSE_VR, expect="on"),
    "sp_qt_vr":    dict(pheno=f"{SPDATA}/small.sp.full.pheno.txt", y="sq1",
                        trait="quantitative", fit=SPARSE_VR, expect="on"),
    # categorical VR: the CV rule walks the pool far past the 30 it starts with
    # (on this data it exhausts all ~950 members), which is the strongest check
    # that the selectively-loaded pool is the same pool.
    "sp_cate_bin": dict(data="vrcate", y="cS1", trait="binary",
                        fit=dict(SPARSE_VR, **CATE), expect="on"),
    # several sample-set groups: the flag has to survive
    # reset_step1_state_for_new_sample_set() between groups.
    "sp_bin_mt":   dict(data="sparse", ys=["sA", "sB", "sC", "sD"], trait="binary",
                        fit=SPARSE_VR, expect="on"),
    "rf_pcg_bin":  dict(data="sparse", y="sA", trait="binary",
                        fit=dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=True),
                        expect="off"),
    "rf_vronly":   dict(data="vrcate", y="cS1", trait="binary",
                        fit=dict(use_sparse_grm_to_fit=False, use_sparse_grm_for_vr=True, **CATE),
                        expect="off"),
    "rf_novr":     dict(data="sparse", y="sA", trait="binary",
                        fit=dict(SPARSE, num_markers_for_vr=0), expect="off"),
}


def config_for(spec, workdir, rundir, selective, scale):
    fit = dict(gate.FIT_DEFAULTS)
    fit.update(trait=spec["trait"], nthreads=1, multi_lockstep=False)
    fit.update(spec.get("fit", {}))
    fit["selective_geno_load"] = bool(selective)
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("binary")
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--outdir", default="/opt/saige/logs/selgeno_gate")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--scale", default="small")
    a = ap.parse_args(argv)
    binary = os.path.abspath(a.binary)
    cases = [c for c in a.cases.split(",") if c]
    for c in cases:
        if c not in CASES:
            raise SystemExit(f"unknown case {c}; known: {list(CASES)}")
    os.makedirs(a.outdir, exist_ok=True)
    gate.log(f"selective_geno_load A/B: {binary} ({gate.md5_file(binary)}); {len(cases)} cases")

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
            used = "[selective] fit.selective_geno_load ON" in logs["on"]
            refused = "requested but NOT used" in logs["on"]
            A, B = tree(dirs["off"]), tree(dirs["on"])
            only_off, only_on = sorted(A - B), sorted(B - A)
            differ = sorted(f for f in A & B
                            if not filecmp.cmp(os.path.join(dirs["off"], f),
                                               os.path.join(dirs["on"], f), shallow=False))
            if spec["expect"] == "on":
                ok = (used and not refused and not only_on and not differ
                      and only_off and all(f.endswith(".grm_diag.txt") for f in only_off))
                why = ("flag not used" if not used else
                       f"only_on={only_on} differ={differ} only_off={only_off}")
            else:
                ok = (refused and not used and not only_on and not only_off and not differ)
                why = ("flag was used but must be refused" if used else
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
