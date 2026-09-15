#!/usr/bin/env python3
"""P=1 byte-identity regression gate (run_p1_byte_gate.sh is the entry point).

A multi-trait-only change must leave single-phenotype runs untouched, so every
output file of every case below must be byte-identical between NEW and BASE
(same file set, same bytes). GPU cases must show "[gpu_matvec] tier=4" in both
logs (a silent CPU fallback would otherwise compare CPU with CPU and pass);
CPU cases must show no gpu_matvec handle. nthreads is 1 unless the case says
otherwise (CPU runs with nthreads > 1 are not reproducible).

  p1_byte_gate.py NEW BASE [--workdir DIR] [--cases a,b] [--label NAME] [--no-cache]

BASE outputs are cached by binary md5 + rendered config; NEW always runs fresh.
Inputs come from make_cases.py (small scale; mid for the *_mid cases).
"""
import argparse
import datetime
import filecmp
import hashlib
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gate  # noqa: E402

SPDATA = "/opt/saige/logs/missing_mt/step1/data"

# name: (device, scale, data_case or external pheno, y_col, trait, fit overrides, extra)
CASES = {
    "bin_cpu":            dict(dev="cpu", data="nomiss", y="b1", trait="binary"),
    "bin_gpu":            dict(dev="gpu", data="nomiss", y="b1", trait="binary"),
    "qt_cpu":             dict(dev="cpu", data="quant", y="qtF1", trait="quantitative"),
    "qt_gpu":             dict(dev="gpu", data="quant", y="qtF1", trait="quantitative"),
    "missing_bin_cpu":    dict(dev="cpu", data="block", y="y9", trait="binary"),
    "missing_bin_gpu":    dict(dev="gpu", data="block", y="y9", trait="binary"),
    "loco_bin_cpu":       dict(dev="cpu", data="loco", y="lF1", trait="binary", plink="loco", fit=dict(loco=True)),
    "loco_bin_gpu":       dict(dev="gpu", data="loco", y="lF1", trait="binary", plink="loco", fit=dict(loco=True)),
    "loco_qt_gpu":        dict(dev="gpu", data="quant", y="qtF1", trait="quantitative", plink="smallloco",
                               fit=dict(loco=True)),
    "sparse_direct_bin_cpu": dict(dev="cpu", data="sparse", y="sA", trait="binary", sparse=True,
                                  fit=dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=False)),
    "sparse_pcg_bin_cpu": dict(dev="cpu", data="sparse", y="sA", trait="binary", sparse=True,
                               fit=dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=True)),
    "sparse_direct_qt_cpu": dict(dev="cpu", pheno=f"{SPDATA}/small.sp.full.pheno.txt", y="sq1",
                                 trait="quantitative", sparse=True,
                                 fit=dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=False)),
    "sparse_pcg_qt_cpu":  dict(dev="cpu", pheno=f"{SPDATA}/small.sp.full.pheno.txt", y="sq1",
                               trait="quantitative", sparse=True,
                               fit=dict(use_sparse_grm_to_fit=True, use_pcg_with_sparse_grm=True)),
    "nrun2_cpu":          dict(dev="cpu", data="nomiss", y="b1", trait="binary", fit=dict(nrun=2)),
    "nrun2_gpu":          dict(dev="gpu", data="nomiss", y="b1", trait="binary", fit=dict(nrun=2)),
    "covoff_off_cpu":     dict(dev="cpu", data="nomiss", y="b1", trait="binary", fit=dict(covariate_offset=False)),
    "covoff_off_gpu":     dict(dev="gpu", data="nomiss", y="b1", trait="binary", fit=dict(covariate_offset=False)),
    "qr_off_cpu":         dict(dev="cpu", data="nomiss", y="b1", trait="binary", fit=dict(covariate_qr=False)),
    "qr_off_gpu":         dict(dev="gpu", data="nomiss", y="b1", trait="binary", fit=dict(covariate_qr=False)),
    "nocov_cpu":          dict(dev="cpu", data="nomiss", y="b1", trait="binary", covars=[]),
    "nocov_gpu":          dict(dev="gpu", data="nomiss", y="b1", trait="binary", covars=[]),
    "novr_cpu":           dict(dev="cpu", data="nomiss", y="b1", trait="binary", fit=dict(num_markers_for_vr=0)),
    "novr_gpu":           dict(dev="gpu", data="nomiss", y="b1", trait="binary", fit=dict(num_markers_for_vr=0)),
    "vrcate_cpu":         dict(dev="cpu", data="vrcate", y="cS1", trait="binary", sparse=True,
                               fit=dict(use_sparse_grm_for_vr=True, isCateVarianceRatio=True,
                                        cateVarRatioMinMACVecExclude=[10.5, 20.5],
                                        cateVarRatioMaxMACVecInclude=[20.5])),
    "mid_bin_gpu":        dict(dev="gpu", scale="mid", data="nomiss", y="b1", trait="binary"),
    "mid_bin_gpu_t8":     dict(dev="gpu", scale="mid", data="nomiss", y="b1", trait="binary", nthreads=8),
    "mid_qt_gpu":         dict(dev="gpu", scale="mid", data="quant", y="qtF1", trait="quantitative"),
}


def config_for(spec, workdir, rundir):
    scale = spec.get("scale", "small")
    fit = dict(gate.FIT_DEFAULTS)
    fit.update(trait=spec["trait"], nthreads=spec.get("nthreads", 1), multi_lockstep=False)
    fit.update(spec.get("fit", {}))
    inputs = {}
    if "data" in spec:
        man = gate.load_manifest(workdir, scale, spec["data"])
        pheno = man["pheno"]
        inputs = dict(man["inputs"], pheno=man["pheno_md5"])
    else:
        man = None
        pheno = spec["pheno"]
        inputs = {pheno: gate.md5_file(pheno)}
    plink = f"/opt/saige/data/{scale}"
    if spec.get("plink") == "loco":
        plink = man["loco_plink"]
    elif spec.get("plink") == "smallloco":
        plink = "/opt/saige/data/smallloco"
    paths = {"plinkFile": plink, "out_prefix": f"{rundir}/m", "out_prefix_vr": f"{rundir}/mvr",
             "overwrite_varratio": True}
    if spec.get("sparse"):
        paths["sparse_grm"] = f"{SPDATA}/{scale}.sgrm2.mtx"
        paths["sparse_grm_ids"] = f"{SPDATA}/{scale}.sgrm2.ids"
        inputs[paths["sparse_grm"]] = gate.md5_file(paths["sparse_grm"])
    design = {"csv": pheno, "iid_col": "IID", "y_col": spec["y"], "covar_cols": spec.get("covars", ["x1", "x2"])}
    return {"paths": paths, "design": design, "fit": fit}, inputs


def tree_files(root):
    out = set()
    for d, _, fs in os.walk(root):
        for f in fs:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("new")
    ap.add_argument("base")
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--cases", default=",".join(CASES))
    ap.add_argument("--label")
    ap.add_argument("--no-cache", action="store_true")
    a = ap.parse_args(argv)
    new, base = os.path.abspath(a.new), os.path.abspath(a.base)
    mn, mb = gate.md5_file(new), gate.md5_file(base)
    label = a.label or f"{mn[:8]}_vs_{mb[:8]}"
    root = os.path.join(a.workdir, "p1", label)
    os.makedirs(root, exist_ok=True)
    cases = [c for c in a.cases.split(",") if c]
    for c in cases:
        if c not in CASES:
            raise SystemExit(f"unknown case {c}; known: {list(CASES)}")
    gate.log(f"P=1 byte gate: new {new} ({mn}) vs base {base} ({mb}); {len(cases)} cases")
    results = []
    for c in cases:
        spec = CASES[c]
        gpu = spec["dev"] == "gpu"
        cfg0, inputs = config_for(spec, a.workdir, "@OUT@")
        key = hashlib.sha256(json.dumps(dict(v=1, bin=mb, gpu=gpu, cfg=cfg0, inputs=inputs),
                                        sort_keys=True).encode()).hexdigest()[:20]
        bcache = os.path.join(a.workdir, "p1_cache", f"{c}-{key}")
        bdir = os.path.join(bcache, "run")
        if a.no_cache or not os.path.exists(os.path.join(bcache, "DONE")):
            if os.path.exists(bcache):
                shutil.rmtree(bcache)
            os.makedirs(bcache)
            cfg, _ = config_for(spec, a.workdir, bdir)
            gate.log(f"{c}: base run ...")
            rc_b, sb = gate.run_saige(base, cfg, bdir, gpu)
            json.dump(dict(rc=rc_b, seconds=sb), open(os.path.join(bcache, "meta.json"), "w"))
            if rc_b == 0:
                open(os.path.join(bcache, "DONE"), "w").write("ok\n")
        else:
            meta = json.load(open(os.path.join(bcache, "meta.json")))
            rc_b, sb = meta["rc"], meta["seconds"]
            gate.log(f"{c}: base cached")
        ndir = os.path.join(root, c, "new")
        cfg, _ = config_for(spec, a.workdir, ndir)
        gate.log(f"{c}: new run ...")
        rc_n, sn = gate.run_saige(new, cfg, ndir, gpu)
        problems = []
        if rc_b != 0 or rc_n != 0:
            problems.append(f"exit codes base={rc_b} new={rc_n}")
        fb, fn = tree_files(bdir), tree_files(ndir)
        if fb != fn:
            problems.append(f"file sets differ: base-only {sorted(fb - fn)[:5]} new-only {sorted(fn - fb)[:5]}")
        diff = sorted(r for r in fb & fn if not filecmp.cmp(os.path.join(bdir, r), os.path.join(ndir, r), shallow=False))
        if diff:
            problems.append(f"{len(diff)} file(s) differ: {diff[:8]}")
        if not fb:
            problems.append("no output files")
        logs = {"base": gate.read_text(bdir + ".log"), "new": gate.read_text(ndir + ".log")}
        for side, txt in logs.items():
            has4 = "[gpu_matvec] tier=4" in txt
            if gpu and not has4:
                problems.append(f"{side}: GPU case without '[gpu_matvec] tier=4'")
            if not gpu and "[gpu_matvec] tier=" in txt:
                problems.append(f"{side}: CPU case created a GPU handle")
            if spec["fit"].get("loco") if "fit" in spec else False:
                if "LOCO: on" not in txt:
                    problems.append(f"{side}: LOCO requested but log lacks 'LOCO: on'")
        ok = not problems
        results.append(dict(case=c, ok=ok, files=len(fb & fn), base_s=sb, new_s=sn, problems=problems))
        gate.log(f"{c}: {'IDENTICAL' if ok else 'DIFFER'} ({len(fb & fn)} files) {'; '.join(problems)}")
    lines = [f"== P=1 byte gate {label}", f"new  {new} ({mn})", f"base {base} ({mb})", ""]
    for r in results:
        lines.append(f"  {'IDENTICAL' if r['ok'] else 'DIFFER   '}  {r['case']:<24} {r['files']:>4} files  "
                     f"base {r['base_s']:6.0f}s new {r['new_s']:6.0f}s  {'; '.join(r['problems'])}")
    allok = all(r["ok"] for r in results)
    lines.append("")
    lines.append(f"P1 BYTE GATE {'PASS' if allok else 'FAIL'}: {sum(r['ok'] for r in results)}/{len(results)} identical")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(root, "summary.txt"), "w").write(txt + "\n")
    json.dump(results, open(os.path.join(root, "summary.json"), "w"), indent=1)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
