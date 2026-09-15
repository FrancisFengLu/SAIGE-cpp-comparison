#!/usr/bin/env python3
"""Multi-trait tolerance gate driver (run_gate.sh is the entry point).

For every case in cases/: run each trait alone (P=1, cached by binary md5 +
config hash) and all traits together (P>1), verify the case really exercised
its scenario (checks), compare every trait's multi output with its solo output
(compare.py), and print a case x trait table.

  gate.py BIN [--gpu|--cpu] [--small] [--cases a,b] [--workdir DIR]
          [--lockstep] [--nthreads N] [--solo-vs-solo] [--label NAME]

Exit 0 iff every check and every comparison passed.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402

CASES_DIR = os.path.join(HERE, "cases")
ALL_CASES = ["nomiss", "block", "indep", "tiny", "qc", "fill", "weak", "quant", "loco",
             "sparse", "sparse_pcg", "vrcate"]
CACHE_VERSION = 1

FIT_DEFAULTS = dict(loco=False, maxiter=20, tol=0.02, tolPCG=1e-5, maxiterPCG=500, nrun=30,
                    trace_seed=200, num_markers_for_vr=30, min_maf_grm=0.01,
                    use_sparse_grm_to_fit=False, use_pcg_with_sparse_grm=False)


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def log(msg):
    print(f"[gate {datetime.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def load_case(name):
    with open(os.path.join(CASES_DIR, f"{name}.yaml")) as f:
        c = yaml.safe_load(f)
    c["name"] = name
    c.setdefault("fit", {})
    c.setdefault("checks", [])
    c.setdefault("sparse", "none")
    c.setdefault("plink", "base")
    c.setdefault("data_case", name)
    return c


def load_manifest(workdir, scale, data_case):
    path = os.path.join(workdir, "data", scale, f"{data_case}.manifest.json")
    if not os.path.exists(path):
        log(f"manifest missing for {scale}/{data_case}: running make_cases.py")
        subprocess.run([sys.executable, os.path.join(HERE, "make_cases.py"), "--workdir", workdir,
                        "--scale", scale, "--cases", data_case], check=True)
    m = json.load(open(path))
    if os.path.exists(m["pheno"]) and md5_file(m["pheno"]) != m["pheno_md5"]:
        raise SystemExit(f"{m['pheno']} changed since its manifest was written; rerun make_cases.py --force")
    return m


def render_config(case, man, scale_plink, nthreads, traits, rundir, lockstep):
    fit = dict(FIT_DEFAULTS)
    fit["trait"] = case["trait"]
    fit["nthreads"] = nthreads
    fit.update(case["fit"])
    fit.update((case.get("fit_scale") or {}).get(man["scale"], {}))
    fit["multi_lockstep"] = bool(lockstep)
    plink = man["loco_plink"] if case["plink"] == "loco" else scale_plink
    paths = {"plinkFile": plink, "out_prefix": f"{rundir}/m", "out_prefix_vr": f"{rundir}/mvr",
             "overwrite_varratio": True}
    if case["sparse"] != "none":
        paths["sparse_grm"] = man["sparse_grm"]
        paths["sparse_grm_ids"] = man["sparse_grm_ids"]
    design = {"csv": man["pheno"], "iid_col": "IID", "covar_cols": ["x1", "x2"]}
    if len(traits) == 1:
        design["y_col"] = traits[0]
    else:
        design["y_cols"] = list(traits)
    return {"paths": paths, "design": design, "fit": fit}


def clean_env():
    env = {k: v for k, v in os.environ.items() if not k.startswith("SAIGE_")}
    return env


def run_saige(binary, cfg, rundir, gpu, timeout=6 * 3600):
    """Fresh rundir; config at rundir.yaml, log at rundir.log. -> (rc, seconds)."""
    if os.path.exists(rundir):
        shutil.rmtree(rundir)
    os.makedirs(rundir)
    cfg_path = rundir + ".yaml"
    with open(cfg_path, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    cmd = [binary, "-c", cfg_path] + (["--gpu"] if gpu else [])
    t0 = time.time()
    with open(rundir + ".log", "w") as lf:
        try:
            rc = subprocess.run(cmd, cwd=rundir, stdout=lf, stderr=subprocess.STDOUT, env=clean_env(),
                                timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = -999
    return rc, time.time() - t0


def solo_cache_key(bin_md5, cfg, gpu, man):
    c = json.loads(json.dumps(cfg))
    c["paths"]["out_prefix"] = "@OUT@/m"
    c["paths"]["out_prefix_vr"] = "@OUT@/mvr"
    c["fit"].pop("multi_lockstep", None)      # irrelevant at P=1 (the binary says so and takes the per-trait path)
    blob = json.dumps(dict(v=CACHE_VERSION, bin=bin_md5, gpu=gpu, cfg=c, pheno_md5=man["pheno_md5"],
                           inputs=man["inputs"]), sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:20]


# ----------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------


def read_text(p):
    try:
        return open(p, errors="replace").read()
    except OSError:
        return ""


def check_results(case, man, solos, multi, mode):
    """solos: {trait: rundir}; multi: rundir or None. -> list of (name, ok|None, detail)."""
    out = []
    traits = case["traits"]
    slog = {t: read_text(solos[t] + ".log") for t in traits}
    mlog = read_text(multi + ".log") if multi else ""
    nm = {}
    for t in traits:
        try:
            nm[t] = json.load(open(os.path.join(solos[t], "m", "nullmodel.json")))
        except Exception:  # noqa: BLE001
            nm[t] = None

    def add(name, ok, detail):
        out.append((name, ok, detail))

    add("coverage(make_cases)", bool(man.get("coverage_ok")), "; ".join(man.get("coverage_problems", [])) or "met")
    bad = [t for t in traits if nm[t] is None]
    add("solo outputs exist", not bad, f"missing nullmodel.json: {bad}" if bad else f"{len(traits)} traits")
    bad = [t for t in traits if "[multi-pheno] P=" in slog[t]]
    add("solo runs are P=1", not bad, f"multi-pheno banner in solo log: {bad}" if bad else "ok")
    if multi:
        m = re.search(r"\[multi-pheno\] P=(\d+)", mlog)
        add("multi run is P>1", bool(m) and int(m.group(1)) == len(traits),
            f"banner P={m.group(1) if m else None}, expected {len(traits)}")
    ns = {t: (nm[t] or {}).get("n") for t in traits}
    bad = {t: (ns[t], man["traits"][t]["n"]) for t in traits if ns[t] != man["traits"][t]["n"]
           or len((nm[t] or {}).get("sampleIDs", [])) != man["traits"][t]["n"]}
    add("n per trait == manifest", not bad, f"mismatch (got, want): {bad}" if bad else str(ns))
    sets = {tuple((nm[t] or {}).get("sampleIDs", [])) for t in traits}
    want_sets = len({man["traits"][t]["set"] for t in traits})
    add("distinct sample sets == manifest", len(sets) == want_sets, f"{len(sets)} (manifest {want_sets})")
    pred = man.get("predicted", {})
    bad = {}
    for t in traits:
        got = re.findall(r"^(\d+) markers with MAF >= ", slog[t], re.M)
        want = pred.get("grm_markers", {}).get(t)
        if want is None or len(got) != 1 or int(got[0]) != want:
            bad[t] = (got, want)
    add("GRM marker count == genoqc prediction", not bad,
        f"(log, predicted): {bad}" if bad else str({t: pred['grm_markers'][t] for t in traits}))
    bad = {}
    for t in traits:
        got = [[int(a), int(b), int(c)] for a, b, c in
               re.findall(r"^Sample \d: 0=(\d+), 1=(\d+), 2=(\d+)$", slog[t], re.M)]
        want = pred.get("first5_counts", {}).get(t)
        if want is None or got != want:
            bad[t] = dict(log=got, predicted=want)
    add("first-5-sample filled genotype counts == prediction", not bad,
        json.dumps(bad) if bad else "all traits")

    dense = case["sparse"] != "fit"
    if mode["gpu"]:
        if dense:
            bad = [t for t in traits if "[gpu_matvec] tier=4" not in slog[t]]
            if multi and "[gpu_matvec] tier=4" not in mlog:
                bad.append("(multi)")
            add("GPU tier=4 used", not bad, f"no tier=4 line: {bad}" if bad else "solo + multi")
        else:
            add("GPU tier=4 used", None, "n/a: sparse fit does not use the dense GPU kernel")
    else:
        bad = [t for t in traits if "[gpu_matvec] tier=" in slog[t]] + (["(multi)"] if "[gpu_matvec] tier=" in mlog else [])
        add("CPU run has no GPU handle", not bad, f"gpu_matvec line in: {bad}" if bad else "ok")

    for ch in case["checks"]:
        name, arg = (ch, {}) if isinstance(ch, str) else next(iter(ch.items()))
        if name == "lockstep_group":
            shared = [ts for ts in man["sets"].values() if len([x for x in ts if x in traits]) >= 2]
            if not mode["lockstep"]:
                continue
            if not multi:
                continue
            n_lock = len(re.findall(r"^### lockstep multi-phenotype fit: P=", mlog, re.M))
            if shared:
                add("lockstep exercised", n_lock >= 1, f"{n_lock} lockstep fit(s); {len(shared)} set(s) with >=2 traits")
            else:
                add("lockstep exercised", None, "n/a: no sample set holds >=2 traits")
        elif name == "loco":
            want = " ".join(str(c) for c in man["loco_chroms"])
            bad = []
            for t in traits:
                bl = compare.parse_logs(solos[t] + ".log")
                if len(bl) != 1 or not bl[0].get("loco_line", "").endswith(f"chroms: {want}") \
                        or "LOCO: on" not in bl[0].get("loco_line", ""):
                    bad.append(f"solo {t}")
                dirs = sorted(d for d in os.listdir(os.path.join(solos[t], "m")) if d.startswith("chr")) \
                    if os.path.isdir(os.path.join(solos[t], "m")) else []
                if dirs != sorted(f"chr{c}" for c in man["loco_chroms"]):
                    bad.append(f"solo {t} dirs {dirs}")
                if multi:
                    mb = [b for b in compare.parse_logs(multi + ".log") if b.get("phenotype") == t]
                    if len(mb) != 1 or "LOCO: on" not in mb[0].get("loco_line", "") \
                            or not mb[0]["loco_line"].endswith(f"chroms: {want}"):
                        bad.append(f"multi {t}")
            add("LOCO on, chr dirs written", not bad, f"failing: {bad}" if bad else f"chroms {want}")
        elif name == "sparse":
            bad = []
            flag = "use_pcg_with_sparse_grm=" + ("true" if arg.get("pcg") else "false")
            for t in traits:
                d = os.path.join(solos[t], "m")
                if not (nm[t] or {}).get("flagSparseGRM"):
                    bad.append(f"{t}: flagSparseGRM false")
                for f_ in ("sparseGRM_locationMat.arma", "sparseGRM_valueVec.arma"):
                    if not os.path.exists(os.path.join(d, f_)):
                        bad.append(f"{t}: no {f_}")
                if "[sparse] Subsetted GRM" not in slog[t] or flag not in slog[t]:
                    bad.append(f"{t}: log lacks subset/{flag}")
            if multi and ("[sparse] Subsetted GRM" not in mlog or flag not in mlog):
                bad.append("multi log lacks sparse lines")
            add(f"sparse fit ({flag})", not bad, "; ".join(bad) if bad else "ok")
        elif name == "nondegenerate":
            th = {t: (nm[t] or {}).get("theta") for t in traits}
            zero = [t for t in traits if not th[t] or len(th[t]) < 2 or th[t][1] <= 0]
            add("tau1 > 0 in every solo fit (GRM enters the fit)", not zero,
                f"tau1 == 0: {zero}; theta: {th}" if zero else str({t: th[t][1] for t in traits}))
        elif name == "theta_zero":
            z = [t for t in traits if nm[t] and len(nm[t].get("theta", [])) > 1 and nm[t]["theta"][1] == 0]
            add("tau1 == 0 boundary reached", len(z) >= arg.get("min_traits", 1),
                f"solo traits with tau1==0: {z}; theta: " + str({t: (nm[t] or {}).get('theta') for t in traits}))
        elif name == "strong_signal":
            r = {t: (nm[t]["theta"][1] / (nm[t]["theta"][0] + nm[t]["theta"][1])) if nm[t] else None for t in traits}
            add("strong genetic signal", all(v is not None and v >= arg.get("min_ratio", 0.3) for v in r.values()),
                str({t: (round(v, 3) if v is not None else None) for t, v in r.items()}))
        elif name == "cate_vr":
            bad, counts = [], {}
            for t in traits:
                p = os.path.join(solos[t], "mvr.varianceRatio.txt")
                rows = [l.split() for l in read_text(p).splitlines() if l.strip()]
                kinds = sorted({(r[1], r[2]) for r in rows}) if rows else []
                want = sorted({(k, b) for k in ("null", "null_noXadj", "sparse") for b in ("1", "2")})
                if kinds != want:
                    bad.append(f"{t}: rows {kinds}")
                counts[t] = {int(b): int(n) for b, n in
                             re.findall(r"^\[VR\] Bin (\d+): null=\S+ null_noXadj=\S+ sparse=\S+ \(n=(\d+)\)", slog[t], re.M)}
            add("categorical VR: 2 bins x {null,null_noXadj,sparse}", not bad, "; ".join(bad) if bad else "ok")
            empty = [t for t in traits if counts[t].get(1, 0) < 1 or counts[t].get(2, 0) < 1]
            add("categorical VR: both bins estimated from markers", not empty,
                f"markers per bin (solo logs): {counts}")
        elif name in ("qc_design", "fill_design"):
            d = man.get("design", {})
            if name == "qc_design":
                s = {k: d.get(k) for k in ("designed_i", "designed_ii", "designed_iii", "designed_iv")}
                add("QC-boundary markers present", all((v or 0) >= 1 for v in s.values()),
                    f"{s}, union n={d.get('union_n')}")
            else:
                add("fill-boundary markers present", d.get("flip_markers", 0) >= 1 and d.get("first5_sensitive", 0) >= 1,
                    f"flip markers {d.get('flip_markers')}, cells {d.get('flip_cells')}, "
                    f"fill-sensitive first-5 samples {d.get('first5_sensitive')}")
        else:
            add(f"unknown check {name}", False, "fix cases/*.yaml")
    return out


# ----------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("binary")
    dev = ap.add_mutually_exclusive_group()
    dev.add_argument("--gpu", action="store_true")
    dev.add_argument("--cpu", action="store_true")
    ap.add_argument("--small", action="store_true", help="small-scale data (N=10,000 x M=20,000)")
    ap.add_argument("--cases", default=",".join(ALL_CASES))
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--lockstep", action="store_true", help="multi run with fit.multi_lockstep: true")
    ap.add_argument("--nthreads", type=int, default=8)
    ap.add_argument("--solo-vs-solo", action="store_true",
                    help="also rerun every solo uncached and compare it with the cached one (run-to-run noise)")
    ap.add_argument("--no-cache", action="store_true", help="ignore cached solo references")
    ap.add_argument("--reuse-runs", action="store_true",
                    help="do not run anything new for the multi side if this label already has a finished "
                         "multi run; recompute checks and comparisons (e.g. after editing thresholds.yaml)")
    ap.add_argument("--label", help="run label (default derived from binary/mode)")
    ap.add_argument("--thresholds", default=compare.DEFAULT_THRESHOLDS)
    a = ap.parse_args(argv)

    binary = os.path.abspath(a.binary)
    if not os.access(binary, os.X_OK):
        raise SystemExit(f"not executable: {binary}")
    gpu = bool(a.gpu)
    scale = "small" if a.small else "mid"
    bin_md5 = md5_file(binary)
    cases = [c for c in a.cases.split(",") if c]
    for c in cases:
        if not os.path.exists(os.path.join(CASES_DIR, f"{c}.yaml")):
            raise SystemExit(f"unknown case {c}")
    mode = dict(gpu=gpu, lockstep=a.lockstep)
    label = a.label or f"{scale}-{'gpu' if gpu else 'cpu'}{'-lock' if a.lockstep else ''}-t{a.nthreads}-{bin_md5[:8]}"
    rundir_root = os.path.join(a.workdir, "runs", label)
    os.makedirs(rundir_root, exist_ok=True)
    scale_plink = f"/opt/saige/data/{scale}"
    log(f"binary {binary} md5 {bin_md5}; {'GPU' if gpu else 'CPU'}; scale {scale}; nthreads {a.nthreads}; "
        f"lockstep {a.lockstep}; label {label}")

    rows, case_rows, all_ok = [], [], True
    for cname in cases:
        case = load_case(cname)
        man = load_manifest(a.workdir, scale, case["data_case"])
        traits = case["traits"]
        cdir = os.path.join(rundir_root, cname)
        os.makedirs(cdir, exist_ok=True)
        solos, solo_secs, solo_cached = {}, 0.0, 0
        fail_runs = []
        for t in traits:
            cfg = render_config(case, man, scale_plink, a.nthreads, [t], "@OUT@", False)
            key = solo_cache_key(bin_md5, cfg, gpu, man)
            cache = os.path.join(a.workdir, "ref_cache", scale, cname, f"{t}-{key}")
            rd = os.path.join(cache, "run")
            meta_p = os.path.join(cache, "meta.json")
            if not a.no_cache and os.path.exists(os.path.join(cache, "DONE")):
                meta = json.load(open(meta_p))
                solo_cached += 1
                solo_secs += meta["seconds"]
                log(f"{cname}/{t}: solo cached ({meta['seconds']:.0f}s originally)")
            else:
                if os.path.exists(cache):
                    shutil.rmtree(cache)
                os.makedirs(cache)
                cfg = render_config(case, man, scale_plink, a.nthreads, [t], rd, False)
                log(f"{cname}/{t}: solo run ...")
                rc, secs = run_saige(binary, cfg, rd, gpu)
                solo_secs += secs
                json.dump(dict(rc=rc, seconds=secs, binary=binary, bin_md5=bin_md5, gpu=gpu, key=key,
                               when=datetime.datetime.now().isoformat()), open(meta_p, "w"), indent=1)
                log(f"{cname}/{t}: solo rc={rc} {secs:.0f}s")
                if rc == 0:
                    open(os.path.join(cache, "DONE"), "w").write("ok\n")
                else:
                    fail_runs.append(f"solo {t} rc={rc}")
            solos[t] = rd
        mdir = os.path.join(cdir, "multi")
        mmeta = os.path.join(cdir, "multi_meta.json")
        cfg = render_config(case, man, scale_plink, a.nthreads, traits, mdir, a.lockstep)
        if a.reuse_runs and os.path.exists(mmeta):
            mm = json.load(open(mmeta))
            rc, msecs = mm["rc"], mm["seconds"]
            log(f"{cname}: reusing multi run (rc={rc}, {msecs:.0f}s)")
        else:
            log(f"{cname}: multi run P={len(traits)} ...")
            rc, msecs = run_saige(binary, cfg, mdir, gpu)
            json.dump(dict(rc=rc, seconds=msecs, bin_md5=bin_md5), open(mmeta, "w"))
            log(f"{cname}: multi rc={rc} {msecs:.0f}s")
        if rc != 0:
            fail_runs.append(f"multi rc={rc}")
        checks = check_results(case, man, solos, mdir, mode)
        orph = compare.orphans(mdir, traits) if os.path.isdir(mdir) else ["(no multi dir)"]
        checks.append(("multi files all claimed by a trait", not orph, f"unclaimed: {orph[:10]}" if orph else "ok"))
        if fail_runs:
            checks.append(("runs exit 0", False, "; ".join(fail_runs)))
        case_ok = all(ok is not False for _, ok, _ in checks)
        for t in traits:
            res = compare.compare_trait(solos[t], mdir, t, "solo", "multi", thresholds=a.thresholds)
            with open(os.path.join(cdir, f"compare_{t}.json"), "w") as f:
                json.dump(res, f, indent=1, default=str)
            with open(os.path.join(cdir, f"compare_{t}.txt"), "w") as f:
                f.write(compare.summary_text(res, verbose=True) + "\n")
            w = res["worst_numeric"]
            fails = [it for it in res["items"] if not it["passed"]]
            # Documented discrete flips (cases/<case>.yaml: known_flips: {trait: {items: [...], reason: ...}}):
            # only the listed exact items may differ; they are reported as FLIP, never as a plain PASS.
            kf = (case.get("known_flips") or {}).get(t)
            flip_items = []
            if kf and fails:
                allowed = set(kf.get("items", []))
                if all(it["item"] in allowed and it["file"] == "(log)" for it in fails):
                    flip_items = [f"{it['item']} {it.get('ref_value')}->{it.get('test_value')}" for it in fails]
                    res["passed"] = True
                    res["known_flip"] = dict(items=flip_items, reason=kf.get("reason"))
                    fails = []
            rows.append(dict(case=cname, trait=t, passed=res["passed"], n_failed=len(fails), flip=flip_items,
                             worst=(f"{w['file']}::{w['item']} {w['dev']:.2g}/{w['tol']:.0g}" if w else "-"),
                             first_fail=(compare.describe(fails[0])[:150] if fails else ""),
                             iters=res["info"].get("iterations")))
            all_ok &= res["passed"]
            if a.solo_vs_solo:
                rd2 = os.path.join(cdir, f"solo2_{t}")
                cfg2 = render_config(case, man, scale_plink, a.nthreads, [t], rd2, False)
                if not (a.reuse_runs and os.path.exists(rd2 + ".log")):
                    log(f"{cname}/{t}: solo rerun for run-to-run noise ...")
                    run_saige(binary, cfg2, rd2, gpu)
                r2 = compare.compare_trait(solos[t], rd2, None, "solo", "solo", thresholds=a.thresholds)
                with open(os.path.join(cdir, f"solo_vs_solo_{t}.json"), "w") as f:
                    json.dump(r2, f, indent=1, default=str)
                w2 = r2["worst_numeric"]
                rows.append(dict(case=cname, trait=t + " (solo rerun)", passed=r2["passed"],
                                 n_failed=r2["n_failed"],
                                 worst=(f"{w2['file']}::{w2['item']} {w2['dev']:.2g}/{w2['tol']:.0g}" if w2 else "-"),
                                 first_fail=(compare.describe([i for i in r2["items"] if not i["passed"]][0])[:150]
                                             if r2["n_failed"] else ""), iters=r2["info"].get("iterations")))
                all_ok &= r2["passed"]
        all_ok &= case_ok
        case_rows.append(dict(case=cname, ok=case_ok, checks=checks, solo_secs=solo_secs, solo_cached=solo_cached,
                              multi_secs=msecs, n_traits=len(traits)))
        with open(os.path.join(cdir, "checks.json"), "w") as f:
            json.dump(checks, f, indent=1)
        for name, ok, detail in checks:
            log(f"{cname}: check {'PASS' if ok else ('n/a ' if ok is None else 'FAIL')} {name}: {detail[:300]}")

    lines = [f"== multi-trait tolerance gate: {label}", f"binary {binary} (md5 {bin_md5})",
             f"device {'GPU' if gpu else 'CPU'}, scale {scale}, nthreads {a.nthreads}, lockstep {a.lockstep}",
             f"thresholds {os.path.abspath(a.thresholds)}", ""]
    lines.append(f"{'case':<11} {'trait':<22} {'result':<6} {'iters':>5}  worst numeric item (dev/tol)    first failure")
    for r in rows:
        verdict = ("FLIP*" if r.get("flip") else "PASS") if r["passed"] else "FAIL"
        lines.append(f"{r['case']:<11} {r['trait']:<22} {verdict:<6} {str(r['iters']):>5}  "
                     f"{r['worst']:<34} {r['first_fail'] or ('known flip: ' + '; '.join(r['flip']) if r.get('flip') else '')}")
    lines.append("")
    lines.append(f"{'case':<11} {'checks':<7} {'solo s':>8} {'(cached)':>8} {'multi s':>8}  failing/na checks")
    for c in case_rows:
        bad = [f"{n}: {d[:120]}" for n, ok, d in c["checks"] if ok is False]
        na = [n for n, ok, _ in c["checks"] if ok is None]
        lines.append(f"{c['case']:<11} {'PASS' if c['ok'] else 'FAIL':<7} {c['solo_secs']:>8.0f} "
                     f"{str(c['solo_cached']) + '/' + str(c['n_traits']):>8} {c['multi_secs']:>8.0f}  "
                     + ("; ".join(bad) if bad else "-") + (f"  [n/a: {', '.join(na)}]" if na else ""))
    lines.append("")
    lines.append(f"GATE {'PASS' if all_ok else 'FAIL'}")
    txt = "\n".join(lines)
    print(txt, flush=True)
    with open(os.path.join(rundir_root, "summary.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(rundir_root, "summary.json"), "w") as f:
        json.dump(dict(label=label, binary=binary, bin_md5=bin_md5, gpu=gpu, scale=scale, nthreads=a.nthreads,
                       lockstep=a.lockstep, passed=all_ok, rows=rows, cases=case_rows), f, indent=1, default=str)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
