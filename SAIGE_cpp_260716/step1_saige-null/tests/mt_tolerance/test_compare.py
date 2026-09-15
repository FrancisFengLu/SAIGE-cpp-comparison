#!/usr/bin/env python3
"""Fault-injection tests for compare.py.

Takes a solo/multi output pair that is known to pass (found in the gate's run
tree, or given explicitly), copies the multi side into a scratch directory,
injects the symptom of a typical bug, and asserts the comparator's verdict
(and, for failures, that the failing item is the injected one).

  test_compare.py [--workdir /opt/saige/logs/mt_gate]
                  [--ref SOLO --test MULTI --trait T] [--loco-ref ... --loco-test ... --loco-trait ...]

Exit 0 iff every injection gives the expected verdict.
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402


def find_pair(workdir, case, trait):
    """A passing, byte-identical pair from an earlier gate run (grouped mode)."""
    best = None
    for p in sorted(glob.glob(os.path.join(workdir, "runs", "*", case, f"compare_{trait}.json"))):
        if "-lock" in p:
            continue
        r = json.load(open(p))
        if r.get("passed") and os.path.isdir(r["ref"]) and os.path.isdir(r["test"]):
            best = (r["ref"], r["test"])
    return best


def write_arma_like(path, arr, header_type):
    r, c = arr.shape
    dt = compare._ARMA_TYPES[header_type]
    with open(path, "wb") as f:
        f.write(f"ARMA_MAT_BIN_{header_type}\n{r} {c}\n".encode())
        f.write(np.asarray(arr, dtype=dt).T.copy().tobytes())


def model_dir(test, trait):
    return os.path.join(test, "m", trait)


def bump_sig_digit(txt, k):
    """Add one unit in the k-th significant digit of the decimal number `txt`."""
    v = float(txt)
    e = math.floor(math.log10(abs(v)))
    return repr(v + 10.0 ** (e - (k - 1)))


# ---- injections: each takes (scratch_test_dir, scratch_ref_dir, trait) ----------------------


def inj_none(t, r, tr):
    return "no change"


def inj_tau_1e4(t, r, tr):
    p = os.path.join(model_dir(t, tr), "nullmodel.json")
    s = open(p).read()
    m = re.search(r'"theta": \[([^\]]*)\]', s)
    vals = [float(x) for x in m.group(1).split(",")]
    k = int(np.argmax(np.abs(vals[1:]))) + 1 if len(vals) > 1 and max(abs(x) for x in vals[1:]) > 0 else 0
    vals[k] *= 1 + 1e-4
    s = s[:m.start(1)] + ",".join(f"{x:.10g}" for x in vals) + s[m.end(1):]
    open(p, "w").write(s)
    return f"theta[{k}] *= 1+1e-4"


def inj_tau_noise(t, r, tr):
    p = os.path.join(model_dir(t, tr), "nullmodel.json")
    s = open(p).read()
    m = re.search(r'"theta": \[([^\]]*)\]', s)
    vals = [float(x) for x in m.group(1).split(",")]
    vals = [x * (1 + 2e-6) for x in vals]
    s = s[:m.start(1)] + ",".join(f"{x:.10g}" for x in vals) + s[m.end(1):]
    open(p, "w").write(s)
    return "theta *= 1+2e-6 (noise level, must pass)"


def inj_mu_frac(t, r, tr):
    p = os.path.join(model_dir(t, tr), "mu.arma")
    ty, a = compare.read_arma(p)
    a = a.copy()
    rng = np.random.default_rng(1)
    idx = rng.choice(a.shape[0], max(1, a.shape[0] // 1000), replace=False)
    a[idx, 0] *= 1 + 1e-4
    write_arma_like(p, a, ty)
    return f"mu: {len(idx)} of {a.shape[0]} elements (0.1%) *= 1+1e-4"


def inj_mu_noise(t, r, tr):
    p = os.path.join(model_dir(t, tr), "mu.arma")
    ty, a = compare.read_arma(p)
    rng = np.random.default_rng(2)
    a = a * (1 + 2e-7 * rng.standard_normal(a.shape))
    write_arma_like(p, a, ty)
    return "mu *= 1 + 2e-7*N(0,1) on every element (fp32-noise level, must pass)"


def _vr_path(t, tr):
    return os.path.join(t, f"mvr_{tr}.varianceRatio.txt")


def _bump_vr(t, tr, k):
    p = _vr_path(t, tr)
    rows = [l.split("\t") for l in open(p).read().splitlines()]
    old = rows[0][0]
    rows[0][0] = bump_sig_digit(old, k)
    open(p, "w").write("\n".join("\t".join(x) for x in rows) + "\n")
    return old, rows[0][0]


def inj_vr_5th(t, r, tr):
    o, n = _bump_vr(t, tr, 5)
    return f"VR null bin1: {o} -> {n} (5th significant digit +1)"


def inj_vr_4th(t, r, tr):
    o, n = _bump_vr(t, tr, 4)
    return f"VR null bin1: {o} -> {n} (4th significant digit +1)"


def inj_vr_6th(t, r, tr):
    o, n = _bump_vr(t, tr, 6)
    return f"VR null bin1: {o} -> {n} (6th = last printed digit +1: inside print rounding, must pass)"


def inj_theta_5e5(t, r, tr):
    p = os.path.join(model_dir(t, tr), "nullmodel.json")
    s = open(p).read()
    m = re.search(r'"theta": \[([^\]]*)\]', s)
    vals = [float(x) for x in m.group(1).split(",")]
    k = 1 if len(vals) > 1 and vals[1] > 0 else 0
    vals[k] *= 1 + 6e-5
    s = s[:m.start(1)] + ",".join(f"{x:.10g}" for x in vals) + s[m.end(1):]
    open(p, "w").write(s)
    return f"theta[{k}] *= 1+6e-5 (just over the tolerance)"


def inj_delete_file(t, r, tr):
    os.remove(os.path.join(model_dir(t, tr), "XV.arma"))
    return "deleted XV.arma"


def inj_extra_file(t, r, tr):
    open(os.path.join(model_dir(t, tr), "stray.txt"), "w").write("x\n")
    return "added stray.txt to the trait's model dir"


def inj_iterations(t, r, tr):
    p = t + ".log"
    lines = open(p).read().split("\n")
    in_block = False
    for i, l in enumerate(lines):
        if l.startswith("== SAIGE Null Fit Completed =="):
            in_block = True
            ph = None
        elif in_block and l.startswith("Phenotype: "):
            ph = l.split(": ", 1)[1]
        elif in_block and l.startswith("Iterations: ") and ph == tr:
            lines[i] = f"Iterations: {int(l.split(': ')[1]) + 1}"
            in_block = False
    open(p, "w").write("\n".join(lines))
    return "Iterations +1 in the multi log"


def inj_converged(t, r, tr):
    p = t + ".log"
    s = open(p).read()
    s2 = re.sub(r"(Phenotype: %s\nConverged: )yes" % re.escape(tr), r"\1NO", s)
    assert s2 != s
    open(p, "w").write(s2)
    return "Converged yes -> NO"


def inj_sampleids_swap(t, r, tr):
    p = os.path.join(model_dir(t, tr), "nullmodel.json")
    j = json.load(open(p))
    ids = j["sampleIDs"]
    ids[0], ids[1] = ids[1], ids[0]
    s = open(p).read()
    m = re.search(r'"sampleIDs": \[(.*)\]', s, re.S)
    s = s[:m.start(1)] + ",".join(f'"{x}"' for x in ids) + s[m.end(1):]
    open(p, "w").write(s)
    return "sampleIDs[0] <-> [1]"


def inj_grm_diag(t, r, tr):
    p = os.path.join(t, "m", f"{tr}.grm_diag.txt")
    L = open(p).read().splitlines()
    L[7] = repr(float(L[7]) * (1 + 1e-4))
    open(p, "w").write("\n".join(L) + "\n")
    return "grm_diag line 8 *= 1+1e-4"


def inj_truncate_arma(t, r, tr):
    p = os.path.join(model_dir(t, tr), "res.arma")
    b = open(p, "rb").read()
    open(p, "wb").write(b[:-8])
    return "res.arma truncated by 8 bytes"


def inj_X_one_cell(t, r, tr):
    p = os.path.join(model_dir(t, tr), "X.arma")
    ty, a = compare.read_arma(p)
    a = a.copy()
    a[3, a.shape[1] - 1] += 1e-6
    write_arma_like(p, a, ty)
    return "X.arma one cell += 1e-6 (input must be exact)"


def inj_S_a_noise(t, r, tr):
    """S_a is a near-zero sum; a fp32-noise-sized change must pass, a real one must fail."""
    p = os.path.join(model_dir(t, tr), "S_a.arma")
    ty, a = compare.read_arma(p)
    a = a.copy()
    a[0, 0] += 1e-3
    write_arma_like(p, a, ty)
    return "S_a[0] += 1e-3 (fp32 summation noise for n>=10k, must pass)"


def inj_S_a_real(t, r, tr):
    p = os.path.join(model_dir(t, tr), "S_a.arma")
    _, X = compare.read_arma(os.path.join(model_dir(t, tr), "X.arma"))
    _, res = compare.read_arma(os.path.join(model_dir(t, tr), "res.arma"))
    ty, a = compare.read_arma(p)
    a = a.copy()
    mass = float(np.abs(X[:, 0] * res[:, 0]).sum())
    a[0, 0] += 1e-4 * mass
    write_arma_like(p, a, ty)
    return f"S_a[0] += 1e-4 * sum|X res| ({1e-4 * mass:.3g})"


def _add_pcol(path, shift_log10=0.0):
    L = open(path).read().splitlines()
    out = [L[0] + "\tp.value"]
    for i, l in enumerate(L[1:]):
        p = 10.0 ** (-(1 + (i % 7))) * 0.37
        out.append(l + "\t" + repr(p * 10.0 ** shift_log10))
    open(path, "w").write("\n".join(out) + "\n")


def _markers_file(d, tr=None):
    pat = f"mvr_{tr}.*markers.SAIGE.results.txt" if tr else "mvr.*markers.SAIGE.results.txt"
    g = glob.glob(os.path.join(d, pat))
    assert len(g) == 1, g
    return g[0]


def inj_pvalue_same(t, r, tr):
    _add_pcol(_markers_file(t, tr))
    _add_pcol(_markers_file(r))
    return "p.value column added to both sides, identical (control, must pass)"


def inj_pvalue_1e3(t, r, tr):
    _add_pcol(_markers_file(t, tr), shift_log10=1e-3)
    _add_pcol(_markers_file(r))
    return "p.value column added; test side log10 p shifted by 1e-3"


def inj_vr_markers_count(t, r, tr):
    f = _markers_file(t, tr)
    os.rename(f, f.replace(".30markers.", ".31markers."))
    return "VR markers file renamed 30markers -> 31markers (tested-marker count changed)"


def inj_loco_missing(t, r, tr):
    d = model_dir(t, tr)
    chrs = sorted(x for x in os.listdir(d) if x.startswith("chr"))
    os.remove(os.path.join(d, chrs[len(chrs) // 2], "mu.arma"))
    return f"deleted {chrs[len(chrs) // 2]}/mu.arma"


def inj_loco_mu(t, r, tr):
    d = model_dir(t, tr)
    p = os.path.join(d, "chr2", "mu.arma")
    ty, a = compare.read_arma(p)
    a = a.copy()
    a[:, 0] *= 1 + 1e-4
    write_arma_like(p, a, ty)
    return "chr2/mu.arma *= 1+1e-4"


# (name, injector, expect_pass, item expected among failures (regex on "file :: item"), which pair)
TESTS = [
    ("unchanged", inj_none, True, None, "main"),
    ("tau +1e-4 rel", inj_tau_1e4, False, r"nullmodel\.json :: theta", "main"),
    ("tau +2e-6 rel (noise)", inj_tau_noise, True, None, "main"),
    ("tau +6e-5 rel (just over tol)", inj_theta_5e5, False, r"nullmodel\.json :: theta", "main"),
    ("mu 0.1% elements +1e-4 rel", inj_mu_frac, False, r"mu\.arma :: mu", "main"),
    ("mu fp32-level noise", inj_mu_noise, True, None, "main"),
    ("VR 4th sig digit +1", inj_vr_4th, False, r"varianceRatio\.txt :: ratio", "main"),
    # 1e-5 relative on the VR is BELOW the vr.ratio tolerance: recorded here as the sensitivity limit
    ("VR 5th sig digit +1 (1e-5 rel)", inj_vr_5th, True, None, "main"),
    ("VR 6th sig digit +1 (rounding)", inj_vr_6th, True, None, "main"),
    ("delete a file", inj_delete_file, False, r"\(file set\) :: files", "main"),
    ("extra file", inj_extra_file, False, r"\(file set\) :: files", "main"),
    ("iterations +1", inj_iterations, False, r"\(log\) :: iterations", "main"),
    ("converged -> NO", inj_converged, False, r"\(log\) :: converged", "main"),
    ("sampleIDs swapped", inj_sampleids_swap, False, r"nullmodel\.json :: sampleIDs", "main"),
    ("grm_diag one value +1e-4", inj_grm_diag, False, r"grm_diag\.txt :: diag", "main"),
    ("truncated arma", inj_truncate_arma, False, r"res\.arma :: parse", "main"),
    ("X one cell +1e-6", inj_X_one_cell, False, r"X\.arma :: X", "main"),
    ("S_a +1e-3 (noise)", inj_S_a_noise, True, None, "main"),
    ("S_a +1e-4 of score mass", inj_S_a_real, False, r"S_a\.arma :: S_a", "main"),
    ("p.value column, identical", inj_pvalue_same, True, None, "main"),
    ("p.value log10 +1e-3", inj_pvalue_1e3, False, r"results\.txt :: p\.value", "main"),
    ("VR tested-marker count", inj_vr_markers_count, False, r"n_markers_tested", "main"),
    ("LOCO: unchanged", inj_none, True, None, "loco"),
    ("LOCO: chr dir missing a file", inj_loco_missing, False, r"\(file set\) :: files", "loco"),
    ("LOCO: chr2 mu +1e-4 rel", inj_loco_mu, False, r"chr2/mu\.arma :: mu", "loco"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--ref")
    ap.add_argument("--test")
    ap.add_argument("--trait")
    ap.add_argument("--loco-ref")
    ap.add_argument("--loco-test")
    ap.add_argument("--loco-trait")
    ap.add_argument("--thresholds", default=compare.DEFAULT_THRESHOLDS)
    a = ap.parse_args()
    pairs = {}
    if a.ref:
        pairs["main"] = (a.ref, a.test, a.trait)
    else:
        fp = find_pair(a.workdir, "nomiss", "b2")
        if not fp:
            raise SystemExit("no passing nomiss/b2 pair found under WORKDIR/runs; run run_gate.sh first or pass --ref/--test/--trait")
        pairs["main"] = (fp[0], fp[1], "b2")
    if a.loco_ref:
        pairs["loco"] = (a.loco_ref, a.loco_test, a.loco_trait)
    else:
        fp = find_pair(a.workdir, "loco", "lS")
        if not fp:
            raise SystemExit("no passing loco/lS pair found under WORKDIR/runs")
        pairs["loco"] = (fp[0], fp[1], "lS")
    scratch = os.path.join(a.workdir, "test_compare")
    if os.path.exists(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    for k, (r, t, tr) in pairs.items():
        print(f"pair {k}: ref {r}\n         test {t} trait {tr}")
    n_ok = 0
    for i, (name, fn, expect_pass, item_re, which) in enumerate(TESTS):
        ref, test, tr = pairs[which]
        d = os.path.join(scratch, f"t{i:02d}")
        tcopy = os.path.join(d, "test")
        rcopy = os.path.join(d, "ref")
        # copy only this trait's files (+ log) of the multi side; the solo side as is
        os.makedirs(os.path.join(tcopy, "m"))
        for rel in compare.walk(test):
            if compare.claimed_by(rel, tr):
                dst = os.path.join(tcopy, rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(os.path.join(test, rel), dst)
        shutil.copy2(test + ".log", tcopy + ".log")
        shutil.copytree(ref, rcopy)
        shutil.copy2(ref + ".log", rcopy + ".log")
        what = fn(tcopy, rcopy, tr)
        res = compare.compare_trait(rcopy, tcopy, tr, "solo", "multi", thresholds=a.thresholds)
        fails = [compare.describe(it) for it in res["items"] if not it["passed"]]
        ok = res["passed"] == expect_pass
        if ok and not expect_pass and item_re:
            ok = any(re.search(item_re, f) for f in fails)
        n_ok += ok
        verdict = "PASS" if res["passed"] else "FAIL"
        print(f"[{'ok ' if ok else 'BAD'}] {name:<34} expect {'PASS' if expect_pass else 'FAIL'}, got {verdict}  ({what})")
        for f in fails[:3]:
            print(f"        {f[:220]}")
        if not ok and res["passed"]:
            w = res["worst_numeric"]
            print(f"        worst numeric: {compare.describe(w) if w else '-'}")
        json.dump(res, open(os.path.join(d, "result.json"), "w"), indent=1, default=str)
    print(f"\ntest_compare: {n_ok}/{len(TESTS)} as expected")
    return 0 if n_ok == len(TESTS) else 1


if __name__ == "__main__":
    sys.exit(main())
