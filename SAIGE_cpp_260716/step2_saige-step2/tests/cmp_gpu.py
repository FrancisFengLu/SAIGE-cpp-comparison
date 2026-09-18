#!/usr/bin/env python3
"""Compare a CPU step-2 output directory against a GPU one (config key useGPU).

Two kinds of column, and they are held to different standards.

EXACT -- any difference is a bug, not a rounding artefact:
    the row set and its order          (the QC decisions)
    CHR POS MarkerID Allele1 Allele2   (the marker and its allele orientation)
    AC_Allele2 AF_Allele2              (allele counts, incl. the imputed cells)
    MissingRate imputationInfo N       (imputation bookkeeping)
Those are all computed on the host from the marker's 2-bit code counts, by the
same code on both paths, so they must agree bit for bit through the printed
form.

APPROXIMATE -- the reduction ran on a different machine in a different
association order:
    BETA SE Tstat var p.value
Reported as max / median of |delta(-log10 p)| and of |delta Tstat| / sd(Tstat).
sd(Tstat) is the empirical spread of the statistic over the compared markers,
which is the scale a p-value actually depends on: an absolute Tstat error is
only meaningful against it.

Exit status 1 if any exact column differs, or if either approximate metric
exceeds --tol-logp / --tol-sd.
"""
import argparse
import math
import os
import sys

EXACT = ["CHR", "POS", "MarkerID", "Allele1", "Allele2",
         "AC_Allele2", "AF_Allele2", "MissingRate", "imputationInfo", "N"]
APPROX = ["BETA", "SE", "Tstat", "var", "p.value"]


def read(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f]
    return head, rows


def logp(s):
    """-log10(p) from the printed field, which is '%.6E' or the '%.1fE%d'
    underflow form that format_score_result falls back to below ~1e-300."""
    try:
        v = float(s)
    except ValueError:
        return float("nan")
    if v > 0:
        return -math.log10(v)
    return float("inf")


def quant(v, q):
    if not v:
        return float("nan")
    v = sorted(v)
    i = min(len(v) - 1, max(0, int(q * (len(v) - 1))))
    return v[i]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", required=True)
    ap.add_argument("--gpu", required=True)
    ap.add_argument("--P", type=int, required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--tol-logp", type=float, default=1e-3,
                    help="max allowed |delta(-log10 p)|")
    ap.add_argument("--tol-sd", type=float, default=1e-4,
                    help="max allowed |delta Tstat| / sd(Tstat)")
    a = ap.parse_args()

    bad = 0
    allLogp, allSd = [], []
    nrows = 0
    for t in range(1, a.P + 1):
        fc = os.path.join(a.cpu, f"y{t}.txt")
        fg = os.path.join(a.gpu, f"y{t}.txt")
        if not (os.path.exists(fc) and os.path.exists(fg)):
            print(f"  MISSING output for y{t}")
            bad += 1
            continue
        hc, rc = read(fc)
        hg, rg = read(fg)
        if hc != hg:
            print(f"  y{t}: HEADER differs")
            bad += 1
            continue
        if len(rc) != len(rg):
            print(f"  y{t}: row count {len(rc)} vs {len(rg)} -- the two paths "
                  f"disagree about which markers pass QC")
            bad += 1
            continue

        idx = {n: i for i, n in enumerate(hc)}
        exact_bad = {}
        dlp, dsd = [], []
        tstats = []
        for x, y in zip(rc, rg):
            for n in EXACT:
                i = idx.get(n)
                if i is not None and x[i] != y[i]:
                    exact_bad.setdefault(n, [0, (x[idx["MarkerID"]], x[i], y[i])])
                    exact_bad[n][0] += 1
            i = idx.get("p.value")
            if i is not None:
                d = abs(logp(x[i]) - logp(y[i]))
                if not math.isnan(d):
                    dlp.append(d)
            i = idx.get("Tstat")
            if i is not None:
                try:
                    cx, cy = float(x[i]), float(y[i])
                    tstats.append(cx)
                    dsd.append(abs(cx - cy))
                except ValueError:
                    pass
        nrows = len(rc)
        if exact_bad:
            bad += 1
            for n, (cnt, ex) in exact_bad.items():
                print(f"  y{t}: EXACT column {n} differs in {cnt} rows, e.g. "
                      f"{ex[0]}: {ex[1]!r} vs {ex[2]!r}")
        if tstats:
            m = sum(tstats) / len(tstats)
            sd = math.sqrt(sum((v - m) ** 2 for v in tstats) / max(1, len(tstats) - 1))
            if sd > 0:
                dsd = [d / sd for d in dsd]
            else:
                dsd = []
        allLogp += dlp
        allSd += dsd

    tag = f"[{a.tag}] " if a.tag else ""
    if allLogp:
        mx, md = max(allLogp), quant(allLogp, 0.5)
        print(f"  {tag}{nrows} markers x {a.P} traits = {len(allLogp)} pairs")
        print(f"  {tag}|d(-log10 p)|   max {mx:.3e}  median {md:.3e}  p99 {quant(allLogp,0.99):.3e}")
        if mx > a.tol_logp:
            print(f"  {tag}FAIL: max |d(-log10 p)| {mx:.3e} > tol {a.tol_logp:.3e}")
            bad += 1
    if allSd:
        mx, md = max(allSd), quant(allSd, 0.5)
        print(f"  {tag}|dS|/sd(S)      max {mx:.3e}  median {md:.3e}  p99 {quant(allSd,0.99):.3e}")
        if mx > a.tol_sd:
            print(f"  {tag}FAIL: max |dS|/sd(S) {mx:.3e} > tol {a.tol_sd:.3e}")
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
