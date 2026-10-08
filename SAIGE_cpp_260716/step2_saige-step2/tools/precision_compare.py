#!/usr/bin/env python3
"""precision_compare.py A B [--json FILE] [--top K]

Compare two step-2 single-variant text outputs, A (reference, e.g. all-fp64)
and B (e.g. one stage in fp32 / int8), field by field. The accuracy gate for
the per-stage GPU precision modes (gpu/gpu_precision.hpp).

A and B are directories holding one <trait>.txt per trait (the outDir of a
`saige-gpu-cpp step2` run, or a directory with an out/ subdirectory holding
them), or two single .txt files. Traits are matched by file name, rows by
(CHR, POS, MarkerID, Allele1, Allele2).

Reported per trait and overall:
  rows only in A / only in B (and traits only in one of them)
  for p.value, p.value.NA, BETA, SE, Tstat, var:
      max |d|, max relative |d| = |a-b| / max(|a|,|b|) (0 when both are 0),
      rows that differ at all, rows NA on one side only;
      for the two p columns also max |d log10 p|, read from the printed
      mantissa and exponent so that underflow strings ("1.9E-350") count
  p.value crossing 5e-8 and 1e-5: rows significant on one side only
  Is.SPA: rows where it differs
Exit status 0 always (it reports; the caller decides what passes).
"""
import argparse
import json
import math
import os
import sys

NUM_FIELDS = ["p.value", "p.value.NA", "BETA", "SE", "Tstat", "var"]
P_FIELDS = ("p.value", "p.value.NA")
THRESH = (5e-8, 1e-5)
KEY = ("CHR", "POS", "MarkerID", "Allele1", "Allele2")


def trait_files(path):
    if os.path.isfile(path):
        return {os.path.basename(path)[:-4] if path.endswith(".txt") else os.path.basename(path): path}
    for d in (path, os.path.join(path, "out")):
        if not os.path.isdir(d):
            continue
        out = {}
        for f in sorted(os.listdir(d)):
            if not f.endswith(".txt"):
                continue
            fp = os.path.join(d, f)
            with open(fp) as h:
                first = h.readline()
            if first.startswith("CHR\t"):
                out[f[:-4]] = fp
        if out:
            return out
    sys.exit(f"precision_compare: no step-2 result files (<trait>.txt with a CHR header) in {path}")


def to_float(s):
    if s in ("NA", "nan", "NaN", "-nan", ""):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def log10p(s):
    """log10 of a printed p-value, exact for strings below the double range."""
    if s in ("NA", "nan", "NaN", "-nan", ""):
        return None
    t = s.upper()
    if "E" in t:
        m, e = t.split("E", 1)
        try:
            mv, ev = float(m), int(e)
        except ValueError:
            return None
        if mv <= 0:
            return -math.inf if mv == 0 else None
        return math.log10(mv) + ev
    try:
        v = float(t)
    except ValueError:
        return None
    if v <= 0:
        return -math.inf if v == 0 else None
    return math.log10(v)


def read(fp):
    rows = {}
    with open(fp) as h:
        hdr = h.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(hdr)}
        miss = [k for k in KEY if k not in col]
        if miss:
            sys.exit(f"precision_compare: {fp}: no column(s) {miss}")
        ki = [col[k] for k in KEY]
        fi = {f: col[f] for f in NUM_FIELDS + ["Is.SPA"] if f in col}
        for line in h:
            t = line.rstrip("\n").split("\t")
            rows[tuple(t[i] for i in ki)] = {f: t[i] for f, i in fi.items()}
    return rows


class Stat:
    def __init__(self):
        self.n = 0
        self.ndiff = 0
        self.nNA1 = 0          # NA on one side only
        self.maxabs = 0.0
        self.maxrel = 0.0
        self.maxdl = 0.0       # p columns: max |d log10 p|
        self.worst = None      # key of the max relative difference

    def add(self, key, sa, sb, isp):
        self.n += 1
        if sa == sb:
            return
        self.ndiff += 1
        a, b = to_float(sa), to_float(sb)
        if (a is None) != (b is None):
            self.nNA1 += 1
            return
        if a is not None:
            d = abs(a - b)
            m = max(abs(a), abs(b))
            r = d / m if m > 0 else 0.0
            if math.isfinite(d) and d > self.maxabs:
                self.maxabs = d
            if r > self.maxrel:
                self.maxrel = r
                self.worst = key
        if isp:
            la, lb = log10p(sa), log10p(sb)
            if la is not None and lb is not None and math.isfinite(la) and math.isfinite(lb):
                self.maxdl = max(self.maxdl, abs(la - lb))

    def merge(self, o):
        self.n += o.n; self.ndiff += o.ndiff; self.nNA1 += o.nNA1
        self.maxabs = max(self.maxabs, o.maxabs)
        if o.maxrel > self.maxrel:
            self.maxrel, self.worst = o.maxrel, o.worst
        self.maxdl = max(self.maxdl, o.maxdl)


def compare_trait(fa, fb):
    A, B = read(fa), read(fb)
    onlyA = [k for k in A if k not in B]
    onlyB = [k for k in B if k not in A]
    st = {f: Stat() for f in NUM_FIELDS}
    cross = {t: [0, 0] for t in THRESH}   # [significant in A only, in B only]
    spa = 0
    common = 0
    for k, ra in A.items():
        rb = B.get(k)
        if rb is None:
            continue
        common += 1
        for f in NUM_FIELDS:
            if f in ra and f in rb:
                st[f].add(k, ra[f], rb[f], f in P_FIELDS)
        if "Is.SPA" in ra and "Is.SPA" in rb and ra["Is.SPA"] != rb["Is.SPA"]:
            spa += 1
        la, lb = log10p(ra.get("p.value", "NA")), log10p(rb.get("p.value", "NA"))
        if la is not None and lb is not None:
            for t in THRESH:
                lt = math.log10(t)
                sa, sb = la < lt, lb < lt
                if sa and not sb:
                    cross[t][0] += 1
                elif sb and not sa:
                    cross[t][1] += 1
    return dict(rowsA=len(A), rowsB=len(B), common=common, onlyA=len(onlyA), onlyB=len(onlyB),
                onlyA_ex=[":".join(k) for k in onlyA[:3]], onlyB_ex=[":".join(k) for k in onlyB[:3]],
                stats=st, cross=cross, spa=spa)


def fmt(x):
    return "0" if x == 0 else f"{x:.3g}"


def print_block(name, r):
    print(f"== {name}: rows A {r['rowsA']}, B {r['rowsB']}, common {r['common']}, "
          f"only in A {r['onlyA']}, only in B {r['onlyB']}")
    if r["onlyA_ex"]:
        print(f"   only in A e.g. {', '.join(r['onlyA_ex'])}")
    if r["onlyB_ex"]:
        print(f"   only in B e.g. {', '.join(r['onlyB_ex'])}")
    print(f"   {'field':<11} {'max|d|':>10} {'max rel':>10} {'max|dlog10p|':>13} {'rows differ':>12} {'NA one side':>12}")
    for f in NUM_FIELDS:
        s = r["stats"][f]
        dl = fmt(s.maxdl) if f in P_FIELDS else "-"
        print(f"   {f:<11} {fmt(s.maxabs):>10} {fmt(s.maxrel):>10} {dl:>13} {s.ndiff:>12} {s.nNA1:>12}")
    cr = r["cross"]
    print("   p.value crossings: " + "; ".join(
        f"{t:g}: {cr[t][0]} significant in A only, {cr[t][1]} in B only" for t in THRESH))
    print(f"   Is.SPA differs: {r['spa']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("A")
    ap.add_argument("B")
    ap.add_argument("--json", help="also write the numbers as JSON to this file")
    ap.add_argument("--quiet", action="store_true", help="overall block only")
    a = ap.parse_args()
    TA, TB = trait_files(a.A), trait_files(a.B)
    if os.path.isfile(a.A) and os.path.isfile(a.B):
        TA, TB = {"trait": a.A}, {"trait": a.B}
    traits = [t for t in TA if t in TB]
    onlyTA = [t for t in TA if t not in TB]
    onlyTB = [t for t in TB if t not in TA]
    print(f"A: {a.A}\nB: {a.B}")
    print(f"traits: {len(traits)} in both" +
          (f"; only in A: {', '.join(onlyTA)}" if onlyTA else "") +
          (f"; only in B: {', '.join(onlyTB)}" if onlyTB else ""))
    tot = dict(rowsA=0, rowsB=0, common=0, onlyA=0, onlyB=0, onlyA_ex=[], onlyB_ex=[],
               stats={f: Stat() for f in NUM_FIELDS}, cross={t: [0, 0] for t in THRESH}, spa=0)
    per = {}
    for t in traits:
        r = compare_trait(TA[t], TB[t])
        per[t] = r
        if not a.quiet:
            print_block(t, r)
        for k in ("rowsA", "rowsB", "common", "onlyA", "onlyB", "spa"):
            tot[k] += r[k]
        for f in NUM_FIELDS:
            tot["stats"][f].merge(r["stats"][f])
        for th in THRESH:
            tot["cross"][th][0] += r["cross"][th][0]
            tot["cross"][th][1] += r["cross"][th][1]
        if len(tot["onlyA_ex"]) < 3:
            tot["onlyA_ex"] += [f"{t}:{x}" for x in r["onlyA_ex"]][:3 - len(tot["onlyA_ex"])]
        if len(tot["onlyB_ex"]) < 3:
            tot["onlyB_ex"] += [f"{t}:{x}" for x in r["onlyB_ex"]][:3 - len(tot["onlyB_ex"])]
    print_block(f"overall ({len(traits)} traits)", tot)
    if a.json:
        def js(r):
            return dict({k: r[k] for k in ("rowsA", "rowsB", "common", "onlyA", "onlyB", "spa")},
                        cross={f"{t:g}": r["cross"][t] for t in THRESH},
                        fields={f: dict(maxabs=s.maxabs, maxrel=s.maxrel, maxdlog10p=s.maxdl,
                                        ndiff=s.ndiff, naOneSide=s.nNA1,
                                        worst=":".join(s.worst) if s.worst else None)
                                for f, s in r["stats"].items()})
        with open(a.json, "w") as h:
            json.dump(dict(A=a.A, B=a.B, traitsOnlyA=onlyTA, traitsOnlyB=onlyTB,
                           overall=js(tot), traits={t: js(r) for t, r in per.items()}), h, indent=1)


if __name__ == "__main__":
    main()
