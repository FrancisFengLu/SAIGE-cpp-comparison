#!/usr/bin/env python3
"""firth_precision_report.py SETDIR [SETDIR ...]

Firth-specific view of a run_precision_gate.sh set (SETDIR = WORKDIR/<set>,
holding fp64/ and test/ runs with out/<trait>.txt and log.txt), for
gpuPrecisionFirth. Prints per set:
  rows whose BETA / SE differ, max relative |d| of BETA and SE (all rows, and
  without MAC = 1 rows, the ones the legacy fit can leave oscillating at maxit),
  BETA sign flips, rows whose p.value differs (Firth must not change p) and
  p.value crossings of 5e-8 / 1e-5;
  from the logs: Firth fits applied / "successfully converged" (maxit counts),
  the SAIGE_FIRTH_STATS line (strict / maxit / singular / mean iterations) when
  the runs had it set, and the device Firth us/pair line.
"""
import os
import re
import sys


def rows(d):
    out = {}
    od = os.path.join(d, "out")
    for f in sorted(os.listdir(od)):
        if not f.endswith(".txt"):
            continue
        with open(os.path.join(od, f)) as h:
            hdr = h.readline().rstrip("\n").split("\t")
            ix = {k: hdr.index(k) for k in ("MarkerID", "BETA", "SE", "p.value", "AC_Allele2", "CHR", "POS", "Allele1", "Allele2")}
            for line in h:
                a = line.rstrip("\n").split("\t")
                key = (f, a[ix["CHR"]], a[ix["POS"]], a[ix["MarkerID"]], a[ix["Allele1"]], a[ix["Allele2"]])
                out[key] = (a[ix["BETA"]], a[ix["SE"]], a[ix["p.value"]], a[ix["AC_Allele2"]])
    return out


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def rel(a, b):
    m = max(abs(a), abs(b))
    return abs(a - b) / m if m > 0 else 0.0


def logstats(path):
    txt = open(path).read()
    app = sum(int(m) for m in re.findall(r"Firth approx was applied to (\d+) markers", txt))
    conv = sum(int(m) for m in re.findall(r"(\d+) successfully converged", txt))
    st = re.findall(r"^firth stats.*$", txt, re.M)
    dev = re.findall(r"device Firth: .*?us/pair\)", txt)
    return app, conv, st, dev


def mac(s):
    v = num(s)
    return v if v is not None else 0.0


for sd in sys.argv[1:]:
    A, B = rows(os.path.join(sd, "fp64")), rows(os.path.join(sd, "test"))
    common = A.keys() & B.keys()
    nb = ns = npd = flips = 0
    mb = ms = mb2 = ms2 = 0.0
    cross = {5e-8: [0, 0], 1e-5: [0, 0]}
    for k in common:
        a, b = A[k], B[k]
        if a[0] != b[0]:
            nb += 1
        if a[1] != b[1]:
            ns += 1
        if a[2] != b[2]:
            npd += 1
        ba, bb, sa, sb = num(a[0]), num(b[0]), num(a[1]), num(b[1])
        if ba is not None and bb is not None:
            r = rel(ba, bb)
            mb = max(mb, r)
            if not (0.5 < mac(a[3]) < 1.5):
                mb2 = max(mb2, r)
            if ba * bb < 0:
                flips += 1
        if sa is not None and sb is not None:
            r = rel(sa, sb)
            ms = max(ms, r)
            if not (0.5 < mac(a[3]) < 1.5):
                ms2 = max(ms2, r)
        pa, pb = num(a[2]), num(b[2])
        if pa is not None and pb is not None:
            for t in cross:
                if pa < t <= pb:
                    cross[t][0] += 1
                if pb < t <= pa:
                    cross[t][1] += 1
    print(f"== {os.path.basename(sd.rstrip('/'))}: {len(common)} common rows (only fp64 {len(A) - len(common)}, only test {len(B) - len(common)})")
    print(f"   BETA differ {nb} rows, max rel {mb:.3g} (MAC != 1: {mb2:.3g}); SE differ {ns} rows, max rel {ms:.3g} (MAC != 1: {ms2:.3g})")
    print(f"   BETA sign flips {flips}; p.value differ {npd} rows; crossings 5e-8 {cross[5e-8]}, 1e-5 {cross[1e-5]} (fp64-only, test-only)")
    for lab in ("fp64", "test"):
        app, conv, st, dev = logstats(os.path.join(sd, lab, "log.txt"))
        print(f"   {lab}: Firth applied {app}, converged {conv} (not converged {app - conv})")
        for s in st:
            print(f"      {s}")
        for s in dev:
            print(f"      {s}")
