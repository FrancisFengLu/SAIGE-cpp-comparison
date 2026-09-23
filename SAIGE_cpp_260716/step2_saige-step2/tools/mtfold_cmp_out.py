#!/usr/bin/env python3
"""Field-by-field comparison of two step-2 output directories.

usage: cmp_out.py DIR_A DIR_B [--prefix q_]
Reports, over every (marker, trait) pair: byte identity, max abs/rel delta per
numeric field, max |d(-log10 p)|, and any pair that crosses 5e-8 in one run but
not the other.
"""
import sys, os, math, glob

NUM = ["AC_Allele2","AF_Allele2","MissingRate","BETA","SE","Tstat","var","p.value","N"]

def read(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {h:i for i,h in enumerate(hdr)}
        rows = [l.rstrip("\n").split("\t") for l in f]
    return hdr, idx, rows

def nlp(s):
    # p.value is printed as %.6E, or "%.1fE%d" for underflow
    v = float(s)
    if v > 0: return -math.log10(v)
    return float('inf')

def main():
    A, B = sys.argv[1], sys.argv[2]
    pref = None
    if len(sys.argv) > 3 and sys.argv[3].startswith("--prefix"):
        pref = sys.argv[4] if len(sys.argv) > 4 else sys.argv[3].split("=")[1]
    files = sorted(os.path.basename(p) for p in glob.glob(A + "/*.txt"))
    if pref: files = [f for f in files if f.startswith(pref)]
    agg = {k: [0.0, 0.0] for k in NUM}     # max abs, max rel
    maxdlp = 0.0; maxdlp_where = ""
    nbyte = 0; nrows = 0; ndiffrows = 0
    cross = []
    for fn in files:
        pa, pb = A + "/" + fn, B + "/" + fn
        ba, bb = open(pa,'rb').read(), open(pb,'rb').read()
        if ba == bb:
            nbyte += 1
        ha, ia, ra = read(pa); hb, ib, rb = read(pb)
        assert ha == hb, fn
        assert len(ra) == len(rb), (fn, len(ra), len(rb))
        nrows += len(ra)
        for r1, r2 in zip(ra, rb):
            assert r1[ia["MarkerID"]] == r2[ib["MarkerID"]]
            if r1 != r2: ndiffrows += 1
            for k in NUM:
                s1, s2 = r1[ia[k]], r2[ib[k]]
                if s1 == s2: continue
                v1, v2 = float(s1), float(s2)
                d = abs(v1 - v2)
                agg[k][0] = max(agg[k][0], d)
                den = max(abs(v1), abs(v2))
                if den > 0: agg[k][1] = max(agg[k][1], d/den)
            p1, p2 = float(r1[ia["p.value"]]), float(r2[ib["p.value"]])
            if p1 != p2:
                d = abs(nlp(r1[ia["p.value"]]) - nlp(r2[ib["p.value"]]))
                if d > maxdlp: maxdlp, maxdlp_where = d, "%s %s p=%g/%g" % (fn, r1[ia["MarkerID"]], p1, p2)
            if (p1 < 5e-8) != (p2 < 5e-8):
                cross.append((fn, r1[ia["MarkerID"]], p1, p2))
    print("files %d, byte-identical %d/%d, rows %d, rows with any field change %d (%.4f%%)"
          % (len(files), nbyte, len(files), nrows, ndiffrows, 100.0*ndiffrows/max(nrows,1)))
    for k in NUM:
        if agg[k][0] or agg[k][1]:
            print("  %-12s max|d| %.3e   max rel %.3e" % (k, agg[k][0], agg[k][1]))
    print("  max |d(-log10 p)| = %.3e   %s" % (maxdlp, maxdlp_where))
    print("  pairs crossing 5e-8: %d" % len(cross))
    for c in cross[:10]: print("    ", c)

main()
