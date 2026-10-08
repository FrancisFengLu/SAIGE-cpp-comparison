#!/usr/bin/env python3
"""redact.py DIR [--fam F]... [--psam F]... [--bim F]... [--pvar F]... [--pheno F:COL]...

Last safety net before results leave the analysis environment. Every text file under DIR is
scanned; a line that contains a sample ID (from .fam / .psam / the phenotype file's ID column)
or a marker ID (from .bim / .pvar) as a token is replaced by
    [line removed by redact.py: contained a sample ID | marker ID]
Writes DIR/REDACTION_REPORT.txt (how many lines were removed per file and why; no IDs).
Sample IDs shorter than --min-len characters (default 5) are not searched for, because short
numbers occur everywhere in logs; the report says how many IDs that skipped.
"""
import argparse, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("dir")
for k in ("fam", "psam", "bim", "pvar", "pheno"):
    ap.add_argument("--" + k, action="append", default=[])
ap.add_argument("--min-len", type=int, default=5)
a = ap.parse_args()

samples, markers = set(), set()
for f in a.fam:
    for l in open(f):
        t = l.split()
        if len(t) >= 2:
            samples.update(t[:2])
for f in a.psam:
    hdr = None
    for l in open(f):
        t = l.split()
        if l.startswith("#"):
            hdr = [x.lstrip("#") for x in t]; continue
        if hdr and "IID" in hdr:
            samples.add(t[hdr.index("IID")])
            if "FID" in hdr: samples.add(t[hdr.index("FID")])
        elif t:
            samples.add(t[0])
for spec in a.pheno:
    f, col = spec.rsplit(":", 1)
    with open(f) as fh:
        hdr = re.split(r"[\t ,]+", fh.readline().strip())
        i = hdr.index(col)
        for l in fh:
            t = re.split(r"[\t ,]+", l.strip())
            if len(t) > i:
                samples.add(t[i])
for f in a.bim:
    for l in open(f):
        t = l.split()
        if len(t) >= 2:
            markers.add(t[1])
for f in a.pvar:
    idc = 2
    for l in open(f):
        if l.startswith("##"):
            continue
        t = l.split()
        if l.startswith("#"):
            h = [x.lstrip("#") for x in t]; idc = h.index("ID") if "ID" in h else 2; continue
        if len(t) > idc:
            markers.add(t[idc])
nshort = sum(1 for s in samples if len(s) < a.min_len)
samples = {s for s in samples if len(s) >= a.min_len}
markers = {m for m in markers if m not in (".", "") and len(m) >= 3}
SPLIT = re.compile(r"[\s,;'\"()\[\]{}=<>|]+")

report = []
tot = 0
for root, _, files in os.walk(a.dir):
    for fn in sorted(files):
        p = os.path.join(root, fn)
        try:
            lines = open(p, encoding="utf-8").read().split("\n")
        except (UnicodeDecodeError, OSError):
            report.append("%s: not text, left out" % os.path.relpath(p, a.dir)); os.remove(p); continue
        ns = nm = 0
        out = []
        for line in lines:
            toks = set()
            for t in SPLIT.split(line):
                if t:
                    toks.add(t); toks.add(t.strip(".:"))
            if toks & samples:
                out.append("[line removed by redact.py: contained a sample ID]"); ns += 1
            elif toks & markers:
                out.append("[line removed by redact.py: contained a marker ID]"); nm += 1
            else:
                out.append(line)
        if ns or nm:
            open(p, "w").write("\n".join(out))
            report.append("%s: %d line(s) with a sample ID, %d with a marker ID removed" % (os.path.relpath(p, a.dir), ns, nm))
        tot += ns + nm
with open(os.path.join(a.dir, "REDACTION_REPORT.txt"), "w") as f:
    f.write("searched for %d sample IDs and %d marker IDs (%d sample IDs shorter than %d characters not searched)\n"
            % (len(samples), len(markers), nshort, a.min_len))
    f.write("lines removed in total: %d\n" % tot)
    f.write("\n".join(report) + ("\n" if report else ""))
print(open(os.path.join(a.dir, "REDACTION_REPORT.txt")).read(), end="")
