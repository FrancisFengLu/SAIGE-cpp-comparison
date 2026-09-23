#!/usr/bin/env python3
"""Full-precision (fp64) comparison of two runs' .sgs trait files."""
import sys, os, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from mtfold_sgsread import read_trait

A, B = sys.argv[1], sys.argv[2]
pref = sys.argv[3] if len(sys.argv) > 3 else ""
files = sorted(os.path.basename(p) for p in glob.glob(A + "/*.sgs"))
files = [f for f in files if f.startswith(pref) and ".markers." not in f]
agg = {}
maxdlp = 0.0; where = ""; ncross = 0; nbyte = 0
npair = 0; ndiff = 0
for fn in files:
    ba, bb = open(A+"/"+fn,'rb').read(), open(B+"/"+fn,'rb').read()
    if ba == bb: nbyte += 1
    da, _ = read_trait(A+"/"+fn); db, _ = read_trait(B+"/"+fn)
    n = len(da["BETA"]); npair += n
    ndiff += int(np.sum(da["BETA"] != db["BETA"]))
    for k in ("BETA","SE","Tstat","var"):
        x, y = da[k], db[k]
        d = np.abs(x-y); den = np.maximum(np.abs(x), np.abs(y))
        rel = np.where(den > 0, d/np.maximum(den,1e-300), 0.0)
        m = agg.setdefault(k, [0.0,0.0])
        m[0] = max(m[0], float(d.max())); m[1] = max(m[1], float(rel.max()))
    # -log10 p recomputed at fp64 from Tstat and var (the printed p is only 7 digits)
    for tag, (x1,v1,x2,v2) in (("",(da["Tstat"],da["var"],db["Tstat"],db["var"])),):
        c1 = x1*x1/v1; c2 = x2*x2/v2
        l1 = neg = None
    def nlp(chi):
        out = np.empty_like(chi)
        for i, c in enumerate(chi):
            if not np.isfinite(c) or c < 0: out[i] = np.nan; continue
            z = math.sqrt(c/2.0); e = math.erfc(z)
            out[i] = -math.log10(e) if e > 1e-300 else (z*z)/math.log(10)+math.log10(z*math.sqrt(math.pi))
        return out
    c1 = da["Tstat"]**2/da["var"]; c2 = db["Tstat"]**2/db["var"]
    sel = np.where(c1 != c2)[0]
    if sel.size:
        l1 = nlp(c1[sel]); l2 = nlp(c2[sel])
        d = np.abs(l1-l2); j = int(np.nanargmax(d))
        if d[j] > maxdlp:
            maxdlp = float(d[j]); where = "%s row %d chi2=%.6f" % (fn, sel[j], c1[sel[j]])
        # 5e-8 crossing: chi2 for p=5e-8 is 29.7168
        thr = 29.716821
        ncross += int(np.sum((c1[sel] > thr) != (c2[sel] > thr)))
print("files %d, byte-identical %d/%d, pairs %d, pairs with BETA change %d (%.4f%%)"
      % (len(files), nbyte, len(files), npair, ndiff, 100.0*ndiff/max(npair,1)))
for k in ("BETA","SE","Tstat","var"):
    print("  %-6s max|d| %.3e   max rel %.3e" % (k, agg[k][0], agg[k][1]))
print("  max |d(-log10 p)| (fp64, from Tstat^2/var) = %.3e   %s" % (maxdlp, where))
print("  pairs crossing p=5e-8: %d" % ncross)
