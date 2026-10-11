#!/usr/bin/env python3
"""firth_table.py [FIRTHB]: the Firth-status gate (firthstress: gen_firthstress.py, r_runs_firth.sh,
cpp_runs_firth.sh) as Markdown --
  1. per trait, the Firth counts (fits / strictly converged / maxit / singular) from the logs of the
     CPU single-trait, CPU multi-trait, GPU fp64 and GPU fp32 runs, and R's "applied to / converged";
  2. the per-row Firth.Status of the CPU single-trait run against the GPU runs (fp64, fp32) and the
     CPU multi-trait run, and the .sgs -> sgs2txt round trip (byte compare);
  3. BETA / SE of the Firth rows against R, separately on the converged and the maxit rows: max abs
     and relative difference, sign flips; p-values identical or not;
  4. GPU fp32 vs fp64 on the converged and the maxit rows;
  5. Is.SPA on the MAC <= 4 rows against R (the exact-test rows), and the default-format run's header.
Paths: FIRTHB (default /opt/saige/logs/integrate/firth)."""
import sys, os, re, math
B = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FIRTHB", "/opt/saige/logs/integrate/firth")
TR = ["b%d" % k for k in range(1, 9)]
def read(path):
    rows = {}
    if not os.path.exists(path): return None, rows
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for l in f:
            x = l.rstrip("\n").split("\t")
            rows[x[ix["MarkerID"]]] = x
    return ix, rows
def counts(logpath, trait=None):
    if not os.path.exists(logpath): return None
    for l in open(logpath):
        m = re.match(r"(\[(\w+)\] )?Firth fits: (\d+); strictly converged (\d+), stopped at maxit \(50\) (\d+), singular (\d+)", l.strip())
        if m and (trait is None or m.group(2) == trait): return tuple(int(m.group(k)) for k in (3, 4, 5, 6))
    return None
def rcounts(logpath):
    if not os.path.exists(logpath): return None
    for l in open(logpath):
        m = re.search(r"Firth approx was applied to (\d+) markers. (\d+) suc", l)
        if m: return int(m.group(1)), int(m.group(2))
    return None
def fl(s):
    try: return float(s)
    except ValueError: return math.nan
print("## firthstress: N = 20,000, 2,000 common + 5,000 rare (MAC 1..30) markers, 8 binary traits with 86-272 cases "
      "(b7 / b8 with their own missing phenotypes), R defaults + is_Firth_beta, pCutoffforFirth 0.05\n")
print("### 1. Firth counts per trait (fits; strictly converged / stopped at maxit / singular), and R's line\n")
print("| trait | cases | CPU single | CPU multi | GPU fp64 | GPU fp32 | R: applied / converged |\n|---|--:|---|---|---|---|---|")
allsame = True
for t in TR:
    ix, R = read(f"{B}/R/s2/defF05/{t}.txt")
    ncase = next(iter(R.values()))[ix["N_case"]] if R else "?"
    cs = counts(f"{B}/cpp/single/log_{t}.txt"); cm = counts(f"{B}/cpp/multi/log.txt", t)
    cg = counts(f"{B}/cpp/gpu/log.txt", t); c32 = counts(f"{B}/cpp/gpu32/log.txt", t)
    rc = rcounts(f"{B}/R/s2/defF05/{t}.log")
    f = lambda c: "%d; %d / %d / %d" % c if c else "(missing)"
    if not (cs == cm == cg == c32): allsame = False
    print(f"| {t} | {ncase} | {f(cs)} | {f(cm)} | {f(cg)} | {f(c32)} | {rc[0]} / {rc[1]} |" if rc else f"| {t} | {ncase} | {f(cs)} | {f(cm)} | {f(cg)} | {f(c32)} | (missing) |")
print(f"\ncounts identical across the four C++ paths: {'yes' if allsame else 'NO'}; R counts every fit as converged (its isfirthconverge is true at maxit).\n")

print("### 2. per-row Firth.Status: CPU single-trait vs the other paths (rows with a status in either), and the .sgs round trip\n")
print("| trait | rows | vs CPU multi: status differs | vs GPU fp64 | vs GPU fp32 | single == gpu (bytes) | gpusgs -> sgs2txt == gpu (bytes) |\n|---|--:|--:|--:|--:|---|---|")
import filecmp
for t in TR:
    ixs, S = read(f"{B}/cpp/single/{t}.txt")
    out = [t, str(len(S))]
    for v in ["multi", "gpu", "gpu32"]:
        ixv, V = read(f"{B}/cpp/{v}/{t}.txt")
        if ixv is None: out.append("(missing)"); continue
        d = sum(1 for k, x in S.items() if k in V and x[ixs["Firth.Status"]] != V[k][ixv["Firth.Status"]])
        out.append(str(d))
    out.append("yes" if filecmp.cmp(f"{B}/cpp/single/{t}.txt", f"{B}/cpp/gpu/{t}.txt", shallow=False) else "no")
    out.append("yes" if os.path.exists(f"{B}/cpp/gpusgs/{t}.txt") and filecmp.cmp(f"{B}/cpp/gpusgs/{t}.txt", f"{B}/cpp/gpu/{t}.txt", shallow=False) else "no")
    print("| " + " | ".join(out) + " |")

def cmp_rows(A, ixa, Bm, ixb, keys, cols=("BETA", "SE")):
    res = {}
    for c in cols:
        mx = 0.0; mr = 0.0; flips = 0; n = 0; mxk = None
        for k in keys:
            a = fl(A[k][ixa[c]]); b = fl(Bm[k][ixb[c]])
            if math.isnan(a) or math.isnan(b): continue
            n += 1
            d = abs(a - b); r = d / max(abs(a), abs(b)) if max(abs(a), abs(b)) > 0 else 0.0
            if d > mx: mx = d; mxk = (k, a, b)
            if r > mr: mr = r
            if c == "BETA" and a * b < 0: flips += 1
        res[c] = (n, mx, mr, flips, mxk)
    return res
def peq(A, ixa, Bm, ixb, keys):
    return sum(1 for k in keys if A[k][ixa["p.value"]] != Bm[k][ixb["p.value"]])

print("\n### 3. Firth rows against R SAIGE 1.5.2 (the CPU single-trait run; the GPU fp64 run is byte-identical to it): converged rows and maxit rows separately\n")
print("| trait | rows | converged: n, BETA max abs / rel, SE max abs / rel, sign flips | maxit: n, BETA max abs / rel, SE max abs / rel, sign flips | p.value differs (all Firth rows) | largest maxit BETA gap (marker: C++, R) |\n|---|--:|---|---|--:|---|")
for t in TR:
    ixs, S = read(f"{B}/cpp/single/{t}.txt"); ixr, R = read(f"{B}/R/s2/defF05/{t}.txt")
    if ixr is None: print(f"| {t} | (R missing) |"); continue
    conv = [k for k, x in S.items() if x[ixs["Firth.Status"]] == "converged" and k in R]
    mx = [k for k, x in S.items() if x[ixs["Firth.Status"]] == "maxit" and k in R]
    rc = cmp_rows(S, ixs, R, ixr, conv); rm = cmp_rows(S, ixs, R, ixr, mx)
    f = lambda r: "%d, %.1e / %.1e, %.1e / %.1e, %d" % (r["BETA"][0], r["BETA"][1], r["BETA"][2], r["SE"][1], r["SE"][2], r["BETA"][3])
    big = rm["BETA"][4]
    bigs = "%s: %.6g, %.6g" % big if big else "-"
    print(f"| {t} | {len(S)} | {f(rc)} | {f(rm)} | {peq(S, ixs, R, ixr, conv + mx)} | {bigs} |")

print("\n### 4. GPU fp32 (gpuPrecisionFirth: fp32) vs GPU fp64, converged rows and maxit rows\n")
print("| trait | converged: n, BETA max abs / rel, SE max abs / rel, sign flips | maxit: n, BETA max abs / rel, SE max abs / rel, sign flips | status differs | p.value differs |\n|---|---|---|--:|--:|")
for t in TR:
    ixg, G = read(f"{B}/cpp/gpu/{t}.txt"); ix3, G3 = read(f"{B}/cpp/gpu32/{t}.txt")
    if ix3 is None: print(f"| {t} | (fp32 missing) |"); continue
    conv = [k for k, x in G.items() if x[ixg["Firth.Status"]] == "converged" and k in G3]
    mx = [k for k, x in G.items() if x[ixg["Firth.Status"]] == "maxit" and k in G3]
    rc = cmp_rows(G, ixg, G3, ix3, conv); rm = cmp_rows(G, ixg, G3, ix3, mx)
    sd = sum(1 for k in conv + mx if G[k][ixg["Firth.Status"]] != G3[k][ix3["Firth.Status"]])
    f = lambda r: "%d, %.1e / %.1e, %.1e / %.1e, %d" % (r["BETA"][0], r["BETA"][1], r["BETA"][2], r["SE"][1], r["SE"][2], r["BETA"][3])
    print(f"| {t} | {f(rc)} | {f(rm)} | {sd} | {peq(G, ixg, G3, ix3, conv + mx)} |")

print("\n### 5. Is.SPA on the MAC <= 4 rows (the exact-test rows) against R, and the default-format header\n")
print("| trait | MAC <= 4 rows | R: Is.SPA true / false | C++ single: true / false | Is.SPA differs | p.value differs on those rows |\n|---|--:|---|---|--:|--:|")
for t in TR:
    ixs, S = read(f"{B}/cpp/single/{t}.txt"); ixr, R = read(f"{B}/R/s2/defF05/{t}.txt")
    if ixr is None: continue
    def mac(x, ix):
        ac = fl(x[ix["AC_Allele2"]]); n = fl(x[ix["N_case"]]) + fl(x[ix["N_ctrl"]]); return min(ac, 2 * n - ac)
    low = [k for k, x in R.items() if mac(x, ixr) <= 4 and k in S]
    rt = sum(1 for k in low if R[k][ixr["Is.SPA"]] == "true"); st = sum(1 for k in low if S[k][ixs["Is.SPA"]] == "true")
    d = sum(1 for k in low if R[k][ixr["Is.SPA"]] != S[k][ixs["Is.SPA"]])
    print(f"| {t} | {len(low)} | {rt} / {len(low) - rt} | {st} / {len(low) - st} | {d} | {peq(S, ixs, R, ixr, low)} |")
h = open(f"{B}/cpp/cpu_def/b1.txt").readline().rstrip("\n").split("\t") if os.path.exists(f"{B}/cpp/cpu_def/b1.txt") else []
hr = open(f"{B}/R/s2/defF05/b1.txt").readline().rstrip("\n").split("\t") if os.path.exists(f"{B}/R/s2/defF05/b1.txt") else []
print(f"\ndefault format (outputFirthStatus off), b1 header == R's header: {'yes' if h == hr else 'NO: ' + str(h)}")
hs = open(f"{B}/cpp/single/b1.txt").readline().rstrip("\n").split("\t")
print(f"outputFirthStatus on, b1 header: {' '.join(hs)}")
