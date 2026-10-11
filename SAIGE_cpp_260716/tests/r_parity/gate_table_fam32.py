#!/usr/bin/env python3
"""gate_table_fam32.py: fam32 (32 traits, every one with its own missing-phenotype pattern, 1000 families
of 4) and fam32sp (8 of them with a sparse GRM): every C++ path x configuration against R SAIGE 1.5.2,
and the multi-trait / GPU paths against the single-trait path (worst over the traits), the per-trait
marker counts, and the device log lines. Paths: GPB (default /opt/saige/logs/gpuprep)."""
import sys, os, json, subprocess, re
B = os.environ.get("GPB", "/opt/saige/logs/gpuprep")
CMP = os.environ.get("CMP_R", "/opt/saige/logs/rdefaults/cmp_r.py")
BIN = ["b%d" % k for k in range(1, 17)]; QNT = ["q%d" % k for k in range(1, 17)]
SP = ["b1", "b2", "b3", "b4", "q1", "q2", "q3", "q4"]
cases = [("fam32", "def", BIN + QNT), ("fam32", "adj", BIN + QNT), ("fam32", "defF", BIN), ("fam32", "adjF", BIN),
         ("fam32sp", "def", SP), ("fam32sp", "adj", SP), ("fam32sp", "defT", SP), ("fam32sp", "adjT", SP),
         ("fam32sp", "defF", SP[:4]), ("fam32sp", "defTF", SP[:4])]
def cmp(a, b, label):
    if not (os.path.exists(a) and os.path.exists(b)): return None
    return json.loads(subprocess.check_output([sys.executable, CMP, a, b, "--label", label, "--json"]).decode())
def fmt(r):
    if r is None: return "| (missing) |||||||||"
    m = r["maxrel"]
    return (f"| {r['label']} | {r['rows']} | {m['p.value']:.1e} | {m['BETA']:.1e} | {m['SE']:.1e} | {m['Tstat']:.1e} | {m['var']:.1e} "
            f"| {r['beyondPrint']['p.value']}/{r['beyondPrint']['SE']}/{r['beyondPrint']['var']} | {r['cross5e8']}/{r['cross1e5']} "
            f"| {r['isSPAmismatch']} ({r['nSPA_A']}) | {r['firthMismatch']} ({r['firthA']}/{r['firthB']}) |")
hdr = ("| run | rows | p rel | BETA rel | SE rel | Tstat rel | var rel | beyond print p/SE/var | cross 5e-8/1e-5 | Is.SPA mism (SPA rows) | Firth mism (A/B rows) |\n"
       "|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
def worst(rs):
    rs = [r for r in rs if r]
    if not rs: return None
    w = dict(rs[0]); w["rows"] = sum(r["rows"] for r in rs)
    for k in ["p.value", "BETA", "SE", "Tstat", "var"]:
        w["maxrel"][k] = max(r["maxrel"][k] for r in rs); w["beyondPrint"][k] = sum(r["beyondPrint"][k] for r in rs)
    for k in ["cross5e8", "cross1e5", "isSPAmismatch", "nSPA_A", "firthMismatch", "firthA", "firthB"]:
        w[k] = sum(r[k] for r in rs)
    return w
def tested(logpath, trait=None):
    if not os.path.exists(logpath): return None
    for l in open(logpath):
        m = re.match(r"(\[(\w+)\] )?(\d+) markers were tested", l.strip())
        if m and (trait is None or m.group(2) == trait): return int(m.group(3))
    return None
print("## fam32 / fam32sp: C++ vs R SAIGE 1.5.2 (worst over the traits)\n")
print(hdr)
cnt = []
for ds, cfg, traits in cases:
    for path in ["single", "multi", "gpu"]:
        rs = []
        for t in traits:
            a = f"{B}/cpp/{ds}/{cfg}/{path}/{t}.txt"
            rs.append(cmp(a, f"{B}/R/{ds}/{cfg}/{t}.txt", f"{ds} {cfg} {path} {t}"))
            if path != "single":
                cnt.append((ds, cfg, path, t, tested(f"{B}/cpp/{ds}/{cfg}/{path}/log.txt", t),
                            tested(f"{B}/cpp/{ds}/{cfg}/single/log_{t}.txt")))
        w = worst(rs)
        if w: w["label"] = f"{ds} {cfg} {path} vs R"
        print(fmt(w))
print("\n## multi-trait and GPU paths vs the single-trait path\n")
print(hdr)
for ds, cfg, traits in cases:
    for path in ["multi", "gpu"]:
        rs = [cmp(f"{B}/cpp/{ds}/{cfg}/{path}/{t}.txt", f"{B}/cpp/{ds}/{cfg}/single/{t}.txt", f"{ds} {cfg} {path} {t}") for t in traits]
        w = worst(rs)
        if w: w["label"] = f"{ds} {cfg} {path} vs single"
        print(fmt(w))
print("\n## markers tested per trait (multi / gpu vs single)\n")
bad = [c for c in cnt if c[4] != c[5]]
print(f"{len(cnt)} cells; differing counts: {len(bad)}" + (": " + str(bad[:10]) if bad else ""))
print("\n## device / host log lines (GPU runs)\n")
for ds, cfg, traits in cases:
    L = f"{B}/cpp/{ds}/{cfg}/gpu/log.txt"
    if not os.path.exists(L): continue
    for l in open(L):
        s = l.strip()
        if (s.startswith("gpuDeviceStats: ") and ("pairs took" in s or "on --" in s)) or s.startswith("gpuPrep:") or s.startswith("GPU coverage") \
           or s.startswith("gpuSparse: the sparse statistic") or s.startswith("gpuSparse: cross-term") or s.startswith("gate:") or s.startswith("useGPU: refused") \
           or s.startswith("gpuOwnSampleSets:"):
            print(f"- {ds} {cfg}: {s[:240]}")
