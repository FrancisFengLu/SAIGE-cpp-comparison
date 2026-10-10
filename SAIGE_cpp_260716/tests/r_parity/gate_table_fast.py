#!/usr/bin/env python3
"""gate_table_fast.py OUT [datasets...]: the fast-test runs under OUT/<ds>/<cfg>/{single,gpu} (cpp_fast_runs.sh):
the GPU multi-trait path against R SAIGE 1.5.2 (is_fastTest TRUE; r_runs_fast.sh / r_runs_fam32.sh) and against
the single-trait scalar path (worst over the dataset's traits), the per-trait marker counts, and the device log lines (dense recompute pairs on the device, host
hand-backs by reason)."""
import sys, os, json, subprocess, re
OUT = sys.argv[1]
B = os.environ.get("RDEF", "/opt/saige/logs/rdefaults")
G = os.environ.get("GPB", "/opt/saige/logs/gpuprep")
CMP = os.path.join(B, "cmp_r.py")
DS = sys.argv[2:] or ["audit", "bvs", "qt12", "fam32", "fam32sp"]
BINT = ["b%d" % k for k in range(1, 17)]; QNT = ["q%d" % k for k in range(1, 17)]
def traits(ds, cfg):
    if ds in ("audit", "bvs"): return ["b1", "b2", "b3", "b4"]
    if ds == "qt12": return ["b1", "b2"] if cfg.endswith("F") else ["q1", "q2", "b1", "b2"]
    if ds == "fam32": return BINT if cfg.endswith("F") else BINT + QNT
    if ds == "fam32sp": return ["b1", "b2", "b3", "b4"] if cfg.endswith("F") else ["b1", "b2", "b3", "b4", "q1", "q2", "q3", "q4"]
def cfgs(ds): return ["defT", "adjT", "defTF"] if ds == "fam32sp" else ["defT", "defTF", "adjT", "adjTF"]
def rdir(ds, cfg): return f"{G}/R/{ds}/{cfg}" if ds.startswith("fam32") else f"{B}/R/{ds}/{cfg}"
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
print("## fast test on: GPU multi-trait path vs R SAIGE 1.5.2 (is_fastTest TRUE) and vs the single-trait scalar path\n")
print(hdr)
cnt = []
for ds in DS:
    for cfg in cfgs(ds):
        rsR, rsS = [], []
        for t in traits(ds, cfg):
            g = f"{OUT}/{ds}/{cfg}/gpu/{t}.txt"
            rsR.append(cmp(g, f"{rdir(ds, cfg)}/{t}.txt", f"{ds} {cfg} gpu {t}"))
            rsS.append(cmp(g, f"{OUT}/{ds}/{cfg}/single/{t}.txt", f"{ds} {cfg} gpu {t}"))
            cnt.append((ds, cfg, t, tested(f"{OUT}/{ds}/{cfg}/gpu/log.txt", t), tested(f"{OUT}/{ds}/{cfg}/single/log_{t}.txt")))
        for rs, lab in [(rsR, "vs R"), (rsS, "vs single")]:
            w = worst(rs)
            if w: w["label"] = f"{ds} {cfg} gpu {lab}"
            print(fmt(w))
print("\n## markers tested per trait (GPU / single)\n")
bad = [c for c in cnt if c[3] != c[4]]
print(f"{len(cnt)} (dataset, cfg, trait) cells; differing counts: {len(bad)}" + (": " + str(bad[:10]) if bad else ""))
print("\n## device / host log lines (GPU runs)\n")
for ds in DS:
    for cfg in cfgs(ds):
        L = f"{OUT}/{ds}/{cfg}/gpu/log.txt"
        if not os.path.exists(L): continue
        for l in open(L):
            s = l.strip()
            if (s.startswith("gpuDeviceStats: ") and ("pairs took" in s or "on --" in s)) or s.startswith("fast-test recompute") \
               or s.startswith("gpuSparse: cross-term") or s.startswith("gate:") or s.startswith("useGPU: refused") or s.startswith("GPU coverage"):
                print(f"- {ds} {cfg}: {s[:300]}")
