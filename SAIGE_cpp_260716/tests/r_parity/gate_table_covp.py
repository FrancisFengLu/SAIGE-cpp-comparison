#!/usr/bin/env python3
"""gate_table_covp.py OUT [variants...]: the covariate-count gate (cpp_runs_covp.sh) as Markdown --
per p (4 9 13 24 40) and configuration, every GPU variant against R SAIGE 1.5.2 and against this
binary's CPU single-trait path (worst over the run's traits), the single path against R, the
markers tested per trait, and the device log lines (gpuSpa / gpuFirth setup, device SPA / Firth
pair counts and us per pair)."""
import sys, os, json, subprocess, re
OUT = sys.argv[1]
VARS = sys.argv[2:] or ["same", "own", "sameown", "same32", "own32"]
B = os.environ.get("COVLIM", "/opt/saige/logs/covlimit")
PS = os.environ.get("PS", "4 9 13 24 40").split()
CMP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cmp_r.py")
CFGS = ["def", "defF", "adj", "adjF"]
TR = {"same": ["b1", "b2"], "own": ["b1", "b2", "b3", "b4"], "sameown": ["b1", "b2"], "same32": ["b1", "b2"], "own32": ["b1", "b2", "b3", "b4"]}
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
    w = json.loads(json.dumps(rs[0])); w["rows"] = sum(r["rows"] for r in rs)
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
print("## CPU single-trait path vs R SAIGE 1.5.2 (worst over b1..b4)\n")
print(hdr)
for p in PS:
    for cfg in CFGS:
        rs = [cmp(f"{OUT}/p{p}/{cfg}/single/{t}.txt", f"{B}/R/covp/s2/p{p}/{cfg}/{t}.txt", f"p={p} {cfg} single {t}") for t in TR["own"]]
        w = worst(rs)
        if w: w["label"] = f"p={p} {cfg} single vs R"
        print(fmt(w))
for v in VARS:
    print(f"\n## GPU variant `{v}` (traits {' '.join(TR[v])}) vs R and vs the single-trait path\n")
    print(hdr)
    for p in PS:
        for cfg in CFGS:
            rsR = [cmp(f"{OUT}/p{p}/{cfg}/{v}/{t}.txt", f"{B}/R/covp/s2/p{p}/{cfg}/{t}.txt", f"p={p} {cfg} {v} {t}") for t in TR[v]]
            rsS = [cmp(f"{OUT}/p{p}/{cfg}/{v}/{t}.txt", f"{OUT}/p{p}/{cfg}/single/{t}.txt", f"p={p} {cfg} {v} {t}") for t in TR[v]]
            for rs, lab in [(rsR, "vs R"), (rsS, "vs single")]:
                w = worst(rs)
                if w: w["label"] = f"p={p} {cfg} {v} {lab}"
                print(fmt(w))
print("\n## markers tested per trait (GPU variant / single)\n")
bad = []; n = 0
for v in VARS:
    for p in PS:
        for cfg in CFGS:
            for t in TR[v]:
                ng = tested(f"{OUT}/p{p}/{cfg}/{v}/log.txt", t); ns = tested(f"{OUT}/p{p}/{cfg}/single/log_{t}.txt")
                n += 1
                if ng != ns: bad.append((v, p, cfg, t, ng, ns))
print(f"{n} (variant, p, cfg, trait) cells; differing counts: {len(bad)}" + (": " + str(bad) if bad else ""))
print("\n## device log lines (setup, pairs, us per pair)\n")
print("| p | cfg | variant | gpuSpa | gpuFirth | device SPA | device Firth |\n|--:|---|---|---|---|---|---|")
def pick(lines, key):
    for s in lines:
        if s.startswith(key): return s[len(key):].strip()
    return "-"
for v in VARS:
    for p in PS:
        for cfg in CFGS:
            L = f"{OUT}/p{p}/{cfg}/{v}/log.txt"
            if not os.path.exists(L): continue
            lines = [l.strip() for l in open(L)]
            spa = pick(lines, "gpuSpa:")
            spa = re.sub(r", up to \d+ pairs per device batch.*", "", spa)
            spa = re.sub(r"^gpu/spa_gpu library \(gpuSpaImpl: lib\), ", "lib, ", spa)
            spa = re.sub(r"^gpu/gpu_spa.cu kernel \(gpuSpaImpl: own\), ", "own, ", spa)
            fi = pick(lines, "gpuFirth:")
            fi = re.sub(r", up to \d+ pairs per device batch.*", "", fi)
            ds = pick(lines, "device SPA:")
            ds = re.sub(r" s of kernel time \(", " s (", ds)
            df = pick(lines, "device Firth:")
            df = re.sub(r" s of kernel time \(", " s (", df)
            print(f"| {p} | {cfg} | {v} | {spa[:90]} | {fi[:80]} | {ds[:120]} | {df[:120]} |")
