#!/usr/bin/env python3
"""Compare a CPU step-2 output directory against a GPU one, binary traits
included (S2_BINARY_GPU.md accuracy gates).

Per trait file (same basename in both directories):

  EXACT   row set and order; CHR POS MarkerID Allele1 Allele2; AC_Allele2
          AF_Allele2 MissingRate; N / N_case / N_ctrl; AF_case / AF_ctrl and the
          four N_*_hom/het columns (integer counts, or an exact replay of the
          sequential sum); Is.SPA.
  ROUTING "took SPA and it converged" (Is.SPA) and "SPA changed the p-value"
          (p.value != p.value.NA) must agree pair for pair. With route dumps
          (SAIGE_STEP2_ROUTE_DUMP, one byte per pair: bit0 needSPA, bit1
          needFirth, bit2 needFast, bit3 ER, bit4 isSPAConverge, bit5 is_Firth,
          bit6 is_FirthConverge) the gate decisions themselves are compared.
  APPROX  BETA SE Tstat var p.value p.value.NA: max |delta| and max relative
          |delta| per field, max |delta(-log10 p)|, and how many pairs cross
          5e-8 and 1e-5 in either direction.

usage: s2bingpu_cmp.py --cpu DIR --gpu DIR [--routes-cpu DIR --routes-gpu DIR]
                       [--sgs-cpu DIR --sgs-gpu DIR] [--tol-logp 1e-8] [--tag NAME]
Exit 1 if any exact column differs, any routing differs, or the p-value
tolerance is exceeded.
"""
import argparse, glob, math, os, sys
import numpy as np
RT = np.dtype([("r", "u1"), ("p", "<f8")])

EXACT = ["CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2",
         "MissingRate", "imputationInfo", "N", "N_case", "N_ctrl", "AF_case", "AF_ctrl",
         "N_case_hom", "N_case_het", "N_ctrl_hom", "N_ctrl_het", "Is.SPA"]
APPROX = ["BETA", "SE", "Tstat", "var", "p.value", "p.value.NA"]


def read(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in f]
    return head, rows


def pnum(s):
    try:
        return float(s)
    except ValueError:
        return float("nan")


def neglog10(s):
    v = pnum(s)
    if math.isnan(v):
        return float("nan")
    if v > 0:
        return -math.log10(v)
    return float("inf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", required=True)
    ap.add_argument("--gpu", required=True)
    ap.add_argument("--routes-cpu")
    ap.add_argument("--routes-gpu")
    ap.add_argument("--tol-logp", type=float, default=1e-8)
    ap.add_argument("--tag", default="")
    ap.add_argument("--sgs-cpu", help="dir with the CPU run's .sgs trait files (full-precision compare)")
    ap.add_argument("--sgs-gpu", help="dir with the GPU run's .sgs trait files")
    a = ap.parse_args()

    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(a.cpu, "*.txt")))
    if not files:
        print("no *.txt in", a.cpu); sys.exit(2)
    bad = 0
    tot = 0
    agg = {f: [0.0, 0.0] for f in APPROX}      # max |d|, max rel |d|
    maxlogp = {"p.value": 0.0, "p.value.NA": 0.0}
    cross = {"p.value": [0, 0], "p.value.NA": [0, 0]}  # 5e-8, 1e-5
    rout_mis = 0
    isspa_mis = 0
    spachg_mis = 0
    nspa_cpu = nspa_gpu = 0
    nchg = 0
    exact_mis = {}
    for fn in files:
        fg = os.path.join(a.gpu, fn)
        if not os.path.exists(fg):
            print(f"  {fn}: missing in --gpu"); bad += 1; continue
        hc, rc = read(os.path.join(a.cpu, fn))
        hg, rg = read(fg)
        if hc != hg:
            print(f"  {fn}: header differs"); bad += 1; continue
        if len(rc) != len(rg):
            print(f"  {fn}: {len(rc)} vs {len(rg)} rows"); bad += 1; continue
        col = {h: k for k, h in enumerate(hc)}
        for r1, r2 in zip(rc, rg):
            tot += 1
            for h in EXACT:
                if h in col and r1[col[h]] != r2[col[h]]:
                    exact_mis[h] = exact_mis.get(h, 0) + 1
                    if exact_mis[h] <= 3:
                        print(f"  EXACT {fn} {r1[col['MarkerID']]} {h}: {r1[col[h]]} vs {r2[col[h]]}")
            if "Is.SPA" in col:
                s1, s2 = r1[col["Is.SPA"]], r2[col["Is.SPA"]]
                nspa_cpu += (s1 == "true"); nspa_gpu += (s2 == "true")
                if s1 != s2: isspa_mis += 1
            if "p.value.NA" in col:
                c1 = r1[col["p.value"]] != r1[col["p.value.NA"]]
                c2 = r2[col["p.value"]] != r2[col["p.value.NA"]]
                nchg += c1
                if c1 != c2:
                    spachg_mis += 1
                    if spachg_mis <= 3:
                        print(f"  SPA-changed-p differs {fn} {r1[col['MarkerID']]}: cpu {r1[col['p.value']]}/{r1[col['p.value.NA']]} gpu {r2[col['p.value']]}/{r2[col['p.value.NA']]}")
            for h in APPROX:
                if h not in col: continue
                v1, v2 = pnum(r1[col[h]]), pnum(r2[col[h]])
                if math.isnan(v1) and math.isnan(v2): continue
                if math.isnan(v1) != math.isnan(v2):
                    agg[h][0] = float("inf"); continue
                d = abs(v1 - v2)
                if d > agg[h][0]: agg[h][0] = d
                sc = max(abs(v1), abs(v2))
                if sc > 0 and d / sc > agg[h][1]: agg[h][1] = d / sc
                if h in maxlogp:
                    l1, l2 = neglog10(r1[col[h]]), neglog10(r2[col[h]])
                    if l1 == l2: continue
                    if math.isinf(l1) or math.isinf(l2):
                        dl = float("inf")
                    else:
                        dl = abs(l1 - l2)
                    if dl > maxlogp[h]: maxlogp[h] = dl
                    for k, thr in enumerate((5e-8, 1e-5)):
                        if (v1 < thr) != (v2 < thr): cross[h][k] += 1
    # route dumps
    rtot = 0
    gate_rel = 0.0; gate_dl = 0.0; gate_n = 0; gate_cross = {5e-8: 0, 1e-5: 0}
    if a.routes_cpu and a.routes_gpu:
        for fn in files:
            base = fn[:-4]
            p1 = os.path.join(a.routes_cpu, base + ".route")
            p2 = os.path.join(a.routes_gpu, base + ".route")
            if not (os.path.exists(p1) and os.path.exists(p2)):
                print(f"  route dump missing for {base}"); bad += 1; continue
            d1 = np.fromfile(p1, dtype=RT); d2 = np.fromfile(p2, dtype=RT)
            if len(d1) != len(d2):
                print(f"  route dump length differs for {base}: {len(d1)} vs {len(d2)}"); bad += 1; continue
            rtot += len(d1)
            b1, b2 = d1["r"], d2["r"]
            if not np.array_equal(b1, b2):
                w = np.nonzero(b1 != b2)[0]
                rout_mis += len(w)
                k = int(w[0])
                print(f"  ROUTE {base}: {len(w)} pairs differ, first at row {k}: cpu {int(b1[k]):#04x} gpu {int(b2[k]):#04x}")
            # full-precision pre-SPA p-value of the batch kernel, where both have it
            g1, g2 = d1["p"], d2["p"]
            m = ~np.isnan(g1) & ~np.isnan(g2)
            if (np.isnan(g1) != np.isnan(g2)).any():
                print(f"  gateP presence differs for {base}"); bad += 1
            if m.any():
                x, y = g1[m], g2[m]
                with np.errstate(divide="ignore", invalid="ignore"):
                    rel = np.abs(x - y) / np.maximum(np.maximum(np.abs(x), np.abs(y)), 1e-300)
                    lx = np.where(x > 0, -np.log10(np.where(x > 0, x, 1)), np.abs(x) / math.log(10))   # log-domain values are negative
                    ly = np.where(y > 0, -np.log10(np.where(y > 0, y, 1)), np.abs(y) / math.log(10))
                    dl = np.abs(lx - ly)
                gate_rel = max(gate_rel, float(rel.max())); gate_dl = max(gate_dl, float(dl.max()))
                gate_n += int(m.sum())
                for thr in (5e-8, 1e-5):
                    gate_cross[thr] += int(((x < thr) != (y < thr)).sum())
    tag = f"[{a.tag}] " if a.tag else ""
    sgs_lines = []
    if a.sgs_cpu and a.sgs_gpu:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from mtfold_sgsread import read_trait
        sg = {}
        for fn in files:
            base = fn[:-4]
            p1 = os.path.join(a.sgs_cpu, base + ".txt.sgs"); p2 = os.path.join(a.sgs_gpu, base + ".txt.sgs")
            if not (os.path.exists(p1) and os.path.exists(p2)):
                sgs_lines.append(f"  sgs missing for {base}"); bad += 1; continue
            d1, n1 = read_trait(p1); d2, n2 = read_trait(p2)
            if n1 != n2: sgs_lines.append(f"  sgs rows differ for {base}"); bad += 1; continue
            for k in ("BETA", "SE", "Tstat", "var", "AF_case", "AF_ctrl"):
                if k not in d1: continue
                x, y = d1[k], d2[k]
                both = np.isnan(x) & np.isnan(y)
                d = np.where(both, 0.0, np.abs(x - y))
                den = np.maximum(np.abs(x), np.abs(y))
                with np.errstate(invalid="ignore", divide="ignore"):
                    rel = np.where(den > 0, d / np.maximum(den, 1e-300), 0.0)
                m = sg.setdefault(k, [0.0, 0.0, 0])
                m[0] = max(m[0], float(np.nanmax(d)) if d.size else 0.0)
                m[1] = max(m[1], float(np.nanmax(rel)) if rel.size else 0.0)
                m[2] += int((d > 0).sum())
        for k, m in sg.items():
            sgs_lines.append(f"{tag}sgs fp64 {k:8s} max|d| {m[0]:.3e}  max rel {m[1]:.3e}  pairs differing {m[2]:,}")
            if k in ("AF_case", "AF_ctrl") and m[2]: bad += 1
    print(f"{tag}{len(files)} files, {tot:,} pairs")
    print(f"{tag}exact columns: " + ("all identical" if not exact_mis else
          ", ".join(f"{h} {n} differ" for h, n in exact_mis.items())))
    print(f"{tag}routing: Is.SPA mismatch {isspa_mis} (TRUE: cpu {nspa_cpu:,} gpu {nspa_gpu:,}); "
          f"SPA-changed-p mismatch {spachg_mis} (changed in cpu: {nchg:,})"
          + (f"; gate route bytes compared {rtot:,}, differ {rout_mis}" if rtot else ""))
    if gate_n:
        print(f"{tag}batch-kernel p (fp64, {gate_n:,} gated pairs): max rel {gate_rel:.3e}  max|d(-log10 p)| {gate_dl:.3e}"
              f"  cross 5e-8: {gate_cross[5e-8]}  cross 1e-5: {gate_cross[1e-5]}")
    for h in APPROX:
        if agg[h][0] == 0.0 and agg[h][1] == 0.0:
            print(f"{tag}{h:11s} max|d| 0 (identical)")
        else:
            extra = f"  max|d(-log10 p)| {maxlogp[h]:.3e}  cross 5e-8: {cross[h][0]}  cross 1e-5: {cross[h][1]}" if h in maxlogp else ""
            print(f"{tag}{h:11s} max|d| {agg[h][0]:.3e}  max rel {agg[h][1]:.3e}{extra}")
    if exact_mis: bad += 1
    if isspa_mis or spachg_mis or rout_mis: bad += 1
    for h in maxlogp:
        if maxlogp[h] > a.tol_logp: bad += 1
    if gate_n and gate_dl > a.tol_logp: bad += 1
    for L in sgs_lines: print(L)
    print(f"{tag}RESULT: " + ("PASS" if bad == 0 else "FAIL"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
