#!/usr/bin/env python3
"""compare_outputs.py DIR_A DIR_B --label L --json OUT.json [--rows]

Compare two sets of step-2 single-variant result files (one <trait>.txt per trait in each
directory; C++ CPU vs C++ GPU, or C++ vs R SAIGE 1.5.2). Files are paired by name.

1. md5 of every pair. If all pairs are byte-identical and --rows is not given, stop there.
2. Otherwise (or with --rows), a row-level comparison per trait. Rows are matched on
   (CHR, POS, MarkerID, Allele1, Allele2). Reported, per trait and in total:
   rows in A / B, rows only in A / only in B, rows whose text differs, and for p.value, BETA, SE:
   number of rows that differ, max |A-B|, max |A-B| / max(|A|,|B|); max |-log10 pA + log10 pB|;
   number of p < 5e-8 (and < 1e-5) in A and in B and how many cross the threshold in one but not the
   other; Is.SPA agreement.

Only aggregate numbers are written -- no marker IDs, positions or per-variant values -- so the
JSON can leave the analysis environment.
"""
import argparse, glob, hashlib, json, math, os

ap = argparse.ArgumentParser()
ap.add_argument("a"); ap.add_argument("b")
ap.add_argument("--label", default="")
ap.add_argument("--json", required=True)
ap.add_argument("--rows", action="store_true", help="row-level comparison even when md5s agree")
a = ap.parse_args()


def md5(fn):
    h = hashlib.md5()
    with open(fn, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def num(s):
    try:
        v = float(s)
        return v if not math.isnan(v) else None
    except ValueError:
        return None


def load(fn):
    with open(fn) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        k = [ix[c] for c in ("CHR", "POS", "MarkerID", "Allele1", "Allele2") if c in ix]
        rows = {}
        dup = 0
        for line in f:
            t = line.rstrip("\n").split("\t")
            key = tuple(t[i] for i in k)
            if key in rows:
                dup += 1
            rows[key] = t
    return hdr, ix, rows, dup


def compare(fa, fb):
    ha, ia, ra, da = load(fa)
    hb, ib, rb, db = load(fb)
    common = [k for k in ra if k in rb]
    s = {"rows_a": len(ra), "rows_b": len(rb), "only_a": len(ra) - len(common), "only_b": len(rb) - len(common),
         "dup_keys_a": da, "dup_keys_b": db, "header_equal": ha == hb, "rows_text_differ": 0}
    for c in ("p.value", "BETA", "SE"):
        s[c] = {"n_differ": 0, "max_abs": 0.0, "max_rel": 0.0, "na_mismatch": 0}
    s["max_abs_dlog10p"] = 0.0
    for thr in ("5e-8", "1e-5"):
        s["n_p_lt_" + thr + "_a"] = s["n_p_lt_" + thr + "_b"] = 0
        s["cross_" + thr + "_a_only"] = s["cross_" + thr + "_b_only"] = 0
    s["is_spa_true_a"] = s["is_spa_true_b"] = s["is_spa_disagree"] = 0
    pa_i, pb_i = ia.get("p.value"), ib.get("p.value")
    for k in common:
        x, y = ra[k], rb[k]
        if x != y:
            s["rows_text_differ"] += 1
        for c in ("p.value", "BETA", "SE"):
            if c not in ia or c not in ib:
                continue
            u, v = x[ia[c]], y[ib[c]]
            if u == v:
                continue
            st = s[c]
            st["n_differ"] += 1
            fu, fv = num(u), num(v)
            if fu is None or fv is None:
                st["na_mismatch"] += (fu is None) != (fv is None)
                continue
            d = abs(fu - fv)
            den = max(abs(fu), abs(fv))
            st["max_abs"] = max(st["max_abs"], d)
            if den > 0:
                st["max_rel"] = max(st["max_rel"], d / den)
        if pa_i is not None and pb_i is not None:
            pa, pb = num(x[pa_i]), num(y[pb_i])
            if pa is not None and pb is not None:
                la, lb = -math.log10(max(pa, 1e-320)), -math.log10(max(pb, 1e-320))
                s["max_abs_dlog10p"] = max(s["max_abs_dlog10p"], abs(la - lb))
                for thr, tv in (("5e-8", 5e-8), ("1e-5", 1e-5)):
                    s["n_p_lt_" + thr + "_a"] += pa < tv
                    s["n_p_lt_" + thr + "_b"] += pb < tv
                    s["cross_" + thr + "_a_only"] += (pa < tv) and not (pb < tv)
                    s["cross_" + thr + "_b_only"] += (pb < tv) and not (pa < tv)
        if "Is.SPA" in ia and "Is.SPA" in ib:
            u, v = x[ia["Is.SPA"]].lower(), y[ib["Is.SPA"]].lower()
            s["is_spa_true_a"] += u == "true"
            s["is_spa_true_b"] += v == "true"
            s["is_spa_disagree"] += u != v
    return s


fa = {os.path.basename(f): f for f in glob.glob(os.path.join(a.a, "*.txt"))}
fb = {os.path.basename(f): f for f in glob.glob(os.path.join(a.b, "*.txt"))}
names = sorted(set(fa) & set(fb))
res = {"label": a.label, "n_files_a": len(fa), "n_files_b": len(fb), "n_paired": len(names),
       "files_only_a": len(set(fa) - set(fb)), "files_only_b": len(set(fb) - set(fa))}
md = {n: md5(fa[n]) == md5(fb[n]) for n in names}
res["n_md5_identical"] = sum(md.values())
res["identical"] = bool(names) and all(md.values()) and res["files_only_a"] == 0 and res["files_only_b"] == 0
if a.rows or not res["identical"]:
    per = {}
    for n in names:
        per[os.path.splitext(n)[0]] = compare(fa[n], fb[n])
    tot = {}
    for s in per.values():
        for k, v in s.items():
            if isinstance(v, dict):
                t = tot.setdefault(k, {})
                for kk, vv in v.items():
                    t[kk] = max(t.get(kk, 0), vv) if kk.startswith("max") else t.get(kk, 0) + vv
            elif isinstance(v, bool):
                tot[k] = tot.get(k, True) and v
            elif k.startswith("max"):
                tot[k] = max(tot.get(k, 0), v)
            else:
                tot[k] = tot.get(k, 0) + v
    res["total"] = tot
    res["per_trait"] = per
json.dump(res, open(a.json, "w"), indent=1, sort_keys=True)
t = res.get("total", {})
print("%s: %d/%d files md5-identical%s" % (
    a.label, res["n_md5_identical"], res["n_paired"],
    "" if not t else "; rows A %d B %d only A %d only B %d, text differs %d; p: n %d max abs %.3g max rel %.3g; "
    "max|dlog10p| %.3g; cross 5e-8 A-only %d B-only %d" % (
        t["rows_a"], t["rows_b"], t["only_a"], t["only_b"], t["rows_text_differ"], t["p.value"]["n_differ"],
        t["p.value"]["max_abs"], t["p.value"]["max_rel"], t["max_abs_dlog10p"],
        t["cross_5e-8_a_only"], t["cross_5e-8_b_only"])))
