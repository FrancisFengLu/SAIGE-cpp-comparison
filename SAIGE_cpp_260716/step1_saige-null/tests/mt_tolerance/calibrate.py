#!/usr/bin/env python3
"""Aggregate observed deviations from gate runs, per threshold rule (for THRESHOLDS.md).

  calibrate.py [--workdir /opt/saige/logs/mt_gate] LABEL [LABEL ...] [--markdown]

Reads WORKDIR/runs/<LABEL>/<case>/compare_<trait>.json and solo_vs_solo_<trait>.json and
prints, for every rule key, the largest deviation seen (with case/trait/file/index), the
current tolerance and the ratio; plus every exact item that differed (iterations,
converged, file sets, ...), which a calibration run must not have unless explained.
"""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("labels", nargs="+")
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--markdown", action="store_true")
    ap.add_argument("--recompute", action="store_true",
                    help="re-run compare.py on the stored run directories first (use after editing "
                         "thresholds.yaml or the metric definitions; no saige-null run is repeated)")
    a = ap.parse_args()
    if a.recompute:
        for lab in a.labels:
            for p in sorted(glob.glob(os.path.join(a.workdir, "runs", lab, "*", "*.json"))):
                b = os.path.basename(p)
                if not (b.startswith("compare_") or b.startswith("solo_vs_solo_")):
                    continue
                r = json.load(open(p))
                if not (os.path.isdir(r["ref"]) and os.path.isdir(r["test"])):
                    print(f"  skip {p}: run directory gone")
                    continue
                layout = "multi" if b.startswith("compare_") else "solo"
                res = compare.compare_trait(r["ref"], r["test"], r.get("trait"), "solo", layout)
                json.dump(res, open(p, "w"), indent=1, default=str)
            print(f"recomputed {lab}")
    worst = {}
    per_label = defaultdict(dict)
    exact_bad = []
    n_files = 0
    for lab in a.labels:
        for p in sorted(glob.glob(os.path.join(a.workdir, "runs", lab, "*", "*.json"))):
            base = os.path.basename(p)
            if not (base.startswith("compare_") or base.startswith("solo_vs_solo_")):
                continue
            kind = "solo-vs-solo" if base.startswith("solo_vs_solo_") else "multi-vs-solo"
            case = os.path.basename(os.path.dirname(p))
            trait = base.split("_", 1)[1][:-5] if kind == "multi-vs-solo" else base[len("solo_vs_solo_"):-5]
            r = json.load(open(p))
            n_files += 1
            for it in r["items"]:
                key = it["rule"]
                if it.get("metric") in ("rel", "log10p"):
                    dev = it["dev"] if it["dev"] is not None else 0.0
                    rec = dict(dev=dev, tol=it["tol"], label=lab, kind=kind, case=case, trait=trait,
                               file=it["file"], item=it["item"], at=it.get("worst_index"),
                               ref=it.get("ref_value"), test=it.get("test_value"))
                    for store in (worst, per_label[lab]):
                        if key not in store or dev > store[key]["dev"]:
                            store[key] = rec
                elif not it["passed"]:
                    exact_bad.append(dict(label=lab, kind=kind, case=case, trait=trait, file=it["file"],
                                          item=it["item"], note=it.get("note")))
    print(f"{n_files} comparison reports from {a.labels}")
    fmt = "| {:<18} | {:>9} | {:>7} | {:>6} | {} |" if a.markdown else "{:<18} {:>9} {:>7} {:>6}  {}"
    print(fmt.format("rule", "max dev", "tol", "tol/dev", "where"))
    if a.markdown:
        print("|---|---:|---:|---:|---|")
    for key in sorted(worst):
        w = worst[key]
        ratio = (w["tol"] / w["dev"]) if w["dev"] else float("inf")
        where = f"{w['label']} {w['kind']} {w['case']}/{w['trait']} {w['file']}::{w['item']}@{w['at']}"
        print(fmt.format(key, f"{w['dev']:.2e}", f"{w['tol']:.0e}", f"{ratio:.1f}" if ratio != float('inf') else "inf", where))
    if len(a.labels) > 1:
        print("\nper label:")
        for lab in a.labels:
            print(f"  {lab}: " + ", ".join(f"{k}={v['dev']:.1e}" for k, v in sorted(per_label[lab].items())))
    print("\ntimings (seconds; solo = sum over the case's traits as first run, multi = one P>1 run):")
    for lab in a.labels:
        sj = os.path.join(a.workdir, "runs", lab, "summary.json")
        if not os.path.exists(sj):
            continue
        S = json.load(open(sj))
        print(f"  {lab}: " + ", ".join(f"{c['case']} {c['solo_secs']:.0f}+{c['multi_secs']:.0f}" for c in S["cases"]))
    print(f"\nexact items that differed: {len(exact_bad)}")
    for e in exact_bad[:60]:
        print(f"  {e['label']} {e['kind']} {e['case']}/{e['trait']} {e['file']}::{e['item']}  {str(e['note'])[:160]}")


if __name__ == "__main__":
    main()
