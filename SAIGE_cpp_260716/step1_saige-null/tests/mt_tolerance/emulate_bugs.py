#!/usr/bin/env python3
"""Size of real bugs, for THRESHOLDS.md: emulate a scheme-C bug at the DATA level with the
current binary (no C++ change), run the affected trait alone, and compare it with its
correct solo run using compare.py. The deviations are what the gate would see if the bug
were in the multi-trait path.

  emulate_bugs.py BIN [--gpu] [--small] [--workdir /opt/saige/logs/mt_gate] [--only e1,e2,...]

  e1_grm_drop   trait qA3 (qc case) loses from its GRM the markers that pass QC on its own
                samples but fail on the union (what "use the union's QC list" does in that
                direction). Emulated by making those markers all-missing in a BED copy.
  e2_fill       trait fS (fill case) gets the union's fill value at every marker where it
                differs from its own. Emulated by writing the union fill into the BED copy for
                fS's missing genotypes there.
  e3_one_more   trait t10 (tiny case) fitted with one extra sample (off-by-one sample set).
  e4_one_less   trait lS (loco case, LOCO on) fitted with one sample fewer.
"""
import argparse
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402
import gate  # noqa: E402
import genoqc as gq  # noqa: E402
import make_cases as mc  # noqa: E402

CODE = {0: 0b11, 1: 0b10, 2: 0b00, 3: 0b01}   # A1 count -> PLINK 2-bit code (3 = missing)


def write_bed_variant(src_prefix, dst_prefix, edits):
    """edits: list of (marker j, sample indices, A1 count or 3). Symlinks .bim/.fam."""
    N, M = gq.bed_dims(src_prefix)
    nb = (N + 3) // 4
    raw = np.fromfile(src_prefix + ".bed", dtype=np.uint8)
    body = raw[3:].reshape(M, nb)
    for j, idx, val in edits:
        idx = np.asarray(idx, np.int64)
        byte, sh = idx // 4, (idx % 4) * 2
        row = body[j]
        row[byte] = (row[byte] & ~(np.uint8(3) << sh).astype(np.uint8)) | (np.uint8(CODE[val]) << sh).astype(np.uint8)
    raw.tofile(dst_prefix + ".bed")
    for ext in (".bim", ".fam"):
        if os.path.lexists(dst_prefix + ext):
            os.remove(dst_prefix + ext)
        os.symlink(os.path.realpath(src_prefix + ext), dst_prefix + ext)


def solo_ref(workdir, scale, case_name, trait, bin_md5, gpu, nthreads):
    case = gate.load_case(case_name)
    man = gate.load_manifest(workdir, scale, case["data_case"])
    cfg = gate.render_config(case, man, f"/opt/saige/data/{scale}", nthreads, [trait], "@OUT@", False)
    key = gate.solo_cache_key(bin_md5, cfg, gpu, man)
    d = os.path.join(workdir, "ref_cache", scale, case_name, f"{trait}-{key}")
    if not os.path.exists(os.path.join(d, "DONE")):
        raise SystemExit(f"no cached solo reference for {case_name}/{trait} ({d}); run run_gate.sh first")
    return case, man, os.path.join(d, "run")


def run_variant(binary, case, man, scale, trait, rundir, gpu, nthreads, plink=None, pheno=None):
    cfg = gate.render_config(case, man, f"/opt/saige/data/{scale}", nthreads, [trait], rundir, False)
    if plink:
        cfg["paths"]["plinkFile"] = plink
    if pheno:
        cfg["design"]["csv"] = pheno
    rc, secs = gate.run_saige(binary, cfg, rundir, gpu)
    return rc, secs


def summarize(res):
    out = {}
    for it in res["items"]:
        if it.get("metric") in ("rel", "log10p") or not it["passed"]:
            k = f"{it['file']}::{it['item']}"
            out[k] = dict(dev=it.get("dev"), tol=it.get("tol"), passed=it["passed"], note=it.get("note"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("binary")
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--small", action="store_true")
    ap.add_argument("--nthreads", type=int, default=8)
    ap.add_argument("--workdir", default="/opt/saige/logs/mt_gate")
    ap.add_argument("--only", default="e1_grm_drop,e2_fill,e3_one_more,e4_one_less")
    a = ap.parse_args()
    binary = os.path.abspath(a.binary)
    bmd5 = gate.md5_file(binary)
    scale = "small" if a.small else "mid"
    out_root = os.path.join(a.workdir, "emulate", f"{scale}-{'gpu' if a.gpu else 'cpu'}-{bmd5[:8]}")
    os.makedirs(out_root, exist_ok=True)
    only = set(a.only.split(","))
    ctx = mc.Ctx(scale, a.workdir)
    report = {}

    def finish(name, ref, rundir, trait, desc, extra):
        res = compare.compare_trait(ref, rundir, None, "solo", "solo", thresholds=compare.DEFAULT_THRESHOLDS)
        json.dump(res, open(os.path.join(out_root, f"{name}.compare.json"), "w"), indent=1, default=str)
        s = summarize(res)
        report[name] = dict(trait=trait, description=desc, gate_verdict="PASS" if res["passed"] else "FAIL",
                            n_failed=res["n_failed"], items=s, **extra)
        print(f"== {name}: {desc}\n   gate verdict {report[name]['gate_verdict']} ({res['n_failed']} items failed)")
        for k, v in s.items():
            print(f"   {k:<60} dev={compare._fmt(v['dev'])} tol={compare._fmt(v['tol'])} {'ok' if v['passed'] else 'FAIL'}")

    if "e1_grm_drop" in only:
        case, man, ref = solo_ref(a.workdir, scale, "qc", "qA3", bmd5, a.gpu, a.nthreads)
        A3 = mc.nonmissing(mc.read_tsv(man["pheno"])[1]["qA3"])
        tr = man["traits"]
        U = np.unique(np.concatenate([mc.nonmissing(mc.read_tsv(man["pheno"])[1][t]) for t in tr]))
        sA, sU = ctx.qc(A3), ctx.qc(U)
        drop = np.nonzero(sA["passQC"] & ~sU["passQC"])[0]
        dst = os.path.join(out_root, "e1_bed")
        write_bed_variant(ctx.cfg["plink"], dst, [(int(j), np.arange(ctx.N), 3) for j in drop])
        rd = os.path.join(out_root, "e1_run")
        rc, _ = run_variant(binary, case, man, scale, "qA3", rd, a.gpu, a.nthreads, plink=dst)
        finish("e1_grm_drop", ref, rd, "qA3",
               f"qA3 GRM without the {len(drop)} markers that pass QC on qA3 but fail on the union "
               f"(of {int(sA['passQC'].sum())} GRM markers)", dict(rc=rc, markers_dropped=int(len(drop))))
        os.remove(dst + ".bed")

    if "e2_fill" in only:
        case, man, ref = solo_ref(a.workdir, scale, "fill", "fS", bmd5, a.gpu, a.nthreads)
        S = mc.nonmissing(mc.read_tsv(man["pheno"])[1]["fS"])
        sS, sU = ctx.qc(S), ctx.full
        flip = np.nonzero((sS["fill"] != sU["fill"]) & sS["passQC"] & (sS["miss"] > 0))[0]
        edits, cells = [], 0
        G = ctx.G
        for j in flip:
            miss_S = S[G[j, S] == 3]
            cells += len(miss_S)
            edits.append((int(j), miss_S, int(sU["fill"][j])))
        dst = os.path.join(out_root, "e2_bed")
        write_bed_variant(ctx.cfg["plink"], dst, edits)
        rd = os.path.join(out_root, "e2_run")
        rc, _ = run_variant(binary, case, man, scale, "fS", rd, a.gpu, a.nthreads, plink=dst)
        finish("e2_fill", ref, rd, "fS",
               f"fS with the union fill value at {len(flip)} markers ({cells} missing genotypes)",
               dict(rc=rc, markers=int(len(flip)), cells=int(cells)))
        os.remove(dst + ".bed")

    for name, case_name, trait, delta in (("e3_one_more", "tiny", "t10", +1), ("e4_one_less", "loco", "lS", -1)):
        if name not in only:
            continue
        case, man, ref = solo_ref(a.workdir, scale, case_name, trait, bmd5, a.gpu, a.nthreads)
        hdr, cols = mc.read_tsv(man["pheno"])
        vals = cols[trait]
        if delta > 0:
            i = next(k for k, v in enumerate(vals) if v == "NA")
            src = cols["tF"] if "tF" in cols else cols[hdr[1]]
            vals[i] = src[i] if src[i] != "NA" else "0"
        else:
            i = next(k for k, v in enumerate(vals) if v != "NA")
            vals[i] = "NA"
        cols[trait] = vals
        p = os.path.join(out_root, f"{name}.pheno.txt")
        with open(p, "w") as f:
            f.write("\t".join(hdr) + "\n")
            for r in range(len(cols["IID"])):
                f.write("\t".join(cols[h][r] for h in hdr) + "\n")
        rd = os.path.join(out_root, f"{name}_run")
        rc, _ = run_variant(binary, case, man, scale, trait, rd, a.gpu, a.nthreads, plink=(man.get("loco_plink") if case["plink"] == "loco" else None), pheno=p)
        finish(name, ref, rd, trait, f"{trait} with one sample {'more' if delta > 0 else 'fewer'} (FAM row {i})",
               dict(rc=rc))
    json.dump(report, open(os.path.join(out_root, "report.json"), "w"), indent=1, default=str)
    print(f"report: {os.path.join(out_root, 'report.json')}")


if __name__ == "__main__":
    sys.exit(main())
