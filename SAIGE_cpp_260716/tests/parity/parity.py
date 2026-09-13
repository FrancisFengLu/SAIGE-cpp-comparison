#!/usr/bin/env python3
"""R-SAIGE 1.5.2 vs standalone C++ Step-2 parity harness.

One case spec -> a semantically equivalent pair of runs -> a column-wise diff.

    ./parity.py --list                       # cases and datasets
    ./parity.py --knobs                      # the path-space inventory
    ./parity.py base_binary_single           # run a registered case
    ./parity.py base_binary_single --set minMAF=0.2
    ./parity.py base_binary_single --only r  # re-run one side only
    ./parity.py base_binary_single --skip-run # just re-compare existing output

WHY THE NULL MODEL COMES FROM R
-------------------------------
Both sides read the SAME step-1 fit: R reads the .rda, C++ reads an .arma
directory converted from that identical .rda by tools/rda_to_arma.R.  Running
step 1 twice would fold step-1's own divergence (tau is a stochastic AI-REML
estimate; PATH_TESTS/RESULTS_202609.md measures run-to-run spread of ~4e-6 on
tau alone) into every step-2 p-value, and a step-2 bug would be indistinguish-
able from step-1 Monte-Carlo noise.  So: convert, never re-fit.

The corollary is that this harness proves nothing about step 1.  Step-1 parity
is a separate exercise with its own tolerance story.

LOCO
----
tools/rda_to_arma.R refuses LOCO models, so the converted-model trick does not
reach the LOCO path, and NO LOCO CASE IS REGISTERED HERE.  What is missing is
only the converter: R's .rda carries modglmm$LOCOResult[[j]]$however
{fitted.values, residuals, obj.noK, offset}, which is exactly the per-chromosome
file set LOCO_FORMAT.md asks for (mu, res, V, offset, XV, XVX, XVX_inv,
XVX_inv_XV, XXVX_inv, S_a).  Until that converter exists, a "LOCO comparison"
would have to fit step 1 twice, which is not a step-2 comparison at all.
Do not fake it.
"""
import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import knobs as K                      # noqa: E402
import compare_out                     # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CPP_BIN = os.path.join(ROOT, "step2_saige-step2", "saige-step2")
RDA2ARMA = os.path.join(ROOT, "step2_saige-step2", "tools", "rda_to_arma.R")
R_STEP2 = "/opt/saige/SAIGE-upstream/extdata/step2_SPAtests.R"
R_LIBS = "/opt/saige/Rlib-upstream"
CONDA_SH = "/home/francisfenglu4/miniforge3/etc/profile.d/conda.sh"
WORKROOT = os.environ.get("PARITY_WORK", "/opt/saige/logs/parity")

# ---------------------------------------------------------------------------
# datasets and null models
# ---------------------------------------------------------------------------
DATASETS = {
    # 50k samples, 2000 markers -- the fast one, same samples as `mid`
    "mid2k": dict(kind="plink", prefix="/opt/saige/data/mid2k"),
    # 50k samples, 40k markers
    "mid":   dict(kind="plink", prefix="/opt/saige/data/mid"),
}

MODELS = {
    # R 1.5.x step-1 fit on mid + mid.mp8.pheno.txt column y4 (binary),
    # non-LOCO, single variance ratio, n=50000 p=3 theta=[1, 0.4566934].
    # `arma` was produced from `rda` by tools/rda_to_arma.R -- same fit.
    "y4_binary": dict(
        trait="binary",
        rda="/opt/saige/logs/mp/single_y4.rda",
        arma="/opt/saige/logs/step2/null_y4_arma",
        vr="/opt/saige/logs/mp/single_y4.varianceRatio.txt"),
}

# ---------------------------------------------------------------------------
# Baseline: every knob that the two sides default differently is PINNED here,
# explicitly, on both sides.  Nothing in a parity run may rely on a default.
# ---------------------------------------------------------------------------
BASE = dict(
    # QC -- pinned so a filter difference can never be mistaken for a maths one
    minMAF=0.0,
    minMAC=0.5,
    maxMissing=0.15,
    minInfo=0.0,
    is_imputed_data=False,
    dosage_zerod_cutoff=0.2,
    dosage_zerod_MAC_cutoff=10.0,
    AlleleOrder="alt-first",
    impute_method="mean",          # R default is best_guess; cpp json is mean
    # model path
    LOCO=False,                    # R default TRUE, cpp default FALSE
    SPAcutoff=2.0,
    max_MAC_for_ER=4.0,
    is_noadjCov=False,             # R CLI default TRUE -- see knobs.py note (b)
    is_fastTest=False,             # R default FALSE, cpp json default TRUE
    pval_cutoff_for_fastTest=0.05,
    is_Firth_beta=False,
    pCutoffforFirth=0.01,
    # output / perf
    is_output_moreDetails=False,
    markers_per_chunk=10000,
    nThreads=1,                    # region row order is nondeterministic above 1
)

REGION_BASE = dict(
    r_corr=0.0,
    maxMAF_in_groupTest=[0.0001, 0.001, 0.01],
    MACCutoff_to_CollapseUltraRare=10.0,
    markers_per_chunk_in_groupTest=100,   # R default; cpp default is 500
    groups_per_chunk=100,
    is_single_in_groupTest=True,          # forced TRUE anyway when r_corr=0
    is_output_markerList_in_groupTest=False,
    **{"weights.beta": [1, 25]},
)

CASES = {
    "base_binary_single": dict(
        desc="binary, single-variant, non-LOCO, is_noadjCov=FALSE, "
             "impute_method=mean -- the configuration that should agree",
        data="mid2k", model="y4_binary", set={}),
    "binary_single_minmaf": dict(
        desc="same as base but minMAF=0.2 on BOTH sides -- a real filter change, "
             "still expected to agree (used to prove the harness is not blind)",
        data="mid2k", model="y4_binary", set=dict(minMAF=0.2)),
    "binary_single_firth": dict(
        desc="Firth effect sizes on, with a loose p cutoff so the branch is "
             "actually entered",
        data="mid2k", model="y4_binary",
        set=dict(is_Firth_beta=True, pCutoffforFirth=0.5)),
    "binary_single_moredetails": dict(
        desc="is_output_moreDetails -- extra columns must match too",
        data="mid2k", model="y4_binary", set=dict(is_output_moreDetails=True)),
    "binary_single_noadjcov": dict(
        desc="is_noadjCov=TRUE on BOTH sides. MEASURED 2026-09-13: bit-identical "
             "over 2000 markers, so the C++ port REPRODUCES R's noadjCov "
             "arithmetic exactly -- including the AF>0.5 defect (R 1.5.x centers "
             "scoreTestFast_noadjCov on 2*altFreq after the genotype has already "
             "been flipped, UTIL.cpp:75-122). The C++/R divergence here is only "
             "the DEFAULT (R CLI TRUE, cpp FALSE), not the maths. Against the "
             "same run with is_noadjCov=FALSE the branch moves p by up to 4 "
             "orders of magnitude, so this case really does enter it.",
        data="mid2k", model="y4_binary", set=dict(is_noadjCov=True),
        vr_add_noXadj=True),
    "binary_single_spa1": dict(
        desc="SPAcutoff=1 -- pushes many more markers onto the SPA branch",
        data="mid2k", model="y4_binary", set=dict(SPAcutoff=1.0)),
    "binary_single_bestguess": dict(
        desc="impute_method=best_guess on both sides (R's own default)",
        data="mid2k", model="y4_binary", set=dict(impute_method="best_guess")),
}

BOOLS = {True: "TRUE", False: "FALSE"}


# ---------------------------------------------------------------------------
# knob -> command line / YAML
# ---------------------------------------------------------------------------
def r_value(name, v):
    if isinstance(v, bool):
        return BOOLS[v]
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)


def make_vr(model, wd, add_noXadj):
    """Return the variance-ratio file both sides will read.

    `add_noXadj` synthesises a `null_noXadj` row by copying the `null` value.
    That value is NOT the real no-covariate-adjustment variance ratio -- step 1
    only emits one with --skipModelFitting -- so an is_noadjCov run built this
    way is a CODE-PATH comparison (same inputs into both implementations), not
    a statistically meaningful analysis.  It is still the only way to reach the
    branch at all: R 1.5.2 aborts with
        Mat::init(): requested size is not compatible with row vector layout
    when is_noadjCov is on and the VR file carries no null_noXadj row
    (readInGLMM.R:418-427 lets a zero-length vector through unguarded).
    """
    src = MODELS[model]["vr"]
    if not add_noXadj:
        return src
    dst = os.path.join(wd, "vr_with_noXadj.txt")
    rows = [l.rstrip("\n") for l in open(src) if l.strip()]
    if any(r.split()[1] == "null_noXadj" for r in rows if len(r.split()) > 1):
        shutil.copyfile(src, dst)
        return dst
    # Normalise every row to tabs. data.table::fread sniffs one separator for
    # the whole file, so appending a tab-separated row to a space-separated
    # file silently yields ncol==1 and R then reports an EMPTY ratioVec_null
    # and dies later inside mainMarkerInCPP with "index out of bounds".
    out = ["\t".join(r.split()) for r in rows]
    for r in rows:
        f = r.split()
        if len(f) >= 2 and f[1] == "null":
            out.append("%s\tnull_noXadj\t%s" % (f[0], f[2] if len(f) > 2 else "1"))
    with open(dst, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return dst


def build_r_args(vals, data, model, out_path, vr):
    args = []
    d = DATASETS[data]
    if d["kind"] == "plink":
        args += ["--bedFile=%s.bed" % d["prefix"],
                 "--bimFile=%s.bim" % d["prefix"],
                 "--famFile=%s.fam" % d["prefix"]]
    else:
        raise SystemExit("dataset kind %r has no R mapping yet" % d["kind"])
    args += ["--GMMATmodelFile=%s" % MODELS[model]["rda"],
             "--varianceRatioFile=%s" % vr,
             "--SAIGEOutputFile=%s" % out_path]
    for name, v in sorted(vals.items()):
        spec = K.BY_NAME.get(name)
        if spec is None:
            raise SystemExit("unknown knob %r" % name)
        if spec["r"] is None:
            continue                       # no R CLI flag (e.g. pval_cutoff_for_fastTest)
        if v is None or v == "":
            continue
        args.append("%s=%s" % (spec["r"], r_value(name, v)))
    return args


# canonical name -> C++ YAML key, for knobs whose cpp slot is yaml/yaml+json
def build_cpp_yaml(vals, data, model, out_path, model_dir, vr):
    d = DATASETS[data]
    y = {}
    lines = ["# generated by tests/parity/parity.py -- do not edit"]
    lines.append("modelFile:         %s" % model_dir)
    lines.append("varianceRatioFile: %s" % vr)
    lines.append("outputFile:        %s" % out_path)
    if d["kind"] == "plink":
        lines.append("genoType:   plink")
        lines.append("plinkFile:  %s" % d["prefix"])
    else:
        raise SystemExit("dataset kind %r has no C++ mapping yet" % d["kind"])
    for name, v in sorted(vals.items()):
        spec = K.BY_NAME.get(name)
        if spec is None:
            raise SystemExit("unknown knob %r" % name)
        if spec["cpp"] is None:
            continue                      # not ported -- reported by check_case()
        where, key = spec["cpp"]
        if where == "json":
            continue                      # handled in patch_model()
        if v is None or v == "":
            continue
        if isinstance(v, bool):
            lines.append("%s: %s" % (key, "true" if v else "false"))
        elif isinstance(v, (list, tuple)):
            lines.append("%s:" % key)
            for x in v:
                lines.append("  - %s" % (('"%s"' % x) if isinstance(x, str) else x))
        elif isinstance(v, str):
            lines.append('%s: "%s"' % (key, v))
        else:
            lines.append("%s: %s" % (key, v))
        y[key] = v
    return "\n".join(lines) + "\n", y


def patch_model(vals, model, dst):
    """Materialise the C++ model dir with the json-resident knobs patched in.

    .arma files are symlinked (they are the identical fit); only
    nullmodel.json is rewritten.
    """
    src = MODELS[model]["arma"]
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(dst)
    for fn in os.listdir(src):
        if fn == "nullmodel.json":
            continue
        os.symlink(os.path.join(src, fn), os.path.join(dst, fn))
    with open(os.path.join(src, "nullmodel.json")) as fh:
        j = json.load(fh)
    patched = {}
    for name, v in vals.items():
        spec = K.BY_NAME.get(name)
        if spec is None or spec["cpp"] is None:
            continue
        where, key = spec["cpp"]
        if where in ("json", "yaml+json"):
            j[key] = v
            patched[key] = v
    with open(os.path.join(dst, "nullmodel.json"), "w") as fh:
        json.dump(j, fh, indent=1)
    return patched


# ---------------------------------------------------------------------------
def sh(cmd, log_path, env_setup):
    script = "set -o pipefail\n" + env_setup + "\nexec " + cmd + "\n"
    with open(log_path, "w") as lf:
        t0 = time.time()
        p = subprocess.run(["bash", "-c", script], stdout=lf,
                           stderr=subprocess.STDOUT)
        dt = time.time() - t0
    return p.returncode, dt


R_ENV = ("source %s\nconda activate RSAIGE_GPU\n"
         "export R_LIBS_USER=%s\n"
         "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}\n"
         % (CONDA_SH, R_LIBS))
CPP_ENV = ("source %s\nconda activate saige-build\n"
           "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}\n"
           % CONDA_SH)


def check_case(vals):
    """Refuse to pretend a knob is comparable when one side cannot express it."""
    problems = []
    for name, v in vals.items():
        spec = K.BY_NAME.get(name)
        if spec is None:
            problems.append("unknown knob: %s" % name)
            continue
        default_like = (v == spec["r_def"]) or (v in (False, 0, 0.0, "", None, []))
        if spec["cpp"] is None and not default_like:
            problems.append(
                "%s=%r is NOT PORTED to C++ (%s) -- the two runs would not be "
                "the same experiment" % (name, v, spec["note"]))
        if spec["r"] is None and not default_like and name != "pval_cutoff_for_fastTest":
            problems.append("%s=%r has no R CLI flag" % (name, v))
    return problems


def run_case(case, overrides, only, skip_run, rtol, work, desync=None):
    spec = CASES[case]
    vals = dict(BASE)
    if "groupFile" in spec.get("set", {}) or spec.get("region"):
        vals.update(REGION_BASE)
    vals.update(spec.get("set", {}))
    vals.update(overrides)
    # --desync deliberately breaks the semantic equivalence on ONE side. It is
    # the harness's own negative control: if a one-sided knob change does not
    # show up in the report, the report is worthless.
    desync = desync or {}
    r_vals = dict(vals, **desync.get("r", {}))
    c_vals = dict(vals, **desync.get("cpp", {}))

    problems = check_case(vals)
    if problems:
        print("REFUSING to run %s:" % case)
        for p in problems:
            print("  - " + p)
        return 2

    wd = os.path.join(work, case)
    os.makedirs(wd, exist_ok=True)
    r_out = os.path.join(wd, "r.txt")
    c_out = os.path.join(wd, "cpp.txt")
    model_dir = os.path.join(wd, "model")

    print("=== case %s" % case)
    print("    %s" % spec["desc"])
    print("    data=%s model=%s (trait=%s)"
          % (spec["data"], spec["model"], MODELS[spec["model"]]["trait"]))
    nd = {k: v for k, v in vals.items() if v != BASE.get(k, object())}
    if nd:
        print("    non-baseline knobs: %s"
              % ", ".join("%s=%r" % kv for kv in sorted(nd.items())))
    for side, dd in sorted(desync.items()):
        print("    !! DESYNCED on the %s side only: %s"
              % (side, ", ".join("%s=%r" % kv for kv in sorted(dd.items()))))
    print("    work: %s" % wd)

    vr = make_vr(spec["model"], wd, spec.get("vr_add_noXadj", False))
    if vr != MODELS[spec["model"]]["vr"]:
        print("    synthesised variance-ratio file: %s (see make_vr docstring)" % vr)
    r_args = build_r_args(r_vals, spec["data"], spec["model"], r_out, vr)
    patched = patch_model(c_vals, spec["model"], model_dir)
    yaml_txt, _ = build_cpp_yaml(c_vals, spec["data"], spec["model"], c_out, model_dir, vr)
    yaml_path = os.path.join(wd, "cpp.yaml")
    with open(yaml_path, "w") as fh:
        fh.write(yaml_txt)
    with open(os.path.join(wd, "r.cmd"), "w") as fh:
        fh.write("Rscript %s \\\n  %s\n" % (R_STEP2, " \\\n  ".join(r_args)))
    if patched:
        print("    nullmodel.json patched: %s"
              % ", ".join("%s=%r" % kv for kv in sorted(patched.items())))

    if not skip_run:
        if only in (None, "r"):
            cmd = "Rscript %s %s" % (R_STEP2, " ".join(shlex.quote(a) for a in r_args))
            rc, dt = sh(cmd, os.path.join(wd, "r.log"), R_ENV)
            print("    R   : rc=%d  %.1fs  (%s)" % (rc, dt, os.path.join(wd, "r.log")))
            if rc != 0:
                print("    R FAILED -- tail of log:")
                print(tail(os.path.join(wd, "r.log")))
                return 3
        if only in (None, "cpp"):
            if not os.path.exists(CPP_BIN):
                print("    cpp binary missing: %s" % CPP_BIN)
                return 3
            cmd = "%s %s" % (shlex.quote(CPP_BIN), shlex.quote(yaml_path))
            rc, dt = sh(cmd, os.path.join(wd, "cpp.log"), CPP_ENV)
            print("    cpp : rc=%d  %.1fs  (%s)" % (rc, dt, os.path.join(wd, "cpp.log")))
            if rc != 0:
                print("    cpp FAILED -- tail of log:")
                print(tail(os.path.join(wd, "cpp.log")))
                return 3

    if only is not None and not skip_run:
        return 0

    rcs = []
    pairs = [(r_out, c_out, "main")]
    for sfx in (".singleAssoc.txt", ".markerList.txt"):
        if os.path.exists(r_out + sfx) and os.path.exists(c_out + sfx):
            pairs.append((r_out + sfx, c_out + sfx, sfx.strip(".")))
    for rp, cp, label in pairs:
        if not (os.path.exists(rp) and os.path.exists(cp)):
            print("    MISSING OUTPUT: %s / %s" % (rp, cp))
            rcs.append(1)
            continue
        print("\n--- %s ---" % label)
        rep = compare_out.compare(rp, cp, rtol)
        compare_out.print_report(rep)
        with open(os.path.join(wd, "diff_%s.json" % label), "w") as fh:
            json.dump(rep, fh, indent=2)
        ok = rep["verdict"] == "IDENTICAL"
        if desync:
            print("    (desync is active; a verdict of IDENTICAL would mean the "
                  "harness is blind)")
            ok = not ok
        elif spec.get("expect") == "differ":
            print("    (case is declared expect=differ; agreement would be the surprise)")
            ok = not ok
        rcs.append(0 if ok else 1)
    return max(rcs) if rcs else 1


def tail(path, n=25):
    try:
        with open(path) as fh:
            return "".join(fh.readlines()[-n:])
    except OSError:
        return "(no log)"


def parse_set(s):
    k, _, v = s.partition("=")
    if v.upper() in ("TRUE", "FALSE"):
        return k, v.upper() == "TRUE"
    if "," in v:
        return k, [float(x) if _isnum(x) else x for x in v.split(",")]
    if _isnum(v):
        return k, float(v) if ("." in v or "e" in v.lower()) else int(v)
    return k, v


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case", nargs="?")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--knobs", action="store_true")
    ap.add_argument("--set", action="append", default=[],
                    help="knob=value override applied to BOTH sides, repeatable")
    ap.add_argument("--desync", action="append", default=[],
                    help="side:knob=value applied to ONE side only "
                         "(side = r | cpp). Negative control: the report must "
                         "flag it. Inverts the exit status.")
    ap.add_argument("--only", choices=["r", "cpp"])
    ap.add_argument("--skip-run", action="store_true",
                    help="re-compare existing outputs without re-running")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--work", default=WORKROOT)
    ap.add_argument("--all", action="store_true", help="run every registered case")
    a = ap.parse_args()

    if a.knobs:
        K.print_table()
        return 0
    if a.list or (not a.case and not a.all):
        print("cases:")
        for n, c in CASES.items():
            print("  %-28s %s" % (n, c["desc"].split("--")[0].strip()))
            if c.get("expect") == "differ":
                print("  %-28s   [expect=differ]" % "")
        print("\ndatasets: %s" % ", ".join(DATASETS))
        print("models:   %s" % ", ".join(MODELS))
        print("\nno LOCO case is registered -- see the module docstring.")
        return 0

    ov = dict(parse_set(s) for s in a.set)
    desync = {}
    for s in a.desync:
        side, _, rest = s.partition(":")
        if side not in ("r", "cpp") or not rest:
            print("--desync takes side:knob=value with side in {r, cpp}")
            return 2
        k, v = parse_set(rest)
        desync.setdefault(side, {})[k] = v
    cases = list(CASES) if a.all else [a.case]
    rc = 0
    for c in cases:
        if c not in CASES:
            print("unknown case %r (see --list)" % c)
            return 2
        rc = max(rc, run_case(c, ov, a.only, a.skip_run, a.rtol, a.work, desync))
        print()
    return rc


if __name__ == "__main__":
    sys.exit(main())
