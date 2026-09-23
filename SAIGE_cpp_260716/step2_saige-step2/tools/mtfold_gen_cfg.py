#!/usr/bin/env python3
"""Write one step-2 config for the mtFoldQuantProj A/B.

usage: gen_cfg.py OUTDIR PLINKPREFIX FOLD(0|1) SPEC...
SPEC is  q:<n>   n quantitative traits from the tg2 P128 step-1 models (y1..yn)
         b:<n>   n binary traits from bin_models/nospa (y1..yn)
Trait names are prefixed q_/b_ so a mixed run has unique names.
"""
import os, sys

QM = "/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out"
BM = "/opt/saige/logs/tg2_step2/bin_models/nospa"

outdir, plink, fold = sys.argv[1], sys.argv[2], sys.argv[3]
specs = sys.argv[4:]
os.makedirs(outdir + "/out", exist_ok=True)

L = ["genoType: plink",
     "plinkFile: " + plink,
     "minMAF: 0", "minMAC: 1", "maxMissRate: 0.15",
     "AlleleOrder: alt-first", "LOCO: false",
     "isnoadjCov: false", "isMoreOutput: false",
     "isFirth: false", "is_Firth_beta: false",
     "MACCutoffforER: 4", "relatednessCutoff: 0",
     "nThreads: %s" % os.environ.get("NTHREADS", "8"),
     "mtFoldQuantProj: " + ("true" if fold == "1" else "false"),
     "models:"]
if os.environ.get("SGS") == "1":
    L.insert(-1, "outputFormat: sgs")
for sp in specs:
    kind, n = sp.split(":"); n = int(n)
    for i in range(1, n + 1):
        y = "y%d" % i
        if kind == "q":
            nm, md, vr = "q_" + y, QM + "/m/" + y, QM + "/mvr_" + y + ".varianceRatio.txt"
        else:
            nm, md, vr = "b_" + y, BM + "/m/" + y, BM + "/mvr_" + y + ".varianceRatio.txt"
        L += ["  - traitName: " + nm,
              "    modelFile: " + md,
              "    varianceRatioFile: " + vr,
              "    outputFile: " + outdir + "/out/" + nm + ".txt"]
open(outdir + "/cfg.yaml", "w").write("\n".join(L) + "\n")
print(outdir + "/cfg.yaml")
