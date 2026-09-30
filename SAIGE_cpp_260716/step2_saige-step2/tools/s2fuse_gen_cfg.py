#!/usr/bin/env python3
"""Write one step-2 config for the mtFuseGemm / mtPopcountAF A/B.

usage: s2fuse_gen_cfg.py OUTDIR PLINKPREFIX SPEC...
SPEC is  q:<n>   n quantitative traits from the tg2 P128 step-1 models (y1..yn)
         b:<n>   n binary traits from bin_models/nospa (y1..y8, reused
                 cyclically past 8 under distinct names -- the same model
                 twice costs the loop exactly what two models cost, nothing
                 is deduplicated on the binary side)
Switches come from the environment, each 0/1, default 0:
  FOLD  mtFoldQuantProj    FUSE  mtFuseGemm    PCAF  mtPopcountAF
  SGS=1 adds outputFormat: sgs (fp64 binary output)
  MORE=1 sets isMoreOutput: true (hom/het counts)
"""
import os, sys

QM = "/opt/saige/logs/tg2_step2/runs/s1_cpp_P128/out"
BM = "/opt/saige/logs/tg2_step2/bin_models/nospa"

outdir, plink = sys.argv[1], sys.argv[2]
specs = sys.argv[3:]
os.makedirs(outdir + "/out", exist_ok=True)
E = os.environ.get
onoff = lambda k: "true" if E(k, "0") == "1" else "false"

L = ["genoType: plink",
     "plinkFile: " + plink,
     "minMAF: 0", "minMAC: 1", "maxMissRate: 0.15",
     "AlleleOrder: alt-first", "LOCO: false",
     "isnoadjCov: false", "isMoreOutput: " + onoff("MORE"),
     "isFirth: false", "is_Firth_beta: false",
     "MACCutoffforER: 4", "relatednessCutoff: 0",
     "nThreads: %s" % E("NTHREADS", "8"),
     "mtFoldQuantProj: " + onoff("FOLD"),
     "mtFuseGemm: " + onoff("FUSE"),
     "mtPopcountAF: " + onoff("PCAF"),
     "models:"]
if E("SGS") == "1":
    L.insert(-1, "outputFormat: sgs")
for sp in specs:
    kind, n = sp.split(":"); n = int(n)
    for i in range(1, n + 1):
        if kind == "q":
            y = "y%d" % i
            nm, md, vr = "q_" + y, QM + "/m/" + y, QM + "/mvr_" + y + ".varianceRatio.txt"
        else:
            y = "y%d" % ((i - 1) % 8 + 1)
            nm, md, vr = "b_y%d" % i, BM + "/m/" + y, BM + "/mvr_" + y + ".varianceRatio.txt"
        L += ["  - traitName: " + nm,
              "    modelFile: " + md,
              "    varianceRatioFile: " + vr,
              "    outputFile: " + outdir + "/out/" + nm + ".txt"]
open(outdir + "/cfg.yaml", "w").write("\n".join(L) + "\n")
print(outdir + "/cfg.yaml")
