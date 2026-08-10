#!/usr/bin/env python3
"""Compare a C++ LOCO run against the R SAIGE reference.

    Rscript tests/loco/r_reference_loco.R <plink> <pheno> <Rout> y_binary binary
    ./saige-null -c <cfg with loco: true>
    python3 tests/loco/compare_with_r.py <cpp_model_dir> <Rout>_LOCO_mu.csv

Reports, per chromosome:
  * max|mu_cpp_chr - mu_R_chr|                (absolute agreement)
  * max|mu_chr - mu_full| on each side        (LOCO effect size)
  * correlation of the two LOCO delta vectors (does the exclusion move the same
    samples in the same direction?)

Reference results recorded 2026-08-10 (gcc 12.4, nthreads=1, trace_seed=10,
covariate_offset on, x1+x2, y_binary):

  A) plinkforGRM_1000samples_10kMarkers (2 autosomes: 9638 + 12 QC'd markers)
     tau[1]: C++ 0.245284 vs R 0.239694 (2.3% — the documented trace-RNG /
     marker-QC gap; the R package installed here writes its random-vector bypass
     to a hard-coded macOS path, so the usual bypass cannot be used on Linux).
       chr2 (small exclusion): LOCO delta 6.19e-3 (C++) vs 6.78e-3 (R),
                               corr(delta) = 0.992
       chr1 (removes 9638/9650 markers, near-degenerate GRM):
                               1.638e-1 vs 1.723e-1, corr(delta) = 0.904

  B) nfam_100_..._poly_22chr (22 balanced autosomes, ~5020 QC'd markers each)
     tau[1]: C++ 0.261039 vs R 0.258116 (1.1%)
       max|mu_cpp_chr - mu_R_chr| ~ 8.5e-4 for every chromosome, whereas
       max|mu_cpp_full - mu_R_full| = 3.4e-3. The per-chromosome fits agree with
       R about 4x BETTER than the full-genome fits do.
       This is diagnostic, not luck: the C++ binary_glmm_solver returns the eta
       from the last inner IRLS at the PREVIOUS tau and never re-solves
       Get_Coef at the final tau, while R does exactly that
       (SAIGE_fitGLMM_fast.R:967). LOCO re-solves the fixed effects at the final
       tau by construction, so it incidentally performs the refresh R does and
       lands closer to R. The leftover ~2.5e-3 "LOCO delta" the C++ shows on
       every chromosome (vs R's chromosome-specific ~8e-4) is that missing
       full-genome refresh, NOT the chromosome exclusion. Tightening tolPCG from
       1e-5 to 1e-9 changes it by <1%, ruling out PCG tolerance.
       => pre-existing non-LOCO bug, deliberately NOT fixed here because it
          would break the byte-identical non-LOCO guarantee.
"""
import os
import re
import sys

import numpy as np


def read_arma(path):
    with open(path, "rb") as f:
        f.readline()
        dims = f.readline().split()
        r, c = int(dims[0]), int(dims[1])
        a = np.fromfile(f, dtype=np.float64, count=r * c)
    return a.reshape((c, r)).T if c > 1 else a


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    d, rcsv = sys.argv[1], sys.argv[2]

    import csv
    with open(rcsv) as f:
        rows = list(csv.DictReader(f))
    R = {k: np.array([float(row[k]) for row in rows]) for k in rows[0]}

    full_c = read_arma(os.path.join(d, "mu.arma"))
    full_r = R["full"]
    theta = re.search(r'"theta":\s*\[([^\]]*)\]',
                      open(os.path.join(d, "nullmodel.json")).read()).group(1)
    print(f"cpp theta = [{theta}]")
    print(f"full  : max|cpp-R| = {np.abs(full_c-full_r).max():.3e}")

    chroms = sorted(int(k[3:]) for k in R if k.startswith("chr"))
    print(f"{'chr':>4} {'max|cpp-R|':>12} {'|d|cpp':>10} {'|d|R':>10} {'corr(d)':>9}")
    for j in chroms:
        p = os.path.join(d, f"chr{j}", "mu.arma")
        if not os.path.exists(p):
            print(f"{j:>4}   (no cpp chr{j}/)")
            continue
        c = read_arma(p)
        r = R[f"chr{j}"]
        dc, dr = c - full_c, r - full_r
        corr = np.corrcoef(dc, dr)[0, 1] if dc.std() and dr.std() else float("nan")
        print(f"{j:>4} {np.abs(c-r).max():12.3e} {np.abs(dc).max():10.3e} "
              f"{np.abs(dr).max():10.3e} {corr:9.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
