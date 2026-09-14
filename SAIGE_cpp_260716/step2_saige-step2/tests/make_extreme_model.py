#!/usr/bin/env python3
"""Scale a null model's residuals so that ordinary markers produce p == 0.

MULTITRAIT_DESIGN.md section 9.4: when the chi-square upper tail underflows,
the score test stops printing "%.6E" and switches to the fixed-point
"%.1fE%d" form built from log_chisq1_uppertail. That branch has to behave
identically in the batch kernel and in the scalar path, and the only pairs that
reach it WITHOUT also being routed to the fallback are quantitative ones
(quantitative never triggers SPA or Firth). Random test phenotypes never get
anywhere near p < 1e-300, so this manufactures the case.

For a quantitative trait the score test reads res only through

    S   = ( g'res - S_a'Z ) / tau0        with  S_a = colSums(X * res)
    var2 = Z'XVX Z * tau0 + g'g - 2 g'X Z   (no res at all)

so multiplying res and S_a by the same k multiplies S by k and leaves var
alone: stat grows by k^2 and every marker's p-value underflows. The model stays
internally consistent, and the golden single-trait run and the multi-trait run
read the identical files, so the comparison between them is still exact.

usage: make_extreme_model.py <src model dir> <dst model dir> <k>
"""

import os
import shutil
import sys

import numpy as np

from make_reduced_p_model import read_arma, write_arma


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    src, dst, k = sys.argv[1], sys.argv[2], float(sys.argv[3])
    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    res = read_arma(os.path.join(src, "res.arma")).ravel() * k
    S_a = read_arma(os.path.join(src, "S_a.arma")).ravel() * k
    write_arma(os.path.join(dst, "res.arma"), res)
    write_arma(os.path.join(dst, "S_a.arma"), S_a)
    print(f"make_extreme_model: {src} -> {dst}, res and S_a scaled by {k:g}")


if __name__ == "__main__":
    main()
