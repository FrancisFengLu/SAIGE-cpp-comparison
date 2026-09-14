#!/usr/bin/env python3
"""Derive a lower-p null model from an existing .arma model directory.

Why this exists: every null model we can fit from the test data happens to have
the same number of covariates, so nothing in the suite exercises the case the
multi-trait design calls out explicitly -- one OpenMP thread alternating
between traits whose p differs, which is what the tl_X1 / tl_A1 grow test in
SAIGEClass::scoreTestFast has to survive (MULTITRAIT_DESIGN.md section 9.1).

The output is NOT a refit: mu / res / y / V are copied unchanged and only the
p-dependent blocks are rebuilt from the retained covariate columns, using the
same formulas tools/rda_to_arma.R uses when obj.noK is absent:

    XV          = (X * V)'                 p x N
    XVX         = X' (X * V)               p x p
    XVX_inv     = inv(XVX)
    XXVX_inv    = X XVX_inv                N x p
    XVX_inv_XV  = XXVX_inv scaled by V     N x p
    S_a         = colSums(X * res)         p

That makes the model internally consistent, which is all a numerical-equality
test needs: the golden single-trait run and the multi-trait run read the very
same files, so any difference between them is the implementation's, not the
model's.

usage: make_reduced_p_model.py <src model dir> <dst model dir> <keep columns>
       keep columns: number of leading columns of X to keep (>= 1, < p)
"""

import json
import os
import shutil
import struct
import sys

import numpy as np

HEADER = b"ARMA_MAT_BIN_FN008"


def read_arma(path):
    with open(path, "rb") as f:
        hdr = f.readline().strip()
        if hdr != HEADER:
            raise ValueError(f"{path}: unexpected header {hdr!r}")
        rows, cols = (int(x) for x in f.readline().split())
        buf = f.read(rows * cols * 8)
    if len(buf) != rows * cols * 8:
        raise ValueError(f"{path}: short read")
    return np.frombuffer(buf, dtype="<f8").reshape((rows, cols), order="F").copy()


def write_arma(path, m):
    m = np.asarray(m, dtype="<f8")
    if m.ndim == 1:
        m = m.reshape((-1, 1))
    with open(path, "wb") as f:
        f.write(HEADER + b"\n")
        f.write(f"{m.shape[0]} {m.shape[1]}\n".encode())
        f.write(np.asfortranarray(m).tobytes(order="F"))


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    src, dst, keep = sys.argv[1], sys.argv[2], int(sys.argv[3])

    X = read_arma(os.path.join(src, "X.arma"))
    V = read_arma(os.path.join(src, "V.arma")).ravel()
    res = read_arma(os.path.join(src, "res.arma")).ravel()
    n, p = X.shape
    if not 1 <= keep < p:
        sys.exit(f"keep must be in [1, {p}), got {keep}")

    Xk = np.ascontiguousarray(X[:, :keep])
    XVt = Xk * V[:, None]              # N x keep, == (X'V)'
    XVX = Xk.T @ XVt                   # keep x keep
    XVX_inv = np.linalg.inv(XVX)
    XXVX_inv = Xk @ XVX_inv            # N x keep
    XVX_inv_XV = XXVX_inv * V[:, None]  # N x keep
    S_a = (Xk * res[:, None]).sum(axis=0)

    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    write_arma(os.path.join(dst, "X.arma"), Xk)
    write_arma(os.path.join(dst, "XV.arma"), XVt.T)
    write_arma(os.path.join(dst, "XVX.arma"), XVX)
    write_arma(os.path.join(dst, "XVX_inv.arma"), XVX_inv)
    write_arma(os.path.join(dst, "XXVX_inv.arma"), XXVX_inv)
    write_arma(os.path.join(dst, "XVX_inv_XV.arma"), XVX_inv_XV)
    write_arma(os.path.join(dst, "S_a.arma"), S_a)

    jpath = os.path.join(dst, "nullmodel.json")
    with open(jpath) as f:
        j = json.load(f)
    j["p"] = keep
    j["alpha"] = list(np.asarray(j.get("alpha", [0.0] * keep), dtype=float)[:keep]) or [0.0]
    with open(jpath, "w") as f:
        json.dump(j, f)
    print(f"make_reduced_p_model: {src} (p={p}) -> {dst} (p={keep}), n={n}")


if __name__ == "__main__":
    main()
