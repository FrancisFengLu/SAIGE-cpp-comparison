#!/usr/bin/env python3
"""Copy a step 1 .arma model directory with its sample order permuted and/or
nullmodel.json keys edited.

Why this exists: the C++ step 1 always writes sampleIDs in .fam order, so two
real models can differ in which samples they hold but never in how those samples
are ordered. The multi-trait union / per-trait index maps (MULTITRAIT_DESIGN.md
section 4.7) must not depend on that, and SPA_Cutoff -- which decides how many
pairs leave the batch kernel -- is only settable in nullmodel.json.

--permute SEED reorders the samples consistently: sampleIDs and every
per-sample row of mu / res / V / offset / y / X / XVX_inv_XV / XXVX_inv, and
every per-sample column of XV. p-level blocks (XVX, XVX_inv, S_a) are sums over
samples and do not change. The result is the same model with its samples
listed in another order, so a single-trait run on it is a valid golden for the
same model inside a multi-trait run.

--json key=value sets a top-level nullmodel.json key (value parsed as JSON).

usage: make_model_variant.py <src model dir> <dst model dir> [--permute SEED] [--json k=v ...]
"""

import json
import os
import shutil
import sys

import numpy as np

HEADER = b"ARMA_MAT_BIN_FN008"
ROW_FILES = ["mu", "res", "V", "offset", "y", "X", "XVX_inv_XV", "XXVX_inv"]
COL_FILES = ["XV"]


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
    with open(path, "wb") as f:
        f.write(HEADER + b"\n")
        f.write(f"{m.shape[0]} {m.shape[1]}\n".encode())
        f.write(np.asfortranarray(m).tobytes(order="F"))


def main():
    args = sys.argv[1:]
    if len(args) < 2:
        sys.exit(__doc__)
    src, dst = args[0], args[1]
    seed = None
    edits = []
    k = 2
    while k < len(args):
        if args[k] == "--permute":
            seed = int(args[k + 1]); k += 2
        elif args[k] == "--json":
            key, val = args[k + 1].split("=", 1)
            edits.append((key, json.loads(val))); k += 2
        else:
            sys.exit(f"unknown argument {args[k]}")

    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    jpath = os.path.join(dst, "nullmodel.json")
    j = json.load(open(jpath))

    if seed is not None:
        n = len(j["sampleIDs"])
        perm = np.random.default_rng(seed).permutation(n)   # new k <- old perm[k]
        j["sampleIDs"] = [j["sampleIDs"][i] for i in perm]
        for name in ROW_FILES:
            p = os.path.join(dst, name + ".arma")
            if not os.path.exists(p):
                continue
            m = read_arma(p)
            if m.shape[0] != n:
                sys.exit(f"{p}: expected {n} rows, got {m.shape}")
            write_arma(p, m[perm, :])
        for name in COL_FILES:
            p = os.path.join(dst, name + ".arma")
            m = read_arma(p)
            if m.shape[1] != n:
                sys.exit(f"{p}: expected {n} columns, got {m.shape}")
            write_arma(p, m[:, perm])
        # obj_noK.json is not read by step 2; drop it rather than leave a
        # file that silently disagrees with the permuted arrays.
        o = os.path.join(dst, "obj_noK.json")
        if os.path.exists(o):
            os.remove(o)

    for key, val in edits:
        j[key] = val
    json.dump(j, open(jpath, "w"))


if __name__ == "__main__":
    main()
