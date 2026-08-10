#!/usr/bin/env python3
"""Build a synthetic LOCO null-model directory for step-2 LOCO tests.

Copies an existing (non-LOCO) step-1 output directory, then writes chr1/ and
chr2/ subdirectories holding the 10 per-chromosome files, each perturbed by a
different, chromosome-specific factor.  The perturbation is what makes the
tests meaningful: if the loader silently reads the top-level file instead of
the chr<j>/ copy, the LOCO results collapse onto the non-LOCO results and the
test fails.

X.arma and y.arma are deliberately NOT copied into chr<j>/ -- they are
chromosome-invariant per LOCO_FORMAT.md.

.arma format: ASCII header "ARMA_MAT_BIN_FN008\\n<rows> <cols>\\n" followed by
raw little-endian float64, column-major.
"""
import os
import shutil
import struct
import sys

PER_CHROM = ["mu", "res", "V", "offset", "XV", "XVX",
             "XVX_inv", "XVX_inv_XV", "XXVX_inv", "S_a"]


def read_arma(path):
    with open(path, "rb") as f:
        raw = f.read()
    nl1 = raw.index(b"\n")
    header = raw[:nl1].decode()
    nl2 = raw.index(b"\n", nl1 + 1)
    rows, cols = (int(x) for x in raw[nl1 + 1:nl2].split())
    body = raw[nl2 + 1:]
    n = rows * cols
    assert len(body) == 8 * n, (path, len(body), n)
    vals = list(struct.unpack("<%dd" % n, body))
    return header, rows, cols, vals


def write_arma(path, header, rows, cols, vals):
    with open(path, "wb") as f:
        f.write(header.encode() + b"\n")
        f.write(("%d %d" % (rows, cols)).encode() + b"\n")
        f.write(struct.pack("<%dd" % len(vals), *vals))


def perturb(name, vals, factor):
    """Chromosome-specific, value-preserving-ish perturbation.

    mu must stay strictly inside (0, 1) for the binary trait, so it is pulled
    toward 0.5 by `factor`.  Everything else is scaled.
    """
    if name == "mu":
        return [0.5 + (v - 0.5) * factor for v in vals]
    return [v * factor for v in vals]


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    # chrom -> perturbation factor.  Distinct per chromosome AND != 1.0.
    factors = {1: 0.90, 2: 1.15}

    if os.path.exists(dst):
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

    for chrom, factor in factors.items():
        cdir = os.path.join(dst, "chr%d" % chrom)
        os.makedirs(cdir, exist_ok=True)
        for name in PER_CHROM:
            p = os.path.join(src, name + ".arma")
            if not os.path.exists(p):
                raise SystemExit("missing source file: " + p)
            header, rows, cols, vals = read_arma(p)
            write_arma(os.path.join(cdir, name + ".arma"),
                       header, rows, cols, perturb(name, vals, factor))
        # Sanity: X / y must NOT be duplicated into chr<j>/
        for name in ("X", "y"):
            assert not os.path.exists(os.path.join(cdir, name + ".arma"))

    # Patch nullmodel.json: add "loco": true and "loco_chroms".
    jpath = os.path.join(dst, "nullmodel.json")
    with open(jpath) as f:
        txt = f.read()
    assert '"loco"' not in txt
    i = txt.rindex("}")
    chroms = ", ".join(str(c) for c in sorted(factors))
    txt = txt[:i].rstrip().rstrip(",") + ',\n  "loco": true,\n  "loco_chroms": [%s]\n}\n' % chroms
    with open(jpath, "w") as f:
        f.write(txt)

    print("wrote LOCO model at %s (chroms %s)" % (dst, chroms))


if __name__ == "__main__":
    main()
