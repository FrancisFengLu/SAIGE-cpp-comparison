#!/usr/bin/env python3
"""Generate a SAIGE benchmark dataset: PLINK bed/bim/fam + pheno + group file.

Two genotype sets are produced:

  <out>/grm.{bed,bim,fam}   common markers (MAF 0.05-0.5) used for the step-1 GRM
  <out>/wes.{bed,bim,fam}   rare markers organised into genes, used for step 2
  <out>/group.txt           SAIGE-GENE group file over the wes markers
  <out>/pheno.txt           IID, y_quantitative, y_binary, x1, x2  (tab separated)
  <out>/design.csv          same content, comma separated (C++ step1 wants CSV)
  <out>/meta.txt            dimensions

Genetic model
  - population structure: K latent ancestry factors perturb per-marker allele
    frequencies per individual, so the GRM is not the identity.
  - polygenic: y = 0.5*x1 - 0.3*x2 + g + e, g built from the standardised
    common markers with h2 = --h2.  This is what makes tau > 0.
  - rare gene markers are null by default; --causal-genes makes a fraction of
    genes carry a burden effect so the region test is not testing a pure null.

Usage:
  gen_bench_data.py --out DIR --n 10000 --grm-markers 20000 \
                    --genes 400 --min-var 10 --max-var 100 --seed 1
"""
import argparse
import os
import sys

import numpy as np

BED_MAGIC = bytes([0x6C, 0x1B, 0x01])

# PLINK .bed 2-bit codes, SNP-major.  With AlleleOrder=alt-first the .bim A1
# column is the ALT allele, so dosage(ALT) maps as:
#   dosage 2 -> 0b00, dosage 1 -> 0b10, dosage 0 -> 0b11, missing -> 0b01
_CODE = np.array([0b11, 0b10, 0b00], dtype=np.uint8)  # index = alt dosage


def pack_marker(dos, nbytes):
    """dos: uint8 array of alt dosages (0/1/2), length n. -> packed bytes."""
    codes = _CODE[dos]
    n = codes.size
    pad = nbytes * 4 - n
    if pad:
        # pad with 0b00 (hom-alt); those slots are never read back
        codes = np.concatenate([codes, np.zeros(pad, dtype=np.uint8)])
    c = codes.reshape(nbytes, 4)
    return (c[:, 0] | (c[:, 1] << 2) | (c[:, 2] << 4) | (c[:, 3] << 6)).tobytes()


def draw_genotypes(rng, freq, n, U=None, loadings=None):
    """Binomial(2, p_i) with optional per-individual frequency perturbation."""
    if U is None:
        p = np.full(n, freq)
    else:
        p = freq + U @ loadings
        np.clip(p, 1e-4, 1.0 - 1e-4, out=p)
    return rng.binomial(2, p).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--grm-markers", type=int, default=20000)
    ap.add_argument("--genes", type=int, default=400)
    ap.add_argument("--min-var", type=int, default=10)
    ap.add_argument("--max-var", type=int, default=100)
    ap.add_argument("--nfactors", type=int, default=4)
    ap.add_argument("--h2", type=float, default=0.4)
    ap.add_argument("--causal-genes", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    n = a.n
    nbytes = (n + 3) // 4
    os.makedirs(a.out, exist_ok=True)

    iids = ["s%d" % i for i in range(n)]

    # ---- latent ancestry (drives GRM off-diagonal structure) --------------
    U = rng.normal(0.0, 1.0, size=(n, a.nfactors)) * 0.05

    # ---- common markers / GRM file ---------------------------------------
    Mg = a.grm_markers
    # spread over 22 autosomes
    grm_chr = np.sort(rng.integers(1, 23, size=Mg))
    grm_freq = rng.uniform(0.05, 0.5, size=Mg)
    beta = rng.normal(0.0, 1.0, size=Mg)
    g = np.zeros(n)

    sys.stderr.write("[gen] writing %d common markers x %d samples\n" % (Mg, n))
    with open(os.path.join(a.out, "grm.bed"), "wb") as bed, \
         open(os.path.join(a.out, "grm.bim"), "w") as bim:
        bed.write(BED_MAGIC)
        pos = {}
        for j in range(Mg):
            loadings = rng.normal(0.0, 0.5, size=a.nfactors) * grm_freq[j]
            dos = draw_genotypes(rng, grm_freq[j], n, U, loadings)
            bed.write(pack_marker(dos, nbytes))
            c = int(grm_chr[j])
            pos[c] = pos.get(c, 0) + 1000
            bim.write("%d\tc%d_%d\t0\t%d\tA\tG\n" % (c, c, pos[c], pos[c]))
            # accumulate polygenic value from standardised genotypes
            d = dos.astype(np.float64)
            p = d.mean() / 2.0
            sd = np.sqrt(2.0 * p * (1.0 - p)) if 0 < p < 1 else 1.0
            g += beta[j] * (d - 2.0 * p) / sd
            if (j + 1) % 5000 == 0:
                sys.stderr.write("       %d/%d\n" % (j + 1, Mg))
    g /= np.sqrt(Mg)
    g = (g - g.mean()) / g.std()

    # ---- rare markers organised into genes -------------------------------
    nvar = rng.integers(a.min_var, a.max_var + 1, size=a.genes)
    Mr = int(nvar.sum())
    gene_chr = np.sort(rng.integers(1, 23, size=a.genes))
    is_causal = rng.random(a.genes) < a.causal_genes
    burden = np.zeros(n)

    sys.stderr.write("[gen] writing %d rare markers in %d genes\n" % (Mr, a.genes))
    annos = np.array(["lof", "missense", "synonymous"])
    group_lines = []
    gpos = {}
    with open(os.path.join(a.out, "wes.bed"), "wb") as bed, \
         open(os.path.join(a.out, "wes.bim"), "w") as bim:
        bed.write(BED_MAGIC)
        for gi in range(a.genes):
            c = int(gene_chr[gi])
            ids, anns = [], []
            for _ in range(int(nvar[gi])):
                # rare spectrum: MAF ~ 10^U(-4,-2), floored so MAC >= 1
                maf = max(10.0 ** rng.uniform(-4.0, -2.0), 1.5 / (2.0 * n))
                dos = draw_genotypes(rng, maf, n)
                if dos.sum() == 0:            # guarantee at least one carrier
                    dos[rng.integers(0, n)] = 1
                bed.write(pack_marker(dos, nbytes))
                gpos[c] = gpos.get(c, 0) + 100
                p = gpos[c]
                bim.write("%d\t%d:%d:A:C\t0\t%d\tC\tA\n" % (c, c, p, p))
                ids.append("%d:%d:A:C" % (c, p))
                anns.append(str(rng.choice(annos, p=[0.25, 0.5, 0.25])))
                if is_causal[gi]:
                    burden += 0.35 * dos
            name = "GENE%04d" % gi
            group_lines.append("%s var %s\n" % (name, " ".join(ids)))
            group_lines.append("%s anno %s\n" % (name, " ".join(anns)))
            if (gi + 1) % 100 == 0:
                sys.stderr.write("       gene %d/%d\n" % (gi + 1, a.genes))

    with open(os.path.join(a.out, "group.txt"), "w") as f:
        f.writelines(group_lines)

    # ---- fam files --------------------------------------------------------
    fam = "".join("%s %s 0 0 0 -9\n" % (i, i) for i in iids)
    for pref in ("grm", "wes"):
        with open(os.path.join(a.out, pref + ".fam"), "w") as f:
            f.write(fam)

    # ---- phenotype --------------------------------------------------------
    x1 = rng.normal(size=n)
    x2 = rng.binomial(1, 0.5, size=n).astype(float)
    fixed = 0.5 * x1 - 0.3 * x2
    ve = np.sqrt((1.0 - a.h2) / a.h2)
    yq = fixed + g + burden + rng.normal(0.0, ve, size=n)
    eta = -1.0 + 0.4 * x1 - 0.2 * x2 + 0.7 * g + burden
    yb = rng.binomial(1, 1.0 / (1.0 + np.exp(-eta)))

    with open(os.path.join(a.out, "pheno.txt"), "w") as f:
        f.write("IID\ty_quantitative\ty_binary\tx1\tx2\n")
        for i in range(n):
            f.write("%s\t%.6f\t%d\t%.6f\t%.1f\n" % (iids[i], yq[i], yb[i], x1[i], x2[i]))
    with open(os.path.join(a.out, "design.csv"), "w") as f:
        f.write("IID,y_quantitative,y_binary,x1,x2\n")
        for i in range(n):
            f.write("%s,%.6f,%d,%.6f,%.1f\n" % (iids[i], yq[i], yb[i], x1[i], x2[i]))

    with open(os.path.join(a.out, "meta.txt"), "w") as f:
        f.write("n_samples          %d\n" % n)
        f.write("grm_markers        %d\n" % Mg)
        f.write("genes              %d\n" % a.genes)
        f.write("wes_markers        %d\n" % Mr)
        f.write("variants_per_gene  min=%d max=%d mean=%.1f\n"
                % (nvar.min(), nvar.max(), nvar.mean()))
        f.write("causal_genes       %d\n" % int(is_causal.sum()))
        f.write("h2                 %.2f\n" % a.h2)
        f.write("case_rate          %.4f\n" % yb.mean())
        f.write("seed               %d\n" % a.seed)
    sys.stderr.write("[gen] done -> %s\n" % a.out)
    sys.stderr.write(open(os.path.join(a.out, "meta.txt")).read())


if __name__ == "__main__":
    main()
