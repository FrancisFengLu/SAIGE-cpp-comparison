#!/usr/bin/env python3
"""sim_cohort.py -- simulated test data for the collaborator test protocol (no real data).
Pure Python standard library; output depends only on the arguments (fixed seeds).

  sim_cohort.py geno  OUT_PREFIX --n N --m M [--chroms C] [--rare-frac F] [--sib-frac S] [--seed K]
      PLINK .bed/.bim/.fam. Markers spread evenly over chromosomes 1..C. A fraction F of the
      markers is rare (expected minor-allele count 1..19), the rest has MAF uniform in [0.01, 0.5].
      A fraction S of the samples comes in full-sibling pairs (kinship 0.25), so a sparse GRM has
      off-diagonal entries. 0.5% missing calls. Writes OUT_PREFIX.pheno_score (a polygenic score per
      sample, used by `pheno --score`).
  sim_cohort.py pheno FAM OUT_TXT --nbin K [--score FILE] [--gscale G] [--seed K]
      Phenotype/covariate table: IID, b1..bK (0/1), x1 (normal), x2 (0/1). Case fractions cycle
      through 5, 10, 20, 30, 50 %. b2 is missing (NA) for 10% of the samples, so it has its own
      sample set. With --score the liability includes GSCALE x the polygenic score (default 0.6).
"""
import argparse, math, random, struct, sys

def geno(a):
    rng = random.Random(a.seed)
    n, m = a.n, a.m
    # samples: sibling pairs first, then unrelated
    nsib = int(n * a.sib_frac) // 2 * 2
    fam = [(f"s{i+1}", f"s{i+1}") if i >= nsib else (f"f{i//2+1}", f"s{i+1}") for i in range(n)]
    with open(a.out + ".fam", "w") as f:
        for fid, iid in fam:
            f.write(f"{fid}\t{iid}\t0\t0\t0\t-9\n")
    per = (m + a.chroms - 1) // a.chroms
    nscore = min(300, m)
    score = [0.0] * n
    bytes_per = (n + 3) // 4
    bim = open(a.out + ".bim", "w")
    bed = open(a.out + ".bed", "wb")
    bed.write(bytes([0x6C, 0x1B, 0x01]))
    code = (0b00, 0b10, 0b11)          # PLINK: 00 hom A1, 10 het, 11 hom A2, 01 missing
    for j in range(m):
        chrom = j // per + 1
        if rng.random() < a.rare_frac:
            maf = rng.randint(1, 19) / (2.0 * n)
        else:
            maf = rng.uniform(0.01, 0.5)
        g = [0] * n
        i = 0
        while i < nsib:                 # two sibs from four parental haplotypes
            h = [rng.random() < maf for _ in range(4)]
            for k in (0, 1):
                g[i + k] = h[rng.randrange(2)] + h[2 + rng.randrange(2)]
            i += 2
        for i in range(nsib, n):
            g[i] = (rng.random() < maf) + (rng.random() < maf)
        if j < nscore and maf > 0.01:
            sd = math.sqrt(2 * maf * (1 - maf))
            for i in range(n):
                score[i] += (g[i] - 2 * maf) / sd
        row = bytearray(bytes_per)
        for i in range(n):
            c = 0b01 if rng.random() < 0.005 else code[g[i]]
            row[i >> 2] |= c << ((i & 3) * 2)
        bed.write(row)
        bim.write(f"{chrom}\tm{j+1}\t0\t{1000 * (j % per + 1)}\tA\tG\n")
    bed.close(); bim.close()
    sd = math.sqrt(sum(s * s for s in score) / n) or 1.0
    with open(a.out + ".pheno_score", "w") as f:
        for (fid, iid), s in zip(fam, score):
            f.write(f"{iid}\t{s / sd:.6f}\n")
    print(f"wrote {a.out}.bed/.bim/.fam: {n} samples ({nsib} in sibling pairs) x {m} markers on {a.chroms} chromosome(s)")

def pheno(a):
    rng = random.Random(a.seed)
    ids = [l.split()[1] for l in open(a.fam)]
    score = {}
    if a.score:
        for l in open(a.score):
            k, v = l.split(); score[k] = float(v)
    prev = [0.05, 0.10, 0.20, 0.30, 0.50]
    gauss = lambda: math.sqrt(-2 * math.log(1 - rng.random())) * math.cos(2 * math.pi * rng.random())
    def logit_inv(z): return 1 / (1 + math.exp(-z))
    with open(a.out, "w") as f:
        f.write("\t".join(["IID"] + [f"b{k+1}" for k in range(a.nbin)] + ["x1", "x2"]) + "\n")
        for iid in ids:
            x1 = gauss(); x2 = 1 if rng.random() < 0.5 else 0
            s = score.get(iid, 0.0)
            ys = []
            for k in range(a.nbin):
                p = prev[k % len(prev)]
                z = math.log(p / (1 - p)) + 0.3 * x1 + 0.2 * x2 + a.gscale * s
                y = 1 if rng.random() < logit_inv(z) else 0
                if k == 1 and rng.random() < 0.10:
                    y = "NA"
                ys.append(str(y))
            f.write("\t".join([iid] + ys + [f"{x1:.4f}", str(x2)]) + "\n")
    print(f"wrote {a.out}: {len(ids)} samples, {a.nbin} binary traits")

ap = argparse.ArgumentParser()
sp = ap.add_subparsers(dest="cmd", required=True)
g = sp.add_parser("geno"); g.add_argument("out"); g.add_argument("--n", type=int, required=True)
g.add_argument("--m", type=int, required=True); g.add_argument("--chroms", type=int, default=2)
g.add_argument("--rare-frac", type=float, default=0.3); g.add_argument("--sib-frac", type=float, default=0.2)
g.add_argument("--seed", type=int, default=1)
p = sp.add_parser("pheno"); p.add_argument("fam"); p.add_argument("out")
p.add_argument("--nbin", type=int, required=True); p.add_argument("--score", default="")
p.add_argument("--seed", type=int, default=7); p.add_argument("--gscale", type=float, default=0.6)
a = ap.parse_args()
geno(a) if a.cmd == "geno" else pheno(a)
