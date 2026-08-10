#!/usr/bin/env python3
"""Build a group file from a PLINK .bim.

Usage:
    make_group_file.py <bim> <out> [mode]

Modes:
    all      (default) one gene per chromosome, named GENE<chr>.
    chr<N>   only the gene for chromosome N (e.g. "chr1").
    mixed    a SINGLE gene named GENEMIX whose variants span every chromosome
             in the .bim.  Used to prove the LOCO region check rejects a
             chromosome-spanning region instead of silently trimming it.
    nonstdchr<N>
             like chr<N>, but the gene also carries variant IDs in spellings
             that do NOT parse as chr:pos:ref:alt (rsIDs and
             "chr1_12345_A_G" underscore form) and that the genotype file does
             not contain.  Real UKB group files use such spellings; the LOCO
             filter must take each variant's chromosome from the genotype file
             and simply drop the unknown IDs, not hard-error on the string.

Marker IDs use the chr:pos:ref:alt convention the region reader expects.
"""
import sys
from collections import defaultdict

bim, out = sys.argv[1], sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else "all"

by_chr = defaultdict(list)
with open(bim) as f:
    for line in f:
        c, rsid, cm, pos, a1, a2 = line.split()
        # .bim column order is A1 (alt) then A2 (ref) -> chr:pos:ref:alt
        by_chr[c].append("%s:%s:%s:%s" % (c, pos, a2, a1))

chroms = sorted(by_chr, key=lambda x: (len(x), x))


def write_gene(g, name, ids):
    g.write("%s\tvar\t%s\n" % (name, "\t".join(ids)))
    g.write("%s\tanno\t%s\n" % (name, "\t".join(["lof"] * len(ids))))


with open(out, "w") as g:
    if mode == "mixed":
        ids = []
        for c in chroms:
            ids.extend(by_chr[c])
        write_gene(g, "GENEMIX", ids)
        print("wrote %s with 1 gene spanning %d chromosomes" % (out, len(chroms)))
    elif mode.startswith("nonstdchr"):
        want = mode[len("nonstdchr"):]
        if want not in by_chr:
            sys.exit("no markers on chromosome %s in %s" % (want, bim))
        ids = list(by_chr[want])
        ids.append("rs9999991")                 # rsID: no colon at all
        ids.append("%s_999999_A_G" % want)      # underscore spelling
        ids.append("SOME_GENE_VARIANT_7")       # free-form label
        write_gene(g, "GENE%s" % want, ids)
        print("wrote %s with 1 gene (chromosome %s + 3 non-standard IDs)"
              % (out, want))
    elif mode.startswith("chr"):
        want = mode[3:]
        if want not in by_chr:
            sys.exit("no markers on chromosome %s in %s" % (want, bim))
        write_gene(g, "GENE%s" % want, by_chr[want])
        print("wrote %s with 1 gene (chromosome %s)" % (out, want))
    elif mode == "all":
        for c in chroms:
            write_gene(g, "GENE%s" % c, by_chr[c])
        print("wrote %s with %d genes" % (out, len(chroms)))
    else:
        sys.exit("unknown mode: %s" % mode)
