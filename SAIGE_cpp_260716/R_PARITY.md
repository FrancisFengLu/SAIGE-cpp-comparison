# Parity with R SAIGE 1.5.2 — catalogue of every known difference

Reference implementation: **R SAIGE 1.5.2** (`/opt/saige/SAIGE-upstream`, installed
into `R_LIBS_USER=/opt/saige/Rlib-upstream`), entry points
`extdata/step1_fitNULLGLMM.R` and `extdata/step2_SPAtests.R`.

Everything below was measured, not inferred, unless a line says otherwise. Both
sides always read **one** step-1 fit: R reads the `.rda`, the C++ reads an
`.arma` directory converted from that same `.rda` by
`step2_saige-step2/tools/rda_to_arma.R` (LOCO: `tools/rda_to_arma_loco.R`).
Re-fitting step 1 twice would fold step 1's own Monte-Carlo spread into every
step-2 p-value and make a step-2 defect indistinguishable from noise.
All comparisons at `nThreads=1` (region row order is nondeterministic above 1).

Classification used throughout:

| | meaning | action |
|---|---|---|
| **(a)** | defect in this C++ port | fixed here |
| **(b)** | defect in R | C++ keeps its own behaviour; argument recorded below |
| **(c)** | acceptable difference | magnitude and source recorded |

The bar for **(b)** is an argument that does not depend on which implementation
you find more agreeable: R contradicting itself, R contradicting the file format,
R contradicting a third tool, or undefined behaviour. "The C++ looks more
sensible" is not enough — anything that weak is filed under (c) with the
uncertainty stated.

**Out of scope this round:** survival traits (not ported —
`step2_saige-step2/spa.hpp:7`), everything involving a sparse GRM.

---

## 1. Summary table

| # | path | difference | magnitude | class | disposition |
|---|---|---|---|---|---|
| A1 | region, SKAT-O | converged Davies integral discarded over qfc round-off, replaced by Liu | ≤1.8e-2 rel in the null range; up to **2.8× anti-conservative** near genome-wide significance | **(a)** | fixed — `skat.cpp` |
| A2 | region, all test types | `isOutputMarkerList` parsed but no file ever written | feature absent | **(a)** | fixed — `main.cpp` |
| A3 | region, SKAT-O | `SE_Burden` left over from an unprinted p-value when `Pvalue_Burden == 1` | 1 row in 300; `Inf` vs `0.0012` | **(a)** | fixed — `skat.cpp` |
| A4 | VCF, HDS field | missing haplotype read as dosage 0, `MissingRate` reported 0 | \|Δlog10 p\| 0.126, BETA sign flip | **(a)** | fixed — `genotype_reader.cpp` |
| A5 | BGEN / PGEN | global `AlleleOrder` default `alt-first` silently swaps REF/ALT | AF 0.694794 → 0.305206, BETA reversed | **(a)** | fixed — `main.cpp` |
| A6 | PGEN mode 0x01 | decoded with the PLINK-2 code table | 2000 rows → 798, \|Δlog10 p\| 2.19 | **(a)** | fixed — `genotype_reader.*` |
| A7 | single-variant + region, `--condition` | shared SPA inputs read uninitialised (inherited from R) | p 1.37e-95 vs 4.38e-02 | **(b)** in R, **(a)** as inherited | fixed — `saige_test.cpp` |
| A8 | all outputs | non-finite doubles spelled `inf`/`nan`, not `Inf`/`NaN` | cells silently lost by `read.table` | **(a)** | fixed — `UTIL.*`, `main.cpp` |
| A9 | region, PLINK group file | variant IDs resolved through disjoint key spaces; both sides fail silently | row sets completely disjoint | interface split | fixed — C++ now accepts both |
| A10 | VCF | `vcfField` defaulted to `GT`; R defaults to `DS` | `MissingRate` 0.80 vs 0.01 on the same file | **(a)** | fixed — `main.cpp` |
| B1 | single-variant + region, PLINK + `--condition` | R misreads the marker at genoIndex 0 | that marker's whole row is another marker's data | **(b)** | C++ keeps the absolute seek |
| B2 | region, non-default `weights.beta` | R applies it in some places and not others | `BETA_Burden` off by 24–31×, p by up to 0.76 rel | **(b)** | C++ applies it everywhere |
| B3 | binary ER, `max_MAC_for_ER ≥ 11` | 32-bit factorial overflow → integer division by zero → SIGFPE | R process dies | **(b)** | C++ uses a saturating `C(n,r)` |
| B4 | `is_noadjCov=TRUE`, AF > 0.5 | centring inconsistent with the flip | var up to 35× too large, p up to 10.5 orders too conservative | **(b)** | C++ copies the arithmetic exactly but **defaults the flag off** |
| B5 | VCF `DS` with a truncated sample field | R reads a missing `./.` as dosage 0 and reports `MissingRate` 0 | \|Δlog10 p\| 0.281, BETA sign flips | **(b)** | C++ treats it as missing |
| C1 | region, BURDEN-only output | R writes 6 significant digits | ≤4.9e-6 rel | **(c)** | printing only |
| C2 | region, `max_MAF` column | `1e-04` (R) vs `0.0001` (C++) | none | **(c)** | same value |
| C3 | `--condition` on the conditioning marker itself | 0/0 after catastrophic cancellation | `BETA_c` O(1) noise; `p.value_c` = 1 both sides | **(c)** | no statistical content |
| C4 | region, SKAT-O / Liu / CCT arithmetic | summation order, LAPACK | ≤2.2e-8 rel on `Pvalue` | **(c)** | float noise |
| C5 | binary ER resampling branch | different RNGs | not yet reachable | **(c)** | unverifiable today |

---

## 2. (a) — defects in this port, fixed

### A1 · SKAT-O threw away a converged Davies integral over qfc round-off

`step2_saige-step2/skat.cpp` — `davies_pvalue_strict` and the SKAT-O integrand.

SKAT-O's outer integral calls Davies once per quadrature node. Deep in the tail
qfc (`acc=1e-6`) returns `Qq = 1 + O(1e-8)` **with `ifault == 0`** — its own
series-truncation error, not a failure. The port rejected `pval <= 0 || pval > 1`
as failure, and the caller responded by discarding the entire, already converged,
Davies integral and recomputing the p-value with the Liu approximation.

R does not: `SKAT:::SKAT_Optimal_Integrate_Func_Davies` `stop()`s only on
`ifault != 0` and otherwise clamps, `if (temp > 1) temp = 1`.

Evidence it is a rule difference and not a qfc porting difference: the same
`q = 2596.8211502474296` handed to R's own `SKAT:::SKAT_davies` returns
`Qq = 1.0000000169605003, ifault = 0` — agreeing with our qfc to 1e-16. All five
triggering rows had raw `Qq ∈ [1+1.7e-8, 1+2.7e-7]`. Reconciling
`region_binary_nocollapse` row by row: 80/90 rows equalled R's Davies result
(rel < 1e-9), 10/90 equalled R's Liu fallback (rel < 1e-8), 0 unexplained. R is
the correct side — recomputing R's integral at `rel.tol=1e-11` leaves R within
~1e-8 of the high-precision value while the C++ was 2e-3 to 7e-3 away.

Damage was worst exactly where it matters. Scaling `Score` into the significant
range and calling each side's SKAT-O on R's own Score/Phi:

```
RGENE026_0.001 ×3   R=1.53e-08  cpp=5.53e-09   2.8× anti-conservative
RGENE000_1e-4  ×4   R=9.86e-09  cpp=1.42e-08   1.4× conservative
RGENE000_0.001 ×3   R=1.94e-08  cpp=2.50e-08
```

Note the deliberate asymmetry that remains: `davies_pvalue` (non-strict, used for
the per-ρ p-values) **keeps** the `p > 1 || p <= 0 → Liu` rule, because that is
what `SKAT:::Get_PValue.Lambda` does at that call site. The two sites have
different rules in R and now have different rules here. Measured: per-ρ
p-values agreed to ≤2.9e-7 across 54 scaled test cases, so that path was already
right.

**This is also the true source of the "gene-based p, max relative difference
8.6e-3" residual** that had previously been attributed to imputation plus
"SKAT-O internal integration".

Verified: `region_binary_nocollapse` 18 differing `Pvalue` rows → IDENTICAL;
`region_binary_skato` 36 → IDENTICAL (max rel 7.0e-12); `region_binary_mid`
12 rows at ≤8.6e-3 → IDENTICAL (5.2e-10); `region_quant_mid` 6 → IDENTICAL
(2.2e-8); `region_binary_collapse20`, `region_binary_weightsfile`,
`region_quant_skato` likewise. All `.singleAssoc.txt` output stayed bit-identical.

### A2 · `isOutputMarkerList` was parsed and then ignored

`main.cpp` read the key, echoed it in the config dump, and that was all —
`grep isOutputMarkerList` found no other use and no file-opening code anywhere.
R writes `<OutputFile>.markerList.txt` (`SAIGE_SPATest_Region.R:1039-1085` builds
it, `:1172-1215` writes it): one row per (region, annotation, max_MAF) mask
holding at least one variant, with the mask's rare variants and its collapsed
ultra-rare variants as two comma-separated fields. It is the only record of
*which* variants a gene-level p-value was computed from, so without it a region
result cannot be audited.

Now implemented, before the BURDEN / SKAT-O split (R builds the list outside its
own `regionTestType != "BURDEN"` block, so burden-only runs get the file too).
Verified: 90 rows (1 annotation × 3 maxMAF) and 270 rows (3 × 3) — same keys,
same order, every variant string byte-identical to R, 0 content differences; on
a burden-only run likewise identical; with the flag off no file is created.

### A3 · `SE_Burden` could come from a p-value that is never printed

`se_Burden` is assigned twice in `SKATO_region_test`: first from the exact
chi²(1) burden p at the top of the function, then again from the final ρ=1 Davies
p that is actually written out. Both assignments were guarded by `if (z != 0.0)`.
When the final burden p is exactly 1, `qnorm(0.5) == 0`, the second assignment
was skipped, and the field kept the value derived from the *first* p — back-solving
`0.00120314264519663` gives an internal p of about 0.9944, which appears nowhere
in the output. R has no such guard (`SAIGE_SPATest_Region_Func.R:343,351`
evaluates `abs(BETA_Burden/qnorm(p/2))` unconditionally, giving `Inf`).

Neither value is statistically meaningful. The point is that the printed SE is
now derived from the printed p, so a reader cannot mistake a stale finite number
for a real standard error.

### A4 · HDS: a missing haplotype was read as dosage 0

The HDS branch started at `dose_val = 0` and *skipped* missing or vector-end
haplotypes before summing the rest, so the `isnan` test below it could never
fire for an HDS sample. `HDS=".,."` therefore became a confident dosage of 0 and
`HDS="0.5,."` a confident 0.5, both with `MissingRate` reported as 0 — those
samples never reached mean imputation. R gets this right by a different route:
savvy's `stride_reduce` adds the haplotypes directly, so a NaN poisons the sum
and `VCF.cpp:250`'s `isnan` catches it.

```
marker                       R                      cpp (before)
h0  HDS=".,."     (both)     MissRate 0.01 AC 1010.1    MissRate 0 AC 1000
h1  HDS="0.5,."   (one)      MissRate 0.01 AC 1010.1    MissRate 0 AC 1250
```

Fixed, with `bcf_float_vector_end` kept distinct from missing: vector-end marks a
sample with *fewer* haplotypes than the record stride (a haploid call in a
diploid record, e.g. chrX outside the PAR), whose dosage is well defined from the
haplotypes present. Reaches anyone whose VCF header declares HDS, including users
who asked for `--vcfField=GT` or `DS`, since both sides fall back to HDS then.

### A5 · BGEN/PGEN `AlleleOrder` defaulted to `alt-first`

One global default, `alt-first`, was applied to every reader. For BGEN that is
wrong by the format's own definition — the BGEN spec makes the first listed
allele the reference — so `genoType: bgen` without an explicit `AlleleOrder` got
REF and ALT swapped, silently:

```
ref-first : snp0 Allele1=A Allele2=G AF_Allele2=0.694794   (agrees with plink + vcf)
alt-first : snp0 Allele1=G Allele2=A AF_Allele2=0.305206   BETA reversed
```

R never even starts in that configuration (`Geno.R:232-235` `stop()`s), so there
was no R run to compare against and the defect could only surface as a wrong
answer. Defaults now follow `Geno.R` per reader — plink `alt-first` (:172), bgen
`ref-first` (:232), pgen `ref-first` (:203) — and bgen/pgen refuse anything else,
as R does. An `AlleleOrder` that is neither spelling is now rejected outright
instead of being passed to the reader.

### A6 · PGEN mode 0x01 decoded with the PLINK-2 code table

`readPgenHeader()` accepts modes 0x01 and 0x02, but both decode loops hardcoded
the mode-0x02 table `00=0, 01=1, 10=2, 11=missing`. PLINK 1's table is different
— `00=hom A1, 01=missing, 10=het, 11=hom A2` — and plink2's `.pvar` for a PLINK-1
dataset carries ALT = A1, REF = A2. The old decode therefore reversed the
homozygotes *and* turned every real missing call into a heterozygote while
reporting `MissingRate` 0. Since `6c 1b 01` is literally the `.bed` magic, a plain
`.bed` file is a valid mode-0x01 pgen and this was reachable without doing
anything unusual. Result: 2000 R rows against 798 C++ rows, AF reversed
(0.998 where R had 0.009), `|Δlog10 p|` up to 2.19 — accepted, plausible-looking,
wrong.

Now a mode-dependent `m_genoCode[4]` lookup set once in the header reader. The
unsupported-mode message also fed a *decimal* `std::to_string` into an `0x`
prefix, so a mode 0x10 file was reported as "0x16"; fixed.

**Remaining capability gap (not fixed):** `plink2 --make-pgen` writes mode 0x10,
which this reader refuses outright while R reads it through real pgenlib.
Refusing is the honest behaviour and is left as is — supporting it needs
pgenlib. In practice the C++ pgen path is usable only for `.bed`-equivalent
files today.

### A7 · Conditional SPA read uninitialised variables (R's defect, inherited)

Listed under (a) even though the root defect is R's, because the port copied the
structure verbatim and so had the identical bug. R's side is
`SAIGE-upstream/src/SAIGE_test.cpp` — `:525,534-537` declare the variables,
`:555-647` is the only place that fills them, and `:840-860` reads them
unconditionally. Nothing about R's version is fixable from here, so what follows
is about this port; the same argument condemns R's.

`m1`, `p_iIndexComVecSize`, `tol1` and the `gNA / gNB / muNA / muNB / NAmu /
NAsigma` set are consumed by two blocks — the unconditional SPA and the
conditional SPA — but were filled only inside the first. A marker reaching the
conditional SPA without having gone through the unconditional one read
uninitialised stack. Two ways in, both observed on real data:

- `|StdStat| ≤ SPAcutoff < |StdStat_c|` (snp1819: 1.99621 vs 2.00515)
- the marker took the ER branch, which skips the whole non-ER block, while
  `stat_c` still exceeds `SPAcutoff²` (1:1271, 1:2588, MAC=4)

Consequences were not rounding-level: chr1:3275000 gave `1.371413E-95` in R
against `4.384698E-02` here; chr1:3550000 `1.913526E-01` against `9.862770E-284`.
Each side was byte-reproducible run to run and they disagreed with each other —
the signature of UB, not of numerical drift.

Fixed with a memoised `prepare_spa_inputs()` holding one copy of the expressions
that used to live inline, called from both blocks. Every expression is verbatim
from the original and all its operands are final by that point, so a marker that
*does* take the unconditional block sees bit-identical values.

Two independent checks that the C++ now produces **the value R itself produces
wherever R is well defined**, not merely a different one:

*ER toggle* (rare set, `condition=1:1500,1:2000`). Disabling ER routes the MAC=4
markers through the unconditional SPA, so R's variables get filled:

```
marker    R ER=off      cpp ER=off    cpp ER=on     R ER=on (UB)
1:947     4.143997E-05  4.143997E-05  4.143997E-05  3.321920E-02
1:1211    1.086367E-02  1.086367E-02  1.086367E-02  3.805660E-02
1:1785    1.261103E-05  1.261103E-05  1.261103E-05  3.181371E-02
1:2042    1.253951E-02  1.253951E-02  1.253951E-02  4.450010E-02
1:1271    2.023312E-02  2.023312E-02  2.023312E-02  2.023312E-02
1:2588    4.185001E-02  4.185001E-02  4.185001E-02  4.185001E-02
```

The ER=off run differs on only 2 of 3000 rows, both unrelated (B1 and C3).

*SPAcutoff sweep* (mid2k, `condition=snp1200,snp1500`). At `SPAcutoff=0.5`
exactly 7 markers satisfy `StdStat ≤ 0.5 < StdStat_c` — R's inputs never filled —
and exactly those 7 are the rows that differ. Every marker whose unconditional
block ran is bit-identical. snp1498 (StdStat 0.6362, so R's block *does* run at
0.5) matches the C++ default-cutoff value to the last digit:
R=`4.200722E-03`, cpp=`4.200722E-03`, cpp at cutoff 2.0=`4.200722E-03`, against
R's UB value `8.401432E-03`.

> **Honest note on row counts.** The conditional cases now show *more* differing
> rows against R, not fewer — `sb_cond_late` 1 → 3, `sb_cond_rare` 4 → 6,
> `sb_loco_chr1_cond` 5 → 6. Before the fix the C++ was reading its own garbage,
> which on some markers happened to coincide with R's garbage. Every added row is
> a marker where R reads uninitialised memory.

**Survival:** the conditional survival branch assigns `qinv`, not `qinv_c`
(`saige_test.cpp:1289`); R has the identical typo at
`SAIGE-upstream/src/SAIGE_test.cpp:854`. Left byte-compatible with R and
commented, since survival is out of scope. Fix both together if survival is ported.

### A8 · Non-finite doubles were spelled the C++ way

Where R writes `Inf`, the C++ wrote `inf`; `NaN` vs `nan`. Not cosmetic: R's own
readers disagree with the C++ spelling — `read.table()` and `data.table::fread()`
parse `Inf` as infinity but turn `inf` into `NA` or a character level, which can
coerce a whole numeric column to character. That is the entire downstream
audience for these files. First seen on the A3 row, whose `SE_Burden` is
legitimately infinite.

Fixed with an `RNumPut` `num_put` facet installed as the global locale in
`main()`, plus an explicit `imbue()` for the namespace-scope `std::ofstream`
objects constructed before `main()`. Only the non-finite branch is overridden;
finite values are delegated to the standard facet with the stream's own
precision flags, so finite output is byte-for-byte unchanged (regression-checked
on `base_binary_single`: 2000 rows × 14 columns, max rel 0.0).

### A9 · PLINK group-file variant IDs: two disjoint key spaces, both failing silently

Not a numerical difference — an interface split where one group file could not
drive both sides, and *neither* side said so.

- **R** resolves PLINK group-file variants only through `markerInfo$ID` = `.bim`
  column 2 (`Geno.R:182-196`; the PGEN branch at `:215` is the one that builds
  `CHROM:POS:REF:ALT` — the PLINK branch never does). Given a coordinate-spelled
  group file it printed "40000 markers in 'RegionFile' are not in 'GenoFile'",
  wrote a **zero-byte** output file, and exited 0.
- **C++** built only `chr:pos:ref:alt` (both allele orders). Given a `.bim`-name
  group file it printed "Skipping region … no matching variants" per gene plus
  "Total regions skipped: 100", and exited 0.

R is already self-inconsistent here (its PGEN path accepts only coordinates and
*overwrites* the `.pvar` ID to build them; its PLINK path accepts only names), so
there is no single "R behaviour" to copy. The C++ now accepts the union, with
positional keys inserted first so they win any collision. Every group file R
accepts works, and so does the `chr:pos:ref:alt` convention that real
SAIGE-GENE+ group files normally use.

Verified: `region_mid_rsid` went from completely disjoint row sets to IDENTICAL
(221 main rows, 767 `.singleAssoc.txt` rows); coordinate-spelled cases unchanged.

**Not changed:** the BGEN and VCF marker maps are still coordinate-only. The
region × bgen/vcf path was not compared this round.

### A10 · `vcfField` defaulted to `GT` where R defaults to `DS`

R defaults `vcfField` to `"DS"` — `extdata/step2_SPAtests.R:18` for the CLI,
`R/SAIGE_Test_main.R:68` and `R/Geno.R:26,133` for the function. This port
defaulted to `"GT"`.

Not cosmetic: on a `plink2 --export vcf vcf-dosage=DS` file the same records
read two entirely different ways — `MissingRate` 0.80 under GT against 0.01
under DS — because plink2 writes the hard call as `./.` wherever the dosage is
not near-integer. A config that omitted the key analysed a different set of
genotypes than R would, silently.

This one is worth noting as a methodological point: it never appeared as a
*measured* difference, because every parity case pinned `vcfField` explicitly on
both sides. It was hiding in the defaults, and only a source-level comparison
found it. The same is true of A5.

R additionally rejects anything other than `DS` or `GT` (`Geno.R:111-112`),
where this port passed the string through to the reader, which would only warn
that the field was absent and then return all-missing. Both behaviours now
match. The existing header check and `DS → HDS` fallback
(`genotype_reader.cpp:1017-1032`) are unchanged.

---

## 3. (b) — defects in R; the C++ keeps its own behaviour

### B1 · R's PLINK reader misreads the marker at genoIndex 0 under `--condition`

`SAIGE-upstream/src/PLINK.cpp:194-206`:

```cpp
uint64_t posSeek;
if (t_gIndex > 0) {                  // t_gIndex == 0 skips the whole block: no seek at all
   if (t_gIndex_prev == 0) { posSeek = 3 + m_numBytesofEachMarker0 * t_gIndex; fseek(..., SEEK_SET); }
   else { posSeek = m_numBytesofEachMarker0*(t_gIndex - t_gIndex_prev - 1); fseek(..., SEEK_CUR); }
}
fread(...);                          // reads from wherever the handle happens to be
```

In an ordinary run the handle is parked at offset 3 by the constructor and the
first marker happens to be the right one. With `--condition`,
`assign_conditionMarkers_factors` reads the conditioning markers first, leaving
the handle after the last of them — so the main loop's first marker
(`i == 0` → `gIndex_prev = 0`, `gIndex = 0`) reads **someone else's genotypes**
and prints them under its own CHR/POS/MarkerID/Allele1/Allele2. The variant
metadata comes from `t_gIndex` (:220-236) and is therefore correct; only the
genotypes are wrong.

**Why (b) and not a judgement call:** R contradicts itself — the same marker,
with and without `--condition`, yields two different genotype vectors — and the
wrong one is bit-for-bit another marker in the file. Pointer arithmetic predicts
*which* one, and the prediction was made in advance and confirmed four times out
of four (predicted value = the marker right after the last conditioning marker):

| conditioning markers | R's row for the first marker | equals |
|---|---|---|
| snp100, snp500 | AC=37605.8, miss=0.00986 | snp501, every field |
| snp1200, snp1500 | AC=59533.7, miss=0.01056 | snp1501 |
| 1:1500, 1:2000 (rare) | AC=63.3051, miss=0.00482 | 1:2001 |
| chr1:100000, chr1:500000 (LOCO) | AC=19018.4, miss=0.0088 | chr1:501000 |

Trait-independent: reproduced on both binary and quantitative with the same
dataset. A second, independent confirmation on `condition=snp0,snp1`: R's snp0
row is bit-identical to snp2's row from the unconditional run, down to snp2's own
`Allele1/Allele2` being mismatched against the snp0 labels that are printed.

**Minimal reproduction.** Run the same PLINK dataset twice, with and without
`--condition`, and diff the row for the marker at index 0 of the `.bim`:

```bash
tests/parity/parity.py base_binary_single                                  # r.txt
tests/parity/parity.py base_binary_single --set condition=snp100,snp500    # r.txt
# R's snp0 row in the second run == R's snp501 row in the first, field for field.
```

The C++ already forces an absolute `SEEK_SET` whenever `t_gIndex_prev == 0`
(`genotype_reader.cpp:832-845`, done originally so OpenMP threads could request
any marker) and gives snp0's own correct values in every case.

Blast radius in R: any PLINK + `--condition` run silently returns a wrong result
for whichever marker sits at genoIndex 0 of the tested set. The same
`gIndex_prev` code is in the region path (`Main.cpp:1232-1245`, `2253-2270`), so
region `--condition` runs are affected too.

### B2 · R applies a non-default `weights.beta` in some places and not others

`SAIGE-upstream/R/SAIGE_SPATest_Region.R:533` and `:842` **hardcode**
`AnnoWeights = dbeta(MAFVec, 1, 25)` when the group file carries no weight line,
ignoring the `weights.beta` passed on the command line — while the same
command-line value *does* reach upstream's own C++ through
`setAssocTest_GlobalVarsInCPP` and is used at `src/Main.cpp:1295` (per-marker
weights) and `:1578` (the weighted sum in the ultra-rare collapse). Changing
`weights.beta` therefore changes part of the computation and not the rest.

Measured (rare set, lof, collapse 10, only `weights.beta` changed from `1,25` to
`1,1`): every R `BETA_Burden` moves by the *same constant* factor 1.2322, while
the C++ values move by 24.4–31.5× — the magnitude you expect when a weight really
goes from ≈24.9 to 1. 116 of 120 rows differ, up to 0.76 relative on p.
Controlled with a mini group file: not caused by the annotation list or by
chunking (changing the annotation list from 1 to 2 entries leaves the `lof` rows
bit-identical).

**Minimal reproduction.**

```bash
tests/parity/parity.py region_binary_flatweights     # weights.beta = 1,1
tests/parity/parity.py region_binary_nocollapse      # same masks, default 1,25
# R's BETA_Burden moves by a constant 1.2322x between the two; the C++ moves by 24-31x.
```

**Why (b):** R is internally inconsistent, not merely different. The C++ applies
`weights.beta` uniformly. At the default `1,25` the two sides agree — which is
what every other region case in this document demonstrates.

### B3 · R's exact-resampling path overflows a 32-bit factorial and dies

`SAIGE-upstream/src/ER_binary_func.cpp:79-91`:

```cpp
int fact(int n) { ... return n * fact(n - 1); }                            // 13! > INT_MAX
int n_choose_r(int n, int r) { return fact(n) / (fact(r) * fact(n - r)); } // denominator wraps to 0
```

The wrapped denominator makes the division an integer divide-by-zero → **SIGFPE,
process killed**.

The trigger is subtler than "MAC too large": `k` is the length of
`indexNonZero`, not the MAC. `UTIL.cpp:98-102,125-131` shows mean imputation
writing `2*altFreq` into the missing positions (non-zero) *before*
`t_indexNonZero` is built, and the intervening
`if (t_MAC <= t_dosage_zerod_MAC_cutoff) t_GVec.clean(t_dosage_zerod_cutoff)` is
what would have zeroed those small dosages. So any marker whose post-imputation
MAC lands in `(dosage_zerod_MAC_cutoff, max_MAC_for_ER]` keeps them and gets
`k = real carriers + number of missing`, which can be in the hundreds.

Bisected to a single marker, `1:489:A:B` (AC=10.0535, MissingRate=0.00532 → 266
missing → k ≈ 276), which kills R **on its own**:

```
marker 489 alone, max_MAC_for_ER=11 -> 136 (SIGFPE)
marker 489 alone, max_MAC_for_ER=4  -> 0
```

Causal controls, one variable at a time:

| configuration | result |
|---|---|
| `impute=mean`, `dosage_zerod_MAC_cutoff=10` (default) | SIGFPE |
| `impute=mean`, `dosage_zerod_MAC_cutoff=20` (small dosages zeroed) | fine |
| `impute=best_guess` (imputeG=0) | fine |
| `impute=minor` (imputeG=0) | fine |
| `impute=mean`, `dosage_zerod_cutoff=0` (no zeroing at all) | SIGFPE |
| **default `max_MAC_for_ER=4`** + `dosage_zerod_MAC_cutoff=1` + `SPAcutoff=0.5` | **SIGFPE** |
| same but `dosage_zerod_MAC_cutoff=10` (default) | fine |

So R's defaults are safe only because the window
`max_MAC_for_ER=4 < dosage_zerod_MAC_cutoff=10` is empty — not because the code
is right. `max_MAC_for_ER ≥ 11` alone is enough to kill it, and
`max_MAC_for_ER=14` dies a different way:
`vector::_M_range_check: __n (which is 238) >= this->size() (which is 238)`.
R itself half-admits the limit: `ER_binary_func.cpp:25` has
`int ngroup1 = 10; //use the default value as ER will only be used for variants with MAC <= 10`
and the run prints `WARNING: Efficient resampling may not work well for MAC > 4!`.

**Minimal reproduction** (one marker, R only, no C++ needed):

```bash
# /opt/saige/logs/parity_sb/ermac20dbg/{bisect,cause}.sh are the recorded runs
Rscript $SAIGE/extdata/step2_SPAtests.R --bedFile=rare.bed ... \
        --impute_method=mean --max_MAC_for_ER=11 --idstoIncludeFile=<just 1:489:A:B>
# -> Floating point exception (core dumped), exit 136
```

**Why (b):** a crash is not a result. The C++ uses an iterative `C(n,r)` that
saturates at `INT_MAX` instead of wrapping (`er_binary.cpp:747-766`) and a
`long long` accumulator for `n_total` (`:789-801`). Saturating is safe here
precisely because it cannot be mistaken for a small number: an over-large
`n_total` trips the `n_total > NResampling` branch, which is the branch such a
`k` belongs in anyway, so the marker is answered by resampling rather than by a
truncated exact enumeration. For `1:489:A:B` that yields BETA=1.80665,
SE=0.727495, p=1.301412E-02, against the non-ER SPA value 9.055404E-03. At
`max_MAC_for_ER=20` on the rare set, R dies 45 s in and the C++ completes all
3000 rows in 234 s.

Note the consequence for C5: on these markers the C++ is already **in** the
resampling branch, so the RNG difference is live there — it simply cannot be
compared against R, because R does not survive to produce a number.

### B4 · `is_noadjCov=TRUE` with AF > 0.5 — centring inconsistent with the flip

Carried over from the previous round and re-confirmed here. R 1.5.x's
`scoreTestFast_noadjCov` centres on `2·altFreq`, but at AF > 0.5 the genotype has
already been flipped to `2−g` and `altFreq` flipped back
(`SAIGE-upstream/src/UTIL.cpp:75-122`), so the two are inconsistent.

Evidence that does not reference either implementation's opinion: a
**flip-invariance** test. Re-coding a marker's alleles must not change its score
test. With `is_noadjCov=TRUE`, 50 of 50 markers violate invariance; with `FALSE`,
0 of 50. The variance is inflated by up to 35× (MAF 0.05) and p-values are up to
10.5 orders of magnitude too conservative. Introduced upstream by commit
`ff650722` (2025-06-20), first shipped in 1.4.7/1.4.9; 1.1.3/1.3.0 do not have
the option.

Confirmed on quantitative traits too: with both sides set TRUE, 2000 rows are
bit-identical; turning the flag on within the C++ alone moves p by up to 2.20
orders of magnitude and var by up to 1038.7× (snp1127, AF=0.99806), and 1022 of
2000 markers sit past the flip point.

**Minimal reproduction.** Take any marker with AF > 0.5, swap its two alleles
in the `.bim` (which must not change the score test), and run both codings with
`--is_noadjCov=TRUE`: the two disagree. Repeat with `--is_noadjCov=FALSE`: they
agree. 50/50 markers violate invariance in the first case, 0/50 in the second.

**Disposition.** The port reproduces R's arithmetic — defect included — exactly
when the flag is set, so this is *not* a divergence in the maths. The divergence
is the **default**: R's CLI defaults `is_noadjCov` on, and `tools/rda_to_arma.R:93`
writes `isnoadjCov: false`. That default stays, and this is the one documented
exception to using 1.5.2 as the reference.

### B5 · R reads a truncated VCF `DS` sub-field as dosage 0, not as missing

`plink2 --export vcf vcf-dosage=DS` writes missing samples as a bare `./.` —
FORMAT is `GT:DS` and the trailing DS is omitted entirely, which VCF 4.3 §1.6.2
explicitly permits.

R `src/VCF.cpp:215-216` fills the dosage vector with zeros and then `:233` walks
only the elements savvy actually stored. savvy's `compressed_vector::assign`
(`include/savvy/compressed_vector.hpp:204-219`, comment: *"Allows NaN values but
not Zero values"*) stores only non-zero values, and a truncated sub-field never
produces a value at all — so the sample keeps its initialised 0 and the
`std::isnan(*dose_it)` test at `:250` never sees it. The C++ goes through
htslib's `bcf_get_format_float`, which returns `bcf_float_vector_end` (a NaN) at
the truncation point, caught at `genotype_reader.cpp:1335`.

Measured on `dos.vcf.gz`: snp0 genuinely lacks DS for 508 of 50000 samples. The
C++ reports `MissingRate = 0.01016` = 508/50000; **R reports 0** and an AC equal
to the sum of the observed dosages, i.e. it counts those 508 as 0. 235 markers
with true missingness 0.71–0.75 are correctly dropped by the C++ at
`maxMissing=0.15`; R reports missingness 0 and tests all of them. Across the 265
shared rows, 12 columns differ, `|Δlog10 p|` up to **0.281**, `var` up to 3.3e-2
relative, and BETA sign flips (Tstat up to 1.93 relative).

**Minimal reproduction.**

```bash
plink2 --dummy 50000 500 dosage-freq=1 --export vcf vcf-dosage=DS --out dos
# dos.vcf.gz now writes missing samples as a bare "./." under FORMAT GT:DS
cd /opt/saige/data/readers && python3 readers_parity.py rd_dos_vcfds_mean
# R: MissingRate 0 on every marker.  cpp: the true rate.  Row sets disagree.
```

Three independent confirmations:

1. **Third tool.** `plink2 --vcf dos.vcf.gz dosage=DS --freq` gives snp4
   `OBS_CT = 29056 = 2 × (50000 − 35472)`, agreeing with the C++ exactly. R's
   50000 is wrong.
2. **R contradicts itself.** A hand-built VCF writing every missing sample as an
   explicit `./.:.` (`rd_probe`) makes both sides report `MissingRate` 0.01 and
   the outputs IDENTICAL. So R's own handling of a spelled-out missing DS is
   correct; only the truncated spelling is mishandled.
3. **Isolation.** Rewriting every bare `./.` in `dos.vcf.gz` to `./.:.`
   (`rd_dosx_vcfds_mean`) restores the row set and leaves a residue of 0–8
   samples per marker (mean 2.9) — plink2 also writes a few bare hard calls
   (`0/1`, `1/1`, no DS), the same truncation mechanism, not a second defect.

**Why (b):** reporting `MissingRate = 0` is a factual misstatement about the
file, it contradicts both plink2 and R's own handling of `./.:.`, and it
substitutes a specific wrong value (0) for an unknown one. Reach is wide —
`plink2 --export vcf vcf-dosage=DS` is a common way to produce imputed dosage
VCFs.

---

## 4. (c) — acceptable differences

### C1 · BURDEN-only region output carries 6 significant digits

`SAIGE_SPATest_Region.R:1142` gates `fwrite` behind
`if (regionTestType != "BURDEN")`, so the burden-only region output is written by
upstream's own C++ `OutFile <<` (`src/Main.cpp:2511+`), where `std::ofstream`'s
default precision is 6 significant digits. Strictly verified: of the 2400 numeric
cells in `region_binary_burden`, **1080 are bit-identical to the C++ value and
1110 equal the C++ value rounded to 6 significant figures — 0 unexplained.**
Bound ≤4.9e-6 relative. This is printing, not arithmetic; the C++ keeps 15
digits. Affects `region_binary_burden`, `region_quant_burden`,
`region_binary_nosingle`, `region_binary_mingroupmac`.

### C2 · `max_MAF` column spelling

R writes `1e-04`, the C++ writes `0.0001`. Same value. Applies to the main region
output and to the new `.markerList.txt`. The comparison harness joins on the
numeric value and reports the spelling separately.

### C3 · Conditioning a marker on itself

The conditioning marker's own row: `Tstat_c ≈ 1e-15`, `var_c ≈ 5e-13` — the
residue of subtracting two equal large numbers, against a `var` of order 1.9e4
(relative 3e-15, i.e. double-precision noise). `BETA_c = Tstat_c/var_c` is then
noise over noise and lands anywhere: R 0.00142 / cpp 0.00469, `SE_c` 1.26e6 /
1.33e6. **`p.value_c` is 1.000000E+00 on both sides**, the formula is the same,
and the difference is summation order alone. Both sides correctly identify the
row as degenerate; there is no statistical content to disagree about.

### C4 · Floating-point noise in the region path

Maximum relative difference per column, measured across the nine
non-burden region cases after the A1 fix (no column exceeds the harness's 1e-6
tolerance anywhere, i.e. every one of these cases is IDENTICAL):

| column | worst relative | where |
|---|---|---|
| `Pvalue` (SKAT-O) | 2.2e-8 | `region_quant_mid`; most cases ≤5e-13 |
| `Pvalue_Burden` | 6.3e-11 | `region_binary_skato` |
| `Pvalue_SKAT` | 6.3e-11 | `region_binary_skato` |
| `SE_Burden` | 1.9e-7 | `region_binary_skato` |
| `BETA_Burden` | 1.9e-10 | `region_binary_skato` |
| `MAC` | 1.6e-14 | `region_binary_mid` |
| `MAC_case`, `MAC_control` | **0** | every case |
| `Number_rare`, `Number_ultra_rare` | **0** | every case |

Summation order and LAPACK-level differences. The two `Number_*` columns being
exactly equal everywhere means the ultra-rare collapse boundary decisions agree
completely — no marker is classified differently by the two sides.

### C5 · ER resampling RNG — a difference that must exist but cannot be measured yet

The resampling branch is entered only when `n_total = 2^k > NResampling (2e6)`,
i.e. `k ≥ 21`. R uses R's own Mersenne-Twister via `GetRNGstate()`
(`Binary_resampling.cpp:49-55`); the C++ uses `std::mt19937` reseeded per marker
from a splitmix64 stream (`er_binary.cpp:60-90`, deliberately so the result does
not depend on which thread the marker landed on). The two will not agree.
The C++ reaches this branch today (see B3 — a large `k` saturates `n_total` and
routes the marker to resampling). R does not: it hits B3's SIGFPE first. So the
difference is real and live on our side and simply has no R counterpart to be
measured against. Filed as (c), and flagged as the difference that will surface
the moment B3 is fixed upstream. If bit-comparability there ever matters, it
would require porting R's Mersenne-Twister stream and its per-marker seeding
order, which would in turn give up the current property that a marker's ER
result does not depend on which thread it landed on
(`er_binary.cpp:60-90`).

---

## 5. Defaults

Three knobs defaulted differently on the two sides. Two were corrected, one was
deliberately not, and the distinction is worth stating because none of them ever
showed up as a *measured* difference — every parity case pins these keys
explicitly on both sides, so they were invisible to the comparison and only a
source-level reading found them.

| knob | R default | old C++ default | now | why |
|---|---|---|---|---|
| `AlleleOrder` (bgen/pgen) | `ref-first`, nothing else accepted | `alt-first` globally | **changed to R's**, and non-`ref-first` refused | the BGEN spec defines the first allele as the reference — `alt-first` is wrong for the format, not merely different (A5) |
| `vcfField` | `DS` | `GT` | **changed to R's**, non-DS/GT refused | reads a different FORMAT field, so it analyses different genotypes: `MissingRate` 0.80 vs 0.01 on the same file (A10) |
| `impute_method` | `best_guess` | `mean` | **left as `mean`** | see below |

**Why `impute_method` is treated differently.** `best_guess` and `mean` are both
legitimate; neither contradicts a format, a spec, or R's own behaviour elsewhere.
R is not doing anything wrong by defaulting to `best_guess`, so there is no
correctness argument for a change — only a compatibility one, and changing it
would silently move the results of every existing user of this tool. It is also
already an explicit, documented choice here: the README's reference protocol
specifies `--impute_method=mean`, and `tools/rda_to_arma.R` writes
`"impute_method": "mean"` into every converted model, so in the normal workflow
the default is never reached. Recorded as a maintainer decision rather than
taken unilaterally.

Measured effect of the choice, for whoever makes that decision: on the rare set,
`mean` vs `best_guess` changes 819 of 3000 single-variant rows (identically on
both sides), and 11 columns differ with `BETA` moving by up to 1.94 relative.
The two sides agree exactly under either setting.

---

## 6. Input-validation asymmetries (no numerical divergence)

Configurations R rejects at startup and the C++ runs. These produce no cpp↔R
numerical difference — R never starts, so there is nothing to compare — but they
mean the C++ accepts inputs R defines as illegal. Whether to add range checks is
a product decision, not a parity one.

| knob | R | C++ |
|---|---|---|
| `markers_per_chunk < 1000` | `checkArgs.R` halts: "should be a numeric value greater than or equal to 1000" | runs (`marker_chunksize: 137` completed) |
| `dosage_zerod_cutoff > 0.5` | halts | runs (1.5 accepted) |
| `dosage_zerod_MAC_cutoff > 100` | halts | runs (1e6 accepted; snp0's AC moved 69479.4 → 90966 and 5 markers dropped out of the filter) |
| `SPAcutoff < 0.5` | halts | runs |

`AlleleOrder` and `vcfField` used to be on this list; both are now enforced
(A5, A10).

---

## 7. Coverage — what was compared, and what was not

### Suite state at the end of this round

Every registered case re-run against the current build (R 1.5.2 vs C++,
`nThreads=1`, one shared step-1 fit per model):

| group | cases | IDENTICAL | differing, and why |
|---|---|---|---|
| `tests/parity/parity.py` (single + region) | 22 | 17 | 4 burden-only (**C1**, ≤4.9e-6 printing), 1 `region_binary_flatweights` (**B2**) |
| readers overlay | 38 | 31 | 4 `rd_dos_vcfds_*` + 1 `rd_dosx_vcfds_mean` (**B5**); 2 not comparable by design — `rd_bgen_altfirst` (both sides now refuse, **A5**) and `rd_pgen_mean` (R reads mode 0x10, the C++ refuses, **A6**) |
| binary overlay incl. LOCO | spot-checked 8 | 8 | — |
| conditional analysis | 4 | 0 | every differing row accounted for by **B1**, **A7** or **C3**; no unexplained row |

Region `.singleAssoc.txt` output is bit-identical across every region case.
The four conditional cases are the only place where differing rows remain by
design, and each row is attributed: the marker at genoIndex 0 (B1), markers on
R's uninitialised path (A7), and the self-conditioned 0/0 row (C3).

**Compared and IDENTICAL** (R 1.5.2 vs C++, `nThreads=1`, one shared step-1 fit):

- *binary single-variant*: base, minMAF, minMAC, maxMissing (incl. the exact-tie
  boundary), moreDetails, noadjCov (both sides TRUE — bit-identical over 2000
  markers), SPAcutoff=1, best_guess / mean / minor imputation, chunk size, ER at
  `max_MAC_for_ER` ∈ {0, 4, 10}, Firth at the default and at a loose cutoff,
  Firth on the rare set, categorical VR (2 bins, both bins exercised),
  `is_fastTest` both ways, LOCO chr1–5 plus LOCO+moreDetails and LOCO+Firth.
- *quantitative single-variant*: 20 cases — filters, all three imputation
  methods, AlleleOrder, moreDetails, `is_imputed_data`, dosage zeroing, chunking,
  SPAcutoff, ER, Firth, fastTest, noadjCov, plus the rare set. Confirmed by
  source inspection on both sides that quantitative does not reach SPA
  (`SAIGE_test.cpp:545,556` ↔ `saige_test.cpp:957,967`), ER
  (`Main.cpp:492,1324` ↔ `main.cpp:493,1085,1279`) or Firth
  (`SAIGE_test.cpp:704,710` ↔ `saige_test.cpp:1081,1088`).
- *region*: SKAT-O and burden, binary and quantitative, collapse cutoffs
  0/10/20, maxMAF 1e-4/1e-3/1e-2, group-file weights, `minGroupMAC`, chunk sizes
  100 and 500, `is_single_in_groupTest=FALSE`, markerList, and both group-file ID
  spellings. `.singleAssoc.txt` output — 12000+ rows across all cases — is
  bit-identical throughout.
- *readers*: plink / bgen / vcf-GT / vcf-DS × mean / best_guess / minor
  imputation; maxMissing including an exact tie; `AlleleOrder`; fractional-dosage
  bgen and VCF; the `R2` INFO path with `minInfo`; dosage zeroing; a hand-built
  probe VCF with controlled missingness patterns; HDS; PGEN mode 0x01.

**Not covered — stated plainly:**

- **survival** — out of scope, not ported.
- **sparse GRM** in any form — out of scope. This includes the genuine
  two-pass `is_fastTest` path, which needs a sparse GRM; what was exercised is
  the second pass triggered by `is_noadjCov=TRUE`, which is a code-path
  comparison rather than a statistically meaningful one.
- **`pval_cutoff_for_fastTest` at a non-default value** — R's CLI has no such
  flag (it exists only in the `SPAGMMATtest()` signature), so from the command
  line it is always 0.05; the C++ reads it only from `nullmodel.json`, default
  also 0.05. The defaults match; a non-default comparison cannot be constructed.
- **ER's resampling branch** — see C5.
- **region + `--condition`** — differences there are expected to be B1 + A7 and
  were not independently re-attributed.
- **region × bgen / vcf / pgen** — only plink was exercised for region tests;
  see A9's closing note on the coordinate-only marker maps.
- **LOCO with a `chrom` outside `loco_chroms`** — R's silent fallback
  (`readInGLMM.R:107-113`) is untested (the test data has only chromosomes 1–5,
  and chrom=6 fails earlier for having no markers).
- **LOCO with per-chromosome offsets** — `tools/rda_to_arma_loco.R` defaults to
  writing the genome-wide offset, which is what R uses in non-Firth runs and for
  every model from SAIGE ≤1.3.3. Its `--loco-offset` option (per-chromosome) is
  untested. Note the asymmetry: R swaps the offset per chromosome only when
  `is_Firth_beta=TRUE` and `LOCOResult[[j]]$offset` exists
  (`readInGLMM.R:97-101`), whereas a C++ model directory holds exactly one.
- **PGEN modes 0x10/0x11** — refused by the C++, read by R. See A6.
- **`vcfFilters`, `idstoIncludeFile`, `rangestoIncludeFile`, `subSampleFile`,
  `maxMAC_in_groupTest`, `is_no_weight_in_groupTest`** — not ported; the harness
  refuses to run cases that set them.
- **multi-allelic sites, chrX/PAR, sav/bcf containers, BGEN layout 1, BGEN with
  zstd, BGEN with embedded sample IDs and no `.sample` file** — none exercised.
  (BGEN with `B != 8` bits is refused identically by both sides —
  `SAIGE-upstream/src/BGEN.cpp:238` ↔ `genotype_reader.cpp:1858/2233` — verified
  at `bits=16`.)
- **`nThreads > 1` on the R side** — only the C++ was run multi-threaded.
  Re-checked against the current build: the C++ at `nThreads=8` reproduces R's
  single-threaded single-variant output with **0 differing cells** over
  2000 rows × 19 columns, and the new `.markerList.txt` at `nThreads=4` is
  row-for-row identical (sorted) to the `nThreads=1` file — 180 rows, no
  malformed rows, no duplicate keys.
- **Scale** — 50000 samples, ≤40000 markers, 30–100 genes. Whether A1's trigger
  rate (11–30% of mask rows in this data) varies with gene or mask size is
  untested.

**Downgraded from the risk list:** the categorical-VR top-bin boundary
(R appends `nsample`, the C++ uses 1e10) is **unreachable**, since
`MAC = min(AC, 2N−AC) ≤ N = nsample`.

---

## 8. Reproducing

```bash
source /home/francisfenglu4/miniforge3/etc/profile.d/conda.sh
conda activate saige-build
cd step2_saige-step2 && make -j4

cd ../tests/parity
./parity.py --list                  # registered cases
./parity.py region_binary_skato     # run one, both sides, and diff
./parity.py base_binary_single --set condition=snp100,snp500
```

R-side runs use `conda activate RSAIGE_GPU` with
`R_LIBS_USER=/opt/saige/Rlib-upstream`; `parity.py` handles that. `knobs.py`
holds the knob inventory and the R-flag mapping; `compare_out.py` does the
column-wise diff (it joins `max_MAF` numerically and counts Inf-vs-finite
mismatches explicitly, which a bare relative test silently passes).

## 9. Addendum 2026-09-15 — three divergences found after the sparse pipe was connected

All three were C++ implementation errors and are fixed; each fix was checked
against R SAIGE 1.5.2 on identical inputs. None of them was reachable by the
comparisons in §7: mid has no marker with MAF < 0.1 (so step 2 never took the
sparse-Sigma path), and every step-1 configuration there had tau1 well above tol.

| # | Where | What C++ did | What R does | Fix |
|---|---|---|---|---|
| S1 | step 2, sparse Sigma (`fdf1c7a`) | used the sparse GRM K itself as Sigma | builds Sigma = tau1*K + diag(1/mu2) (binary) or tau1*K + tau0*I (quant) in step 2 (`setSparseSigma_new`) | build Sigma in `null_model_loader.cpp`; per-marker identical to R after the fix |
| S2 | step 1, binary AI-REML (`6bfb935`) | kept a tau1 step that landed below tol and halved negative steps | `fitglmmaiRPCG` (`SAIGE_fitGLMM_fast.cpp:5435`) sets tau < tol to 0 before its `while(tau<0)` halving, then `SAIGE_fitGLMM_fast.R:385` breaks on tau[2]==0 | zero tau1 below tol, max(0,.) on the conservative first step, no halving; single-trait and lockstep paths |
| S3 | step 2, sparse GRM (`9535a41`) | kept every K entry | drops K entries with x < relatednessCutoff after subsetting, before scaling (`SAIGE_SPATest_Region.R:65`); step-2 CLI default 0, so negative kinships are dropped | same filter, new step-2 YAML key `relatednessCutoff` (default 0), applied per model |

**Impact.**
- S1: every marker on the sparse-Sigma path had its score variance inflated
  (4.7x binary, 1.2-1.9x quant), pushing p-values toward 1: all markers when
  isFastTest is false; with isFastTest true, markers with first-pass p < 0.05 and
  4 < MAC <= 20.5; in region tests every marker with MAC <= 20.5.
- S2: **not specific to sparse**. Any binary fit (dense CPU/GPU, sparse,
  lockstep) whose AI-REML step lands below tol=0.02 -- traits with weak or no
  genetic signal -- ended at a small positive tau1 where R returns 0. Examples:
  y1 0.00327 vs 0; a weak simulated trait 0.0172 vs 0; mu differed by up to
  5.5e-3 before the fix, 1.8e-7 after. Traits with clearly positive tau1 were
  already identical and are byte-identical before/after the fix.
- S3: with negative kinships in the sparse GRM, variance differed by up to 0.10%
  and p-values by up to 0.2 log10 units; after the fix out.txt is byte-identical
  to R (cutoff 0 and 0.3).

**Verification.** Step-1 P=1 byte gate vs `fdf1c7a` (10 configs: dense binary
GPU, dense quant CPU, LOCO binary CPU+GPU, LOCO quant CPU, sparse direct x3,
lockstep CPU+GPU) all IDENTICAL -- none reach tau1 < tol. Step 2:
run_p1_regression 11/0, run_mt_correctness 147/0, run_mt_subset_tests 184/0.

**Left as is, deliberately.**
- R treats the GRM cutoff three different ways: the step-1 fit uses
  `drop0(tol=cutoff)` (drops |x| <= cutoff, keeps large negatives), step-1 VR drops
  nothing, step 2 drops x < cutoff. C++ step 1 ignores the cutoff when reading a
  GRM file (its `relatedness_cutoff`, default 0.05, only applies when C++ builds
  the GRM). At R's default cutoff 0 the numbers agree; they differ only when a
  user passes a non-zero cutoff to step 1 with a GRM containing small or negative
  entries. Because R is internally inconsistent here, C++ does not copy the
  step-1 behaviour; revisit if a user needs that exact combination.
- If the step-2 cutoff exceeds a diagonal value, R (binary) recycles 1/W over the
  surviving diagonal entries out of position (read from code, not run). C++
  warns instead of reproducing the misalignment.
