---
title: "Step 2: test genetic variants"
nav_order: 3
---

# Step 2: test genetic variants

`saige-step2` tests every marker of a genotype file against one or more step-1
models (single-variant score test, with SPA, Firth and exact-test corrections for
binary traits). All settings are in a YAML config; the only argument is the
config file.

```bash
saige-step2 step2.yaml > step2.log 2>&1
saige-step2                # no argument: prints every config key with its default
```

## Important inputs

**Models.** A `models:` list, one entry per trait:

- `traitName` — label used in the log (default: the model directory name).
- `modelFile` — step-1 model directory (contains `nullmodel.json`).
- `varianceRatioFile` — step-1 `.varianceRatio.txt`.
- `outputFile` — result file for this trait. **The directory must already exist.**

For one trait the three keys can also be written at the top level instead of
`models:` (`modelFile`, `varianceRatioFile`, `outputFile`); the two forms are
mutually exclusive. Binary and quantitative models can be mixed in one run.
Models fitted on different sample sets (e.g. different missing phenotypes) can
be mixed with PLINK or hard-call PGEN input; each trait is tested on its own
samples.

**Genotypes:**

- `genoType` — `plink` (default), `pgen`, `bgen` or `vcf`.
- `plinkFile` — PLINK prefix (`.bed/.bim/.fam`).
- `pgenFile`, `pvarFile`, `psamFile` — PLINK 2 files (hard calls or dosages).
- `bgenFile`, `bgenSampleFile` — BGEN and its `.sample` file. A `.bgi` index next to the `.bgen` is needed with `LOCO: true`.
- `vcfFile` — VCF / VCF.GZ / BCF; `vcfField` — `DS` (default) or `GT`.
- `AlleleOrder` — `alt-first` (default for PLINK: Allele2 in the output is the `.bim` A1 allele) or `ref-first`. PGEN and BGEN require `ref-first` (the default for them); any other value stops the run.

Genotype samples are matched to the model's samples by ID; genotype samples
not in the model are ignored.

**Marker filters:**

- `minMAF` — default 0.
- `minMAC` — default 0.5.
- `maxMissRate` — default 0.15.
- `minINFO` — imputation INFO filter (default 0).
- `dosage_zerod_cutoff` (0.2), `dosage_zerod_MAC_cutoff` (10) — dosage inputs: dosages ≤ cutoff are set to 0 for markers with MAC ≤ the MAC cutoff.

**Tests (binary traits):**

- `is_Firth_beta` — Firth-corrected effect for markers with p < `pCutoffforFirth`. Default: the value stored by step 1 (`fit.firth_beta`). Can be set per model inside a `models:` entry.
- `pCutoffforFirth` — default: the value stored by step 1 (`fit.p_cutoff_for_firth`, 0.01).
- `MACCutoffforER` — markers with MAC ≤ this use the exact test (default 4).
- `isFirth` — only controls the Firth summary lines in the log; whether Firth is applied is decided by `is_Firth_beta`.

The SPA cutoff (`fit.spa_cutoff`) and the fast test (`fit.fast_test`) are set in
step 1 and read from the model.

**Run:**

- `nThreads` — threads (default 1).
- `useGPU` — use the GPU (default `false`). See [GPU](gpu.md) for the sub-switches and what is accepted.
- `outputFormat` — `text` (default) or `sgs` (binary, see [below](#binary-output-sgs)).
- `isMoreOutput` — extra columns (binary: hom/het counts in cases and controls). Default `false`.
- `LOCO`, `chrom` — see [LOCO](#loco).
- `mtRequireSameSamples` — `true` stops the run when the models' sample lists differ (default `false`).

## Example config: four binary traits, GPU

From [`examples/04_step2_binary.sh`](examples/04_step2_binary.sh); the models are
those of [Step 1](step1.md#example-config-four-binary-traits-full-grm-gpu):

```yaml
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first               # PLINK .bed: Allele2 in the output = the .bim A1 column
minMAF: 0
minMAC: 1
maxMissRate: 0.15
LOCO: false
isFirth: true                        # binary traits: Firth correction for p < pCutoffforFirth
pCutoffforFirth: 0.01
MACCutoffforER: 4                    # binary traits: exact test for MAC <= 4
nThreads: 8
useGPU: true                         # every GPU sub-switch defaults to on
outputFormat: text
models:
  - traitName: b1
    modelFile: $M/models/b1
    varianceRatioFile: $M/vr_b1.varianceRatio.txt
    outputFile: $O/out/b1.txt
  - traitName: b2
    modelFile: $M/models/b2
    varianceRatioFile: $M/vr_b2.varianceRatio.txt
    outputFile: $O/out/b2.txt
  # ... b3, b4 the same way
```

```bash
mkdir -p $O/out
saige-step2 step2.yaml > step2.log 2>&1
```

With `useGPU: false` (or the CPU build) the output files are byte-identical.

## Example: other genotype formats

From [`examples/09_step2_formats.sh`](examples/09_step2_formats.sh) — the same
models, only the genotype keys change:

```yaml
# hard-call PGEN (runs on the GPU)
genoType: pgen
pgenFile: $D/geno.pgen
pvarFile: $D/geno.pvar
psamFile: $D/geno.psam
AlleleOrder: ref-first
```

```yaml
# BGEN (CPU)
genoType: bgen
bgenFile: $D/geno.bgen
bgenSampleFile: $D/geno.sample
AlleleOrder: ref-first
```

```yaml
# VCF (CPU)
genoType: vcf
vcfFile: $D/geno.vcf.gz
vcfField: DS
```

Results from the hard-call PGEN are byte-identical to the PLINK results. BGEN gives
the same p-values with Allele1/Allele2 (and the sign of BETA) swapped, because BGEN
is read ref-first. With VCF input and `nThreads` > 1 the output rows are not in
file order (they are with `nThreads: 1`); sort by CHR and POS if you need order.

BGEN, VCF and dosage PGEN need all models to share one sample list.

## LOCO

Models fitted with `fit.loco: true` are tested one chromosome per run.
From [`examples/06_step2_quant_loco.sh`](examples/06_step2_quant_loco.sh):

```yaml
genoType: plink
plinkFile: $D/geno
AlleleOrder: alt-first
minMAF: 0.01
minMAC: 1
LOCO: true
chrom: "1"                           # only this chromosome's markers are tested
nThreads: 8
useGPU: true
models:
  - traitName: q1
    modelFile: $M/models/q1
    varianceRatioFile: $M/vr_q1.varianceRatio.txt
    outputFile: $O/chr1/q1.txt
  - traitName: q2
    modelFile: $M/models/q2
    varianceRatioFile: $M/vr_q2.varianceRatio.txt
    outputFile: $O/chr1/q2.txt
```

The log shows `LOCO: restricting to chromosome 1: 2500 of 5000 markers retained.`
`chrom` must match the genotype file's chromosome codes. `LOCO: false` (default)
uses the whole-genome model and tests every marker. See the
[HPC example](hpc_example.md) for running all chromosomes and merging.

## Binary output (sgs)

`outputFormat: sgs` writes, instead of text, `<outputFile>.sgs` for each trait plus
one shared `<first outputFile>.markers.sgs`. It needs at least two models, or
one model with `useGPU: true` (with or without a device); a one-model run with
`useGPU: false` and region tests stop with an error. `sgsPrecision: fp64` (default) converts back to exactly the text
output; `fp32` halves the size but is not exact.

Convert with `tools/sgs2txt`, on any machine (no GPU needed). From
[`examples/08_step2_sgs.sh`](examples/08_step2_sgs.sh), after copying the files
to another directory:

```bash
for t in b1 b2 b3 b4; do
  sgs2txt -m moved/b1.txt.markers.sgs -o text/$t.txt moved/$t.txt.sgs
done
```

| Option | Meaning |
|---|---|
| `-o FILE` | text file to write (one input `.sgs`) |
| `-m FILE` | the `.markers.sgs` file (needed when the files were moved) |
| `-j N` | convert N traits at a time |

Without `-o` the text is written to the path recorded in the `.sgs` header (the
original `outputFile`).

## Results

One tab-separated file per trait, one row per tested marker.

Binary trait:

| CHR | POS | MarkerID | Allele1 | Allele2 | AC_Allele2 | AF_Allele2 | MissingRate | BETA | SE | Tstat | var | p.value | p.value.NA | Is.SPA | AF_case | AF_ctrl | N_case | N_ctrl |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1000 | snp1 | A | C | 7921.8 | 0.79218 | 0.0102 | -0.0482665 | 0.0828339 | -7.03443 | 145.742 | 5.601024E-01 | 5.601024E-01 | false | 0.784387 | 0.793063 | 509 | 4491 |

Quantitative trait:

| CHR | POS | MarkerID | Allele1 | Allele2 | AC_Allele2 | AF_Allele2 | MissingRate | BETA | SE | Tstat | var | p.value | N |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1000 | snp1 | A | C | 7921.8 | 0.79218 | 0.0102 | 0.0431276 | 0.0242652 | 73.2466 | 1698.37 | 7.551185E-02 | 5000 |

- `Allele2` is the tested allele; `BETA` is its effect (log-odds for binary traits).
- `p.value` is the p-value to use. For binary traits it is SPA-corrected when `Is.SPA` is `true`, and `BETA`/`SE` are Firth-corrected when p < `pCutoffforFirth`; `p.value.NA` is the uncorrected p-value.
- `isMoreOutput: true` adds `N_case_hom N_case_het N_ctrl_hom N_ctrl_het` (binary).

## Not covered here

Region / gene-based tests (`groupFile`, `annotationList`, `maxMAFList`; one
trait per run, CPU only), conditional analysis (`condition`) and LD matrices
(`isLDMatrix`) exist in the same binary; `saige-step2` without arguments lists
their keys. They were not run for this guide.
