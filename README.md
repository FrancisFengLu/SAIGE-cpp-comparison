# SAIGE standalone C++ port

A standalone C++ port of SAIGE Step 1 (null model fitting) and Step 2
(association testing: single-variant and SAIGE-GENE+ region tests), validated
against the R SAIGE package.

```
SAIGE_cpp_260716/
  step1_saige-null/     null model fit        -> binary: saige-null
  step2_saige-step2/    association testing   -> binary: saige-step2
  LOCO_FORMAT.md        step1 -> step2 on-disk contract for LOCO
  tests/                end-to-end test suite
code_copy/SAIGE_isolated/   R SAIGE, kept as the reference implementation
```

## Build

Requires a C++17 compiler plus Armadillo, OpenBLAS, LAPACK, SuperLU, yaml-cpp,
htslib, zstd, zlib, sqlite3, and R (Step 1 embeds an R runtime).

A mamba/conda environment is the easiest way to get all of them:

```bash
mamba create -n saige-build -c conda-forge gxx_linux-64=12 armadillo openblas \
    superlu yaml-cpp htslib zstd zlib sqlite r-base
conda activate saige-build

# Step 1 embeds R -- without R_HOME the build fails on "R.h: No such file"
export R_HOME=$(Rscript -e 'cat(R.home())')

cd SAIGE_cpp_260716/step1_saige-null  && make -j4     # -> saige-null
cd ../step2_saige-step2               && make -j4     # -> saige-step2
```

Notes:
- `-DARMA_64BIT_WORD=1` is required and is already in both Makefiles: a
  whole-cohort sparse GRM is N x N, and with N > 46k the index exceeds 2^32.
- The binaries are deliberately **not** committed. A stale committed binary
  makes `make` believe it is up to date, so a regression run silently compares
  against the wrong build. Always build from source.

## Run

The two binaries take their config differently:

```bash
./saige-null  -c config.yaml     # Step 1: flag
./saige-step2    config.yaml     # Step 2: positional
```

`./saige-null --help` and `./saige-step2` with no arguments list every option
and config key. See the `config_*.yaml` examples in each directory.

Step 1 writes a **directory** of Armadillo binary matrices plus
`nullmodel.json`. Step 2 reads that directory through its `modelFile` key.

## LOCO

Leave-one-chromosome-out is supported end to end.

Step 1 (`loco: true`) estimates the variance component `tau` once on the
full-genome GRM, then holds `tau` fixed and re-solves only the fixed effects
per chromosome, writing a `chr<N>/` subdirectory per autosome. Exclusion is
subtractive rather than a GRM rebuild. Chromosomes are processed serially
because each one warm-starts from the previous one's coefficients, exactly as
in R.

Step 2 takes `LOCO: true` and `chrom: "<N>"`, loads the model from `chr<N>/`,
and restricts the markers it tests to that chromosome.

The on-disk contract is specified in `SAIGE_cpp_260716/LOCO_FORMAT.md`. Read it
before changing either side.

Supported: quantitative and binary traits; single-variant and region tests on
plink, pgen, bgen and vcf.
Not supported: survival traits. Step 1 logs and skips LOCO rather than emitting
a wrong per-chromosome model.

## Tests

```bash
export R_HOME=$(Rscript -e 'cat(R.home())')

# end to end (step1 + step2), including agreement with a stored R reference
bash SAIGE_cpp_260716/tests/run_loco_e2e.sh

cd SAIGE_cpp_260716/step2_saige-step2
bash tests/run_loco_tests.sh          ./saige-step2
bash tests/run_loco_bgen_vcf_tests.sh ./saige-step2
```

73 assertions total. Each script takes an optional second argument: the path to
a baseline `saige-step2` binary, which enables byte-identity checks against a
previous build. Build one from whichever commit you want to compare against
(`git worktree add`); it is deliberately not committed.

Pin `nthreads: 1` whenever you compare numbers. At 4 threads the float32
reduction order makes `tau` nondeterministic -- 0.2796 versus 0.2336 on the
test data -- and you will chase ghosts.

## Agreement with R

- Step 2 single-variant: bit-exact.
- Step 2 region (SAIGE-GENE+): burden and SKAT agree to ~1e-12. The SKAT-O
  omnibus differs on 4 of 2,853 genes, none near significance.
- Step 1: `mu` agrees to ~9e-3, `tau` to ~11%. **This is not a porting bug.**
  The stochastic trace estimator carries that much noise by itself: the same
  binary yields `tau` 0.2796 at one thread and 0.2336 at four, and R adds its
  own unseeded random vectors. Closing the gap further would mean making the
  trace estimator deterministic and identical on both sides.

## Performance

Measured on 4 cores / 16 GB against R SAIGE 1.5.1. Dataset: N=10,000 samples,
300 genes, 16,422 rare markers, quantitative trait. Both sides were given the
same null model and verified to process identical marker and gene counts.
Median of 3 runs with `OPENBLAS_NUM_THREADS=1`.

| stage | C++ | R | ratio |
|---|---|---|---|
| Step 2 region, 1 thread | 270.8 s | 383.8 s | 1.42x |
| Step 2 region, 4 threads | 116.0 s | 383.8 s | 3.31x |
| Step 2 single-variant, compute only | 2.2 s | 7.0 s | 3.1x |
| **Step 1, 4 threads** | **159 s** | **118 s** | **0.74x (slower)** |

**Step 1 is currently slower than R**, and also burns more CPU to get there
(530 s user versus 418 s). Step 2 region tests spend 40-45% of wall clock in
system time versus R's 23-26%. Both are open items, not settled results.

Full methodology, the BLAS-threading variant of these tables, and the raw
per-run numbers are in
`SAIGE_cpp_260716/step2_saige-step2/tools/BENCH_CPP_VS_R.md`.

## Known issues

- Step 1 is slower than R; see above.
- Survival traits do not support LOCO.
- Region tests force the off-diagonal variance blocks to be symmetric, but `V`
  is genuinely asymmetric because `P2Vec` switches operator at MAC 20.5. This
  reproduces identically in R SAIGE 1.5.1, so the port is faithful. It is an
  upstream design question and has deliberately not been patched.
- A gene whose group-file variants span several chromosomes is a hard error
  under LOCO. R instead trims such a gene silently and reports a burden
  statistic computed over the remaining subset.

## Memory

An earlier version had a large transient peak in region tests. The cause was
`max_markers_region` (default 100,000) being passed to the working-matrix
allocation in place of `markers_per_chunk_in_groupTest` (default 500), so
`P1Mat`/`P2Mat` were sized 100,000 x N per thread. The peak is linear in N, not
N x N:

```
per-thread peak ~ 2 * markers_per_chunk_in_groupTest * N * 8 bytes
```

Fixed, with results byte-identical and a 28-48x reduction. At N=165,000 a
region test now peaks around 1.4 GB per thread instead of being unrunnable.

## Earlier fixes retained from the 2026-05-24 snapshot

**Step 1** -- sparse-GRM sample-ordering fix (`preprocess_engine.cpp`): under
`use_sparse_grm_to_fit` the working vectors stay in the design's FAM-ascending
order rather than being reordered to the sparse-GRM sampleID-file order, so the
GLMM working vector `PY` stays aligned with `K`. The misalignment used to
scramble the genetic quadratic form and collapse `tau[1]` to 0. Also the GMMAT
zero-the-boundary trick in `quant_glmm_solver` (`glmm.cpp`), ported from R's
`fitglmmaiRPCG_q`, so a boundary variance component cannot freeze the joint
AI-REML step.

**Step 2** -- `.bgi` byte-offset region-read fix (`genotype_reader.cpp`):
`populateFromBgi` reads `file_start_position`, the marker-ID lookups return the
byte offset rather than a sequential row index, and `getOneMarker` uses an
absolute `fseek`. Region/by-ID reads previously sought to the wrong file
position and crashed with `inflate failed` or a segfault.
