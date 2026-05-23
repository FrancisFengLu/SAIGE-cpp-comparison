# SAIGE standalone C++ — source snapshot (2026-05-24)

Validated standalone C++ port of SAIGE Step-1 (null model) and Step-2 (association
testing), developed to reproduce the R SAIGE package. Source-only snapshot copied
2026-05-24 from the two working trees:

- `step1_saige-null/`  ← `seokho92-SAIGE-cpp-comparison/code_copy/cpp_standalone`  (binary: `saige-null`)
- `step2_saige-step2/` ← `saige_step2_cpp_repo/SAIGE-step2-cpp/code_copy/cpp_standalone` (binary: `saige-step2`)

Build artifacts (`*.o`, `*.d`), compiled binaries, frozen snapshots, and run
outputs are intentionally excluded. The R reference package is not included.

## Key fixes baked into this snapshot

**Step-1 (`saige-null`)**
- **Sparse-GRM sample-ordering fix** (`preprocess_engine.cpp`): for `use_sparse_grm_to_fit`,
  the working vectors are kept in the design's (FAM-ascending) order instead of being
  reordered to the sparse-GRM sampleID-file order, so the GLMM working vector `PY` stays
  aligned with the sparse GRM `K`. Previously this misalignment scrambled the genetic
  quadratic form and collapsed `tau[1]` to 0.
- **GMMAT zero-the-boundary trick** in `quant_glmm_solver` (`glmm.cpp`), ported from R's
  `fitglmmaiRPCG_q`, so a boundary variance component does not freeze the joint AI-REML step.

**Step-2 (`saige-step2`)**
- **`.bgi` byte-offset region-read fix** (`genotype_reader.cpp`): `populateFromBgi` now reads
  `file_start_position`; `getMarkerIDToIndex`/`getMarkerNameToIndex` return the byte offset
  (not the sequential row index); `getOneMarker` uses an absolute `fseek`. This fixes the
  region/by-ID (SAIGE-GENE+) reads, which previously sought to the wrong file position and
  crashed with `inflate failed` / segfault. (Single-variant mode streams and was unaffected.)

## Validation status (vs R SAIGE)
- Step-1 nulls match R: LDL (quant) `tau` matches to ~5 decimals; T2D (binary) `tau[1]`
  matches within MC-trace noise.
- Step-2 single-variant: bit-exact cpp↔R.
- Step-2 SAIGE-GENE+: region p-values match R for both LDL (LDLR) and T2D (TCF7L2)
  across all annotation × maxMAF groups and the Cauchy omnibus.

## Known limitation
SAIGE-GENE+ region tests have a large transient peak RSS (~N×N scale; ~240 GB at ~165k
samples). It is a memory-efficiency issue, not a correctness one — runs complete on a
high-memory node (~250 GB+). Not yet optimized.

## Build (Linux, glibc; mamba env `saige-build`)
```bash
ENV=/data/home/seokhojeong/.local/share/mamba/envs/saige-build
export CONDA_PREFIX=$ENV PATH=$ENV/bin:$PATH
export LD_LIBRARY_PATH=$ENV/lib:$ENV/lib/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$ENV/lib CPATH=$ENV/include PKG_CONFIG_PATH=$ENV/lib/pkgconfig
export CXX=g++ CXXFLAGS="-std=c++17 -O3 -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP"
cd step1_saige-null  && make clean && make    # -> saige-null
cd ../step2_saige-step2 && make clean && make  # -> saige-step2
```
Always do a full `make clean && make` (incremental builds can leave stale `.o` → ABI bugs).
`-DARMA_64BIT_WORD=1` is required (whole-cohort sparse GRM is N×N with N>46k → index >2^32).

## Run
Both binaries take a YAML config: `./saige-null -c config.yaml` and `./saige-step2 config.yaml`.
See the `config_*.yaml` examples in each directory.
