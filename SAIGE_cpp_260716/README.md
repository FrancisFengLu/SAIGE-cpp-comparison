# SAIGE C++ standalone — snapshot `SAIGE_cpp_260716`

Source-only snapshot of the standalone C++ port of SAIGE Step-1 (null-model
fitting) and Step-2 (association testing, single-variant + SAIGE-GENE+ region).
Prepared for collaborative work on **SAIGE-GENE+** — start with
[`SAIGE_GENE_STATUS.md`](SAIGE_GENE_STATUS.md).

Snapshot date: 2026-07-16. Copied from the active dev tree `SAIGE_cpp_260617`
(source only, no binaries/`.o`/`.d`/artifacts), then two fixes applied and
propagated back to the dev tree + `260617` so all three stay in sync.

## Changes finalized in this snapshot (2026-07-16)
1. **Step-1 convergence flag at the tau boundary** (`step1_saige-null/glmm.cpp`,
   both quant + binary solvers). When the genetic variance component hits the 0
   boundary (`tau[1] ≤ 0`, residual `tau[0] > 0`), the solver now **clamps
   `tau[1]=0` and finalizes as `converged=TRUE`** (matching R's
   `if(tau[1]<=0|tau[2]<=0) break` boundary semantics), instead of running to
   maxiter and reporting `converged=false`. **Values are unchanged** — only the
   flag/iteration count. Verified on a 10k sparse-LDL null: cpp tau=`[1.17462, 0]`
   converged in 2 iters, matching R 1.3.0 tau=`[1.174652, 0]` converged=TRUE
   (previously cpp printed "Converged: NO" / 20 iters on the identical result).
2. **Step-1 build: `-lpcre2-8`** added to `step1_saige-null/Makefile` so `libR`'s
   versioned pcre2 refs (`PCRE2_10.47`) resolve against the env's libpcre2.
Both engines verified building + running the full 10k pipeline from this snapshot.

## Layout
```
SAIGE_cpp_260716/
├── README.md                 # this file
├── SAIGE_GENE_STATUS.md      # SAIGE-GENE+ status, threading, known issues  <- read first
├── step1_saige-null/         # Step-1: builds `saige-null`  (+ tools/, examples/)
└── step2_saige-step2/        # Step-2: builds `saige-step2` (single-variant + region; + examples/)
```

## Build

Both engines build with the `saige-build` mamba env
(`/data/home/seokhojeong/.local/share/mamba/envs/saige-build`). **Always
`make clean && make`** (incremental make leaves stale `.o` → ABI segfaults).

```bash
ENV=/data/home/seokhojeong/.local/share/mamba/envs/saige-build
export PATH=$ENV/bin:$PATH
export LD_LIBRARY_PATH=$ENV/lib:$ENV/lib/R/lib
export LIBRARY_PATH=$ENV/lib
export CPATH=$ENV/include
export PKG_CONFIG_PATH=$ENV/lib/pkgconfig
export CXX=g++
```

**Step-2 (`saige-step2`)** — verified building clean in this snapshot (binary ~1.72 MB):
```bash
cd step2_saige-step2
export CXXFLAGS="-std=c++17 -O3 -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP"
make clean && make -j8
```
(The Makefile uses `?=`, so conda's exported `CXXFLAGS` shadows it unless you set it
as above.)

**Step-1 (`saige-null`)** — do NOT set `CXXFLAGS` (its Makefile force-appends the
ARMA defs). Links CUDA (nvcc) + embedded R; `nvlink … libdl.a` warnings are benign.
Binary ~2.1 MB. **Verified building + running clean in this snapshot** (2026-07-16).
```bash
cd step1_saige-null
make clean && make -j8
```
> **Note (fixed here):** `libR.so` has versioned pcre2 refs (`pcre2_*@PCRE2_10.47`).
> The env's `libpcre2-8.so` provides them, but the linker didn't pull it in, so a
> fresh Step-1 link failed with `undefined reference to pcre2_*@PCRE2_10.47`. This
> snapshot's `step1_saige-null/Makefile` adds `-lpcre2-8` after `-lR` to resolve it
> (needs `LIBRARY_PATH=$ENV/lib` set, as above). If your env's `libpcre2-8` predates
> PCRE2_10.47, install a matching one (`conda install 'pcre2>=10.47'`). Step-2 is
> unaffected (no R dependency).

## Run (pipeline)

1. **Step-1** → null model. Two null flavors matter for Step-2:
   - **dense-fit / single-VR** (`use_sparse_grm_to_fit: false`) → for **single-variant** Step-2.
   - **sparse-fit / cate-VR** → for **SAIGE-GENE+ region** Step-2.
   (Using the wrong VR type fails: single-variant needs a single VR; region needs cate-VR.)

2. **Step-2 single-variant:** `saige-step2 <config.yaml>` with `modelFile` = Step-1
   output dir, `varianceRatioFile`, `plinkFile`/`bgenFile`/`pgenFile`. See
   `step2_saige-step2/examples/step2_single.yaml`.

3. **Step-2 region (SAIGE-GENE+):** add `groupFile`, `annotationList`, `maxMAFList`,
   `r_corr`, etc. See `step2_saige-step2/examples/step2_region.yaml` and
   **`SAIGE_GENE_STATUS.md` §3** for the threading/memory tradeoffs.

## Reference / provenance
- R comparison reference: **SAIGE 1.3.0** (Docker `wzhou88/saige:1.3.0`) with
  `--impute_method=mean`. Do NOT use 1.5.x as reference — it has an AF>0.5
  score-variance regression (`scoreTestFast_noadjCov` centers on `2·altFreq`,
  inconsistent with the flip; cpp is correct). 1.1.3/1.3.0 agree with cpp to r≈1.0.
- Validated state: Step-1 nulls match R (LDL tau ~5 dp; T2D tau within MC noise);
  Step-2 single-variant bit-exact cpp↔R (incl. WES ultra-rare, r=1.00000);
  SAIGE-GENE+ region matches R for LDL(LDLR) + T2D(TCF7L2).
