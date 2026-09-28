# step1bench

A standalone SAIGE step-1 null-model benchmark whose **only** variable is how
`Sigma` is inverted.

```
Sigma = tau0 * diag(1/w) + tau1 * Psi
```

`Psi` is the sparse GRM. Four interchangeable solvers sit behind one interface;
the AI-REML driver above them makes the same calls in the same order whichever
one is plugged in. Everything that could differ between solvers and must not —
how Sigma is formed, which right-hand sides are solved, which probe vectors are
drawn, how the clock is read — lives in shared code.

This is an instrument, not a product. Its value is that the solvers are measured
on identical work, so fairness matters more than raw speed.

---

## Build

```bash
source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh   # conda env saige-build
cd SAIGE_cpp_260716/step1bench
make -j8
```

Self-contained: its own `Makefile`, no `-lR`, no Rcpp, no linkage against the
`step1_saige-null/` or `step2_saige-step2/` trees. Dependencies are armadillo
(with its bundled SuperLU support), SuperLU, CHOLMOD and CXSparse, all from the
conda env.

`libcholmod 5.3.1` and `libcxsparse 4.4.1` were installed into `saige-build`
with

```bash
conda install -n saige-build -c conda-forge --freeze-installed \
      libcholmod=5.3.1 libcxsparse=4.4.1
```

pinned to exactly the versions TorchGWAS2 links in `torchgwas2`, so the
comparison is against the same CHOLMOD, not a different one. Nothing already in
`saige-build` was upgraded except `openssl`.

### The one copied file

`block_sigma.hpp` and `block_sigma.cpp` are copied **byte for byte** from
`../step1_saige-null/` at commit **`23290c08`**
(`md5 2deec426775e927c1fb17b9dda54b74b` and
`208129b76bc489356f0a7d52d73c2d00`). Nothing was changed to decouple them — the
class already takes `(loc, val, n)` as arguments and includes nothing outside
`<armadillo>` and the standard library, so it linked as-is. What is benchmarked
is therefore the production class, not a reimplementation of it.

---

## CLI

```
./step1bench --grm <prefix> --pheno <file> --pheno-col <name> [options]

  --grm <prefix>        sparse GRM: <prefix>.mtx (MatrixMarket coordinate real
                        symmetric, 1-based) + <prefix>.ids (one id per line)
  --pheno <file>        phenotype table with a header line
  --pheno-col <name>    the trait column
  --trait quant|binary  default quant
  --covar-cols a,b,c    covariates; an intercept column is always added
  --pheno-id-col <name> id column in --pheno (default IID)
  --solver block|cholmod|cholmod_reuse|superlu    default block
  --trace exact|hutchinson                        default hutchinson
  --probes N            Hutchinson probes per iteration (default 30)
  --maxiter N           AI-REML iterations (default 20)
  --tol X               relative change on tau (default 0.02 — SAIGE's)
  --seed N              probe RNG seed (default 1)
  --threads N           OpenMP threads (default 1)
  --report <file.json>  also write a machine-readable report
  --self-check-only     run the Sigma self-check and exit
```

Example — the SAIGE baseline, and the two things the block path changes:

```bash
D=/opt/saige/data
A="--grm $D/mid.fam10.sgrm --pheno $D/mid.fam10.pheno.txt --pheno-col q1 --covar-cols x1,x2"

./step1bench $A --solver superlu --trace hutchinson   # SAIGE baseline
./step1bench $A --solver block   --trace hutchinson   # inverse only
./step1bench $A --solver block   --trace exact        # inverse + exact traces
```

---

## Two independent axes

They are deliberately separate switches because they are separate claims:

| axis | claim it isolates |
|---|---|
| `--solver` | "an explicit inverse is cheaper per solve than a factorisation" |
| `--trace`  | "exact traces delete `--probes` solves per AI-REML iteration" |

`--trace exact` requires a solver with an explicit inverse, so
`superlu + exact` is refused with an explanation (exit 2) rather than silently
substituted.

---

## The four solvers, honestly

### `block`

`blocksigma::BlockSigma`, unchanged. `build()` runs union-find over Psi's
off-diagonal entries once per run; `refresh(w,tau)` forms each connected
component's dense sub-Sigma and inverts it, skipping the whole thing when
`(w,tau)` repeat bitwise; `solve()` is a per-block matrix–vector product with
fp64 accumulation and one rounding on store.

* **Provides**: `Sigma^-1` at every position where Psi is non-zero, so
  `tr(Sigma^-1 Psi)` and `tr(Sigma^-1)` are exact contractions with no probes.
* **Cannot provide**: anything on a GRM that is not block-diagonal, or whose
  blocks are too large to invert densely. `build()` keeps its production
  wall-clock gate; if it refuses, the harness stops and says so rather than
  quietly falling back.

### `cholmod`

TorchGWAS2's shape, ported from
`/opt/saige/logs/priorart_mt/TorchGWAS2/src/SparseInverse.cpp:388`
(`inv_spamat`) and `src/GMMAT.cpp:1112-1135`, with their settings reproduced
exactly:

```
cm.supernodal         = CHOLMOD_SIMPLICIAL
cm.method[0].ordering = CHOLMOD_NATURAL
cm.postorder          = 0
cm.nmethods           = 1
cm.final_ll           = 1
cm.final_pack         = 1
cm.final_asis         = 0,  cm.final_monotonic = 1   (set after analyze)
```

Per refresh it rebuilds the sparse Sigma, runs `cholmod_start`,
`cholmod_analyze`, `cholmod_factorize`, converts the factor to CXSparse, runs
one `cs_di_spsolve` per column against a **sparse identity** to get the columns
of `L^-1`, forms `Sigma^-1 = (L^-1)' (L^-1)`, and `cholmod_finish`es.
Traces are `sum(Sigma^-1 .* Psi)` and `tr(Sigma^-1)`, their
`GMMAT.cpp:1161-1164`.

`cholmod_analyze` is **inside** the per-refresh function in their code, so the
symbolic factorisation is redone every time. That is reproduced faithfully.

* **Provides**: the full explicit `Sigma^-1` as a sparse matrix — strictly more
  than the block path holds, and exact traces from it.
* **Cannot provide**: a useful answer on a GRM whose inverse is not sparse. Here
  Sigma is block-diagonal with max block 10, so `Sigma^-1` has 139,068
  non-zeros; on a GRM with one large component the explicit inverse would be
  dense and this approach would not fit in memory.

### `cholmod_reuse`

Identical to `cholmod` except that `cholmod_start` and `cholmod_analyze` are
hoisted into `build()`. It exists to answer one question — *how much of their
cost is the redundant symbolic step* — and it costs about twenty lines. It
re-analyses (and counts it in the report) if Sigma's sparsity pattern ever
changes; on this data it never does, but `tau1 == 0` would collapse the pattern
to the diagonal, so the check is not decorative.

### `superlu`

SAIGE's shape, and the baseline. Mirrors `gen_spsolve_v4`
(`../step1_saige-null/SAIGE_step1_fast.cpp:6581`) with both of its caches off,
which is what stock SAIGE and SAIGE-GPU 1.3.3 do: every `solve()` rebuilds the
whole n×n `sp_mat` from every nnz via `gen_sp_Sigma` and runs a fresh
`arma::spsolve` — a complete SuperLU symbolic *and* numeric factorisation — for
one right-hand side. `refresh()` only records `(w, tau)`; there is no state.

* **Provides**: `Sigma^-1 b`.
* **Cannot provide**: `Sigma^-1`. `traces()` prints why and aborts rather than
  returning a number it does not have; the CLI refuses the combination up front.

---

## Fairness: all four must form the same Sigma

`gen_sp_Sigma()` in `sigma_solver.cpp` reproduces SAIGE's
`gen_sp_Sigma` (`SAIGE_step1_fast.cpp:6256`) exactly, including the two details
that are easy to get wrong:

* `tau0 / w_i` is computed in **fp32** and added into an fp64 slot;
* the `1e-4` floor is applied to a diagonal entry **after** `tau0/w_i` has been
  added to it, not before.

`block_sigma.cpp::refresh()` already mirrors both (`block_sigma.cpp:428-462`);
`cholmod`, `cholmod_reuse` and `superlu` call the shared function.

Every run starts with a self-check that builds Sigma all three ways — the block
partition's own copy of Psi, and the shared `gen_sp_Sigma` — at two `(w,tau)`
points (`w == 1`, and a small uneven `w` where `tau0/w_i` is large) and reports
`max|diff|`. **It must be 0**, and the harness exits 3 if it is not. It then
does a numeric residual check, `max|Sigma x - b| / max|b|` for a random `b`,
on each solver, which covers the other half: that each solver's machinery
actually inverts the matrix it formed.

Other things held identical:

* **Timing**: one `TimedSolver` wrapper times every solver. A solver cannot
  time itself.
* **Call pattern**: `refresh(w,tau)` runs before **every** `solve(b)`, exactly
  as `gen_spsolve_v4` calls `blocksigma::refresh` on each of its 199–261 calls.
  Solvers that cache return immediately and report the hit; SuperLU does not
  cache, by design, and pays the rebuild inside `solve()`.
* **Probe vectors**: drawn from a `mt19937_64` seeded `seed + 1000*iteration`,
  so the probe stream at iteration *k* is the same for every solver even if two
  solvers took different numbers of iterations.
* **Psi·v**: one shared fp64 sparse product.

---

## The AI-REML driver

Small and legible on purpose. It is **not** a reimplementation of SAIGE's
null-model fit: no LOCO, no variance-ratio stage, no trace-CV growth loop, and
no reproduction of SAIGE's stale-W nesting. It computes the same quantities in
the same order — `Sigma^-1 y`, `Sigma^-1 X`, `tr(P Psi)`, `tr(P A0)`, the AI
matrix, the same covariate projection — and stops on the same relative-change
rule.

The covariate correction is written out rather than hidden:

```
P v = Sigma^-1 v  -  Sigma_iX * ( cov * ( Sigma_iX' * v ) ),
cov = (X' Sigma^-1 X)^-1
```

**quantitative**: `W === 1` throughout (for a gaussian link `W` and the working
response do not depend on `eta`, so the inner loop is exact in one pass). Both
`tau0` and `tau1` are estimated. Initial `tau` splits the OLS residual variance
evenly. SAIGE instead starts at `[1, 0]` and spends its first iteration on a
`tau^2 * score / n` step to leave the boundary; starting off the boundary
reaches the same optimum with less ceremony, which is what the gates check.

**binary**: `tau0` fixed to 1. The **inner** loop is IRLS, updating `W` with
`tau` held fixed; the **outer** loop is AI-REML, updating `tau1`. That nesting
direction is the one stated in `../step1_saige-null/saige_ai.cpp:52` ("the inner
IRLS loop calls it on every iteration of every outer AI-REML iteration") and
`loco_engine.cpp:243` ("PQL Newton loop with tau held FIXED").

One deliberate divergence from upstream: after the inner IRLS converges this
harness **re-solves** `Sigma_iX` and `cov` at the converged `W`, so `P` is the
textbook projection. Upstream leaves them on the previous `W`, making its `M` a
mixed operator (documented at `SAIGE_step1_fast.cpp:8459`). Reproducing that
would make the AI step depend on iteration history, which a benchmark does not
want. It is the main reason this harness's `tau` is near, not equal to,
production's.

### Trace modes

`exact` (needs an explicit inverse):

```
tr(P Psi) = tr(Sigma^-1 Psi) - tr( cov * Sigma_iX' * (Psi Sigma_iX) )
tr(P A0)  = tr(Sigma^-1)     - tr( cov * Sigma_iX' * Sigma_iX )
```

The second identity is `tr(P A0)` with `A0 = diag(1/w)` **only when `w == 1`**,
which is the quantitative case. `SigmaSolver::traces()` returns `tr(Sigma^-1)`,
not `sum_i (Sigma^-1)_ii / w_i`, so the driver only asks for it on the
quantitative path; binary never needs it because `tau0` is fixed.

`hutchinson`: `--probes` Rademacher probes, fixed count, no CV growth loop — so
the work is the same for every solver.

---

## Report

stdout prints the per-iteration trajectory, the final `tau` and `alpha`, the
connected-component summary (`nblocks`, max block, full size histogram), the
`(w,tau)` cache hit rate, a per-solver extra line, and:

```
  stage      calls        seconds
  build      1            0.0038
  refresh    220          0.0151
  solve      40           0.0114
  trace      180          0.0557
```

`solve` counts the coefficient and AI-step solves; `trace` counts probe solves
plus the trace arithmetic — the split is by the stage the driver is in, not by
the solver, so it means the same thing in all four columns. For `superlu`,
`refresh` is near-zero by construction and the per-call rebuild + factorisation
shows up inside `solve`/`trace`; the extra line splits that out.

`--report x.json` writes the same numbers machine-readably.

---

## Data

`/opt/saige/data/mid.fam10.sgrm.{mtx,ids}` + `.pheno.txt` — simulated,
N = 50,000, 34,687 connected components, max block 10, 128 quantitative and 128
binary traits, h2 = 0.30, covariates `x1,x2`. Do not point this at UK Biobank,
MVP or All of Us individual-level data.

Measured results are in
`/opt/saige/SAIGE-work/optimization/blocksigma/STEP1BENCH.md`.
