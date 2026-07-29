# Region-test memory: root cause found, fixed, and what still needs validation

Investigated 2026-07-29 on the `SAIGE_cpp_260716` snapshot (local box: 16 GB
RAM, 4 cores, N=1000 test data from `test_data/Step_2_Feb_11/test/data`).

---

## TL;DR

The region-test "sudden memory burst" is **not** the `~N×N` accumulation that
`SAIGE_GENE_STATUS.md` §4 describes. It is a **single, findable, linear-in-N
allocation pair** caused by a config-key divergence from R SAIGE:

    per-thread peak ≈ 2 × max_markers_region × N × 8 bytes
                        └─ defaults to 100000, never tuned by anyone

At N=165582 that is **246 GB** — matching the ~242 GB Seokho measured for T2D.
A **one-line fix** (pass `markers_per_chunk_in_groupTest`, as R does) cuts local
peak RSS **28–48×** with bit-identical results.

---

## 1. How it was found

`massif` peak snapshot: **97.6% of 1.64 GB** came from two adjacent
`arma::op_resize::apply_mat_inplace<double>` calls inside `mainRegionInCPP`,
800,000,000 B each = 1e8 doubles. Those are `main.cpp:2801-2802`:

```cpp
int m1new = std::max(m1, q_anno_maf);       // m1 = g_region_maxMarkers_cutoff
P1Mat.resize(m1new, P1Mat.n_cols);          // 100000 × 1000  → 800 MB
P2Mat.resize(P2Mat.n_rows, m1new);          // 1000 × 100000  → 800 MB
```

An RSS timeline (`tools/rss_timeline.sh`) independently localized the growth to
the 2.6 s window *after* ultra-rare classification, consistent with this site.

## 2. Root cause: two variables that must be equal, aren't

R SAIGE uses **one** variable for both the allocation and the chunk-flush
threshold:

```r
# SAIGE_SPATest_Region.R:273
P1Mat = matrix(0, markers_per_chunk_in_groupTest, n)
# SAIGE_Test_main.R:223-228
setRegion_GlobalVarsInCPP(maxMAF_in_groupTest,
                          markers_per_chunk_in_groupTest,  # -> m1
                          ...)
```

The C++ port allocates with `markers_per_chunk_in_groupTest` (default 500) but
feeds `setRegion_GlobalVarsInCPP` a **separate** key `max_markers_region`
(default **100000**). So:

- `m1` (flush threshold, `main.cpp:2763`) = 100000
- `P1Mat` rows actually allocated (`main.cpp:4455`) = 500

### Consequence A — memory
`m1new = max(100000, q_anno_maf)` inflates both matrices to `100000 × N`,
**per OMP thread**. This is why T=1→T=2 cost +1.44 GB at N=1000.

### Consequence B — latent crash (verified)
Because the flush at `i1InChunk == m1` can never fire before
`P1Mat.row(i1InChunk)` runs past row 500, any gene with more passing markers
than `markers_per_chunk_in_groupTest` aborts:

```
terminate called after throwing an instance of 'std::out_of_range'
  what():  Mat::row(): index out of bounds
```

Reproduced with a 600-variant gene (262 passing) at
`markers_per_chunk_in_groupTest: 50`. **The chunk-spill mechanism has been dead
code** — it never executed, so it has never been validated.

## 3. The fix

`main.cpp` — two call sites, matching R:

```cpp
setRegion_GlobalVarsInCPP(maxMAFList,
    (unsigned int)markers_per_chunk_in_groupTest,   // was: max_markers_region
    MACCutoff_to_CollapseUltraRare, min_gourpmac_for_burdenonly);

setGlobalVarsInCPP_LDmat(..., minINFO,
    (unsigned int)markers_per_chunk_in_groupTest,   // was: max_markers_region
    outputFile);
```

### Measured effect (N=1000, 2 genes, `markers_per_chunk_in_groupTest: 500`)

| nThreads | peak RSS before | peak RSS after | ratio |
|---:|---:|---:|---:|
| 1 | 1.54 GB | 0.054 GB | 28× |
| 2 | 2.98 GB | 0.063 GB | 47× |
| 4 | 2.99 GB | 0.062 GB | 48× |

- 600-variant gene at chunk=50: **abort(134) → exit 0**, 6 chunk flushes.
- Region + singleAssoc outputs **bit-identical** before vs after (modulo row
  order, which varies with thread scheduling).
- T=1 / T=2 / T=4 outputs identical after the fix.

## 4. What the fix exposes — RESOLVED: an upstream SAIGE bug

**Status 2026-07-29: investigated and settled. The port is faithful; the bug is
in upstream SAIGE (verified against R SAIGE 1.5.1 built from
`code_copy/SAIGE_isolated/`). The chunking fix in §3 is correct and safe.**

### Root cause: `VarMat` is force-symmetrized but isn't symmetric

`main.cpp:3000` (upstream: `SAIGE_isolated/src/Main.cpp:1789`) assembles the
off-diagonal blocks as:

```cpp
VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();  // assumes symmetry
```

But `P1Mat * P2Mat` is genuinely asymmetric. `P2Vec` is computed by **two
different operators** depending on each marker's MAC (`saige_test.cpp:220-229`):

```cpp
if (!ctx.flagSparseGRM_cur)  t_P2Vec = t_gtilde % m_mu2 * m_tauvec[0];   // diagonal approx
else                         t_P2Vec = getPCG1ofSigmaAndGtilde(...);      // sparse-GRM solve
```

`flagSparseGRM_cur` flips when `MAC > m_cateVarRatioMinMACVecExclude.back()`
(=20.5). So `V[i,j] = gtilde_i · P2Vec_j ≠ gtilde_j · P2Vec_i` whenever markers
i and j straddle that MAC boundary.

In the 600-variant test gene (262 passing markers, 42 sparse / 220 diagonal):
- measured asymmetry `max|V − V'| = 24.4` vs `max|V| = 118` (~20%)
- R-side: off-diagonal entries differ by up to **33%**, and **92% of asymmetric
  entries are marker pairs straddling MAC=20.5**
- `sum(Phi)` shifts 0.65% between chunk sizes; `SE_Burden` ratio equals
  `sqrt(sum(Phi))` ratio exactly

Diagonal blocks keep their true (asymmetric) values; off-diagonal blocks get
mirrored. Changing the chunk size changes *which* entries land off-diagonal,
hence which of the two unequal values survives — non-monotonic because it
depends on the exact partition, not its coarseness.

### Evidence (two independent lines)

1. **C++ / NumPy:** all four chunked `VarMat`s (50/100/150/200) were
   reconstructed **bit-exactly (maxdiff 0.0)** from the chunk=300 matrix by
   applying only the "upper := lower transpose" rule at the respective chunk
   boundaries. `trace(VarMat)` is identical across every chunk size.
2. **R SAIGE 1.5.1:** varies with chunk size the same way, non-monotonically,
   with `SE_Burden` matching the C++ port to ~5 significant figures at every
   chunk size:

   | chunk | R SE_Burden | C++ SE_Burden |
   |---:|---|---|
   | 50  | 0.000600545692839663 | 0.000600545308988055 |
   | 100 | 0.000598377104729481 | 0.000598368797877467 |
   | 150 | 0.000597539630160186 | 0.000597532175215812 |
   | 200 | 0.000597644951355938 | 0.000597634598465226 |
   | 300 | 0.000602515406037800 | 0.000602496890772734 |

   `Pvalue_SKAT` is invariant here only because it is diagonal-dominated;
   SKAT-O and Burden both move.

### Scope: this affects every region run, not just chunked ones

The ultra-rare pseudo-markers **always** form their own trailing chunk, so
`nchunks >= 2` even for a "single-chunk" gene. Applying an explicit-transpose
fix makes all chunk sizes converge to `Pvalue = 0.493677584200227`
(`SE_Burden = 0.000602504863626627`), vs `0.493675709324556` from the
nominally-single-chunk run — a 1.9e-6 residual from the mirrored 262×4 UR band.
Small here; magnitude depends on how many marker pairs straddle MAC=20.5.

### Recommendation

- **Ship the §3 chunking fix** — it is correct and independent of this issue.
- Keep `markers_per_chunk_in_groupTest` above the largest gene anyway; that
  minimizes (does not eliminate) the mirrored band.
- **Report upstream to the SAIGE authors.** The open design question is whether
  the two-operator `P2Vec` split is intended. If `V` is *meant* to be symmetric,
  the correct fix is one consistent operator per region, not a faithful
  asymmetric assembly. That call belongs upstream, not in this port.
- A candidate port-side patch (explicit transposed block instead of mirroring)
  exists but costs 2 extra chunk loads + 1 GEMM per block pair, and is only
  worth applying once upstream settles the design question. **Not applied.**

### Reproduction assets
R harness: `/tmp/saige_build/{run_chunks.R,probe.R,group_uniq600.txt,rout/summary.csv}`
(note: R rejects duplicate SNP IDs in the group file; use `group_uniq600.txt`).

---

## 4b. Original observation (superseded by §4 above)

Activating the chunk-spill path revealed that **chunked and unchunked runs do
not agree**. Same 600-variant gene, sweeping `markers_per_chunk_in_groupTest`:

| chunk | nchunks | Pvalue (lof, maxMAF 0.05) | SE_Burden |
|---:|---:|---|---|
| 50  | 6 | 0.517091568376162 | 0.000600545308988055 |
| 100 | 3 | 0.505618532717921 | 0.000598368797877467 |
| 150 | 2 | 0.516218338432968 | 0.000597532175215812 |
| 200 | 2 | 0.506033089276641 | 0.000597634598465226 |
| 262 | 1 | 0.493675709324556 | 0.000602496890772734 |
| 300 | 1 | 0.493675709324556 | 0.000602496890772734 |
| 500 | 1 | 0.493675709324556 | 0.000602496890772734 |

`BETA_Burden`, `MAC`, `Number_rare`, `Number_ultra_rare` are identical across
all rows — only variance-derived quantities move. VarMat is `G'PG` assembled
blockwise, so it should be **invariant** to chunking; the non-monotonic spread
(150 → 0.5162 but 200 → 0.5060) looks like an indexing bug, not float noise.

Note the C++ VarMat block-assembly loop (`main.cpp:2977-3005`) is a faithful
line-by-line copy of R's (`SAIGE_isolated/src/Main.cpp`), including the
`if (P1Mat.n_cols == 0) continue;` early-outs that skip `first_row`/`first_col`
advancement. So this is **inherited from upstream SAIGE**, not introduced by the
port — but upstream may never hit it either if real runs stay single-chunk.

### Practical recommendation
**Set `markers_per_chunk_in_groupTest` above the largest gene's passing-marker
count** so `nchunks == 1` always. This keeps the validated code path *and*
still captures the full memory win, because peak now scales with the chunk knob
instead of the hardcoded 100000:

| config | per-thread peak @ N=165582 |
|---|---|
| before fix (effectively 100000) | 246 GB |
| after fix, chunk = 2000 | 4.9 GB |
| after fix, chunk = 5000 | 12.3 GB |

A 2000-marker chunk covers essentially every gene in a WES gene-based analysis
(only titin-class outliers exceed it), so single-chunk behavior is the norm.

## 5. Corrected memory model

The earlier `base(N) ≈ 8.76e-9·N²` fit in this file was **wrong** — two points
happened to lie on a parabola. The true relation is linear in N:

    peak_GB(N, T, C) ≈ base_small + T · (2 · C · N · 8 / 2^30)

| check | formula | observed |
|---|---|---|
| N=1000,  C=100000, T=1 | 1.60 GB | 1.54 GB |
| N=1000,  C=5000,   T=1 | 0.08 GB + base | 0.19 GB |
| N=1000,  C=500,    T=1 | 0.008 GB + base | 0.054 GB |
| N=165582, C=100000, T=1 | 246.7 GB | ~242 GB (T2D, status doc §4) |
| N=158598, C=100000, T=1 | 236.3 GB | ~230 GB (LDL, status doc §4) |

`tools/plan_budget.py` implements this. Status doc §4's "accumulated across many
sub-10 GB allocations, no single huge allocation" conclusion should be revised:
there *are* two huge allocations, they just weren't attributed because the
resize happens inside Armadillo's `op_resize`.

## 6. Suggested next steps

1. **Validate the chunked path against R** — run the same gene through R SAIGE
   1.3.0 with `markers_per_chunk_in_groupTest` forced small, and see whether R
   reproduces the C++ chunked numbers or the single-chunk numbers. That tells us
   whether §4 above is an upstream bug or a port bug.
2. **Until then, ship the fix with a guard**: warn (or clamp) when a region's
   passing-marker count would force `nchunks > 1`.
3. **Re-run the LDL/T2D R-vs-C++ validation** with the fix in place — it should
   now fit comfortably on a 250 GB node, and even on much smaller ones.

## Tools in this directory

| file | purpose |
|---|---|
| `bench_region_mem.sh` | sweep `(nThreads, markers_per_chunk)`, capture peak RSS/wall/CPU → CSV |
| `rss_timeline.sh` | sample `/proc/<pid>` RSS while running, timestamp-aligned with program stdout, to localize growth to a phase |
| `plan_budget.py` | `--N <n> --budget-gb <g>` → feasible `(nThreads, chunk)` settings |
| `bench_N1000.csv` | raw 12-cell sweep, pre-fix |
