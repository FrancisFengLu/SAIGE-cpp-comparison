# SAIGE step2 C++ — Parallelism Design

**Branch:** `step2-parallel` · **Author:** (drafted by Francis/Claude) · **For sign-off by:** Seokho Jeong

Goal: speed up step2 wall time by parallelizing the marker loops while keeping
numerical results bit-comparable to the current serial implementation and to R
SAIGE. This document does **not** propose code changes yet — it asks for
sign-off on a phased plan before any commits.

---

## 1. Canonical step2 pipeline

Step2 has two entry-loop shapes, both dispatched from `main.cpp` after the null
model and genotype handles are loaded.

### 1a. Single-variant path (`mainMarkerInCPP`, main.cpp:637–1077)
1. **Allocate per-marker output vectors** sized `q = #markers` (main.cpp:649–688). Single shared variance ratio is pre-assigned if `m_varRatio_null.n_elem == 1` (main.cpp:697–703).
2. **For each marker i** (main.cpp:708): parse byte offset → `Unified_getOneMarker` (main.cpp:747, dispatched at genotype_reader.cpp:2158) reads/decompresses one variant.
3. **QC filter** on missing rate / MAF / MAC / INFO (main.cpp:779–783); else `imputeGenoAndFlip` (788).
4. **Variance-ratio + sparse-GRM flag selection** mutates SAIGEClass scalars `m_flagSparseGRM_cur`, `m_varRatioVal` (main.cpp:837–851).
5. **Score test (cheap path)** via `Unified_getMarkerPval` (main.cpp:872 ER branch, 889 normal). Inside `getMarkerPval` (saige_test.cpp:419) one of `scoreTestFast_noadjCov` (saige_test.cpp:317), `scoreTest` (153, sparse GRM), or `scoreTestFast` (233) runs.
6. **SPA** if binary/survival and `|Tstat|/sqrt(var1) > m_SPA_Cutoff` (saige_test.cpp:508).
7. **ER** if binary and `MAC ≤ g_MACCutoffforER` (saige_test.cpp:498–501 + main.cpp:871).
8. **Firth re-fit** if binary, `pval ≤ m_pCutoffforFirth`, and `m_is_Firth_beta` (saige_test.cpp:622, 629).
9. **Fast-test re-evaluation:** if `m_isFastTest` and the cheap path produced `pval < m_pval_cutoff_for_fastTest`, the marker is recomputed once with `m_isnoadjCov_cur=false` and possibly sparse GRM on (main.cpp:919–957).
10. **Append per-marker outputs** into the pre-allocated vectors (main.cpp:972–1023). After the loop, `writeOutfile_single` flushes one file (main.cpp:1061).

### 1b. Region/group path (`mainRegionInCPP`, main.cpp:1645–2400; driver loop main.cpp:3570–3645)
The outer driver reads regions in chunks of `groups_per_chunk` (main.cpp:3573–3577) and calls `mainRegionInCPP` once per region (main.cpp:3620). Inside:
1. **Per-region setup** allocates output vectors of length `q = q0 + q_anno*q_maf` (main.cpp:1685–1750).
2. **Marker loop** over `q0` variants in the region (main.cpp:1783). For each: read → QC → split into "non-URV" (MAC > `g_region_minMAC_cutoff`) vs URV.
3. **Non-URV branch** calls the same `Unified_getMarkerPval` as the single path (main.cpp:1887/1894) and *additionally* writes a column of P1/P2 into in-memory matrices `P1Mat`, `P2Mat` (main.cpp:1913–1915).
4. **URV branch** accumulates a collapsed pseudo-marker per (annotation × MAF) bin (main.cpp:1972–2005).
5. **Chunk flush:** every `m1 = g_region_maxMarkers_cutoff` non-URV markers, `P1Mat`/`P2Mat` are written to `*_P1Mat_Chunk_N.bin` (main.cpp:2011–2023) and the in-memory chunk is reused.
6. **Post-loop:** URV pseudo-markers go through `getMarkerPval` (main.cpp:2055–2207). Then `VarMat = P1*P2` is rebuilt by loading the chunked files back (main.cpp:2219–2252) and BURDEN / SKAT / SKAT-O p-values are computed per (annotation × MAF) group (main.cpp:2474–2566) and CCT-combined.

**Mutable per-marker SAIGEClass state** (the threading blockers): `m_flagSparseGRM_cur` (saige_test.hpp:55), `m_isnoadjCov_cur` (saige_test.hpp:58), `m_varRatioVal` (saige_test.hpp:36). Setters at saige_test.cpp:838–851, 883, 929. Everything else on the class (XV, XVX, mu, mu2, res, sparse Σ, `m_varRatio_*` tables) is read-only after init.

---

## 2. Parallelism opportunities

| # | Step | Parallelism option | Granularity | Expected speedup | Difficulty | Risk |
|---|---|---|---|---|---|---|
| 1 | Single-variant cheap score test (loop at main.cpp:708) | OpenMP `parallel for` over markers, after lifting `m_*_cur` and `m_varRatioVal` into per-iteration locals | per marker | 4–10× on the score-test phase; less on whole step2 because BGEN decode often dominates | Medium | Mutating shared `SAIGEClass` scalars (saige_test.hpp:36,55,58) ⇒ false sharing / wrong-varRatio bugs; output-vector ordering preserved since we index by `i` |
| 2 | Matrix-level first pass: replace inner `scoreTestFast` with a single GEMM over a block of B markers stacked into `G ∈ ℝ^{n×B}` | block of markers | 2–4× on cheap path only (it's already mostly BLAS-bound at small B) | Hard | scoreTestFast uses sparse-G shortcuts (`getadjGFast`, `square(g1_tilde)·mu2`); naive GEMM gives up that sparsity for dense markers, may *slow down* PLINK/PGEN streaming |
| 3 | Candidate recompute (SPA/Firth/ER) after first pass | OpenMP over candidate list | per candidate marker | 3–6× on the candidate phase (SPA/Firth/ER are 10–100× cheap path each) | Medium | Same SAIGEClass mutable scalars; SPA also writes `t_gtilde` etc. on the stack — already safe |
| 4 | BGEN block read + decode pipeline | producer thread reads `m_zBuf` chunks; worker pool runs `Parse2` (genotype_reader.cpp:1106) | per variant | 1.5–3× on BGEN with zstd; near-zero on PLINK | Medium | `m_fin`/`m_buf`/`m_zBuf` are *shared instance state* (genotype_reader.cpp:1305, 1359, 1365, 1371) — needs per-worker scratch buffers |
| 5 | VCF parallel decode | n/a | — | 0× | n/a | VCF is single-stream via htslib here; no random access ⇒ not worth it. Recommend explicitly *not* parallelizing VCF |
| 6 | Region-mode P1/P2 chunk reconstruction | parallel block-product loop at main.cpp:2219–2243 | per off-diagonal block | 1.2–2× on `VarMat` build only, which is rarely the bottleneck | Easy | Race on `VarMat.submat` writes; need per-block tile assignment |
| 7 | Outer region loop (main.cpp:3570–3645) | OpenMP over regions, one SAIGEClass per worker | per region | 4–8× on group-mode wall time if region count ≫ threads | Hard | `ptr_gSAIGEobj` is a singleton; needs either deep-copy per thread or moving `m_*_cur` into a per-call context struct; output-file ordering needs a serializer thread |
| 8 | Quantitative-trait fast path | trivial OpenMP over markers (no SPA/Firth/ER ever fires for quantitative — see saige_test.cpp:508) | per marker | 4–10× on quantitative-only runs | Easy | Same as #1 minus the candidate path; isolate this as a first vertical slice |

Excluded: PCG solve inside `scoreTest` (saige_test.cpp:182) is already a hot Armadillo call; threading it across markers would oversubscribe BLAS — better left as the sparse-GRM fallback. Per-marker BLAS is already multi-threaded via Armadillo's OpenBLAS link.

---

## 3. Recommended design (3 phases, mapped to Seokho's bullets)

### Bullet 1 — "Matrix-level first pass for most markers"

**Concrete approach.** Inside `mainMarkerInCPP`, read markers in blocks of B (default 256) into `G ∈ ℝ^{n×B}`. Run a *vectorized cheap path* that computes `Tstat[1..B]` and `var1[1..B]` jointly: for the no-sparse-GRM, no-noadjCov common case this is `S = (Gᵀ res) / tau`, `var2` per-marker via column-wise dot products with `m_mu2`. This is essentially `scoreTestFast` (saige_test.cpp:233–314) re-expressed as one GEMM + a few column reductions. Output for each marker is identical to the serial path *if* we keep the same dense-vs-sparse heuristic.

**Why it's correct.** `m_X`, `m_XV`, `m_XVX`, `m_XVX_inv`, `m_res`, `m_mu2`, `m_tauvec[0]`, `m_S_a` are all set in the SAIGEClass constructor and never written during the loop (saige_test.hpp:16–25). The cheap path takes only `m_varRatioVal` from per-marker state; we can lift it to a `varRatioVec[B]` precomputed once per block from `m_varRatio_null` using each marker's MAC. The only output is read-then-write into pre-allocated `BetaVec`, `seBetaVec`, `pvalVec` etc. (main.cpp:660–665) at index `i`, which is race-free if each thread owns disjoint i.

**Implementation steps.**
1. Add `void scoreTestFast_block(const arma::mat& G, ..., arma::vec& Beta, ..., arma::vec& pval_num)` next to `scoreTestFast` in saige_test.cpp.
2. In `mainMarkerInCPP`, restructure the loop into outer block iteration + inner per-marker bookkeeping (QC, output write).
3. Keep the existing single-marker `scoreTestFast` path as a fallback for markers that need the sparse-G shortcut (`getadjGFast`-style accumulation).
4. Compute `StdStat = |Tstat|/sqrt(var1)` for the whole block; produce a *candidate mask* for markers needing SPA/Firth/ER.
5. Gate behind a CLI flag `--block-size N` (default 256, set N=1 to fall back to the current serial loop bit-for-bit).

**Verification plan.** Run all existing configs (`config_compare*.yaml`) with `--block-size 1` and confirm zero diff vs main. Then run with `--block-size 256` and require `|Δ pval| < 1e-12` for the cheap path on the `extdata/input/` 1k-sample dataset and the binary/region/sparse configs. Regression risk: `var2_c = arma::accu(m_mu2) * (2*altFreq)^2` term in `scoreTestFast_noadjCov` (saige_test.cpp:336) reuses a precomputed scalar — must be hoisted, not recomputed per column.

### Bullet 2 — "Separate only SPA/Firth/ER-needed markers as candidates"

**Concrete approach.** After the block first pass, every marker has a `(Tstat, var1, pval_noadj, MAC)` tuple. Build three candidate lists:
- `spa_list`: binary/survival markers where `StdStat > m_SPA_Cutoff` and MAC > `g_MACCutoffforER` (saige_test.cpp:498–501, 508).
- `er_list`: binary markers where `MAC ≤ g_MACCutoffforER` (main.cpp:871).
- `firth_list`: built lazily from SPA output where `pval ≤ m_pCutoffforFirth` (saige_test.cpp:622).

Process each list with the existing single-marker code, parallelized over candidates with OpenMP.

**Why it's correct.** SPA/Firth/ER are pure functions of `(GVec_i, m_mu, m_res, m_y, m_n_case)` plus `Tstat`/`var1` already computed. They write only to local stack variables and to output vectors at index `i`. The candidate list is built deterministically from first-pass results, so output ordering is preserved.

**Implementation steps.**
1. Hoist the SPA-only / ER-only / Firth-only parts of `getMarkerPval` (saige_test.cpp:505–706) into standalone functions taking precomputed `Tstat, var1` and returning `(Beta, seBeta, pval, isConverge, …)`.
2. Run `#pragma omp parallel for schedule(dynamic, 8)` over each candidate list.
3. Make sure the lifted functions take *no* `SAIGEClass&` mutable scalars — pass `varRatioVal` and `flagSparseGRM_cur` as args.
4. For Firth, retain the dependency on SPA's output (Firth needs the SPA `pval`), so SPA must finish before Firth list is built.
5. Note for Seokho: **this optimization does nothing for quantitative traits** — `m_traitType != "quantitative"` gates SPA (saige_test.cpp:508) and Firth (binary-only block at 620), and ER is binary-only. Quantitative benefits only from Bullet 1 + Bullet 4.

**Verification plan.** Diff per-marker outputs against `--block-size 1` baseline on `config_compare_binary.yaml` (forces SPA path) and `config_compare_binary_region.yaml`. Acceptance: byte-identical `pvalVec` strings (since both paths print via the same `sprintf("%.6E", …)` at saige_test.cpp:602/611).

### Bullet 3 — "Recompute candidates with the existing exact path"

**Concrete approach.** This is the conservative variant of Bullet 2: instead of refactoring SPA/Firth/ER, just call the existing `Unified_getMarkerPval` (`getMarkerPval`) on each candidate, in parallel. The block first pass produces "preliminary" results that get overwritten for any candidate marker.

**Why it's correct.** No code path changes — we are reusing the audited serial code, just calling it on a subset. Per-call SAIGEClass mutations (`m_flagSparseGRM_cur`, `m_varRatioVal`) become per-thread races; we fix this by introducing a small `struct PerMarkerCtx { bool flagSparseGRM_cur; bool isnoadjCov_cur; double varRatioVal; }` passed as an argument, and a `getMarkerPval_ctx(...)` overload that reads only from the ctx, not the SAIGEClass scalars.

**Implementation steps.**
1. Define `PerMarkerCtx` in saige_test.hpp.
2. Add `getMarkerPval_ctx(...)` that reads `t_isSparseGRM`, `t_isnoadjCov`, `t_varRatioVal` from the ctx and a local `is_region` flag, never touching `m_*_cur`.
3. Replace the three `m_*_cur` field reads inside `scoreTest`/`scoreTestFast`/`scoreTestFast_noadjCov` with passed-in arguments (small surface — ~6 occurrences across saige_test.cpp).
4. Keep the old API as a wrapper around the ctx version so the existing region path doesn't have to change.
5. Parallelize the candidate-recompute over markers with OpenMP.

**Verification plan.** Same diff harness as Bullet 2 plus the conditional analysis config (`config_compare_conditional.yaml`) — this exercises the `m_P2Mat_cond`/`m_VarInvMat_cond` read paths (saige_test.cpp:715–845) which must remain read-only across threads.

### Bullet 4 — "I/O parallelism only for block read/decode support"

**Concrete approach.** Producer-consumer for BGEN only. One I/O thread does `fseek + fread` of the `(snpID, rsID, chr, alleles, C, D, zBuf[C-4])` block per variant (genotype_reader.cpp:1315–1365). A pool of decoder threads runs `Parse2` (genotype_reader.cpp:1106) on the compressed buffer. Output is a thread-safe queue of `(metadata, arma::vec dosages)`. The main thread consumes in original order.

**Why it's correct.** `Parse2` reads only `m_posSampleInModel` (set once at init, genotype_reader.cpp:1089–1097) and writes only to its `dosages` output arg — no per-instance mutable state. The unsafe state is the file pointer `m_fin` and the scratch buffers `m_buf`/`m_zBuf`, which are confined to the single I/O thread (genotype_reader.cpp:1305, 1359, 1365). Decoder workers each get their own `buf` allocated per call.

**Implementation steps.**
1. Add `BgenClass::readRawBlock(uint64_t offset, RawBlock& out)` that does only I/O — no Parse2.
2. Add `BgenClass::parseBlock(const RawBlock&, arma::vec& dosages, ...)` that does only decode (callable from any thread).
3. Add a `BgenStreamer` class with one reader thread and a `std::thread::hardware_concurrency()-2`-sized decoder pool, bounded queue (e.g. 64 raw blocks).
4. Plumb through `mainMarkerInCPP` / `mainRegionInCPP` so they pull from the streamer instead of calling `Unified_getOneMarker` directly when `t_genoType == "bgen"`.
5. For PLINK/PGEN, do *not* introduce a streamer — single-thread read is already fast (no decompression). For VCF, explicitly skip (htslib stream, no random access here).

**Verification plan.** Compare `output/` against serial on `config_test_bgen.yaml`. Measure wall-time on a real BGEN with zstd compression (this is where it pays off). If the BGEN streamer doesn't improve wall time by >20% on the test data, drop the feature — it's not worth the complexity. Regression risk: order desynchronization in the queue ⇒ markers tagged with wrong byte offsets. Mitigation: each `RawBlock` carries its source marker index.

---

## 4. Open questions for Seokho

1. **Target thread count.** Should we expose `--threads N` as a CLI flag, or read `OMP_NUM_THREADS`? UKB-scale runs typically have 8–32 cores; what's the deployment target?
2. **Block size B for the matrix-level first pass.** Auto-tune based on `n` and L2 cache size, or fix it (e.g. 256) and expose as a flag for benchmarking?
3. **VCF parallelism.** I recommend *not* parallelizing VCF (no random access via htslib in our reader). Confirm or push back?
4. **Region-loop parallelism (table row #7).** This is the biggest potential win for group-mode runs but requires per-thread SAIGEClass copies. Is the memory cost acceptable (one full SAIGEClass + `m_SigmaMat_sp` per thread can be hundreds of MB)?
5. **Bit-exact reproducibility.** Are we OK with reordering reductions inside the block GEMM (which can change low bits via FP non-associativity), as long as `|Δ pval| < 1e-10`? Or do we need bit-identical output?
6. **Firth + SPA on log-scale p-values.** The current code has special-cased log-p arithmetic (saige_test.cpp:567–582, 629–640). Should the block first pass also support log-p, or fall back to serial for any marker where the cheap path returned a log-scale string?

---

## 5. Sequencing — smallest end-to-end vertical slice first

**Phase A (1–2 weeks). Foundation + Bullet 4 dry run.**
- A1. Add a benchmark harness: run all `config_compare*` configs serially, record per-config wall time and per-stage breakdown. Establishes baseline.
- A2. Refactor `m_flagSparseGRM_cur` / `m_isnoadjCov_cur` / `m_varRatioVal` into a `PerMarkerCtx` passed by argument (Bullet 3 step 1–3). This is *non-functional* — same code path, same output, just removes mutable shared state. Land this PR first; it unblocks everything else.

**Phase B (2–3 weeks). Single-variant OpenMP — the first slice with measurable wall-time improvement.**
- B1. OpenMP `parallel for` over the marker loop in `mainMarkerInCPP` (Bullet 3 step 5). Output vectors are pre-sized; thread-safe by index. Per-thread buffers for `t_GVec`, `gtildeVec`, `t_P2Vec`.
- B2. Show ≥3× speedup on `config_compare.yaml` (PLINK, quantitative) with 8 threads. **This is the smallest convincing vertical slice.**
- B3. Run all comparison configs; require bit-identical output (with a single thread) and `|Δ pval| < 1e-10` (with N threads).

**Phase C (3–4 weeks). Bullet 1 — matrix-level block path.**
- C1. Implement `scoreTestFast_block`.
- C2. Wire into single-variant loop; measure on quantitative traits first (no SPA/Firth/ER means full block path benefit).
- C3. Extend to binary traits with the candidate split (Bullet 2 — most of the work is already done by phase B).

**Phase D (2 weeks). Bullet 4 — BGEN I/O+decode pipeline.**
- D1. Implement only if profiling in Phase B shows BGEN decode is >20% of wall time.
- D2. PLINK/PGEN/VCF unchanged.

**Phase E (optional, after sign-off).** Region-loop parallelism (table row #7) — only if step2 is bottlenecked on group-mode runs and the SAIGEClass memory cost is acceptable.

**Decision gate after Phase B:** if we don't see a 3× wall-time speedup on a real test config with 8 threads after Phase B, stop and re-plan — Phases C/D add complexity that needs to be justified by Phase B results.

## Stress test caveats

Three local stress tests (N=1k×128k, N=1k×2M, N=10k×200k) showed peak speedup 1.4-2.2× @ 8 threads; N=50k Step 1 failed to converge due to synthetic data having no LD (sparse GRM became pathological for PCG). The synthetic data lacks LD, rare-variant tail, and BGEN decompression — so candidate dispatch and BgenStreamer were not stressed. Real UKB-scale performance must be measured by Seokho on actual data; see CHANGES_SUMMARY.html for full analysis.
