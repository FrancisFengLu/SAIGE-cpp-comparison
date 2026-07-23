# SAIGE-GENE+ (region/gene-based test) — status & developer guide

Snapshot `SAIGE_cpp_260716` (source-only). This document is the entry point for
collaborative work on the C++ SAIGE-GENE+ region test. It records what is
implemented, what is validated, how the multithreading works, and the open
issues to work on.

---

## 1. What's implemented (complete)

Region/gene-based association testing lives in `step2_saige-step2/`:

| Component | File | Notes |
|---|---|---|
| Region orchestration | `main.cpp` → `mainRegionInCPP` (~2376) | reads group file in chunks, dispatches SKAT/Burden/SKAT-O, combines with CCT |
| SKAT / Burden / SKAT-O | `skat.cpp` (71 KB) | Davies method (Kronrod quadrature), Liu fallback, per-ρ eigenvalues, `SKATO_optimal_pvalue` |
| Cauchy combination | `cct.cpp` | `CCT_cpp`, `get_CCT_pvalue` |
| Group-file parsing | `group_file.cpp` | `checkGroupFile`, `readRegionChunk`; annotation + weight columns |
| Exact test (rare/binary) | `er_binary.cpp` | `SKATExactBin_Work` |
| SPA / Firth | `spa.cpp`, `spa_binary.cpp` | SPA Φ-adjustment for region kernel |
| LD-matrix output | `ldmat.cpp` | conditional-analysis / meta LD path |

Supported group-test options (config keys):
- `annotationList` (e.g. `lof`, `missense;lof`) — annotation masks
- `maxMAFList` (e.g. 0.0001 / 0.001 / 0.01) — multiple max-MAF cutoffs
- `r_corr` — SKAT-O ρ grid
- `MACCutoff_to_CollapseUltraRare` — ultra-rare collapsing
- Beta(a,b) weights, conditional analysis, single-variant-in-group output

---

## 2. Validation status

**Validated against R (SAIGE 1.3.0) for LDL (LDLR) and T2D (TCF7L2)** at the
`SAIGE_cpp_260524` snapshot — cpp→cpp GENE+ matched R→R within ~1–3× p
(Cauchy 8.6e-21 vs R 2.4e-20; missense/lof masks reproduced). Residual diff =
null/VR/MC tolerance + a MAC-boundary rare/ultra-rare split.

**Critical fix baked in:** the `.bgi` byte-offset region-read bug (region/by-ID
reads seek by byte offset, not row index) is **fixed** here (`m_byteOffset` in
`genotype_reader.{hpp,cpp}`; absolute `fseek(SEEK_SET)`). Single-variant mode was
never affected (it streams sequentially).

### ⚠️ Validation GAP to close first
The parallelization below (Phase C/E, PLINK/PGEN streamer) was added **after** the
260524 R-validation. The **threaded region path has not been re-diffed against R**,
and single-variant has two known (harmless-so-far) OpenMP races
(`isSPAConvergeVec`, `firstEndIdx`). **First task for any GENE+ work: run a
threaded-vs-T=1-vs-R region diff** to confirm the parallel gene loop is bit-stable.

---

## 3. Multithreading model (Phase E) — READ THIS

The region test parallelizes **one gene per thread**:

```cpp
// main.cpp:4425  (regions read in chunks of `groups_per_chunk`, then:)
#pragma omp parallel for schedule(dynamic, 1)
for (int r = 0; r < chunkSize; r++) { ... mainRegionInCPP(region r) ... }
```
- Same `nThreads` config key as single-variant. `schedule(dynamic,1)` because gene
  sizes vary 5–500 markers.
- Shared state handled: per-thread scratch allocated inside the loop, output writes
  under `#pragma omp critical(outwrite)`, `regionsProcessed` via `atomic capture`,
  per-region tempfile tags to avoid chunk-spill collisions.
- **I/O:** PLINK/PGEN reads are thread-safe (`Unified_getOneMarker_ts`, per-thread
  `FILE*`, absolute seek — no lock). **VCF/BGEN still serialize** through
  `critical(genoread)` (main.cpp:2541); **BGEN region streaming is flagged TBD** —
  not fully wired for parallel region reads.

### The defining constraint: memory scales LINEARLY with threads
Unlike single-variant (flat ~0.7 GB regardless of T), each region thread allocates
its own `P1Mat`/`P2Mat` = `markers_per_chunk_in_groupTest × N`:

| nThreads | P1/P2 only (markers_per_chunk=500, N=165k) |
|---|---|
| 1 | 1.3 GB |
| 8 | 10.6 GB |
| 16 | 21.2 GB |
| 32 | **42.4 GB** |

…plus each thread's per-gene working set and the accepted large per-gene peak. So
the region **peak-RSS OOM issue and threading compound** — high `nThreads`
multiplies an already-large footprint. On a 251 GB node this (not CPU) is the
ceiling on region-test threads.

**Tuning knobs:**
- `nThreads` — gene-loop parallelism (speed ↔ memory).
- `markers_per_chunk_in_groupTest` (default **500**) — sets per-thread P1/P2 size;
  **lower it to cut per-thread memory** (main lever).
- `groups_per_chunk` (default 100) — regions per I/O batch.

---

## 4. Known issue: region peak-RSS ~N×N (ACCEPTED, not fixed)

Region tests have a huge transient peak RSS (~N×N scale: T2D ~242 GB at n=165582,
LDL ~230 GB at n=158598) that OOM-kills (rc=137) on nodes with <~250 GB free. They
**complete and match R on a high-mem node.** Decision (2026-05-24): accept as a
memory-efficiency (not correctness) issue.

Diagnosis so far (to avoid repeating dead ends): it is **NOT** the ER exact-test
`m_total` (NResampling/ExactMax cap + saturating n_choose_r already bound it), and
**NOT** glibc malloc arenas (`mallopt(M_ARENA_MAX,1)` made no difference). Both a
>10 GB operator-new trap and a >10 GB malloc trap stayed silent → there is **no
single huge allocation**; the ~N×N is **accumulated** across many sub-10 GB live
allocations. To revisit: use heaptrack/massif, not single-alloc traps. Note this
interacts with §3 — reducing per-thread footprint helps both.

---

## 5. Suggested work items
1. **Close the validation gap** (§2): threaded region diff vs R + vs T=1.
2. **Memory** (§3/§4): profile with heaptrack; the accumulated N×N + the
   per-thread P1/P2 growth are the two levers. A thread_local grow-only scratch for
   the gene working set (same trick that fixed the single-variant Firth fault storm)
   is the natural first attempt.
3. **BGEN region streaming** (§3): wire a thread-safe BGEN region reader so BGEN
   region tests parallelize (currently serialized / TBD).
4. **Fix the two single-variant OpenMP races** (`isSPAConvergeVec` → `uint8_t`;
   `firstEndIdx` → `std::atomic<int>`) — harmless so far but real.

---

## 6. Build & run
See `README.md`. TL;DR: `step2_saige-step2` builds clean with the `saige-build`
mamba env (`make clean && make`, CXXFLAGS with the ARMA defines). Region run needs
`modelFile` (a sparse-fit/cate-VR null for GENE+), `varianceRatioFile`, `groupFile`,
`annotationList`, `maxMAFList`. Example config in `step2_saige-step2/examples/`.
