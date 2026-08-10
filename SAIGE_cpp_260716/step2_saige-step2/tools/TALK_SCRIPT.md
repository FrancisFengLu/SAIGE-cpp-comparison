# Talk script — region-test memory investigation

Target length: ~5 minutes. English lines are what you say. `[中文]` lines are
stage directions for you only — don't read them out.

---

## 0. Framing (20 s)

> "Seokho reported a sudden memory burst running gene-based region tests at UKB
> scale — around 242 GB at 165,000 samples, OOM-killed. I found the cause. It's
> two lines. I want to spend most of this on *how* I found it, because the method
> generalizes, and because one thing in the status doc needs correcting."

`[中文] 一开口就说清楚结论已经有了，让大家放心听过程，而不是等答案。`

---

## 1. The starting clue: the knob that did nothing (45 s)

> "The code already had a parameter whose documented purpose was to limit region
> memory — `markers_per_chunk_in_groupTest`. So my first question wasn't 'why is
> this so memory-hungry.' It was 'why isn't that knob working?'
>
> I swept it at N=1,000: chunk 100, 500, 1000 — a 10× range — crossed with 1, 2,
> 4, 8 threads. Peak RSS moved by 0.06%."

**EVIDENCE — Table 1** (`bench_N1000.csv`)

| threads | chunk=100 | chunk=500 | chunk=1000 |
|--------:|----------:|----------:|-----------:|
| 1 | 1.536 | 1.537 | 1.536 |
| 2 | 2.977 | 2.989 | 3.029 |
| 4 | 2.949 | 2.993 | 2.993 |
| 8 | 2.993 | 2.970 | 2.942 |

> "Memory tracked thread count and nothing else. That reframed the whole problem
> — from 'the algorithm needs this much' to 'a parameter isn't wired up.'
> Everything after this was just confirming that."

`[中文] 这是第一个转折点。强调它是纯黑盒实验：只改配置、只量 RSS,不看一行源码。`

---

## 2. Killing a plausible wrong answer (30 s)

> "Before going further I had to rule out a boring explanation: the extra 1.4 GB
> per thread might just be BLAS scratch space — OpenBLAS reserving working memory
> per thread, not our data.
>
> I pinned `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` to 1 so only our own
> parallelism ran. Memory changed by about 10%. So it's ours, not the library's."

`[中文] 这一段的价值是展示你会主动否掉自己的假设。开会时"我猜X,测了,不是X"
比"我一猜就中"可信得多。有人问的话:10% 是线程本地缓存的正常波动。`

---

## 3. Locating it in time (40 s)

> "Next, when does it happen? I sampled `/proc/<pid>/status` VmRSS every 20
> milliseconds, and timestamped the program's own stdout against the same clock,
> so I could line the two streams up.
>
> The entire burst lands in a 2.6-second window, immediately after the program
> prints `In chunks 0-0, 20 markers are ultra-rare`. That's one specific line of
> output to anchor on.
>
> And the shape matters: it climbs, then drops straight to zero. Not a staircase
> that never comes down. So this is not a leak — it's a live allocation, held,
> and correctly freed. That saved me from going down the leak-hunting path."

`[中文] 工具是 rss_timeline.sh。"曲线形状排除泄漏"这一点如果有人追问:泄漏的特征
是不 release,RSS 阶梯式只升不降;我们看到的是升上去然后归零,说明是正常的大块申请。`

---

## 4. Locating it in space (40 s)

> "Then, where in the code? I ran it under valgrind's massif heap profiler and
> looked at the snapshot marked peak — not the last snapshot, the peak one.
>
> 97.57% of the 1.52 GB peak came from exactly two calls, both
> `arma::op_resize::apply_mat_inplace`. Two `.resize()` calls. That's it."

**EVIDENCE — massif peak snapshot**: useful-heap 1,635,657,228 B, matching the
T=1 RSS measurement exactly. Two `op_resize` frames = 97.57%.

`[中文] 如果有人不知道 massif:它记录"每一块堆内存是被哪条调用栈申请的",所以能
从"用了多少"直接跳到"哪一行申请的"。`

---

## 5. Naming the variable by arithmetic (40 s) ← the punchline

> "Here's the step I'd highlight. Each of those two blocks was 800 megabytes.
> 800 MB divided by 8 bytes is 10^8 doubles. By design the matrix should be
> chunk × N — 500 × 1,000 — which is 5 × 10^5. So it's 200× too large.
>
> And 10^8 divided by N, which is 1,000, gives exactly 100,000. That is the
> default value of a different config parameter: `max_markers_region`.
>
> So before reading any more source, the arithmetic had already named the
> variable."

**EVIDENCE — the arithmetic**

```
observed block      800 MB / 8 B      = 1e8 doubles
expected (design)   500 × 1000        = 5e5 doubles      → 200× too big
1e8 / N (=1000)                       = 100,000
                                        = default of max_markers_region
```

`[中文] 这是全场最有说服力的一步,慢一点讲。维度分析:不用调试器、不用读代码,
纯靠单位和数量级就锁定了变量名。`

---

## 6. Confirmation, three ways (40 s)

> "Then three independent confirmations.
>
> One — source: `main.cpp` passes `max_markers_region`, default 100,000, to
> `setRegion_GlobalVarsInCPP`, while the code that actually allocates the matrix
> uses `markers_per_chunk_in_groupTest`, default 500. Same place, two different
> numbers.
>
> Two — cross-reference against R SAIGE: R passes *one* variable to both the
> allocation and the flush threshold. So this is a wiring mistake in the port,
> not a design decision I was about to break.
>
> Three — causality: I changed only `max_markers_region` in the config, nothing
> else, and memory moved proportionally. Correlation to causation."

**EVIDENCE — the fix, 2 lines**

```diff
-                max_markers_region,
+                (unsigned int)markers_per_chunk_in_groupTest,
```

`[中文] 两处,但其中一处在 LDmat 路径上,那条路用稀疏存储、实际不吃内存,所以
"生效的只有一行"。有人问就这么说,别夸大。`

---

## 7. Result (50 s)

> "Same config, before and after, at N=1,000: 28 to 48× reduction depending on
> thread count. And post-fix, the matrices are 8 MB of the 54 MB total — they've
> gone from dominant term to rounding error."

**EVIDENCE — Table 2**

| threads | before | after | reduction |
|--------:|-------:|------:|----------:|
| 1 | 1.54 GB | 0.054 GB | 28× |
| 2 | 2.98 GB | 0.063 GB | 47× |
| 4 | 2.99 GB | 0.062 GB | 48× |

> "Then the decisive test. I only have a 16 GB laptop, so I generated 165,000
> synthetic samples and ran both builds under an 8 GB hard cgroup cap."

**EVIDENCE — Table 3**

| build | chunk | threads | exit | peak RSS | model predicts |
|:------|------:|--------:|:-----|---------:|---------------:|
| pre-fix | 100000 | 1 | **137 SIGKILL** | hit cap | 245.9 GB |
| post-fix | 500 | 1 | 0 | 1.405 GB | 1.229 GB |
| post-fix | 500 | 2 | 0 | 2.672 GB | 2.459 GB |
| post-fix | 500 | 4 | 0 | 5.180 GB | 4.917 GB |
| post-fix | 2000 | 1 | 0 | 5.093 GB | 4.917 GB |

> "The pre-fix kill is my favourite piece of evidence. The kernel log reports
> `total-vm: 129857452 kB` — that's 123.8 GiB, which is exactly 100,000 × 165,000
> × 8 bytes. One P1Mat. It died while zero-filling the first matrix, before it
> even asked for the second."

`[中文] 这条是端到端自证的,不依赖 Seokho 的 242。如果只能留一条证据,留这条。`

---

## 8. Correcting the model (30 s)

> "One correction to the status doc. It describes the peak as `~N × N`. It isn't
> — it's linear in N:
>
> `peak ≈ base + min(threads, genes) × 2 × chunk × N × 8 bytes`
>
> Accurate to +14% worst case in the table above. This matters for planning: at
> half a million samples, an N-squared model gives a number that says 'this
> feature is unusable,' and a linear model says 128 GB on 32 threads."

**EVIDENCE — Table 4**

| | per sample, per thread | N=165k, 1 thread | N=500k, 32 threads |
|:--|---:|---:|---:|
| before fix | 1.6 MB | 246 GB | 23.8 TB |
| after fix (chunk=500) | 8 KB | 1.4 GB | 128 GB |

`[中文] 语气上别踩人。说"the doc's observation was right, the attribution wasn't"。`

---

## 9. Correctness + the free bug fix (30 s)

> "This is a pure memory change — no numerical change. 34 output files
> byte-identical, zero differences: quantitative, binary including the Firth and
> exact-test branches, survival, and PGEN, PLINK and BGEN inputs, at 1, 2 and 4
> threads. And I verified the Firth and exact-test branches actually executed,
> with instrumented counters, rather than assuming they did.
>
> And one bug fixed for free. The flush threshold was 100,000 while the matrix
> had 500 rows — so any gene with more than 500 passing markers aborted with an
> out-of-bounds. Which means the entire chunk-spill mechanism was dead code and
> had never once executed. The fix activated it."

`[中文] "verified with counters rather than assuming" 这句值得说,显得严谨。`

---

## 10. Open items (30 s)

> "Three things to flag.
>
> First, activating chunking exposed a pre-existing upstream SAIGE bug — the
> off-diagonal variance blocks are force-symmetrized, but the matrix genuinely
> isn't symmetric because `P2Vec` switches operator at MAC 20.5. I confirmed it
> reproduces in R SAIGE 1.5.1 with the same magnitude, so our port is faithful.
> I deliberately did not patch it — whether V should be symmetric is an upstream
> design question. But note it affects every region run, because the ultra-rare
> pseudo-markers always form a trailing chunk.
>
> Second, on provenance: the 242 GB figure has no measurement artifact in the
> repo — no time log, no OOM log, no scheduler record. It's labelled peak RSS and
> my model agrees to 1.7%, so it's very likely right, but I'd call that
> corroboration rather than a validated match. Seokho, if you still have that
> log, I'd like it.
>
> Third, and this is the honest caveat: I have not measured wall-clock speedup.
> Locally I have four cores and four genes — I can't see thread scaling. All I can
> claim today is that memory is no longer the constraint on parallelism. Getting
> an actual speed number is my next step, and it needs a few hundred genes.
>
> Also — which genotype format is UKB using? The VCF region path is currently
> non-functional and BGEN needs a `.bgi` index. Both pre-existing, but they'll
> block you."

`[中文] 最后这段是加分项。主动说"速度我还没测"比被人问出来强得多。`

---

## Backup — likely questions

| Question | Answer |
|:--|:--|
| Only 200 variants in the 165k test — representative? | Memory depends on `max_markers_region × N`, not variant count. Fewer variants makes it cleaner by removing "big data" as a confounder. |
| Why 500 and not something tuned? | It's R's default. Matching R keeps us on the validated path. It's still over-allocated — genes here have ~50 markers. |
| Is 1.4 GB at 165k still too high? | It's `2 × 500 × 165000 × 8` = 1.32 GB plus base. Expected, not residual. |
| Threads vs cores? | The model counts threads — each thread allocates its own P1Mat/P2Mat. Independent of core count. |
| Could we drop chunk below 500? | Yes, but it forces multi-chunk, which trips the upstream symmetry bug. Keep chunk above the largest expected gene until that's resolved. |
