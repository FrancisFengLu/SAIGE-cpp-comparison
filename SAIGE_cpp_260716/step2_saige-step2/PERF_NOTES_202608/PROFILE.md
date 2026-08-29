# SAIGE step 2 时间构成实测（2026-08-28）

机器：GCP n1-standard-8（8 vCPU，29 GB RAM），本任务只用 CPU。
数据：`/opt/saige/data/mid`（N=50,000 × M=40,000，plink2 --dummy，~1% 缺失，全 chr1）；
`/opt/saige/data/big2`（N=200,000 × M=113,000）。
null model：`/opt/saige/logs/mp/single_y4.rda`（binary，x1,x2 协变量，单 VR=0.8115，非 LOCO），
经 `tools/rda_to_arma.R` 转成 C++ 端格式 `/opt/saige/logs/step2/null_y4_arma/`。
C++ 版：`SAIGE_cpp_260716/step2_saige-step2` 拷贝到 `/opt/saige/logs/step2/cpp_step2/`
编译（saige-build env，-O3 -fopenmp），**加了相位计时插桩**（[PHASEPROF]，omp atomic 累加，
改动只在拷贝里，原仓库未动）。

## 命令

```
# R 版（对照，1 次）
Rscript extdata/step2_SPAtests.R --bedFile=mid.bed ... --minMAF=0 --minMAC=1 --LOCO=FALSE
# C++ 版
./saige-step2 mid_cpp_single.yaml     # nThreads 8 / 1 两档
./saige-step2 mid_cpp_region.yaml     # gene-based，100 genes × 400 markers
./saige-step2 big2_cpp_single.yaml    # N 扩到 200k
```

配置文件都在本目录。差异说明：C++ 版缺失填充固定 mean，R 用的 best_guess；
p 值一致性 corr(−log10 p)=0.9989，max|Δlog10p|=0.72（全部来自 1% 缺失位点的填充方式差异，
与移植 README 的说明一致——用 --impute_method=mean 的 R 才能位级对齐）。

## 1. 单变异扫描总量（mid，40k markers，binary）

| 运行 | 墙钟 | markers/s | 峰值 RSS | CPU 秒 |
|---|---|---|---|---|
| R 版（1 线程） | 91.9 s | 435 | 490 MB | 85.2 |
| C++ T=1 | 75.2 s | 532 | 73 MB | 71.5 |
| C++ T=8 | 16.7 s | 2,390 | 148 MB | 125.6 |

C++ T=8 对 R 加速 5.5×；T=1→T=8 只有 4.5×（8 vCPU 下 CPU 秒从 71.5 涨到 125.6，
内存带宽/超线程竞争，不是锁——单变异路径已无 critical）。

## 2. 分段拆解（关键交付物，C++ T=1，插桩直接测量）

| 相位 | 秒 | 占墙钟 | 每次调用 | 说明 |
|---|---|---|---|---|
| (a) 基因型读取+解码 `getOneMarker_ts` | 19.0 | 25% | 0.48 ms/marker | fread 是 page cache（≈0），全是 2-bit→double 逐样本解码循环 |
| (a') 缺失填充/翻转/索引 `imputeGenoAndFlip` | 13.0 | 17% | 0.33 ms/marker | 又 1–2 遍全长 N 向量扫描（flip 时 `2-G` 整向量重写） |
| (b) score 检验（41,422 次含重算） | 16.4 | 22% | 0.40 ms/次 | scoreTestFast：数遍 N 长 dot |
| (c) SPA 校正（2,624 次） | 18.4 | 24% | **7.0 ms/次** | 根查找 + 每步 K(t) 全 N 遍 exp/log |
| (c') ER 精确检验 | 0 | 0% | — | MAC≤4 的 marker 数为 0（本数据无超稀有位点） |
| (d) 其他（case/ctrl AF 统计、输出、启动） | ≈8.4 | 11% | — | 启动+模型加载 <0.3 s |

计数：nSPA=2,624 = **1,312 个独立 marker（3.3%）× 2 遍**；nRecompute=1,422。
SPA 触发率 3.3% 与理论一致（SPAcutoff=2，|Z|>2 → ≈4.6%，本数据 null 无信号）。

**发现（免费的 12%）**：C++ 移植默认 `isFastTest=true`，p<0.05 的 marker 走第二遍
`Unified_getMarkerPval` "重算"——但本配置（无 sparse GRM、单 VR）下 ctx 与第一遍
**完全相同**，纯重复计算：T=1 下 recompute 9.8 s（13% 墙钟），其中 SPA 重复 ~9 s。
一行判断（ctx 不变则跳过）即可拿回。

## 3. IO vs 算力（比值判断）

- 每 marker 磁盘字节：N/4 = 12,500 B。裸读 bed（page cache）5.2 GB/s → 40k marker 全文件 0.1 s。
- C++ T=1 有效处理速度 6.6 MB/s，T=8 30 MB/s，都远低于磁盘（冷读 GCP pd ~200 MB/s，500 MB ≈ 2.5 s，
  仍只占 75 s 的 3%）。
- 解码相位 0.48 ms/marker = 每字节 ~38 ns = 每样本 ~9.6 ns（逐样本位抽取 + 间接寻址 + push_back）。

**结论：彻底算力/内存带宽绑定，不是 IO 绑定。** 每 marker 全流程要过 5–7 遍 N 长 double 向量
（解码 1、flip 0.5、索引扫描 1、score 2–3、统计 1），N=50k 时 0.4 MB/遍，L2 放不下。
降精度/融合遍数（解码+填充+score 一遍完成）比并行更根本。

## 4. gene-based（SAIGE-GENE+ 区域检验，C++）

配置：mid 全部 40k markers 切成 100 个"基因"×400 markers（`mid_group_400.txt`，anno=lof），
maxMAF ∈ {0.001, 0.01}，r_corr=0（SKAT-O），markers_per_chunk_in_groupTest=500，
ultra-rare collapse MAC≤10。**注意**：用的是单 VR dense null（无 sparse GRM，
nullmodel flagSparseGRM=false）——时间结构有效，p 值不是规范 SAIGE-GENE+ 输出
（规范要 cate-VR sparse null，本机没有对应 step1 产物）。通过 maxMAF 过滤后每 gene
只剩 ~7.6 个稀有 marker（dummy 数据 MAF 均匀分布，MAF<0.01 仅 765/40,000）。

| 运行 | 墙钟 | 峰值 RSS |
|---|---|---|
| T=1 | 89.5 s | — |
| T=8 | 19.0 s | 3.2 GB（P1/P2 = 2×500×N×8B×8 线程，与 BENCH_REPORT 模型一致） |

perf（T=1，35k 样本）按符号/模块：

| 归属 | 占比 | 是什么 |
|---|---|---|
| `getOneMarker_ts` | 18.5% | 解码（每 gene 读全部 400 个 marker，MAF 过滤在解码之后） |
| `imputeGenoAndFlip` | 13.9% | 同上第二遍 |
| kernel（page fault/清零） | ~25–30% | 每 gene 的 P1Mat/P2Mat（2×200 MB）resize+touch |
| libm exp/log | 15.6% | SPA（group 内 single assoc）+ Davies 被积函数 |
| libc（memset/memcpy） | 11.0% | arma 临时向量 |
| `qfc do_integrate`（Davies） | 7.8% | SKAT p 值积分 |
| **OpenBLAS（P1·P2 GEMM）** | **0.3%** | 每 gene 稀有 marker 太少，线代可忽略 |

**每 gene 的时间构成 = 读取/解码 ≫ 内存管理 ≫ SPA/Davies ≫ 线代。**
本数据下 gene-based 本质上是"把整个 bed 再解码一遍 + 每 gene 触碰 400 MB 内存"，
GEMM 不到 1%。真实 WES（每 gene 数百稀有变异）线代占比会升，但解码仍在关键路径
（marker 先解码后按 MAF 丢弃）。

## 5. 规模扩展（big2：N=200k，M=113k，C++ T=8）

null：`trend_big2_single.rda`（binary，单 VR=0.7558）→ `null_big2_arma/`。

| | mid T=8（N=50k） | big2 T=8（N=200k） | 比值（N×4） |
|---|---|---|---|
| 墙钟 | 16.7 s | 193 s | 11.6×（marker 也 ×2.8） |
| 每 marker 墙钟 | 0.42 ms | 1.71 ms | **4.1×，随 N 线性** |
| markers/s | 2,390 | 585 | |
| 峰值 RSS | 148 MB | 505 MB | 3.4× |
| CPU 利用率 | 7.5/8 | 7.4/8 | |

分相（thread-CPU 秒，占循环 CPU）：

| 相位 | mid | big2 | big2 占比 | 每次调用 big2（vs mid） |
|---|---|---|---|---|
| 解码 | 29.4 | 318 | 21% | 2.8 ms（3.8×，线性） |
| impute/flip | 18.1 | 213 | 14% | 1.9 ms |
| score（含重算） | 26.0 | 495 | 33% | 4.3 ms（6.8×，**超线性**：0.4 MB→1.6 MB/向量，cache 溢出+8 线程带宽竞争） |
| SPA | 35.5 | 335 | 23% | 53 ms/次（3.9×，线性）；触发 3,156/113,000=2.8%，仍 ×2 遍重复 |
| recompute（其中） | 19.2 | 184 | 12% | 纯浪费部分同 mid |

结论不随规模改变：解码+impute ≈ 35%，score+SPA ≈ 55%，全体随 N 线性（score 略超线性），
IO 占比始终可忽略（big2 bed 5.65 GB，冷读 ~30 s 也只占墙钟 15%，页缓存热读 ~1 s）。

## 6. 额外发现（顺手记录）

- **fastTest 重算全等浪费**：见 §2；big2 上 recompute=184 CPU s（12%）。
- **region 临时文件**：每个 gene 把 P1Mat/P2Mat chunk 无条件写盘
  （`*_regionK_P{1,2}Mat_Chunk_0.bin`，每 gene ~7 MB）且跑完不删——本次 2 个 region 运行
  留下 398 个文件 ~1.4 GB（已清理）。真实 WES 2 万 gene 会写 ~140 GB 垃圾。
- T=1 与 T=8 单变异输出逐字节一致（排序后 md5 相同）。
- C++ 输出与 R 输出 corr(−log10 p)=0.9989；差异全部来自 impute mean vs best_guess。

## 原始日志

- `mid_single_y4.log`（R + /usr/bin/time -v）
- `mid_cpp_single_t{1,8}.log` / `.time.log`（[PHASEPROF] 在 log 尾部）
- `mid_cpp_region_t{1,8}.log`，`perf_region_t1.data`
- `big2_cpp_single_t8.log`
