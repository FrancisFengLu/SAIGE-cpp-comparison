# step 2 C++ 移植优化 — 第二波（解码+填充/flip 单遍融合）2026-08-29

对象：`SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/`（工作区改动，未 commit）。
基线 = 第一波之后的工作区状态（WAVE1_LOG.md），同一二进制以
`SAIGE_STEP2_SCALAR_DECODE=1` 走旧链路即为基线（decode/impute 代码在回退模式下未动）。
改动文件：`genotype_reader.{hpp,cpp}`、`main.cpp`。

## 设计

把 `getOneMarker_ts`（逐样本 2-bit 解码，1 遍 N）+ `imputeGenoAndFlip`（flip 重写、
缺失填充、clean、sum、零/非零索引，再 2–3 遍 N）换成三段式，**dense N 向量只写一次**：

- **Stage A `getOneMarkerFusedStats_ts`**：读打包字节进 thread_local 缓存，只算 4 种
  2-bit 码的计数。fam 顺序 == 模型顺序（identity 映射，启动时检测并打印）时按 64-bit 字
  popcount：`hi=(w>>1)&0x5555…`，`lo=w&0x5555…`，missing=popcnt(lo&~hi)、
  het=popcnt(hi&~lo)、homref=popcnt(hi&lo)，homalt 用总数减；12.5 KB/marker 只需
  ~1.5k 次 popcount。尾部不满字节逐码处理（bed 填充位不计入）。非 identity 时逐样本
  gather 一遍并把码缓存给 Stage C（仍是单遍）。pre-impute 统计量
  （altFreq/altCounts/missingRate）沿用原函数的逐表达式——原来就是 counts 推的，
  **位级一致**，故 read 后的 QC 判定与基线完全相同。
- **Stage B `finalizeFusedStats`**（纯计算，无 N 遍）：把 imputeGenoAndFlip 的语义
  （altFreq>0.5 翻转、imputeG=best_guess/mean/minor、MAC+=imputeG·nMissing、
  MAC 门控的 `.clean(cutoff)`，arma 语义 |x|≤cutoff→0）整体重放在 **4 元素
  码→剂量表 fd[4]** 上，而不是 N 向量上。post-impute altCounts 用
  Σ fd[c]·counts[c]（固定顺序点积）——与原 `arma::sum` 的逐元素累加只在
  imputeG 非整数（mean 填充）时差最后 ulp，跨运行/线程完全确定。
- **Stage C `fillOneMarkerFusedDense_ts`**：每 marker 由 fd[4] 展开一张
  **256 项 byte→4×double LUT（8 KB，L1 常驻）+ 4-bit 非零掩码表**，对打包缓存单遍扫描：
  每字节一次 32 B 拷贝写出 4 个最终剂量（flip/填充/clean 全部已折进表），零/非零索引
  用掩码**无分支**发射进 thread_local 暂存，最后一次 memcpy 进输出
  `arma::uvec`（顺带取代了原 `copy_index_uvec_reuse`/`conv_to` 的二次拷贝）。

接线：

- 单变异标量路径（默认 blockSize=1 的热路径）：plink 时 Stage A → 原 pre-QC（表达式
  不变）→ Stage B → 原 post-QC → **只有过 QC 的 marker 才付 Stage C**（基线对 QC 失败
  marker 也要先付整遍解码）。
- region 路径：同上；且 **maxMAF/minMAF/minMAC QC 在 Stage B 的 counts 统计上判定，
  ~98% 被 maxMAF 筛掉的 marker 完全跳过 N 长解码**（本数据 100 gene × 400 marker 只有
  ~7.6/gene 存活）——这是 gene-based 2.2× 的主要来源。`GVec.zeros()` 在融合路径下省去。
- block 路径（blockSize>1，默认关）、BGEN/VCF/PGEN、UR 伪 marker、ldmat：未动，走原链。
- 回退：`SAIGE_STEP2_SCALAR_DECODE=1` 完整恢复旧两遍链路。

语义保真点（逐一核对过源码）：flip 判定用 pre-impute altFreq（严格 >）；缺失在 flip
之后写 imputeG（不再被翻转）；clean 的门 MAC 用**调用方自己的表达式**传入（单变异
`MAF*n*(1-mr)*2` 与 region `MAF*2*n*(1-mr)` 乘法顺序不同，各自保持）；clean 作用于
所有码含整数剂量；索引升序；ref-first 的 altFreq/altCounts 变换逐式复刻；
`genoMaps[MISSING]=-1` 的 sentinel（基线里是 (size_t)-1 转 double 的 1.8e19，随后必被
填充覆盖）不再产生。

## 验收（全部通过）

| 检查 | 结果 |
|---|---|
| mid 单变异 T=1，fused vs 基线 | **逐字节相同**（未排序 md5 相同） |
| mid 单变异 T=8，fused vs 基线 | **逐字节相同**；且 T=8 未排序输出 == T=1 输出 |
| gene 100-gene T=1/T=8 singleAssoc | **逐字节相同** |
| gene 100-gene 主输出 | 行/marker 集合相同，302 个数值格差异全为 fp 尾差：max rel Pvalue 1.5e-12，MAC 1.6e-14（counts 点积 vs arma::sum 的 ulp，≪1e-6 容差） |
| 确定性 | fused 两次运行（单变异 T=8、region T=1）逐字节相同；region T=1 与 T=8 排序后相同（基线亦然） |
| ER 路径（自造 300 个 MAC 1–4 位点 + 缺失，N=50k） | fused vs 基线**逐字节相同** |
| best_guess 填充（mid 与 ER 数据，nullmodel.json 改 impute_method） | 均**逐字节相同** |
| 非 identity 样本映射（fam 乱序，走 gather 回退，启动打印 "no (per-sample gather)"） | fused vs 基线**逐字节相同** |
| 回退开关 | `SAIGE_STEP2_SCALAR_DECODE=1` 输出即基线（上表所有"基线"都由它产生） |

## 计时（同机成对，本次实测；基线=同二进制回退模式）

| 运行 | 基线（回退） | fused | 加速 |
|---|---|---|---|
| mid 单变异 T=1（warm，两次） | 65.18 / 65.35 s | 34.98 / 34.81 s | **1.87×** |
| mid 单变异 T=8（warm，两次） | 14.33 / 14.32 s | 9.67 / 9.69 s | **1.48×** |
| mid 单变异 T=8（冷盘 drop_caches） | 14.83 s | 10.17 s | 1.46× |
| gene 100-gene T=1（warm） | 54.56 s | 25.12 s | **2.17×** |
| gene T=8（warm，两次） | 11.97 / 11.72 s | 6.53 / 6.44 s | **1.83×** |

RSS 不变（单变异 ~150 MB；gene T=8 ~1.74 GB）。对照 WAVE1_LOG 的第一波后数字
（T=1 70.4 s、T=8 14.8 s、gene T=8 13.2 s），第二波在其上再拿 ~1.5–2.2×；
目标 ≥1.3× 全配置达成。

## 融合后的瓶颈（perf，fused 单变异 T=1）

| 归属 | 占比 |
|---|---|
| `scoreTestFast`（N 长 dot 数遍） | 32.0% |
| SPA（exp/K1/K2/getadjGFast/dgemv 等合计） | ~15–20% |
| Stage C `fillOneMarkerFusedDense_ts` | 8.7%（Stage A popcount <1%，不进榜） |
| 主循环其余（case/ctrl AF 扫描、输出等） | 8.5% |
| OpenBLAS dgemv/ddot 杂项 | ~10% |

解码+填充相位从第一波后的 ~42%（25%+17%）压到 ~9%。下一个目标按序是
scoreTestFast（fast 路径可改用载体稀疏表示，MAF 低时 dot 遍数×N → ×nnz）和
SPA 的 K(t) 全 N 遍 exp。

## 备注

- 输出列 AF/MAC 在 mean 填充下理论上有 ulp 级差，但打印精度下单变异输出仍逐字节同；
  差异只在 gene 主输出的高精度列可见（≤1.6e-14）。
- region 的 maxMAF 早退用 post-impute 统计（counts 推出）判定，与基线判定值只差 ulp；
  本数据 marker 集合逐一相同（singleAssoc 逐字节同为证）。理论上 MAF 恰落在 cutoff
  ±1ulp 内才可能翻转，dummy 与真实数据均可忽略。
- b2 遗留的 `getOneMarker_carriers_ts`（未接线）其 clean 用了 `<`，与 arma `.clean`
  的 `<=` 不符——本波新代码用 `<=`；carriers 函数未修未用，供后续清理。
- 验证产物（yaml/输出/perf.data/ER 数据集）在本 session scratchpad `w2/`。
