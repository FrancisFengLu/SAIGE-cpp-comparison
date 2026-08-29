# SAIGE step 2 优化机会评估（只读代码分析，2026-08-28）

对象：用户 C++ 移植 `SAIGE_cpp_260716/step2_saige-step2/`（对照上游 fork `SAIGE-work/src/`）。
实测参考：R-fork `mid_single_y4.log`（N=5e4、M=4e4、单线程 **1.9–2.6 ms/marker**）；C++ 移植 `mid_cpp_single_t8.time.log`（8 线程主循环 16.7 s / 4e4 marker = **0.42 ms/marker**，CPU 利用 741%/800%，RSS 145 MB）——移植 8 线程对 R-fork 单线程已是 ~5–6× 墙钟。但注意 CPU 时间 114.9 s / 4e4 ≈ **2.9 ms·core/marker**，与 R-fork 单线程同量级：并行有了，单核成本还没降，第 2 项（解码）就是要砍这块。分段占比待 PROFILE.md 校准，下述收益为估算。

## 1. 结构事实（a–e）

| 项 | 事实 | 出处 |
|---|---|---|
| (a) G 每 marker 用几次 | 单变异：读+解码 1 次、用 1 次即弃（流式，无驻留复用）；fast 路径只碰非零载体 O(k·p)。region：g̃ 存入 P1Mat/P2Mat，gene 内复用 ~m 次 | main.cpp:1143-; saige_test.cpp:305-,866- |
| (b) 多表型可共享量 | G/解码/QC/impute/非零索引可共享（上游 fork 的 saigeObjects 循环已这样做，Main.cpp:542）；per-pheno 私有：res、mu、mu2、XV、varRatio、SPA。移植目前单表型 | SAIGE-work/src/Main.cpp:542-685 |
| (c) gene-based 矩阵规模 | VarMat=P1(m×N)·P2(N×m)，m≤chunk(500)，gene 内通常<100。FLOP/byte=m/8≈6–60，但总量 m=100、N=5e5 才 1e10 FLOP/gene——CPU 亚秒级，不构成 GPU 目标 | Main.cpp:1327,1900 |
| (d) SPA 可否批量 | 仅 \|Z\|>2 触发（~5% marker）；每次 2 个根查找×5–20 次牛顿迭代，fast 版每迭代 O(k)。分支重、长度不齐，批量收益 ∝ 触发率——不值得 | spa_binary.cpp:63-,242- |
| (e) LOCO 切换 | 换 chr 只换 res/mu 向量，O(N·p)×22，可忽略。多表型×LOCO 时 22·P·N·8B 内存才是约束（需按 chr 惰性加载） | LOCO_FORMAT.md |

量纲：T=Gᵀr 每字节基因型（fp64 dosage）只做 2 FLOP → **0.25 FLOP/byte**，读一遍用一遍。整条流水线（读→解码→QC→O(kp) 打分）层层带宽/IO 绑定，ALU 大量空闲：纯打分 FLOP 只值 ~0.05 ms，实测却 ~2 ms/marker。

## 2. step 2 移植的优化清单（按利润排序，前 5 项）

先对照 step 1 清单的体检结果：**调试 IO** 热循环已干净（fusedMode 默认 0，A/B 比对及其 omp critical 不生效；checkpoint 只写 marker0）；**分配** 已治（thread_local 全套 + mmap 阈值修复，注释记录了 1.44B 页错误的来历）；**并行** 已有 `omp parallel for schedule(dynamic,64)` + per-thread FILE*；**确定性** 无 random_device，但见第 5 项。

| # | 改什么 | 预期收益 | 依据 | 工作量/风险 |
|---|---|---|---|---|
| 1 | **Makefile:9 `CXXFLAGS ?=` 换成 `+=`/覆盖保护，并补 `-march=native`** —— 与 step 1 完全相同的 conda 环境变量覆盖坑；且当前无任何 -m 开关，armadillo 全程走标量 SSE2 | 白捡 ~8%（step 1 实测）；若 conda 环境真把 -O3 顶掉则远不止 | Makefile:9 原文 `CXXFLAGS ?= -std=c++17 -O3 …` | 半小时，零风险 |
| 2 | **解码：两遍合一遍 + AVX2 vpshufb 移植** —— PLINK 路径每 marker 扫两遍 packed buffer（counts 遍 + emit 遍，genotype_reader.cpp:454/489），逐样本标量抽位还带 `m_posSampleInPlink[i]` 间接寻址；默认 identity 映射（:800）可走 step 1 的 `tools/avx2_kernel/avx2_kernel.hpp` 快路径 | 解码估占每 marker 30–50%（未实测，待 PROFILE.md）；两遍→一遍立得 ~2× 解码，AVX2 再 5–10×；端到端估 1.3–2× | step 1 ψv 融合解码 10× 的同款模式 | 2–3 天；风险=missing/flip/子集三条分支要进 SIMD 尾部处理，用 A/B 位一致验证（框架已有 fusedMode=1） |
| 3 | **生产者-消费者流水线：读线程预取 + 算线程池** —— 已写好的 `scoreTestFast_block` BLAS-3 路径因 reader 串行 IO 反而慢 2.2×（main.cpp:108-117 尸检注释），VCF 路径整个读取在 `critical(genoread)` 里串行；双缓冲（读第 i+1 块与算第 i 块重叠）同时救活 block-GEMM 和冷缓存流式 | 把已回归的 blockSize>1 路径从 -2.2× 翻成正收益；VCF 多线程扩展性从被 critical 封顶解放出来 | main.cpp:1421, 108-117 | ~1 周；风险=marker 顺序与输出顺序解耦，需保序写出（现有 pvalVec[i] 按索引存，已天然支持） |
| 4 | **SPA 触发路径去重** —— 触发 SPA/Firth/条件分析时 `getadjGFast` 重建 g̃（scoreTestFast 已算过等价量）；p 值 sprintf→string→stod 往返每 marker 一次（fastTest 判据用）。均为小头 | 每项 ~1–5%，仅在二元 trait 高触发率（不平衡表型 cutoff 放宽）时值得 | saige_test.cpp:1155-;主循环 stod | 1 天；低风险 |
| 5 | **ER 分支确定性** —— `er_rng()` 是 thread_local 默认种子 mt19937，重采样序列取决于 marker→线程分配（dynamic 调度），ER 触发的 p 值**跑两次可能不同**。改为按 marker 索引播种 | 不是速度，是可复现性缺陷修复 | er_binary.cpp:57-64 | 半天，零风险 |

另：region 检验的 `max_markers_region` 配置分裂（per-thread 246 GB 峰值内存，已修，28–48×）和上游 VarMat 强制对称化 bug 见 BENCH_REPORT.md，属已完成/待上游反哺。

**未来方向（降级）**：多表型 PheWAS 共享——上游 fork 已实现"G 读一遍服务 t 个 trait"（但同时把 `omp_set_num_threads(1)` 写死、循环无 pragma，Main.cpp:341，等于自废并行）；移植若加多表型，正确形态是 `scoreTestFast_block` 加一维成 (B marker × P pheno) GEMM，量化 trait 无 SPA 可全批量。

## 3. GPU 判断

**结论：常规单表型 step 2 不值得上 GPU（no-go）；仅当"量化 trait 的 PheWAS、P≳32、基因型块一次 H2D 服务全部 P 个表型"时才 go。**

论证：① 单表型可卸载的只有 0.25 FLOP/byte 的流式点积，数据必经 PCIe 一次（~12 GB/s）≈ 单路 DRAM 带宽，GPU 无算力优势可兑现，而墙在解码/IO（CPU 侧）——按 20% 门槛直接出局。② 批 P 个表型后强度 ∝P/4 FLOP/byte：P=100 → 25，进入算力绑定；P=1000（蛋白组）、N=5e5、M=8.6M → ~10¹⁶ FLOP，V100 fp64 ~30 min vs CPU ~10×之外，GEMM 占比 >80%，Amdahl 撑得住。③ 二元 trait 的 SPA/Firth/ER 长尾只能留 CPU，触发次数 ∝5%·P，P 大时重新成为瓶颈——GPU 方案实际只适用量化 trait。M~10⁷ 基因型也放不进显存（2-bit 仍 ~1 TB），"驻留复用"在 step 2 无对象。

## 4. 与 step 1 的关键差异

step 1 全部收益建立在同一份 GRM 基因型被复用 O(PCG 迭代×探针×表型) 次上——2-bit 常驻、片上解码、批量 RHS 都是"驻留+高复用"的红利。step 2 结构相反：基因型是**流**，每 marker 用一次即弃，复用次数=表型数 P。所以 step 1 的招里能搬的是两条：**AVX2 融合解码**（对象从 ψv 换成 bed 解码+载体收集，本表第 2 项）和**多表型批量**（把复用从 1 提到 P，唯一能把算术强度推过带宽墙的杠杆，现为未来方向）；搬不动的是 2-bit 驻留和 GPU 常驻 GEMV。step 2 的"隐藏大头"不在数值核，而在解码两遍、串行 IO、被写死的单线程（上游 fork）——先修 host、后谈 kernel，与 step 1 的教训完全同构。
