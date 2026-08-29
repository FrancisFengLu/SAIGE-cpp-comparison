# step 2 C++ 移植优化 — 第一波（修 bug + 白捡）2026-08-29

对象：`SAIGE-cpp-comparison/SAIGE_cpp_260716/step2_saige-step2/`（未 commit，工作区改动）。
基线 = 改前源码 + 预期 flags（`-std=c++17 -O3 -fopenmp`，即 PROFILE.md 测量所用配置）。
验证均用 `/opt/saige/data/mid` + `null_y4_arma`（binary，单 VR，无 sparse GRM）。

## 改动清单

### ① fastTest 全等重算去除（W1-1）

- `main.cpp` 标量路径：ctx_first 构建后**先验**计算重算会用的 ctx（只依赖 MAC 和模型常量），
  三元组 (flagSparseGRM_cur, isnoadjCov_cur, varRatioVal) 与第一遍完全相同 → 跳过第二遍
  `Unified_getMarkerPval`；同时 `g_firthDefer` 加 `!fastRecomputeSameCtx`（跳过重算时 Firth
  改回第一遍内联执行，A3 语义保持）。block 回退路径（Phase C）加同样的 sameCtx 跳过。
- 上游语义依据：R/SAIGE 在无 sparse GRM 时直接禁用 fastTest
  （SAIGE_parallel_Test_main.R:346 "No sparse GRM is specified, so is_fastTest is not working"）；
  有 sparse GRM / isnoadjCov / cate-VR 时 ctx 不同，重算照常执行——不是无脑砍。
- 验收：mid 单变异 T=8 输出与基线**逐字节相同**；墙钟 17.45→14.81 s（同 flags，−15%）；
  T=1 warm 78.4→70.4 s（−10%，重算相位 PROFILE 实测 9.8 s ≈ 全部拿回）。

### ② Makefile flags（W1-2）

- `CXXFLAGS ?=` → `:=`（conda activate 导出的 CXXFLAGS=-march=nocona -O2 **无 -fopenmp**
  曾把整个并行编译成 no-op），补 `-march=native -funroll-loops -DARMA_NO_DEBUG`（step1 同款）；
  `REQUIRED_DEFS` override 追加保护必需 define；INCLUDES/LIBS 加 `$(CONDA_PREFIX)`
  include/lib+rpath（原来靠 env CXXFLAGS 里的 -isystem 偶然编过）。
- 验收：`make` 输出编译行确认 `-O3 -march=native -fopenmp` 生效；mid 单变异 40,000 marker
  集合相同，p 值 max rel diff = **0**（逐字节相同）；region 输出 104 个数值格有差，
  max rel 1.6e-13（fp 重排，< 1e-8 容差）。

### ③ region 路径（W1-3）

- (a) 结论：P1/P2 chunk 文件**不是纯调试残留**——gene 内 marker 数超过
  markers_per_chunk_in_groupTest（或 常规+ultra-rare 各成一 chunk）时，多 chunk VarMat
  组装会 load 回来。修法：只在"确定会被读回"时才写（末 chunk 写盘条件 = 已有其他 chunk
  或 UR 可能新增 chunk），并记录所有写过的文件、组装后**无条件清理**；
  `SAIGE_STEP2_DUMP_KERNELS=1` 恢复旧的全写+保留行为（已验证：设 1 时留 200 个文件）。
- (b) P1Mat/P2Mat 改为 caller 侧 `thread_local` 持久缓冲（`set_size` 复用，首个 gene 后
  零分配零清零），region 函数内去掉截断重赋值（改 `head_rows/head_cols` 视图），
  多 chunk 组装用局部矩阵 load，不再破坏持久缓冲容量。
- 验收：gene-based（100 genes）T=1 主输出与 singleAssoc 输出均与基线**逐字节相同**
  （T=8 行序本身不定序——基线两次 T=8 也互不相同，排序后 0 diff）；运行目录 0 个遗留
  chunk 文件（基线 188 个）；T=1 91.1→64.1 s（−30%），T=8 19.2→13.6 s，RSS 3.24→1.74 GB。

### ④ ER 可复现性（W1-4）

- `er_binary.cpp`：废除 `base ^ omp_get_thread_num()` 播种；新增 `SL_set_stream(uint64)`，
  引擎按 splitmix64(base ⊕ splitmix64(stream)) 确定性重播种。
- stream id 走 `PerMarkerCtx.erSeedStream`（新字段），调用点：单变异标量/block 路径 = 全局
  marker 索引 i+1；region 常规 marker = hash(regionName) ⊕ (i+1)；UR 伪 marker =
  hash(regionName) ⊕ (i_ur+1)。`getMarkerPval` ER 分支进入 `SKATExactBin_Work` 前
  `ER::SL_set_stream(ctx.erSeedStream)`。与线程/调度完全解耦。
- 验证数据（scratch 生成，mid.fam 的 5 万样本）：
  - `er`（400 个 MAC 1–4 位点，MACCutoffforER=4）：两次 T=8 运行**逐字节相同**。
    注意：MAC≤4 时 2^k ≤ 16 ≪ 2e6，ER 走精确枚举，RNG 实际不参与。
  - `er2`（200 个 25–35 载体位点，MACCutoffforER=40）：确认走重采样（临时计数：单次运行
    RNG 抽样 3.52 亿次），13 个位点 p.value≠p.value.NA（ER 生效）。新二进制两次运行、
    以及 T=1 vs T=8 输出均**逐字节相同**（播种与调度无关的直接证据）。
  - 旁注：基线在该数据上 T∈{1,2,3,5,8} 也碰巧一致（旧代码每次 ER 前 `SL_setseed(1)` 把
    种子重置为 1⊕tid，且该单变异 ER 构型下重采样贡献恰好对种子不敏感）；设计层面的
    调度依赖（seed 含 tid、epoch 全局竞态）确实存在，本修复将其消除。
- mid 单变异/region 输出不受影响（无 ER 触发位点，改后逐字节相同）。

## 计时收尾（全部改动后，saige-step2.final）

冷盘 = `sync; echo 3 > /proc/sys/vm/drop_caches`。

| 运行 | 改前基线 | 改后 run1 | 改后 run2 | Δ |
|---|---|---|---|---|
| mid 单变异 T=1（冷） | 75.2 s | 72.76 s | 72.75 s | −3%（冷读 ~2.5 s 抵掉部分收益） |
| mid 单变异 T=1（warm 对照对） | 78.4 s（基线同机复测） | 70.4 s | — | **−10.2%** |
| mid 单变异 T=8（冷） | 16.7 s | 15.96 s | 15.80 s | −5% |
| mid 单变异 T=8（warm 对照对） | 17.45 s（基线同机复测） | 14.81–14.82 s | — | **−15.1%** |
| gene-based T=8（冷） | 89.5 s（PROFILE 记录为 T=1；T=8 基线 19.2 s） | 13.22 s | — | **−31%**（vs T=8 基线） |
| gene-based T=1（warm） | 91.1 s（基线同机复测） | 64.07 s | — | **−29.7%** |

RSS：单变异 T=8 149 MB（不变）；gene T=8 3.24 GB → 1.74 GB。
基线数字含同机复测值（PROFILE 的 75.2/16.7/89.5 与本次复测有 ±4% 机器噪声，
同机成对比较以复测值为准）。

## 备注

- 输出正确性验证顺序：①③④在旧 flags 下逐字节比对通过后，才切 ② 的新 flags
  （新 flags 下单变异仍逐字节同，region 仅 ≤1.6e-13 fp 重排差异）。
- region T=8 输出**行序**本身不可复现（omp dynamic + critical 写出，改前已如此），
  内容排序后一致；如需行序确定性属后续项。
- `tests/run_loco_tests.sh` 因 test_data 目录不存在而 SKIP。
- ldmat.cpp 也有同款 chunk 文件写盘（LD 矩阵模式），未动（本波范围外）。
- 验证产物在 scratchpad（session 级）：基线/各阶段二进制、er/er2 数据集与 yaml。
