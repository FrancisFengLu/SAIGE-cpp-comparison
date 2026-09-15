# step 1 多表型容差门：阈值与标定（2026-09-15）

验收口径从「逐字节相同」放宽为「数学一致 + 数值容差」之后，这份文档记录
`thresholds.yaml` 里每个阈值是怎么定出来的、用什么数据标定的、以及它必须能抓到什么。
**改阈值就要同步改这里。**

被测二进制：`saige-null` @ `ef455ff`（`/opt/saige/logs/mt_gate/bin/saige-null-ef455ff-new`，
md5 74b6ddd2）。机器 saige-v100-dev（V100-16GB，8 vCPU）。所有标定运行 `nthreads: 8`。
数据：mid（N=50,000 × M=40,000）、small（N=10,000 × M=20,000），见 `make_cases.py`。

## 1. 阈值表

| 规则 (thresholds.yaml) | 度量 | 阈值 | 实测噪声上界 | 余量 | 真实 bug 偏差（§4） |
|---|---|---|---|---|---|
| `nullmodel.theta` (tau) | rel, abs_floor 1e-3 | 5e-5 | 1.9e-5 | 2.6× | 6e-4 … 7e-3（12–130×）|
| `nullmodel.alpha` | rel, abs_floor 1e-2 | 5e-5 | 7.0e-6 | 7.1× | 3e-3 … 6e-3 |
| `arma.mu` | rel, floor_frac 1.0 | 5e-5 | 7.4e-6 | 6.8× | 1.3e-3 … 9.5e-3（25–190×）|
| `arma.res` | rel, floor_frac 1.0 | 5e-5 | 1.1e-5 | 4.6× | 1.2e-3 … 9.7e-3 |
| `arma.S_a` | rel, sa_floor_frac 1.0 | 5e-5 | 6.2e-6 | 8.1× | 见 §4 说明（不敏感）|
| `arma.V` / `arma.XV` | rel | 2e-5 | 1.6e-6 | 12.7× | 5.5e-4 … 4.2e-3 |
| `arma.XVX` / `XVX_inv` | rel | 2e-5 | 1.2e-6 | 17× | 5.4e-5 … 2.5e-3 |
| `arma.XXVX_inv` / `XVX_inv_XV` | rel, floor_frac 1.0 | 2e-5 | 1.5e-6 | 13× | XVX_inv_XV 5.5e-4 … 4.2e-3；XXVX_inv 在 e1/e2 只有 1.1e-5 / 1.6e-5，**没到阈值** |
| `grm_diag` | rel, print_sig 8 | 5e-5 | 0（同实现）/ 1.2e-5（换实现，见 §3.3）| 4× | 3.6e-4 … 8.3e-3（7–165×）|
| `vr.ratio` | rel, print_sig 6 | 2e-5 | 0（打印值完全相同）| — | 4.6e-5 … 6.5e-4（2.3–33×）|
| `vr_markers.var1/var2/ratio` | rel, print_sig 6 | 2e-5 | 2.5e-6 | 8× | 5.6e-5 … 1.3e-3 |
| `vr_markers.AF` | rel, tol 0 + 打印位 | 0 | 0 | — | 3.4e-4（整数计数变了）|
| `vr_markers.pvalue` | log10p | 1e-4 | 未实测（step 1 不输出 p 值列）| — | — |
| `nullmodel` 其余字段、`sampleIDs`、`X`/`y`/`offset`、`m_cov_*.csv`、`sparseGRM_*`、`loco_chroms`、文件清单、VR 测试 marker 数、`SNPIdx`/`MAC` | exact | 完全相等 | 0 | — | 结构性失败 |
| `log.iterations` / `log.converged` / LOCO 行 | exact | 完全相等 | 0（1,272 次比较无一次不同）| — | — |

「实测噪声上界」= §3 所有标定运行里该规则出现过的最大偏差；「真实 bug 偏差」= §4 数据级
bug 模拟里该规则的偏差范围。阈值取在两者之间：噪声之上留 3–13 倍余量，真实 bug 之下
至少 1 个数量级（最弱的通道 `vr.ratio` 是 2.3 倍，见 §4 的说明）。

## 2. 度量定义（compare.py）

- `exact`：解析成数值后必须完全相等（`1e-5` 与 `1.0e-05` 视为相等）。
- `rel`：`dev_i = max(|a_i − b_i| − slack_i, 0) / max(|a_i|, |b_i|, floor_i)`，取 `max_i dev_i`。
  - `floor_i = max(abs_floor, floor_frac × RMS(参考侧该列))`。近 0 的元素不应该把相对差放大成
    无意义的大数：`floor_frac: 1.0` 表示「跟这一列自己的尺度比」，`0.1` 表示「只压住最小的
    10%」，`0` 表示纯相对（值本身量级固定时用，如 grm_diag ≈ 1、VR ≈ 1、tau）。
  - `slack_i`：来自文本的数值，C++ 默认流精度是 6 位有效数字（grm_diag 是 8 位）。两个正确
    舍入的打印值之差最多是真值之差再加 1 个末位单位，所以先减掉 1 个末位单位再做相对比较。
    **副作用：VR 末位（第 6 位有效数字）改 1 是检不出来的**，这是打印精度的下限，不是容差
    松了（§5 的注入测试把这条记录成用例）。
  - `sa_floor_frac`：`S_a = colSums(X ⊙ res)` 在收敛点接近 0，用它自己做分母全是噪声
    （实测两次 CPU 运行 S_a[0] 可以是 0.00225 vs 0.00182）。改用「score 的绝对质量」
    `sum_i |X_ij res_i|` 做 floor，于是 dev 读作「score 相对于它的总量变了多少」。
- `log10p`：`max |log10 a − log10 b|`，只用于名字像 p 值的列。step 1 现在不输出 p 值列，
  这条规则只有 §5 的合成测试走到过。
- 稀疏 GRM（`sparseGRM_locationMat/valueVec.arma`）按 (row, col, value) 三元组集合比，
  顺序不同不算错（另外报告顺序是否一致），值必须完全相等。

## 3. 标定数据

「数学一致但非逐字节」的真实样本用两种来源：

### 3.1 lockstep（`fit.multi_lockstep: true`）

同一个二进制、同一份数据，多表型锁步把一组表型的 AI-REML 迭代合并成一次批量多 Sigma
PCG 求解，psi·B 的归约顺序与逐表型路径不同 —— 正是方案 C 会引入的那类差异。

### 3.2 CPU nthreads=8 的运行间不确定性

`--solo-vs-solo` 把每个表型的单独跑重跑一次再比。实测：

- GPU（tier=4）非 LOCO：两次运行**逐字节相同**；分组多表型跑与单独跑也逐字节相同。
- GPU LOCO：**不确定**，solo vs solo 就有 ~1e-6（mu/res）、~3e-6（XVX_inv_XV）的差。
  LOCO 的 per-chr 重算走 CPU 多线程归约。
- CPU nthreads=8：不确定，量级同上，quant 的 res 最大到 1.0e-5。

### 3.3 换实现算同一个量（grm_diag）

grm_diag 在分组多表型里跟单独跑逐字节相同（噪声 0），所以它的阈值不能从上面两种噪声里
标。改用独立实现对照：用 numpy float64 按 `genoqc.py` 的 QC/填补规则重算
`(1/M)·Σ_j z_ij²`，与 C++（fp32 累加 38,262 个 marker）比：

```
mid nomiss/b1 (n=50,000)  max rel 1.10e-5   median 1.9e-6
mid fill/fS   (n=5,000)   max rel 1.22e-5   median 4.1e-6
```

即「同样的数学、不同的算术」在 grm_diag 上就有 1.2e-5；阈值 5e-5 = 4 倍余量。

### 3.4 标定运行清单

| label | 设备/规模 | 模式 | 比较份数 |
|---|---|---|---|
| `mid-gpu-t8-74b6ddd2` | GPU / mid | 分组 | 46 |
| `mid-gpu-lock-t8-74b6ddd2` | GPU / mid | lockstep | 46 |
| `small-gpu-t8-74b6ddd2` | GPU / small | 分组 + solo-vs-solo | 92 |
| `small-gpu-lock-t8-74b6ddd2` | GPU / small | lockstep | 46 |
| `small-cpu-t8-74b6ddd2` | CPU / small | 分组 + solo-vs-solo | 92 |
| `small-cpu-lock-t8-74b6ddd2` | CPU / small | lockstep | 46 |
| `mid-cpu-t8-74b6ddd2` | CPU / mid | 分组 | 46 |
| `mid-cpu-lock-t8-74b6ddd2` | CPU / mid | lockstep | 46 |

每份比较 60–120 个 item。汇总命令：

```
tests/mt_tolerance/calibrate.py --recompute <label> [<label> ...]
```

（`--recompute` 只重新跑比较器，不重跑 saige-null；改完 `thresholds.yaml` 用它复核。）

## 4. 真实 bug 的偏差（emulate_bugs.py）

方案 C 的典型 bug 不能靠改输出来模拟，所以在**数据层**制造同样的症状，用当前二进制单跑
受影响的表型，再与它正确的单跑比。这给出的是「门要抓的东西有多大」：

| 模拟 | 做法 | theta | mu | grm_diag | VR (bin 均值) | 门的判定 |
|---|---|---|---|---|---|---|
| `e1_grm_drop` | qA3 的 GRM 少掉 76 个「在自己样本上过 QC、在并集上不过」的 marker（共 38,281 个）| 2.3e-3 | 1.3e-3 | 8.3e-3 | 4.6e-5 | FAIL（13 项）|
| `e2_fill` | fS 在 950 个 marker（64,507 个缺失格）上用并集的填补值 | 6.2e-4 | 9.5e-3 | 3.6e-4 | 8.0e-5 | FAIL（14 项）|
| `e3_one_more` | t10 多算 1 个样本 | 6.7e-3 | 结构 | 结构 | 6.5e-4 | FAIL（30 项）|
| `e4_one_less` | lS（LOCO 开）少算 1 个样本 | 3.7e-3 | 结构 | 结构 | 4.3e-4 | FAIL（75 项）|

说明：

- e1/e2 是方案 C 最该防的两类：「用并集的 QC 名单」和「没按表型修正填补值」。两者在
  theta / mu / grm_diag 上都超阈值 1–2 个数量级。
- e3/e4 是样本集算错（掩码或分组 bug）：`n`、`sampleIDs`、各矩阵形状直接不一致，
  比较器先在结构层失败。
- **`S_a` 对这四类 bug 都不敏感**（e1 5.5e-6、e2 3.2e-5，都在容差内或刚好压线）：
  收敛点上 score 本来就接近 0，bug 把它推动的量级和 fp32 求和误差同级。S_a 的容差按噪声
  定，但**不要指望它单独发现问题**。
- **最弱的通道是 `vr.ratio`**：e1 只有 4.6e-5（阈值 2e-5，2.3 倍）。VR 是 30 个 marker 的比值
  均值，本身对 tau 不敏感，而且只有 6 位有效数字。VR 单独看抓不住小 bug；它在门里的价值
  是 e3/e4 这种大偏差和 marker 名单变化（测试 marker 数、SNPIdx、MAC、AF 都是 exact）。

复现：`emulate_bugs.py <bin> --gpu`（结果在 `/opt/saige/logs/mt_gate/emulate/`）。

## 5. 比较器自身的故障注入（test_compare.py）

拿一份已知通过的输出，人为制造症状，确认比较器给出预期判定。25 项全部符合预期：

| 注入 | 期望 | 实测 dev / tol |
|---|---|---|
| 完全不改 | PASS | — |
| tau 改 1e-4 相对量 | FAIL | 9.9e-5 / 5e-5 |
| tau 改 6e-5（刚过阈值）| FAIL | 5.9e-5 / 5e-5 |
| tau 改 2e-6（噪声级）| PASS | 1.9e-6 |
| mu 0.1% 元素改 1e-4 | FAIL | 1.0e-4 / 5e-5 |
| mu 全元素 ±2e-7（fp32 噪声）| PASS | 2.0e-7 |
| VR 第 4 位有效数字 +1（1e-4）| FAIL | 1.0e-4 / 2e-5 |
| VR 第 5 位有效数字 +1（1e-5）| **PASS** | 9.3e-6 —— 记录灵敏度下限，见 §2 |
| VR 第 6 位（末位）+1 | PASS | 0（在打印舍入以内）|
| 删一个文件 / 多一个文件 | FAIL | 文件清单 |
| iterations +1 / converged → NO | FAIL | exact |
| sampleIDs 交换两个 | FAIL | exact |
| grm_diag 某个值改 1e-4 | FAIL | 1.0e-4 / 5e-5 |
| .arma 截断 8 字节 | FAIL | 解析失败 |
| X 某一格改 1e-6 | FAIL | exact（输入不允许变）|
| S_a[0] += 1e-3（fp32 求和噪声）| PASS | 2.1e-7 |
| S_a[0] += 1e-4 × score 质量 | FAIL | 1.0e-4 / 5e-5 |
| 加一列相同的 p.value（对照）| PASS | 0 |
| p 值 log10 改 1e-3 | FAIL | 1.0e-3 / 1e-4 |
| VR 测试 marker 数 30 → 31 | FAIL | exact |
| LOCO：不改 / chr 目录缺一个文件 / chr2 mu 改 1e-4 | PASS / FAIL / FAIL | 文件清单、1.0e-4 / 5e-5 |

## 6. iterations / converged 的规则

默认**严格相等**。标定的 1,272 次比较里 iterations、converged、LOCO 行、文件清单
**没有一次不同**，包括 lockstep 和 CPU 多线程。

如果以后真出现阈值附近的翻转（AI-REML 的 `rc_tau: X (converge if < 0.02)` 停在 tol 附近，
或 VR 的 `CV=... <= 0.001` 压线）：

1. 先确认是翻转不是 bug：看两边日志里最后一步的 `rc_tau` 与 tol 的距离、以及 tau 本身的差
   是否还在 5e-5 以内；
2. 确认后在 `cases/<case>.yaml` 里写
   `known_flips: {<trait>: {items: [iterations], reason: "..."}}`；
3. 门会把它显示成 `FLIP*` 并在行尾写明翻转内容，**不会**变成普通的 PASS，也不会放宽任何
   数值容差。

数据侧已经先做了规避：`make_cases.py` 给小样本子集用强遗传信号的二值表型，并且避开了
tau1 刚好卡在 tol（0.02）附近的抽样 —— `qc` 用例的 qA3 就是因为 tau1 = 0.0205 太贴边而
换了 trait 种子（`QC_TRAIT_SEED_BUMP`）。`nondegenerate` 检查保证每个用例的每个表型
tau1 > 0，否则 GRM 根本不进入拟合，bug 会被 tau1 = 0 掩盖。

## 7. 没覆盖 / 已知限制

- **p 值列**：step 1 的 VR marker 结果文件只有 `SNPIdx MAC AF var1 var2 ratio`，没有 p 值。
  `log10p` 规则写好了但只有合成数据测过，阈值 1e-4 是按 VR 噪声推的，**未实测**。
- **VR 的灵敏度**：1e-5 相对量级的 VR 变化在容差内（6 位有效数字 + 2e-5 阈值）。
- **`S_a`**：对本文档模拟的所有 bug 都不敏感（§4）。
- **survival trait、条件分析、`--blocked-gemv`、`lowmem_loco`、`inv_normalize`、
  `whitelist_ids`、sex-specific fit** 没有进用例集。
- CPU 的 `nthreads: 1` 没有单独标定（那是逐字节路径，由 `run_p1_byte_gate.sh` 覆盖 P=1）。
- 阈值是在 mid / small 两个规模上标的。UKB 规模（N ~ 4e5）下 fp32 累加误差会更大，
  tau 的噪声大概率也更大，到时要重标。
