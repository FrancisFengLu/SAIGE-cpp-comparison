# step 1 路径覆盖验收（2026-09-13）

在 GPU kernel 搬入（tier 4）+ 并行 loader 用于 VR 运行之后，逐条路径核对数值。
机器 saige-v100-dev（V100-16GB，8 vCPU / 4 物理核）。比较数字时 nthreads=1。

## 覆盖矩阵

| 路径 | CPU | GPU | 数值判定 |
|---|---|---|---|
| binary，非 LOCO | ✅ | ✅ | tau/VR **完全相同**，VR 结果文件 md5 逐字节相同（4 次运行） |
| binary，LOCO（smallloco，5 条染色体） | ✅ | ✅ | tau 0.408079 vs 0.408074，rel **1.2e-5**；VR rel 1.1e-6；两侧都产出 5 个 chr 目录 |
| quantitative，非 LOCO（mid，q1） | ✅ | ✅ | tau[0] 完全相同；tau[1] rel **2.7e-6**；VR rel 1.1e-6 |

CPU 侧 LOCO tau=0.408079 与 `OPT_LOG_202608.md` 的历史参照**完全吻合**。

## 差异的性质

GPU 与 CPU 的差在 1e-5 ~ 1e-6 量级，**与 CPU 自身的 run-to-run 噪声同量级**：
- `README.md` 已记录：4 线程下 fp32 归约顺序让 tau 在 0.2796 / 0.2336 之间跳（该测试数据）
- `OPT_LOG_202608.md:177` 记录：同一路径标量 vs AVX2 的 tau 差 1.1e-5
- 本轮实测：CPU 三次跑 mid 得 0.696524 / 0.696527 / 0.696527（spread 4.3e-6），而 GPU 三次**完全相同**

即 GPU 路径比 CPU 路径**更确定**（kernel 内归约顺序固定），差异来自 CPU 侧的不确定性而非 GPU 引入误差。

PCG 迭代序列在 CPU 与 GPU tier-4 之间**逐次相同**（18 次调用：1 1 1 1 1 1 7 4 6 8 7 7 9 7 8 9 7 7，合计 92）——
求解器收敛判定层面完全等价。

## 尚未覆盖（下一轮）

- sparse GRM 路径（本轮修了一个会在并行 loader 下抛异常的既有 bug，但没做数值验收）
- survival trait（按项目决定跳过）
- 分类 variance ratio（categorical VR）
- 条件分析
- step 2 的路径（本轮未改动 step 2）

## 复现

```
bash /opt/saige/logs/cpp/pathtest.sh          # quantitative CPU/GPU + binary LOCO
# LOCO 多染色体对照：
saige-null --config pathtest_loco_cpu.yaml
saige-null --config pathtest_loco_gpu.yaml --gpu
```
注意二进制当前仍需 `R_HOME`（嵌入式 R 运行时只用于 AI-REML 的随机数流，
计划做成 `--rng=r|std` 后即可去掉）。
