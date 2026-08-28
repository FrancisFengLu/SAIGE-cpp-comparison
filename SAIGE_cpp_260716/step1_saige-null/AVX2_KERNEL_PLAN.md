# AVX2 融合解码内核设计（ψv = A(Aᵀx) 热路径）

状态：设计稿，未实现。目标文件 `SAIGE_step1_fast.cpp`。

## 1. 现状解剖（已读代码核实）

- `CorssProd::operator()` (`SAIGE_step1_fast.cpp:1792-1800`)：对每个 marker i，
  `Get_OneSNP_StdGeno(i,&vec)` 解码出 N=50000 的标准化 `arma::fvec`，然后
  `val1 = dot(vec, m_bVec)`（趟1，沿样本归约），`m_bout += val1*vec`（趟2，沿 marker 累加）。
  两趟在**同一 marker 内背靠背**完成，packed 行（12.5 KB）天然驻留 L1。
- TBB 结构：`parallelCrossProd` (`:1911`) 构造 worker 后 `parallelReduce(0, Msub, CorssProd)` (`:2024-2025`)；
  Split 构造各线程私有 `m_bout`（N 个 fp32，零初始化），`join` (`:1803`) 向量相加；结果 `/Msub` (`:2069`)。
- `bVec`：`arma::fcolvec&` 共享只读；`m_bout`：每 worker 一份 N-fp32 连续数组。
- 解码 (`:609-689`)：`packed_byte(SNPIdx, i)` → `packed_flat_.raw()[SNPIdx*nbyte + i]` (`:89-90`，
  本移植 `numMarkersofEachArray=1` `:907`，行指针即 `raw()+SNPIdx*nbyte`，`nbyte=(N+3)/4=12500` `:292`)。
  每字节逐位拆 `b=bit0, a=bit1`，查 3 项表 `stdGenoLookUpArr(2-(a+b))` (`:651`)，
  表 = `(g − 2·freq)·invStd` (`:230-236`)。逐样本数据依赖索引 + arma 堆间接寻址 →
  `-O3 -march=native` 下编译器无法有效向量化，实测 ~90% 时间在此。
- 反例 `parallelCrossProd_blocked` (`:2098`)：把 N×B 解码块物化出 cache 再 GEMV，重造 15.3 GB/趟流量。
  **本设计不物化任何跨 marker 的解码结果**，替换的只是每 marker 的内层。

## 2. 数学重写（秩一分解，GPU 版已验证）

设 g∈{0,1,2} 为原始基因型（编码 g = 2 − popc(2bit)，missing 已填），f=freq[j]，d=invStd[j]。
- 趟1：val1 = d·(Σᵢ g·x[i] − 2f·Sx)，Sx = Σx **每次 ψv 只算一次**。
- 趟2：bout[i] += val1·(g−2f)·d = s·g − 2f·s，s = val1·d。均匀项累进标量 Coffset += 2f·s，
  reduce 结束后 `bout -= Coffset` 一次减掉。
→ 内核只需算 **Σ g·x** 和 **bout += s·g**，标准化完全移出内层。

## 3. 解码方案：vpshufb 路线（选定）

一次 `vmovdqu` 载 32 字节 = 128 个基因型。vpshufb 是 128-bit lane 内 16 项查表，
2-bit code∈{0,1,2,3} 只用低 4 索引，16 项表 `{2,1,1,0,2,1,1,0,...}` 用
`_mm256_setr_epi8` 复制到两个 lane 即可，无需 256 字节表。提取 4 个位平面用移位+掩码
（`vpsrlw` 右移 2k 后 `vpand 0x03`，跨字节窜入的位落在 bit6-7 被掩掉），每平面一次 pshufb：

```c
__m256i B  = _mm256_loadu_si256(p);            // 32B = 128 genos
__m256i m3 = _mm256_set1_epi8(0x03);
__m256i LUT= _mm256_setr_epi8(2,1,1,0, ...×2); // code→g
// k=0..3: gk 字节 i = 样本 4i+k 的 g 值（步长4交错）
__m256i gk = _mm256_shuffle_epi8(LUT,
             _mm256_and_si256(_mm256_srli_epi16(B, 2*k), m3));
// 展宽: gk 拆 4 段 8 字节 → vpmovzxbd → vcvtdq2ps → 8×fp32
```

**交错处理**：不在内核里 unpack 复原顺序（每 128 样本要 12 个 shuffle），而是
**每次 ψv 调用前把 x 重排一次**成 x_perm（128 样本超块内 k-major：
`x_perm[128s+32k+8j+m] = x[128s+4(8j+m)+k]`），bout 也按此序累加，reduce 后逆排一次。
两次 O(N) 重排代价对 M 个 marker 摊销为零；dot 与逐元素 FMA 对同构排列不变。

**cycle 估算**（每 128 基因型）：解码 1 load + 4×(srlw+and+pshufb)=13 uop；
展宽+FMA 每 32 字节 3 shuffle + 4 pmovzx + 4 cvt + 4 fma = 15 uop ×4 = 60 uop。
合计 ~0.57 uop/样本；瓶颈 port5（pshufb/pmovzx/extract 共 32 uop）→ **~4 样本/cycle/核**，
对比现标量 ~0.25 样本/cycle，单核解码 ~16×。N=50000 一趟 ~12.5k cycles ≈ 5.4 µs。

**否决的路线——位平面+popcnt/maskload**：popcnt 只能算无权和，趟1 是 Σg·x（x 任意 fp32），
必须把位扩成 dword 掩码（vpand+vpcmpeqd+vandps+vaddps），2 个位平面 ×4 相位
≈1.0 uop/样本，且多占 load 口；uop 数比 pshufb 路线高 ~75%，无精度收益。弃。

## 4. 两趟内核伪代码

趟1（4 个独立累加器掩 FMA 延迟；g 复解码两次而非存 200KB 中间量——packed 行在 L1，重解码 ~5µs 远低于 L2 往返）：

```c
__m256 acc[4] = {0};                       // k 相位各一
for (s = 0; s < N/128; ++s) {              // 每超块
  B = loadu(row + 32*s);
  for (k = 0; k < 4; ++k) {
    gk = shuffle_epi8(LUT, and(srli16(B,2k), m3));
    for (j = 0; j < 4; ++j)                // 8字节段
      acc[k] = fmadd(cvt_g8(gk,j), loadu_ps(xp+128*s+32*k+8*j), acc[k]);
  }
}
raw = hsum(acc[0]+acc[1]+acc[2]+acc[3]) + scalar_tail(...);
val1 = invStd[i] * (raw - 2*freq[i]*Sx);
```

趟2（s_val 广播；载入-FMA-回存，无水平归约）：

```c
__m256 sv = set1(val1 * invStd[i]);
for (每 8-fp32 段 t 同上遍历) {
  bp = boutp + off;                        // bout_perm 与 x_perm 同布局
  storeu(bp, fmadd(cvt_g8(gk,j), sv, loadu_ps(bp)));
}
Coffset += 2*freq[i]*val1*invStd[i];       // worker 标量，join 相加
```

尾部（N mod 128 = 80）：x_perm/bout 尾段保持自然序，纯标量循环（编译器可自动向量化），
天然跳过末字节 padding 位。

## 5. TBB 接法

`parallelReduce(0, Msub, worker)` 与粒度**完全不动**。改动：
① `parallelCrossProd` 在构造 worker 前算 `Sx` 与 `x_perm`（存于 worker 的 const 指针成员）；
② `operator()` 内层换成上面两个内核，行指针 `packed_flat_.raw() + i*nbyte`（断言 `use_packed_flat_`，否则走原标量路径）；
③ `join` 加 `Coffset`；reduce 后 `bout = unperm(m_bout) − Coffset`，再 `/Msub`。
`Get_OneSNP_StdGeno` 原样保留（VR / Diag / 调试路径还在用）。GPU 分支 (`:1916-1975`) 在前面早返回，无交互。

## 6. 数值验收

求和顺序变为 32 路树 + 秩一重排。fp32 树归约有效串行深度 ~N/32 → 单 marker val1
相对尾差 ~√(N/32)·ε ≈ 5e-6；秩一式在 Σgx ≈ 2f·Sx 时有相消，用绝对容差
`atol = 1e-5·invStd·(|Σgx| + 2f·|Sx|)` 兜底。验收三层：
① microbench 内 fp64 参考逐 marker 比 val1（rtol 1e-5 或 atol）；
② 加 `SAIGE_AVX2_VERIFY=1`：每次 ψv 同时跑旧标量路径，报 max rel-L2 diff，要求 < 1e-5；
③ 端到端：tau、VR 对基线 rtol 1e-4；PCG 迭代数允许 ±2。

## 7. 风险

- **行首对齐**：`nbyte=12500` 非 32 倍数 → 行首任意对齐，全部用 `loadu`（Skylake 上跨线惩罚 ~几%）。
- **尾字节 padding 位**可能非零 → 只由标量尾段覆盖，向量段止于 ⌊N/128⌋·128。
- **基线可能没那么慢**：`-O3 -march=native -funroll-loops` 已开，先 `perf` 确认标量解码占比
  （arma 3 项表的间接寻址大概率阻止了自动向量化，但要实测），收益预期以 microbench 为准。
- **带宽墙**：packed 总量 M×12.5KB 每次 ψv 从 DRAM 流一遍；4 核 n1 ~10-13 GB/s →
  解码提速 16× 后 ψv 上限约受带宽限制，**端到端预期 3-8×**，不是 16×。8 超线程未必优于 4 线程，两档都测。
- GCP 报 Skylake 但屏蔽 AVX-512 → 只依赖 AVX2+FMA，intrinsics 不用任何 512/VL 指令。

## 8. 实施顺序

1. **独立 microbench**（`tools/avx2_kernel/bench.cpp`，仿 GPU 版做法，不连主程序）：
   ```c
   // g++ -O3 -march=native bench.cpp; 无外部依赖
   // 1) 造 N=50000, M=2000 随机 packed 行(含尾部脏位) + freq/invStd + 随机 x
   // 2) ref: 原样标量解码+fp64 dot/axpy
   // 3) kernel_pass1 / kernel_pass2 (上文)
   // 4) 校验 ①val1 ②bout 逐元素; 计时 rdtsc → 样本/cycle, GB/s
   ```
2. 内核以头文件放 `tools/avx2_kernel/`，接进 `CorssProd`，env `SAIGE_USE_AVX2=1` 门控 + VERIFY 开关。
3. A/B：ψv 单趟耗时、tau/VR 端到端比对（examples 数据）。
4. 通过后再移植 `CorssProd_LOCO` (`:1853`，结构相同) 与 VR 路径。
