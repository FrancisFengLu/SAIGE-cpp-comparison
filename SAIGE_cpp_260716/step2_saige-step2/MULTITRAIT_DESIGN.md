# step 2 多表型批量 — 设计文档

状态：设计（未实现）。读者：实现者 + 后续维护者。
基底：`step2_saige-step2`（上游 1.5.1 的 R-free C++ 移植，二进制 `saige-step2`）。
本文所有行号指本目录下的文件，对应 2026-09-13 的 HEAD。

本设计**不参考**上游 R 包的 multitrait 分支。那条分支对每个 (marker × trait) 重建子矩阵，
比上游单 trait 主线还慢 2.2×，只作为反面教材引用。
可以复用的是我们在 `step2-mt`（R 包侧）实测验证过的**思想**——per-trait 缓存、同样本集短路、
拼接 GEMM、并行回落——以及那四层改动踩出来的坑。数值口径以 R 最新版为准。

---

## 0. 范围、术语、总的形状

**术语**

| 记号 | 含义 |
|---|---|
| `N` | 样本数（样本集相同时所有 trait 共享；不同时是并集大小，见 §4.7） |
| `P` | trait 数 |
| `p_t` | trait t 的协变量数（含截距），**允许各 trait 不同** |
| `Σp` | `Σ_t p_t` |
| `B` | 一个 marker 块里的列数（自适应，§6.3） |
| `q` | 本次运行要测的 marker 总数 |
| 「对」 | 一个 (marker, trait) 组合，是 step 2 输出的最小单元 |

**总的形状。** 主循环从「对每个 marker 做一次标量打分」变成：

```
读 B 个 marker → impute+QC（marker 级，全 trait 共享）
  → 建 N×B 基因型块 Gb
  → 6 个 GEMM（trait 维拼接）一次算出 B×P 个 normal-近似打分结果
  → 挑出需要 SPA / Firth / fastTest 重算的对，逐对回落到现有标量路径
  → 写 P 个输出文件
```

**为什么这样能省。** 唯一与 marker 相关、与 trait 无关的输入是 `Gb`。所有 per-trait 的量
（`X_t / A_t / XVX_t / res_t / mu2_t / S_a,t`）在一次运行内是常量。
不拼接就是 P 个独立 GEMM，每个都把 `Gb` 完整流读一遍：N=50k、B=256 时 `Gb` 是 100 MB（double），
P=64 → 6.4 GB 内存流量；拼成一个 GEMM 后 `Gb` 只读一遍，加上 `Astack` 的 Σp·N·8 = 256 MB，
合计 ~0.36 GB。这是 R 包侧 commit `f9374537` 单独做的一层，收益来自**带宽不是 FLOP**。

**本设计相对 R 包侧那版的三处实质改进**（都是因为 C++ 主干结构更好）：

1. 不需要重写 SPA。主干的 `spa.hpp` 已经是纯 arma 传参、`spa.cpp` 里 `grep static` 零命中、
   已经 R-free。R 包侧为了脱 `Rcpp::List` 花掉的那个 commit 这里**不需要**。
2. 不需要「全或无」门控。R 侧的 gate 是运行时整体开关；这里改成**静态 per-trait 划分**（§3.1），
   一个 trait 是否走批量在 marker 循环开始前就定死。
3. 不需要复刻 altFreq/altCounts 的翻转侧 quirk。样本集相同时（§4.3）
   flip / altFreq / altCounts / MAC / MAF / missingRate 全是 marker 级共享量，
   直接沿用主干 `imputeGenoAndFlip` 翻回后的 ALT 侧值，一份输出一套 flip 约定。
   样本集不同时（§4.7）这些量逐表型、用单表型的同一组表达式重算，仍是一套 flip 约定。

---

## 1. 数据布局决策

### 1.1 结论

**混合 SoA：**

- **per-trait 的"小"量（标量、p 级向量、p×p 矩阵）→ AoS**，装在 `std::vector<TraitMeta>` /
  `std::vector<arma::mat>` 里，按 trait 下标访问。
- **per-trait 的"大"量（N 级向量、N×p 矩阵）→ SoA，横向拼接成宽矩阵**：
  `Xstack [N × Σp]`、`Astack [N × Σp]`、`WXstack [N × Σp_bin]`、`RES [N × P]`、`MU2bin [N × n_bin]`。
- **P 个 `SAIGEClass` 实例**（`std::vector<std::unique_ptr<SAIGE::SAIGEClass>>`），
  每个就是今天那个类，一字不改，用来跑回落路径。

三条理由，按重要性排：

**(a) GEMM 的内存布局要求决定了大矩阵必须拼接。** armadillo 是列主序；
`Astack.cols(off_t, off_t+p_t-1)` 是**连续**的 N×p_t 子块，可以直接喂 BLAS。
拼接之后 `Zall = Astackᵀ·Gb` 是一个 `(Σp)×N · N×B` 的 GEMM，`Gb` 只流读一遍（见 §0）。
AoS（P 个独立 `arma::mat`）做不到这一点。

**(b) P 个 `SAIGEClass` 实例是 P=1 零回归的唯一廉价保证。** 回落路径必须与今天的标量路径
**逐运算相同**（I4）。如果把 `SAIGEClass` 改成持有 P 份成员、所有方法加一个 `int trait` 参数，
那是一次覆盖 `saige_test.cpp` 全文的改写，任何一处漏改都会静默改变数值。
保持类本身不变、开 P 个实例，回落就是**字面意义上的同一段代码**，
P=1 时 `objs[0]` 就是今天的 `ptr_gSAIGEobj`（§5）。

**(c) 缓存局部性。** 批量段每个 trait 的收尾只碰 `Z_t [p_t×B]`、`GW_t [p_t×B]`、`GR.col(t) [B]`，
都是小而连续的块；`Gb` 在一个块内被 6 个 GEMM 复用，块大小按 L3 选（§6.3）。
per-trait 的标量（`tau0 / SPA_Cutoff / pCutoffforFirth / ...`）打包在一个 `TraitMeta` 里，
一次 cache line 拿全，避免 P 个 `SAIGEClass` 之间跳指针。

### 1.2 内部 trait 顺序：binary 优先

内部把 trait 重排成 `[binary..., quantitative...]`，保存 `internal → output` 的置换。理由：

- `WXstack` 的第 t 块是 `mu2_t ∘ X_t`，**只有 binary 需要**（quantitative 那一块就是 `X_t` 本身，
  可以直接用 `Xstack` 的对应列区间）。binary 优先 ⇒ `WXstack` 只要 `N × Σp_bin`，省一整块 N×Σp。
- `MU2bin [N × n_bin]` 和 `CCM [N × 2·n_bin]` 的列区间也因此连续。
- quantitative 的列区间 `[qOff, Σp)` 连续 ⇒ 第二个 GW GEMM 也是单次 BLAS 调用。

置换只影响内部索引，输出文件顺序按 config 里 `models:` 的书写顺序。

### 1.3 结构定义（照着写）

新文件 `saige_mt.hpp` / `saige_mt.cpp`。**不改 `saige_test.hpp` 里 `SAIGEClass` 的任何成员**
（唯一例外见 §9.1 的 `tl_X1` 越界修复）。

```cpp
namespace SAIGE {

enum class TraitKind { Binary, Quantitative, Survival };

// 每 trait 一份；marker 循环内只读。
struct TraitMeta {
    std::string name;          // 输出/日志用；默认取 modelFile 目录名
    std::string modelDir, vrFile, outFile;
    TraitKind   kind;
    int   p        = 0;        // 协变量数
    int   colOff   = 0;        // 在 Xstack / Astack 里的起始列
    int   binOff   = -1;       // 在 WXstack 里的起始列；非 binary 为 -1
    int   binIdx   = -1;       // 在 MU2bin / CCM 里的 trait 序号；非 binary 为 -1
    double tau0    = 1.0;      // = tauvec[0]
    double SPA_Cutoff = 2.0;
    bool   is_Firth_beta = false;
    double pCutoffforFirth = 0.0;
    bool   isFastTest = false;
    double pval_cutoff_for_fastTest = 0.0;
    bool   isnoadjCov = false;
    bool   flagSparseGRM = false;
    bool   isCondition = false;
    bool   isMoreOutput = false;
    bool   locoApplied = false;   // 该 trait 是否真的换成了 chr<N>/ 的文件
    int    nCase = 0, nCtrl = 0;
    bool   batchable = false;     // §3.1 的静态门控结果
};

// 全局唯一，构造后只读。P=1 时**不构造**（§5）。
struct MTContext {
    int N = 0, P = 0, nBin = 0;
    int sumP = 0, sumPbin = 0, sumPqnt = 0, qOff = 0;  // qOff = sumPbin

    arma::mat Xstack;    // N × sumP     内部序（binary 优先）
    arma::mat Astack;    // N × sumP     A_t = XVX_inv_XV_t
    arma::mat WXstack;   // N × sumPbin  第 t 块 = mu2_t % X_t
    arma::mat RES;       // N × P        第 t 列 = res_t
    arma::mat MU2bin;    // N × nBin
    arma::mat CCM;       // N × 2*nBin   [case_0, ctrl_0, case_1, ctrl_1, ...]，0/1 指示列

    std::vector<arma::mat> XVX;   // P 个 p_t × p_t
    std::vector<arma::vec> S_a;   // P 个 p_t
    std::vector<TraitMeta> meta;  // 内部序
    std::vector<int> outOrder;    // internal → config 书写序

    std::vector<int> batchTraits;      // batchable 的内部下标（binary 在前）
    std::vector<int> batchQuantTraits; // batchable 且 quantitative
    std::vector<int> scalarTraits;     // !batchable
};

// 每线程一份，跨 block 复用（grow-only，不是每块 new）。
struct MTScratch {
    arma::mat Gb;      // N × B
    arma::mat Gb2;     // N × B   = Gb % Gb
    arma::mat Zall;    // sumP × B
    arma::mat GWbin;   // sumPbin × B
    arma::mat GWqnt;   // sumPqnt × B
    arma::mat GR;      // B × P
    arma::mat G2Mu2;   // B × nBin
    arma::vec Gsq;     // B
    arma::mat AC;      // B × 2*nBin
    arma::mat VR;      // B × P     每对的 variance ratio
    std::vector<std::pair<int,int>> fb;  // 回落队列：(块内列号, 内部 trait 下标)
};

}  // namespace SAIGE
```

`Xstack`/`Astack`/`WXstack` 构建时**直接按 `[N × Σp]` 分配再原地填**，不要先建 per-trait
临时矩阵再拷进去（R 包侧遗留 9 就是这么多占了一倍内存）。填法：

```cpp
ctx.Xstack.set_size(N, sumP);
for (int t : internalOrder) {
    ctx.Xstack.cols(meta[t].colOff, meta[t].colOff + meta[t].p - 1) = nm[t].X;
    // nm[t] 是 NullModelData，填完即可释放它的 X/XVX_inv_XV
}
```

### 1.4 P=1 时的零开销

`MTContext` 在 `P == 1` 时**不构造**。`main()` 里：

```cpp
if (models.size() == 1) {
    // 完全走今天的路径：一个 SAIGEClass、一个 ptr_gSAIGEobj、mainMarkerInCPP 原样
} else {
    // MT 路径：mainMarkerMT()
}
```

分叉点在 `main()`（§4.4），不在 marker 循环里。P=1 不多分配一个字节、不多一次分支判断。
细节见 §5。

---

## 2. 批量核：从 `scoreTestFast_block` 到 B×P

### 2.1 现状

`saige_test.cpp:549 scoreTestFast_block(const arma::mat& G, ...)` 已经是分块 GEMM 打分，
`const`，无成员写入，签名里输出都是长 B 的 `arma::vec`。它显式物化了两个 N×B 中间矩阵：

```cpp
arma::mat Z  = m_XVX_inv_XV.t() * G;   // p×B   ← GEMM 1  (:588)
arma::mat BX = m_X * Z;                // N×B   ← GEMM 2  (:589)
arma::mat tildeG = G - BX;             // N×B            (:594)
```

**扩到 B×P 时这两块必须消失**（不变量 I8）：`BX` 和 `tildeG` 各是 N×B，
乘上 P 个 trait 就是 P 倍的 N×B 写带宽，直接吃掉拼接省下来的那一遍流读。

### 2.2 代数：把 B 和 g̃ 消掉

设 `Z_t = A_tᵀ g`、`B_t = X_t Z_t`、`g̃_t = g − B_t`（A_t = `m_XVX_inv_XV`，X_t = `m_X`）。

**打分统计量**（`scoreTestFast` 的 `S = (S1 + S2)/τ₀`，其中
`S1 = res₁ᵀg̃₁`、`S2 = −(S_a − X₁ᵀres₁)ᵀZ`）：

```
S1 = gᵀres − res₁ᵀX₁Z
S2 = −S_aᵀZ + res₁ᵀX₁Z
S  = (gᵀres_t − S_a,tᵀZ_t) / τ₀_t          ← res₁ᵀX₁Z 对消
```

**方差，binary / survival**（标量式 `ZᵀXVXZ − Σmu2·B² + Σmu2·g̃²`）：

```
Σ mu2 g̃² = Σ mu2 g² − 2 Σ mu2 g B + Σ mu2 B²
var2 = Z_tᵀXVX_t Z_t + Σᵢ mu2_{it} gᵢ² − 2 gᵀ(mu2_t ∘ X_t) Z_t     ← Σmu2·B² 完全对消
```

标量路径里两个 `Σ` 都只在**载体样本**（g≠0）上求和；对消在同一索引集上成立。
又因为 g 在非载体处为 0，两项都可以免费扩到全 N。

**方差，quantitative**（标量式 `ZᵀXVXZ·τ₀ + g₁ᵀg₁ − 2 g₁ᵀB`）：

```
var2 = Z_tᵀXVX_t Z_t · τ₀_t + gᵀg − 2 gᵀX_t Z_t
```

**其余**（与标量完全一致）：

```
var1 = var2 · vr(marker, t)
stat = S²/var1
pval = chisq₁ 上尾；pval==0 时走 log_chisq1_uppertail 的 "%.1fEd" 定点格式
Beta = S/var1;  seBeta = |Beta|/√|stat|;  StdStat = |S|/√var1
```

对消之后，**B 和 g̃ 都不需要落地**，只剩 `(Σp)×B` 和 `B×P` 的小结果。
唯一保留的 N×B 缓冲是 `Gb` 本身和 `Gb2 = Gb∘Gb`（`Gb2` 同时喂 binary 的 `G2Mu2` GEMM
和 quantitative 的 `Gsq`；没有 binary trait 时用手写循环算 `Gsq`，不建 `Gb2`）。

### 2.3 GEMM 组织：拼接，不是循环 trait

**选择：拼接。** 理由在 §0 已经算过——`Gb` 流读次数从 P 降到 1，
N=50k/B=256/P=64 时内存流量 6.4 GB → 0.36 GB。
副作用是需要 `Astack`/`Xstack`/`WXstack` 三块 N×Σp 的常驻内存，这是可以接受的代价：
N=50k、Σp=640（P=64、p=10）时每块 256 MB。真正大的场景（N=400k、p=20、P=64）
每块 4 GB，此时按 §9.2 的别名方案去重。

块内 6 个 GEMM / 规约：

| 量 | 表达式 | 形状 | 跨 trait |
|---|---|---|---|
| `Zall` | `Astack.t() * Gb` | (Σp)×B | **拼接** |
| `GWbin` | `WXstack.t() * Gb` | (Σp_bin)×B | **拼接**（仅 binary） |
| `GWqnt` | `Xstack.cols(qOff, sumP-1).t() * Gb` | (Σp_qnt)×B | **拼接**（仅 quant） |
| `GR` | `Gb.t() * RES` | B×P | 跨 |
| `G2Mu2` | `Gb2.t() * MU2bin` | B×n_bin | 跨 |
| `AC` | `Gb.t() * CCM` | B×2·n_bin | 跨（case/ctrl 等位计数） |
| `Gsq` | `colsum(Gb2)` | B | trait 无关 |

`GWbin`/`GWqnt` 刻意取成 `(Σp)×B` 而不是 `B×(Σp)`，好让它和 `Zall` 同形同步长，
per-trait 收缩写成 `sum(GW_t % Z_t, 0)` 时两个操作数布局一致。
`GR`/`G2Mu2` 取成 `B×P`，因为收尾是「固定 trait、沿 B 向量化」，`GR.col(t)` 要连续。

每 trait 再做 O(p_t·B) 的小收缩（**必须全部向量化成 `arma::rowvec` 级运算**，不变量 I9——
R 包侧实测逐对写标量时每对都在分配 arma 临时对象，曾是批量段的耗时大头）：

```cpp
const auto& M = ctx.meta[t];
arma::subview<double> Z_t = scr.Zall.rows(M.colOff, M.colOff + M.p - 1);   // p_t × B
arma::rowvec zxz = arma::sum(Z_t % (ctx.XVX[t] * Z_t), 0);                 // 1 × B
arma::rowvec saz = ctx.S_a[t].t() * Z_t;                                   // 1 × B
arma::rowvec gwz = (M.kind == TraitKind::Binary)
        ? arma::sum(scr.GWbin.rows(M.binOff, M.binOff + M.p - 1) % Z_t, 0)
        : arma::sum(scr.GWqnt.rows(M.colOff - ctx.qOff, M.colOff - ctx.qOff + M.p - 1) % Z_t, 0);

arma::vec S    = (scr.GR.col(t) - saz.t()) / M.tau0;
arma::vec var2 = (M.kind == TraitKind::Binary)
        ? (zxz.t() + scr.G2Mu2.col(M.binIdx) - 2.0 * gwz.t())
        : (zxz.t() * M.tau0 + scr.Gsq          - 2.0 * gwz.t());
arma::vec var1 = var2 % scr.VR.col(t);
```

`ctx.XVX[t] * Z_t` 是 p_t×p_t · p_t×B 的小 GEMM，P 次合计 `Σp_t²·B` flop，可忽略。

最后是每对一次的 p 值格式化循环（`sprintf("%.6E")`，不可向量化）。
**这段必须与 `scoreTestFast` 的格式化代码逐字相同**——包括 `var1 <= DBL_MIN → pval=1`、
`!isfinite(stat) → pval=1, stat=0`、`pval==0` 时的 `log_chisq1_uppertail` 分支和
`fraction >= 9.95` 的进位。建议把这 30 行从 `scoreTestFast` 里提成一个
`static inline void format_score_result(double S, double v1, double v2, ...)`，
标量路径和批量路径**共用同一份**，杜绝两边漂移。这是本设计里唯一允许动
`saige_test.cpp` 既有代码的地方，且是纯提取（提取后原路径必须逐字节回归通过）。

### 2.4 新签名

不改 `scoreTestFast_block`（P=1 时 `blockSize>1` 仍然走它，§5）。新增自由函数：

```cpp
// saige_mt.hpp
namespace SAIGE {

// 一个 marker 块的 normal 近似批量打分。
// Gb        : N × B，已 impute / QC / flip，列 j 对应 cols[j] 这个 marker
// traitSet  : 参与本次 GEMM 的内部 trait 下标（低 MAC 块只传 quantitative，见 §3.2）
// VR        : B × P，VR(j,t) = 该对的 variance ratio（traitSet 之外的列不读）
// 输出      : 全部 B × P，traitSet 之外的列不写
struct MTBlockResult {
    arma::mat Beta, seBeta, Tstat, var1, var2, StdStat, pvalNum;  // B × P
    std::vector<std::vector<std::string>> pvalStr;                // [P][B]
    std::vector<std::vector<char>>        pvalIsLog;              // [P][B]
};

void scoreTestBatchMT(const MTContext& ctx,
                      const std::vector<int>& traitSet,
                      const arma::mat& Gb,
                      const arma::mat& VR,
                      MTScratch& scr,
                      MTBlockResult& out);
}
```

自由函数而不是 `SAIGEClass` 成员：它不属于任何单个 trait，
放成员会诱使实现者去读 `m_*`，那正是要避免的（I10）。

---

## 3. 门控与回落判据

### 3.1 静态 per-trait 门控（marker 循环开始前定死）

R 包侧的不变量 I1 是「全或无」，本设计**有意收紧成更强的版本**：

> **批量 / 标量的划分必须是静态 per-trait 的，绝不允许 per-marker 地改变 trait 集合。**

理由：R 侧「全或无」是为了防止 per-marker 的混合判定；
静态 per-trait 划分同样排除了那个风险（划分在加载模型后、读第一个 marker 前就固定），
但避免了「一个 trait 配置特殊就把另外 63 个 trait 一起拖回标量路径」的浪费。

```cpp
// 加载完 P 个 null model 之后、进 marker 循环之前，对每个 trait 求一次。
bool isBatchable(const TraitMeta& M) {
    if (M.kind == TraitKind::Survival) return false;   // §10
    if (M.isCondition)                 return false;   // §10
    if (M.isnoadjCov)                  return false;   // 另一套公式，批量段没写
    // 第一趟的 flagSparseGRM_cur：isFastTest 为真时第一趟强制走稠密路径
    const bool sparseFirstPass = M.isFastTest ? false : M.flagSparseGRM;
    if (sparseFirstPass)               return false;   // 稀疏 GRM 要 per-marker 解 PCG
    return true;
}
```

与 R 侧 gate 的三处**放宽**（都是主干结构支持的）：

- **类别 varRatio 不再是障碍。** `VR(j,t)` 是 per-(marker,trait) 标量，
  Phase 2 里用 `computeVarianceRatio(MAC, sparse, noadj, has)`（`const` 纯函数，
  `saige_test.hpp:395`）按该 trait 自己的分档表查出来填进 `VR`，
  `var1 = var2 % VR.col(t)` 即可。R 侧要求单 VR 只是没做。
- **`isFastTest` 不再是障碍。** 第一趟 `flagSparseGRM_cur = false` 正好是批量核算的那条路；
  `needFastRecompute` 变成一个**回落触发条件**（§3.3），而不是整体开关。
  `main.cpp:1616-1645` 的 `fastRecomputeSameCtx` 预判在 MT 下照旧适用（per-trait 求一次）。
- **`isMoreOutput` 不再是障碍。** 它只影响 hom/het 计数，与打分无关；
  Phase 2/3 先让带 `isMoreOutput` 的 trait 的**计数部分**走 per-(marker,trait) 内联循环
  （主干 `main.cpp:1810-1858` 那段），打分仍然批量。
  指示矩阵 GEMM 的做法见 §9.3，属于 Phase 5。

### 3.2 per-marker 的唯一一处分流：ER / 低 MAC

`MAC ≤ g_MACCutoffforER && binary` 的对必须走 ER（`saige_test.cpp:1206-1257`），永远不进批量。
MAC 是 marker 级共享量，所以「哪些 trait 退出」在一个 marker 上只有两种情形。
**不要因此重建 stack**——改成**按列分组**：

```
块内 QC 通过的列按 MAC 切成两组：
  hi = { MAC >  g_MACCutoffforER }  → scoreTestBatchMT(ctx, ctx.batchTraits,      Gb_hi, ...)
  lo = { MAC <= g_MACCutoffforER }  → scoreTestBatchMT(ctx, ctx.batchQuantTraits, Gb_lo, ...)
                                       且 lo × binary 的对全部进回落队列（走 ER）
```

`batchTraits` 和 `batchQuantTraits` 都是启动时算好的常量列表，两次调用共用同一个 stack，
只是 `traitSet` 不同。这比 R 侧「低 MAC marker 整个退出批量、连 quantitative 也一起退」严格更好，
而且没有引入任何 per-marker 的动态结构。

> 注意：`hi`/`lo` 要**物理分成两个连续的 `Gb`**（或者在填 `Gb` 时先按 MAC 分区、
> 用 `Gb.cols(0, nhi-1)` / `Gb.cols(nhi, nB-1)` 两个 subview），别用掩码——
> BLAS 吃不了掩码。分区会打乱块内列序，用 `colToIdx[]` 映射回 marker 下标即可
> （主干 `main.cpp:1073` 已经有这个数组）。

### 3.3 回落判据（不变量 I3 / I4）

对批量算出来的每一对，判：

```cpp
const TraitMeta& M = ctx.meta[t];
const double stdStat = out.StdStat(j, t);
const double pnum    = out.pvalNum(j, t);      // std::stod(pvalStr) 的结果，与标量路径同源

bool needSPA = (!std::isnan(stdStat)) && (stdStat > M.SPA_Cutoff)
               && (M.kind != TraitKind::Quantitative);

bool needFirth = false;
if (M.kind == TraitKind::Binary && M.is_Firth_beta) {
    needFirth = out.pvalIsLog[t][j] ? (out.pvalNum(j,t) <= std::log(M.pCutoffforFirth))
                                    : (pnum <= M.pCutoffforFirth);
}

bool needFastRecompute = false;
if (M.isFastTest &&
    ((M.kind == TraitKind::Binary && MAC > g_MACCutoffforER) || M.kind != TraitKind::Binary)) {
    needFastRecompute = (pnum < M.pval_cutoff_for_fastTest) && !fastRecomputeSameCtx(M, MAC);
}

bool useBatch = !needSPA && !needFirth && !needFastRecompute;
```

这**就是主干 `main.cpp:1102-1147` 的 `pe.useBlock` 判据**，逐条平移到 per-trait。
不要重新发明，照抄并把 `ptr_gSAIGEobj->m_*` 换成 `ctx.meta[t].*`。

**为什么这个并集正好等于标量路径的触发集合**（R 侧的推断，在主干的代码上同样可验证）：
`getMarkerPval` 里 Firth 的门用的是 SPA 之后的 p；但只有 `StdStat > SPA_Cutoff` 才会做 SPA，
而那一支已经被 `needSPA` 收进来了；`StdStat ≤ SPA_Cutoff` 时 `t_pval = t_pval_noSPA`
（`saige_test.cpp:1205-1208`），即 Firth 的门判的就是 normal 近似那个数。既不漏也不多。

`stdStat` 为 NaN 时（`var1 ≤ DBL_MIN` 导致）两边都不做 SPA：
标量路径 `saige_test.cpp:1061` 是 `if(!std::isnan(StdStat) && (StdStat > m_SPA_Cutoff) && ...)`。
一致，不要额外特判。

### 3.4 回落怎么走

**直接调该 trait 的 `SAIGEClass` 实例的现有 ctx 重载**，不要转写：

```cpp
SAIGE::PerMarkerCtx ctx_first;
ctx_first.flagSparseGRM_cur = M.isFastTest ? false : M.flagSparseGRM;
ctx_first.isnoadjCov_cur    = M.isnoadjCov;
ctx_first.varRatioVal       = scr.VR(j, t);
ctx_first.erSeedStream      = (uint64_t)globalMarkerIdx + 1;   // 与 trait 无关，见下

g_firthDefer = (M.isFastTest && M.kind == TraitKind::Binary
                && MAC > g_MACCutoffforER && !sameCtx);        // 每对重设！

objs[t]->getMarkerPval(/* ... 与 main.cpp:1692 一字不差 ... */, ctx_first);
```

这就是 I4（「回落路径是 legacy 标量路径的逐运算转写」）的**最强形式**：
它不是转写，它是同一段代码。R 包侧不得不手写 `mtb_fb_reentrant` 是因为 R 的 `SAIGEClass`
靠 `assign_for_itrait` 写成员切换当前 trait；主干没有这个问题。

三条必须遵守的细节：

1. **`g_firthDefer` 是 `thread_local`**（`saige_test.cpp:14`）。同一个线程会连续处理
   多个 (marker, trait) 回落，**每对调用前都要重设**，否则串味。
2. **`erSeedStream` 用全局 marker 下标，不要掺 trait 下标。** ER 的重采样种子必须
   只由 marker 在输入中的位置决定（W1-4 的原意：结果不依赖调度）。
   同一个 marker 的不同 binary trait 用同一个流是可以的——每次 `SKATExactBin_Work`
   调用前都会 `SL_set_stream(id)` 重置（`saige_test.cpp:1214`），所以是确定性的。
   *测试时必须固定种子并验证这一点；真实运行不带额外随机源。*
3. **回落走的是「整条 `getMarkerPval`」，不是「只补 SPA」。** 它会重算打分，
   结果覆盖批量算出来的 Beta/Tstat/var1/pval。这是「显著位点逐运算等同」的来源，
   不要为了省那一次重算而只拿批量的 S/var 去喂 SPA。

---

## 4. Null model 加载

### 4.1 配置格式

`modelFile` / `varianceRatioFile` / `outputFile` 的标量形式**保持向后兼容**（P=1）。
多表型走新增的 `models` 序列：

```yaml
plinkFile: /opt/saige/data/mid
genoType: plink
nThreads: 8
LOCO: false
# chrom: "1"          # LOCO 时必填，全 trait 共用

models:
  - traitName: y1
    modelFile: /opt/saige/logs/step2/null_y1_arma
    varianceRatioFile: /opt/saige/logs/step2/null_y1.varianceRatio.txt
    outputFile: /tmp/out.y1.txt
  - traitName: y2
    modelFile: ...
    varianceRatioFile: ...
    outputFile: ...

# 可选，MT 专用
mtBlockSize: 0          # 0 = 自适应（§6.3）
mtMemBudgetGB: 8        # 块缓冲的总预算（所有线程合计）
mtRequireSameSamples: false  # 默认 false：允许不同样本集（§4.7）；true = 不同就报错
```

**校验：** `models` 与三个标量 key 互斥（同时出现 → 报错，别猜）。
`models` 长度 ≥ 1；长度 == 1 时走 P=1 路径（§5）。

YAML 的覆盖项（`is_Firth_beta` / `pCutoffforFirth` / `isnoadjCov` /
`cateVarRatioMinMACVecExclude` / `cateVarRatioMaxMACVecInclude`，`main.cpp:4113-4160`）
在 MT 下的规则：**顶层写 = 对所有 trait 生效；`models[i]` 里写 = 只对该 trait 生效，
且覆盖顶层**。两处都没写就用该 trait 的 `nullmodel.json`。

### 4.2 加载顺序

```
for i in 0..P-1:
    nm[i] = loadNullModel(models[i].modelFile, models[i].varianceRatioFile, useLOCO, chrom)
    apply YAML overrides (顶层 → models[i])
一致性校验（§4.3）
决定内部顺序（binary 优先），算 colOff / binOff / binIdx
建 MTContext 的 5 个 stack（原地填，填完释放 nm[i] 的大矩阵）
建 P 个 SAIGEClass（用 nm[i] 的原始矩阵，构造完 nm[i] 可整体释放）
对每个 trait 求 isBatchable()，划分 batchTraits / batchQuantTraits / scalarTraits
按 §4.5 打印一张表（哪个 trait 走批量、哪个不走、为什么）
```

> 顺序上有一个取舍：先建 stack 再建 `SAIGEClass`，两者都从 `nm[i]` 拷，
> 峰值内存 = stack + P 个实例 + 当前那一个 `nm[i]`。
> 如果 P 大到扛不住，见 §9.2 的别名方案（Phase 5）。

### 4.3 一致性校验

| 检查 | 不通过怎么办 |
|---|---|
| `sampleIDs` 完全相同（**内容和顺序都要**） | 不同 ⇒ 走 §4.7（默认）；`mtRequireSameSamples: true` 时硬报错 |
| `n` 相同 | sampleIDs 相同而 n 不同 ⇒ 硬报错；sampleIDs 不同时逐模型查 `n == sampleIDs 长度`（§4.7.1） |
| `traitType` | 允许混合 binary / quantitative；`survival` → 该 trait `batchable=false`（§10） |
| `p` | **允许不同**，`Σp` 拼接天然支持 |
| `flagSparseGRM` / `isFastTest` / `isnoadjCov` / `isCondition` | 允许不同，进 `isBatchable()` |
| `SPA_Cutoff` / `is_Firth_beta` / `pCutoffforFirth` / `pval_cutoff_for_fastTest` | 允许不同，都是 per-trait |
| `impute_method` | **必须相同**。它是全局 `g_impute_method` 的来源；样本集不同时技术上可逐表型，但未放开 → 硬报错 |
| `dimNum > 0`（稀疏 GRM） | 允许；该 trait 若 `isBatchable()` 为假就走标量。多个 trait 各带一份 `m_spSigmaMat` 会很占内存，加一条日志警告 |
| `loco_chroms` 不一致 | **警告，不报错**；见 §4.6 |

「sampleIDs 顺序恒等也要查」是 R 包侧踩过的坑（不变量 I2）：
只查集合相等会让 `GVec.elem(si)` 置换元素。这里直接要求逐位置字符串相等，一次 `==` 比较即可。

同样本集换来的红利（对应 R 侧那 10.3× 里的第 2 层）在这里是**结构性**的，不是一个短路开关：
`m_posSampleInPlink`（`genotype_reader.hpp:55`）只建一份，
`imputeGenoAndFlip` 只跑一次，`indexZero/indexNonZero`、`altFreq/altCounts/MAC/MAF/missingRate/flip`
全是 marker 级共享量，一份 `Gb` 直接喂所有 trait。
**不要**复刻 R 侧 `imputeGenoAndFlip_sub` 那套翻转侧语义（不变量 I7）——
统一用主干 `imputeGenoAndFlip`（`UTIL.cpp:58`）翻回后的 ALT 侧值，
Beta/Tstat 乘 `(1-2*flip)`，AF_case/AF_ctrl 在 flip 时取补。一份输出一套约定。

### 4.4 `main()` 的改动点

| 行 | 现状 | 改成 |
|---|---|---|
| `main.cpp:4081` | 单次 `loadNullModel` | 循环 P 次 |
| `main.cpp:4186` | 单次 `setSAIGEobjInCPP`（32 参） | 循环 P 次，结果进 `std::vector<std::unique_ptr<SAIGEClass>>`；`ptr_gSAIGEobj = objs[0].get()` 保留给 region 路径和 `openOutfile_single` 用 |
| `main.cpp:4235` | `setPLINKobjInCPP(..., nullModel.sampleIDs, ...)` | 传 reader 样本表：样本集相同时即 `nm[0].sampleIDs`，不同时是并集（§4.7.2） |
| `main.cpp:4113-4160` | 标量 YAML 覆盖 | 顶层 + per-model 两级（§4.1） |
| `main.cpp:4687` | 单次 `openOutfile_single` | P 次，每个 trait 一个 `std::ofstream`；见下 |
| `main.cpp:4699` | `mainMarkerInCPP(...)` | `P==1` 时不变；否则 `mainMarkerMT(...)` |
| `main.cpp:154-157` | 4 个全局 `std::ofstream` | 单变异那个换成 `std::vector<std::ofstream> g_OutFiles_single`；region 的三个不动 |

`openOutfile_single`（`main.cpp:555`）的表头依赖 `t_traitType` / `t_isImputation` /
`ptr_gSAIGEobj->m_isCondition` / `t_isMoreOutput`——全是 per-trait 的。
改成接受 `(std::ofstream&, const TraitMeta&, bool isImputation)`。
`writeOutfile_single`（`main.cpp:612`，37 个参数）同理，加一个流参数和 `TraitMeta`。

**输出策略：每个 trait 一个文件**，列和表头与今天单 trait 跑出来的**完全一致**。
这样回归门就是最朴素的 `cmp`（§8）。可选的合并长表（多一列 `TraitID`）留到 Phase 5，
不要在 Phase 1-4 里引入。

样本集相同时 QC 是 marker 级的、全 trait 共享，P 个文件的**行集合逐行相同**，
比对时不用对齐。样本集不同时 QC 逐表型（§4.7.3），行集合可以不同。

### 4.5 启动时打印的门控表

必须打印，否则出了数值问题没法定位：

```
===== Multi-trait: 64 models =====
  idx  name   type          p   batch  reason
    0  y1     binary       11   yes    -
    1  y2     binary       11   yes    -
   ...
   62  z1     quantitative  9   yes    -
   63  w1     binary       11   no     isnoadjCov=true
  batch traits: 63 (binary 48 / quantitative 15), scalar traits: 1
  Sigma p = 701, block size B = 256 (auto, 8 threads, 1.6 GB scratch)
```

### 4.6 LOCO

**一次运行 = 一条染色体**，这条在 MT 下不变：`chrom` 是全局的（marker 集合由它决定），
但每个 trait 从**自己的** `<modelDir>/chr<N>/` 读那 10 个 stem
（`mu res V offset XV XVX XVX_inv XVX_inv_XV XXVX_inv S_a`；`X` 和 `y` 永远读顶层，
`null_model_loader.cpp:552-565`）。

**好消息：LOCO 不影响批量结构。** 同一次运行内 chrom 固定 ⇒
每个 trait 的 `X_t / A_t / XVX_t / res_t / mu2_t / S_a,t` 在整个 marker 循环里仍然是常量，
§2 的数学一字不改，stack 照建。LOCO 只决定「这一族矩阵从哪个目录装载」。
**但这也意味着不变量 I5（per-trait 常量、无失效逻辑）依赖于「运行中途不换染色体」。**
如果将来要支持一次运行跑多条染色体，必须显式加 stack 失效与重建——
在 `MTContext` 里留一个 `int locoChrom;` 字段并在重建处 assert，不要默默依赖。

**混合状态必须显式处理。** `loadNullModel` 的第 3 道闸
（`null_model_loader.cpp:519-526`：`chrom` 不在该模型的 `loco_chroms` 里就**静默回落**
到全基因组拟合，只打一行日志）会让 P 个 trait 出现「有的用 chr 目录、有的用顶层」的混合。
这不是错误（对齐 R 的 `readInGLMM.R:107-113`），但必须：

- 把 `loco_applied` 记进 `TraitMeta.locoApplied`；
- 启动表里多一列显示；
- 如果 P 个 trait 的 `locoApplied` 不全相同，**打印一条醒目的 WARNING**
  （"N of P models fell back to the genome-wide fit on chrom X"），并列出是哪几个。

**测试数据缺口：** `tools/rda_to_arma.R` 只支持 non-LOCO。
多表型 LOCO 的端到端准确度测试**现在造不出数据**。
Phase 1 的先决条件之一是补它的 LOCO 支持（写出 `chr<N>/` 子目录的那 10 个 stem）；
在补上之前，LOCO 只做「加载路径」的单元检查（P 个模型各自读对了目录、混合状态告警正确），
不做数值对比。

### 4.7 不同样本集（每个模型拟合在自己的样本上）

**状态：已实现（PLINK），2026-09-14。** 真实 biobank 里每个表型缺失的人不同，
§4.3 原来的「sampleIDs 必须逐位相同」让多表型在真实数据上基本用不上；取交集不可接受
（`mid.indep16`：16 个表型各缺 5%，交集只剩 44.1%）。

**不变量：每个表型的输出 == 该模型单独跑（P=1）的输出。** 验收口径是逐字节 `cmp`，
与 §8 的其余门一致。

#### 4.7.1 配置与校验

- 不同样本集**默认允许**。`mtRequireSameSamples: true` 把「任一模型的 sampleIDs 与
  models[0] 不逐位相同」变回硬报错（给期望样本集一致、不一致就说明上游出错的流水线用）。
- `validateMTModels` 返回 `differ`。`differ` 时额外要求：每个模型 `sampleIDs` 非空、
  长度 == `n`、无重复 ID；`mainMarkerMT` 再要求 `m_n`（= `y.arma` 行数）== `sampleIDs` 长度。
  `impute_method` 仍须一致（与样本集无关，保持原检查）。
- `differ` 时**明确报错**的组合：
  - `genoType` 不是 `plink`（BGEN/VCF/PGEN 读取器在 reader 样本序上做浮点累加，
    逐表型复刻没做，见 4.7.6）；
  - 条件分析（`assign_conditionMarkers_factors` 按 reader 样本数读条件位点）；
  - region / group / LD 矩阵：P>1 本来就报错（§10），不变。
- 样本集全部相同时 `differ = false`，走的是原来那条路，一个分支都不多（4.7.5）。

#### 4.7.2 并集与下标

- 并集 = R `ReadModel_multiTrait` 的 `union_vector`：按 config 顺序把各模型的 `sampleIDs`
  拼起来、首次出现者保留（`mtUnionSampleIDs`）。PLINK reader 只按并集建一次。
- 每个 trait 一个 `MTTraitSamples`：`sameAsUnion`（sampleIDs 与并集**逐位**相同）、
  `pos[k]`（该 trait 第 k 个样本在并集里的下标，trait 自己的顺序）、
  `comp`（并集里不属于该 trait 的下标）。
- 所有 stack（`Xstack/Astack/WXstack/RES/MU2bin`）按并集长度建，trait 的第 k 行放在
  `pos[k]`，**其余行是精确的 0**。于是任何「stack 列 × 并集长度的基因型列」内积
  只会收集到该 trait 自己的样本。`sameAsUnion` 的 trait 直接整列赋值，与原来逐字节相同。

#### 4.7.3 逐表型的 marker 统计（全部复用单表型的表达式）

共享的只有「读盘 + 并集样本的 2-bit code」（`copyFusedCodes_ts`，每样本 1 字节）。
对每个 trait：

| 量 | 算法 | 为什么与单独跑逐位相同 |
|---|---|---|
| code 计数 `counts[4]` | 并集计数减去 `comp` 上的 code；trait 小于并集一半时直接在 `pos` 上数 | 整数 |
| 缺失前 altFreq / altCounts / missingRate / info | `PlinkClass::fusedPreStatsFromCounts(fs, n_t)`——从 Stage A 里**原样搬出**的函数，Stage A 自己也改为调它 | 同一函数、同一输入 |
| 前置 QC（maxMissRate/minMAF/minMAC/minINFO） | `mainMarkerInCPP` 的表达式，`n = m_n` | 同一表达式 |
| flip、填充值、dosage-zeroing 闸、code→dosage 表 `fd[4]`、缺失后 altFreq/altCounts | `finalizeFusedStats(fs_t, ...)` | 同一纯函数 |
| `SAIGE_STEP2_SCALAR_DECODE=1` 时的 altCounts | 物化 trait 自己的向量后 `arma::sum`（与 `imputeGenoAndFlip` 相同） | 同值同序 |
| 后置 QC、MAC、VR、ER 判定、`fastRecomputeSameCtx`、`g_firthDefer` | 用该 trait 的 MAC | 同一表达式 |
| trait 看到的基因型 `g_t[k]` | `fd_t[code[pos[k]]]` | 单表型的 Stage C 写的就是 `fd[code]` |
| `indexZero/indexNonZero` | 由物化的 `g_t` 按 `==0` 升序重建 | 两个生产者的定义 |
| AF_case / AF_ctrl / hom/het 计数 | 沿 `m_case_indices` 顺序累加 `fd_t[code[pos[case_idx[k]]]]` | 同值同序 |
| 输出的 AC/AF/MissingRate | 逐表型写（`MTTraitChunk` 里改为 per-trait 列） | — |

QC 逐表型 ⇒ **不同表型文件的行集合可以不同**（§4.4 里「P 个文件行集合逐行相同」只在
样本集相同时成立）。

#### 4.7.4 共享解码下的批量核：并集列 + 精确仿射修正

块里的基因型列 `g` 是**并集自己的**列：`fd_u[code]`，`fd_u` 是「单表型跑在恰好并集这些样本上」
会用的表（同一 `finalizeFusedStats`）。对 trait t，在它的样本上逐元素精确成立

```
g_t = a·g + b + d·[该格是缺失基因型]
  a = +1, b = 0   （t 与并集的 flip 相同）
  a = -1, b = 2   （flip 相反，2−G）
  d = fd_t[MISS] − (a·fd_u[MISS] + b)      （两边的填充值之差）
```

前提是三个非缺失 code 上 `fd_t[c] == a·fd_u[c] + b` 精确成立（0/1/2 是整数，默认
`dosage_zerod_cutoff=0.2` 下只有填充值会被清零，所以恒成立）。不成立的 (marker, trait)
（例如 `dosage_zerod_cutoff ≥ 1` 把 1 清零）**不进批量**，走逐对标量路径。

批量核的全部输入都是 g 的线性或二次式，于是逐项精确映射（`scoreTestBatchMT` 的 adj 分支）：

```
L(g_t)   = a·L(g) + b·L(1_t) + d·L(e_miss)            L ∈ {Aᵀ, (mu2∘X)ᵀ 或 Xᵀ, resᵀ}
Q_v(g_t) = Q_v(g) + 2ab·L_v(g) + b²·Σ_t v + q·Σ_miss v   v = mu2_t（binary）/ 1_t（quantitative）
q        = fd_t[MISS]² − (a·fd_u[MISS] + b)²
```

- `L(1_t)`、`Σ_t v`：per-trait 常量，建 context 时算一次（`sumA/sumW/sumR/sumM`）。
- `L(e_miss)`、`Σ_miss v`：对每个块列，把 stack 在该列缺失格那几行上求和（`MissA/MissW/MissR/MissMu2/MissMask`）。
  stack 在 trait 之外为 0，所以对并集的全部缺失格求和即得该 trait 自己的缺失格之和。
- `L_v(g)`：只有块里出现 flip 相反的对时才做的额外 GEMM（binary `Gbᵀ·MU2bin`，
  quantitative `Gbᵀ·MASKq`）。
- quantitative 的 `Σ_t g²` 不能再用 `colsum(Gb2)`：样本集不同的 quantitative trait 用
  `Gb2ᵀ·MASKq`（`MASKq` 是这些 trait 的 0/1 样本指示列）。
- 每一项只在不恒为 0 时才加（`flip / shift / miss` 三个开关），所以同块里有哪些别的对
  不会改变任何一对的结果（G4.1 仍成立，由 `mtBlockSize: 7` 用例验证）。
- `sameAsUnion` 的 trait 不读任何修正，算术与样本集相同时逐字节一致。

回落对（SPA / Firth / fastTest / ER / 不可批量的 trait / 仿射不成立）：用该 trait 的表从
code 物化 `g_t`（长度 `n_t`、trait 自己的样本序），重建下标，调**该 trait 的**
`SAIGEClass::getMarkerPval`——与单表型是同一段代码、同样的输入。

**块列的 hi/lo 分区**：列为 hi ⇔ 存在一个 QC 通过、可批量的 binary trait 其 MAC_t > MACCutoffforER。
一对 binary 取批量结果仍要求自己的 MAC_t > MACCutoffforER，故 hi 列覆盖所有这类对；
quantitative 在两类列上都打分。

#### 4.7.5 样本集相同的运行不受影响

`differ == false` 时：读路径、QC、`Gb` 的写入、批量核（`t_adj` 读都不读）、回落与原实现
逐行相同；唯一的变化是 AC/AF/MissingRate 从共享列拷进 per-trait 列再写出（值相同）。
`run_mt_correctness.sh` 与 `run_mt_config_tests.sh` 全过即是这一条的验收。

#### 4.7.6 每表型额外开销（相对样本集相同的运行）

记 `N_u` 并集大小，`n_t` trait 样本数，`m_j` 列 j 的缺失格数，`Σp` 各 trait 协变量数之和。

| 项 | 量级 | 何时发生 |
|---|---|---|
| code 拷贝 + 并集列填充 + 缺失格扫描 | O(N_u) / marker，与 P 无关 | 每个 marker |
| per-trait code 计数 | O(min(N_u − n_t, n_t)) / (marker, trait) | 每对 |
| per-trait 统计、QC、仿射参数、VR | O(1) / 对 | 每对 |
| 缺失格 stack 行求和 | O(m_j · (2Σp + P + n_bin + n_maskq)) / 列，所有 trait 共享 | 有缺失格的列 |
| 修正本身 | O(p_t) / 对 | 样本集不同的 trait |
| quantitative 的 `Gb2ᵀ·MASKq` | N_u·B·n_maskq flop / 块 | 有样本集不同的 quantitative trait |
| flip 修正的 GEMM | N_u·B·(n_bin 或 n_maskq) flop / 块 | 块内出现 flip 相反的对 |
| 回落对物化 `g_t` 与下标 | O(n_t) / 回落对（标量路径本身就是 O(n_t)） | 回落对 |
| AF_case/AF_ctrl | O(n_t) / binary 对，多一层下标间接（原来也是 O(n)） | binary 对 |
| `SAIGE_STEP2_SCALAR_DECODE=1` | O(n_t) / 对（为了 `arma::sum`） | 仅该回滚开关 |
| 内存 | 每线程 N_u·B 字节的 code；stack 按 N_u；`MASKq` N_u·n_maskq | — |

即：稳态下每个样本集不同的 trait 的额外开销与「并集中它**没有**的样本数」和
「该 marker 的缺失基因型数 × p」成正比，不与 n_t 成正比（binary 的 AF_case/AF_ctrl 除外，
那一项样本集相同时也是 O(n)）。本轮未做墙钟测量。

#### 4.7.7 没做的 / 已知限制

- **BGEN / VCF / PGEN**：报错。三个读取器的 altFreq/altCounts/info 是 reader 样本序上的浮点累加
  （BGEN 还是按文件样本序），`imputeGenoAndFlip` 的 altCount 是 `arma::sum`；
  逐表型逐位复刻需要保留原始剂量（BGEN 存的是 `2−dosage`，反算不精确）并按各自顺序累加，未做。
- 条件分析：报错（见 4.7.1）。
- `impute_method` 仍须各模型一致。
- 上游 R 多表型分支 `missingRate` 的问题（§9.6）在这里不存在：MissingRate 与 QC 都是逐表型的。

---

## 5. P=1 零回归保证

**硬门槛：P=1 时输出与今天的 `saige-step2` 逐字节相同。**

设计上用三层保证，从粗到细：

**(1) 路径分叉在 `main()`，不在循环里。**
`models.size() == 1`（或者用的是老的标量 key）时：
- 不构造 `MTContext`（0 字节额外分配）；
- `objs` 里只有一个实例，`ptr_gSAIGEobj = objs[0].get()`；
- 调的是**今天的 `mainMarkerInCPP`，函数体一行不改**。

**(2) `saige_test.cpp` 的既有函数不许改语义。**
允许的改动只有两处，都要单独 commit 并各自过一遍完整回归：
- §2.3 提取 `format_score_result`（纯提取，不改一个算子）；
- §9.1 `tl_X1`/`tl_A1` 的越界修复（P=1 时 p 恒定，修复前后行为完全相同）。

`scoreTestFast_block`（`saige_test.cpp:549`）**不动**。它的 `BX`/`tildeG` 物化在 P=1 时
只是内存开销不是正确性问题，而它是 `blockSize>1` 唯一的实现。
MT 的批量核是**新函数**（`scoreTestBatchMT`），不复用它。
代价是两份数学并存；收益是 P=1 的字节级不变性不依赖任何论证。
（`scoreTestFast_block` 与 `scoreTestBatchMT` 的一致性由 §8 的 Phase 2 门单独验，
两者算的是同一组数、求和序不同，差异应在 1e-12 量级。）

**(3) 新增的全局状态必须在 P=1 时是惰性的。**
`g_OutFiles_single` 在 P=1 时长度为 1，写出的字节序列与今天的 `OutFile_single` 相同
（`openOutfile_single` / `writeOutfile_single` 的重构只是把全局流换成引用参数）。
不许在 P=1 路径上引入任何新的 `if (mtCtx) ...` 判断到 marker 循环体内部。

**验收：** `tests/` 下现有的全部回归脚本 + §8 的 P=1 门，逐字节 `cmp`。
只要有一个字节不同就是设计失败，不许用「相对误差很小」搪塞。

---

## 6. 并行策略

### 6.1 现状

`main.cpp` + `saige_test.cpp` 共 24 处 pragma，真正的 `parallel for` 只有 4 个：

| 行 | 指令 | 维度 |
|---|---|---|
| `main.cpp:971` | `parallel for schedule(dynamic,16)` | 块内 marker（块预取 Phase 1 读盘，非 bgen） |
| `main.cpp:1000` | `parallel for schedule(static)` | 块内 marker（块预取 Phase 2 impute+QC） |
| **`main.cpp:1152`** | **`parallel for schedule(dynamic,64)`** | **marker（单变异主循环，主力）** |
| `main.cpp:4873` | `parallel for schedule(dynamic,1)` | region |

`main.cpp:3829-3831`：`omp_set_num_threads(g_nThreads)` + **`openblas_set_num_threads(1)`**
（防止 OMP 线程内再调多线程 BLAS 过订阅）。

### 6.2 选择：**block 维并行，BLAS 保持单线程**

```
#pragma omp parallel for schedule(dynamic, 1)
for (int blk = 0; blk < nBlocks; ++blk) {
    // 每个线程完整拥有一个 block：
    //   读 B 个 marker → impute+QC → 建 Gb（hi/lo 分区）
    //   → scoreTestBatchMT ×2 → per-trait 收尾 → 回落队列就地处理
    //   → 结果写进 chunk 级的 [B × P] 输出缓冲（无锁，各线程列区间不重叠）
}
```

**为什么不是 trait 维并行。** trait 维只有 P 路，而且拼接 GEMM 的全部收益就来自
「所有 trait 共用一次 `Gb` 流读」——按 trait 切开就是把这个收益扔掉。

**为什么不是 (marker, trait) 对维并行。** 对维是 B×P，粒度太细，
而且批量段本来就是一个 GEMM，没有可切的循环。

**为什么不是「单线程流水 + 多线程 BLAS」。** 解码、impute/QC、p 值格式化、回落
（SPA/Firth，全是标量代码）加起来不是小头，多线程 BLAS 对它们无能为力。
block 维并行让这四段一起吃满核。

**因此 `openblas_set_num_threads(1)` 保持不变。** 每个线程自己那个
`(Σp)×N · N×B` 的 GEMM 在 Σp=700、B=256 时已经是 ~180 GFLOP 级别的单次调用，
单线程 BLAS 的效率完全够（不是「小 GEMM 被调度开销吃掉」的场景）。
*这一条是设计假设，Phase 3 必须实测：跑一次 `nThreads=1, OPENBLAS_NUM_THREADS=T`
对 `nThreads=T, OPENBLAS_NUM_THREADS=1` 的墙钟对比，把数字记进 `PERF_NOTES_*`。*

### 6.3 块大小 B 与内存

每线程的块缓冲主要是 `Gb` + `Gb2`，各 N×B 个 double：

```
bytes_per_thread ≈ 2 · N · B · 8
B = clamp( floor( mtMemBudgetGB·2^30 / (2·N·8·nThreads) ) 向下取到 2 的幂, 32, 512 )
```

N=50k、8 线程、预算 8 GB → B = 512（实际用 1.6 GB @ B=256，2^幂 取 512 时 3.3 GB）。
N=400k、8 线程、预算 8 GB → B = 128。
下限 32 是 GEMM 效率下限；上限 512 是因为 B 再大对 `Gb` 的复用度已无增益。
`mtBlockSize` 显式给值时跳过自适应（便于复现实验）。

**`Gb`/`Gb2` 必须是 `MTScratch` 里的 grow-only 缓冲，跨 block 复用**，
不能每块 `new`——这正是主干 `main.cpp:1155-1180` 和 `saige_test.cpp:308-321`
用 `thread_local` 解决的那个 mmap/munmap 页错误问题。
MT 下**不要用 `thread_local`**（见 §9.1），改成
`std::vector<MTScratch> scratches(omp_get_max_threads())`，
在并行区里 `MTScratch& scr = scratches[omp_get_thread_num()];`。
理由：`thread_local` 无法按 trait/尺寸做正确的失效判断，而 `MTScratch` 的尺寸
由 `(N, B, Σp, P)` 唯一确定，可以在进并行区前一次性 resize 好，循环里零判断。

### 6.4 同步原语

- 回落在**拥有该 block 的线程里就地处理**，不进全局队列。
  回落对占比小（~1e-3 量级），`schedule(dynamic,1)` 足以吸收不均衡。
  好处：零共享状态，天然满足 I10。
  *如果 Phase 3 实测出显著的尾部不均衡，再改成「全局回落队列 + 第二个 `parallel for` 按对并行」，
  届时队列只在 block 循环结束后才 drain，仍然无锁写入。*
- Firth 计数 `mFirth` / `mFirthConverge` 变成 per-trait 数组，用 `#pragma omp atomic`
  更新（主干 `main.cpp:1770-1780` 的做法平移）。
- 输出写盘仍然在并行区之外（chunk 结束后串行写 P 个文件），不需要 `critical(outwrite)`。
- EOF 哨兵 `firstEndIdx` + `critical(endflag)` 照旧。
- VCF 的 `critical(genoread)` 照旧（htslib 非线程安全）。

### 6.5 `thread_local` 审计清单（必须逐个过）

| 位置 | 是什么 | MT 下的处置 |
|---|---|---|
| `saige_test.cpp:14` `g_firthDefer` | `thread_local bool` | **每 (marker,trait) 回落前重设**（§3.4）。不改类型。 |
| `saige_test.cpp:308-321` `tl_g1/res1/mu21/B/g1t` | `thread_local arma::vec`，按 `n_elem < nnz` 扩 | 安全（只和 nnz 有关） |
| `saige_test.cpp:308-321` `tl_X1/tl_A1` | `thread_local arma::mat`，按 **`n_rows < nnz`** 扩 | **有越界隐患，必须修**，见 §9.1 |
| `saige_test.cpp:1261` Firth 的 `x` | `thread_local arma::mat`，恒 N×2 | 安全 |
| `main.cpp:1159-1163` 的 6 个缓冲 | 尺寸只依赖 N | 安全；MT 路径用自己的 `MTScratch`，不共用 |
| genotype reader 的 per-thread `FILE*` | — | 安全（与 trait 无关） |

**规矩：想按 trait 缓存任何东西，`thread_local` 都不是正确的载体**——
同一个线程会交替处理不同 trait。per-trait 的东西放 `MTContext`，per-线程的放 `MTScratch`。

---

## 7. 输出与内存：marker 分块

### 7.1 现在的问题

`mainMarkerInCPP`（`main.cpp:770`）一次性为**全部 q 个 marker** 分配 21 个长度 q 的输出向量，
跑完再 `writeOutfile_single` 一次写出。`marker_chunksize` 现在只用来打进度
（`main.cpp:1183-1190`），没有真的分块。

P 个 trait 会把其中 16 个结果字段乘上 P。q=10⁷、P=64 时这是不可行的。

### 7.2 改法

**把 marker 循环按 `marker_chunksize` 真的分块**（这也正是 R step 2 的形态，
R 侧在 R 层循环 chunk 再调 `mainMarkerInCPP`）：

```
for (chunkStart = 0; chunkStart < q; chunkStart += g_marker_chunksize) {
    分配/复用 chunk 级缓冲：
      marker 级共享（×1）: chr pos ref alt marker info altFreq altCounts missingRate imputeInfo
      per-trait（×P）    : Beta seBeta Tstat varT pval pvalNA isSPAConverge
                           AF_case AF_ctrl N_case N_ctrl [N_case_hom ...] N
    block 维并行跑完这个 chunk
    串行写 P 个文件（append）
}
```

**关键点：marker 元数据只存一份。** chr/pos/ref/alt/marker 不 ×P。
altFreq/altCounts/missingRate/imputeInfo 在实现里是 per-trait 的 4 列（样本集不同时它们逐表型，
§4.7.3；样本集相同时从共享值拷入），chunk 缓冲 = 10×chunk + 20×chunk×P。

`writeOutfile_single` 末尾的 `numtest` 和 Firth 计数现在是**每次调用打印一行**
（`main.cpp:737-748`）。分块后要改成累加到 per-trait 的计数器，
在全部 chunk 跑完后统一打印 P 行。

**输出文件用 append 模式**（`openOutfile_single` 的 `isappend` 分支已经有了），
第一个 chunk 前写表头，之后只追加。行序 = marker 输入序，天然保持。

### 7.3 BGEN streamer 与分块

`BgenStreamer` 是按整个 `t_genoIndex` 建的（`main.cpp:851-866`），
`getMarker(i)` 按全局下标取。分块只改 OMP 循环的范围和输出缓冲，
**streamer 仍然跨 chunk 存活**，不要每 chunk 重建。
队列容量 64 对 block 维并行可能偏小，Phase 4 实测后调（做成 `bgenQueueCap` 配置项）。

---

## 8. 分阶段实施计划与回归门

每一阶段都必须 commit（不 push），commit message 写清「改了什么、为什么、验证了没有」。
上游源码（`saige_test.*` 的两处允许改动）与 MT 新增代码**分开 commit**。

### 前置：测试资产

- `tools/rda_to_arma.R` 补 LOCO 支持（写 `chr<N>/` 的 10 个 stem）。没有它测不了 §4.6。
- 造 P 个 null model：用 `/opt/saige/data/mid.{bed,bim,fam}`（N=50,000 × M=40,000），
  R 最新版跑 step 1 出 P 个表型的模型，再转 arma。
  表型集合要覆盖：binary（不同 case 比例，含极端不平衡以触发 SPA/Firth）、
  quantitative、以及一个 survival（验证它被正确排除）。
- **黄金基线**：对每个 trait 单独跑一次今天的 `saige-step2`（`nthreads: 1`），
  存成 `golden/<trait>.txt`。后面每一阶段都拿它比。
  同时对每个 trait 跑一次 R 最新版，确认 golden 与 R 一致（数值以 R 为准）。

### Phase 0 — 骨架 + P=1 直通

**产出**：`models:` 配置解析；`P==1` 时完全走今天的路径；
`openOutfile_single` / `writeOutfile_single` 改成接受流引用 + `TraitMeta`；
`saige_test.cpp` 的 `format_score_result` 提取（§2.3）。

**回归门（必须全过）**
- G0.1 `tests/` 下现有全部脚本通过。
- G0.2 P=1（用 `models:` 写法）的输出与 `golden/<trait>.txt` **逐字节相同**（`cmp`）。
- G0.3 P=1（用老的标量 key）的输出同样逐字节相同。
- G0.4 `blockSize: 8` + P=1 的输出与改动前的 `blockSize: 8` 输出逐字节相同
      （证明 `format_score_result` 提取没动到块路径）。

### Phase 1 — MT 加载 + per-pair 标量（无批量）

**产出**：P 个 `NullModelData` → 一致性校验 → P 个 `SAIGEClass` → `MTContext` 的
`meta` / `XVX` / `S_a`（stack 还不建）→ `mainMarkerMT` 的最朴素形态：
marker 维并行，内层 `for t in 0..P-1` 逐对调 `objs[t]->getMarkerPval(...)` →
P 个输出文件 + §7 的 chunk 分块。

**这是整个设计的正确性基线**，后面所有批量都拿它当对照物（它比 golden 更方便：
一次运行出 P 个文件）。

**回归门**
- G1.1 每个 trait 的输出与 `golden/<trait>.txt` **逐字节相同**（P=2/8/64 各跑一次）。
- G1.2 `nthreads: 1` 与 `nthreads: 8` 输出逐字节相同。
- G1.3 一致性校验的负面测试：样本集顺序不同 / `impute_method` 不同 / `models` 与标量 key
      同时出现 → 都必须**报错退出**，错误信息指明是哪个模型哪一项。
- G1.4 LOCO：P 个模型 + `chrom: "1"`，每个 trait 读对了自己的 `chr1/`；
      构造一个 `loco_chroms` 缺 chr1 的模型，确认混合状态 WARNING 打出来且该 trait 用顶层文件。
- G1.5 内存峰值 ≤ 预期（P=64、N=50k 时给出实测数字，记进 `PERF_NOTES_*`）。

### Phase 2 — 批量核（串行，先只求对）

**产出**：`MTContext` 的 5 个 stack；`scoreTestBatchMT`；hi/lo 列分区（§3.2）；
静态门控（§3.1）；回落判据（§3.3）+ 回落走 Phase 1 的 per-pair 路径（§3.4）。
**这一阶段 block 循环保持串行**（`nthreads: 1`），把数值问题和并发问题分开。

**回归门**
- G2.1 与 Phase 1 的输出比：**p 值字符串差异数 = 0**；
      Beta / seBeta / Tstat / var 的相对误差 ≤ 1e-10。
      > 验收口径刻意**不是**「逐字节」。批量对与逐对路径的求和序不同，
      > 差异在 ~1e-12 相对量级，靠 `%.6E`（7 位有效数字）吃掉。
      > R 包侧在 mid 数据 64 binary 表型上实测确实逐字节一致，但按 1e-12 估算，
      > 落在末位舍入边界的概率约 1e-5/值——更大规模上会出现零星末位差。
      > **出现差异时必须逐个核到确实是末位舍入边界，不许统计性地放过。**
- G2.2 回落对（SPA / Firth / fastTest 触发的）与 Phase 1 **逐字节相同**——
      这些走的是同一段代码，任何差异都是 bug。加一个诊断开关打印回落对的清单，
      确认两边触发集合完全一致（数量和 (marker,trait) 身份都要比）。
- G2.3 §11 的验收矩阵铺满。
- G2.4 `scoreTestBatchMT` 与 `scoreTestFast_block`（P=1 强行走 MT 路径的 debug 开关）
      对同一批数据的 Tstat/var1 相对误差 ≤ 1e-10。
- G2.5 P=1 的四条 G0 门仍然全过（防止改回归）。

### Phase 3 — block 维并行

**产出**：§6.2 的并行结构；`MTScratch` 数组；Firth 计数改 per-trait atomic；
`g_firthDefer` 的 per-pair 重设；§9.1 的 `tl_X1` 修复。

**回归门**
- G3.1 `nthreads: 1/2/8/16` 输出**两两逐字节相同**。
- G3.2 与 Phase 2 输出逐字节相同。
- G3.3 ER 的确定性：固定 `erSeedStream`，同一配置跑 3 次输出逐字节相同；
      改 `nthreads` 输出不变。
- G3.4 ThreadSanitizer 跑一遍小规模配置，零 data race 报告。
- G3.5 性能数字：P=64、mid 数据、8 线程，报「s/trait」并与 Phase 1 对比。
      同时记 §6.2 要求的 BLAS 线程切法对比。

### Phase 4 — 内存与块大小自适应

**产出**：§6.3 的 B 自适应；§7 的 chunk 缓冲收敛；`bgenQueueCap`；
启动表（§4.5）打印实际用的 B 和 scratch 大小。

**回归门**
- G4.1 输出与 Phase 3 逐字节相同（B 变化不改数值——**这一条要特别验**：
      B 只影响分块边界，不影响任何一对的算式；如果出现差异说明有跨列的状态泄漏）。
- G4.2 `mtBlockSize` 取 32 / 128 / 512 输出两两逐字节相同。
- G4.3 内存峰值 ≤ `mtMemBudgetGB` + 常驻（stack + P 个实例），实测记录。
- G4.4 大 q 场景（q ≥ 10⁶）下内存不随 q 增长。

### Phase 5 — 可选优化（各自独立，各自带门）

每一项都是「输出逐字节不变」的纯优化，单独 commit、单独验：

- 5a `Xstack`/`Astack` 与 `SAIGEClass::m_X`/`m_XVX_inv_XV` 的别名去重（§9.2）
- 5b `CCM` 减半：`AC_ctrl = colsum(Gb) − AC_case`（§9.3）
- 5c `isMoreOutput` 的 hom/het 指示矩阵 GEMM（§9.3）
- 5d 合并长表输出（多一列 `TraitID`）
- 5e `m_Xt` / `m_At` 仅在 `fusedMode != 0` 时构造（P=64 时省 2 块 N×Σp）

---

## 9. 已发现的坑与专门的技术点

### 9.1 `tl_X1` / `tl_A1` 的越界隐患（**MT 下会真的踩到，必须先修**）

`saige_test.cpp:312-317`：

```cpp
const arma::uword p = m_X.n_cols;
thread_local arma::mat tl_X1, tl_A1;
if (tl_X1.n_rows < nnz) { tl_X1.set_size(nnz, p); tl_A1.set_size(nnz, p); }
arma::mat X1(tl_X1.memptr(), nnz, p, false, false);   // ← alias
```

扩容条件只看 `n_rows`，**没有看 `p`**。单 trait 下 p 恒定，永远正确。
MT 下同一个线程会交替处理 p 不同的 trait：先用 p=5 的 trait 触发
`set_size(nnz=1000, 5)`（缓冲 5000 个 double），再处理 p=20 的 trait 且 nnz=900
（`n_rows=1000 >= 900`，不扩容），于是 alias 是 900×20 = 18000 个 double
盖在 5000 个 double 的缓冲上 → **堆越界写**。

修法（不改任何算式，P=1 时行为完全相同）：

```cpp
if (tl_X1.n_elem < nnz * p) { tl_X1.set_size(nnz, p); tl_A1.set_size(nnz, p); }
```

`tl_g1` 等 5 个向量用的是 `n_elem < nnz`，本来就是对的，不用动。

### 9.2 stack 与 `SAIGEClass` 成员的内存去重（Phase 5a）

`Astack.cols(off_t, off_t+p_t-1)` 与 `objs[t]->m_XVX_inv_XV` 存的是同一份数据，
`Xstack` 与 `m_X` 同理。P=64、N=400k、p=20 时这是 2 × 4 GB 的重复。

去重办法：`Xstack`/`Astack` 作为**唯一**存储，让 `SAIGEClass` 的两个成员变成
指向 stack 列区间的 armadillo 别名：

```cpp
// 列主序 ⇒ 列区间连续 ⇒ colptr 就是 N×p_t 子矩阵的起始地址
m_X          = arma::mat(ctx.Xstack.colptr(off_t), N, p_t, /*copy_aux_mem=*/false, /*strict=*/true);
m_XVX_inv_XV = arma::mat(ctx.Astack.colptr(off_t), N, p_t, false, true);
```

**风险**：`SAIGEClass` 里任何对这两个成员的赋值（`m_X = ...`）都会在 strict 别名上抛异常。
已核实（2026-09-13）：`grep -n "m_X *=\|m_XVX_inv_XV *=" saige_test.cpp` 只有
`:101 m_XVX_inv_XV = t_XVX_inv_XV;` 和 `:109 m_X = t_X;` 两处，都在构造函数里。
实现时重新 grep 一遍确认没有新增赋值。
并且 `m_Xt = m_X.t()` / `m_At = m_XVX_inv_XV.t()` 会各自再建一份 p×N 的副本——
这正是 5e 要处理的（只在 `fusedMode != 0` 时建）。
**正因为有这些副作用，这一项放 Phase 5，不要在 Phase 2 就做。**

### 9.3 case/ctrl 统计的两个优化

- **`CCM` 减半**：binary 的 y 只取 0/1，所以 case ∪ ctrl = 全样本，
  `AC_ctrl[j] = colsum(Gb)[j] − AC_case[j]`。`CCM` 从 `N×2n_bin` 降到 `N×n_bin`，
  那个 GEMM 的 flop 减半。
  **但求和序变了**（标量路径是沿 ctrl 索引显式累加，`main.cpp:1837-1846`）。
  AF 以 6 位有效数字打印、值是 O(0.1)、N 量级 10⁵ 时累加误差 ~1e-12 相对——
  不该可见，但**必须在验收矩阵里专门查 AF_ctrl 列**。所以放 Phase 5b，默认先用 2 列。
- **`isMoreOutput` 的 hom/het**：需要 `Ghom = 1[1.5≤g≤2]`、`Ghet = 1[0.5≤g<1.5]`
  两个 N×B 指示矩阵，各与 `CCM` 做一次 GEMM 得 `B×2n_bin` 的计数。
  计数是整数，GEMM 结果精确（double 能精确表示 ≤2⁵³ 的整数），**不存在求和序问题**。
  代价是两块额外的 N×B scratch。Phase 5c。

### 9.4 `log10p` 输出

主干 `saige_test.hpp:125` 有 `m_islog10p` 成员但**全仓 0 处使用**（死成员）。
现在只有 `pval == 0` 时才走 `log_chisq1_uppertail` 的 `"%.1fE%d"` 定点格式
（`saige_test.cpp:384-395` 与 `scoreTestFast_block` 里的同一段）。
R 包侧的 step2-mt 同样没实现 log10p。

**本设计不新增 log10p 功能**，但批量核必须**完整复刻 `pval == 0` 的那条分支**
（§2.3 的共享 `format_score_result` 保证了这一点）。
验收矩阵里必须有一格专门构造 p < 1e-300 的对，确认批量路径和标量路径打出同样的
`"%.1fE%d"` 字符串。注意这类对通常也满足 `StdStat > SPA_Cutoff` ⇒ 会回落，
所以要专门构造一个 quantitative 的超显著对（quantitative 永不回落）来测批量路径本身。

### 9.5 Firth 的 offset：不要向 R multitrait 对齐

R 包侧的 `mtb_fb_reentrant` 给 `fast_logistf_fit_simple` 传的是**空 offset 向量**，
因为 R 包那边 `m_offset` 从未赋值——那是刻意的 bug-compatible。
**C++ 主干不是这样**：`loadNullModel` 读 `offset.arma`，构造函数 `saige_test.cpp:152`
在 binary 分支里 `m_offset = t_offset`，`saige_test.cpp:1265` 传的是真实 offset。
回落路径直接调 `objs[t]->getMarkerPval`，自然就用对了。
**不要为了"和 R 包 multitrait 一致"去把它改成空向量。** 基准是 R 最新版的单 trait 路径。

### 9.6 `missingRate` 的上游 bug

R 侧子集 trait 的 `missingRate` 是错的（`Main.cpp:1310` 那行被注释掉，
打印和 QC 用的都是 marker 级值）。
样本集相同时 marker 级值就是正确值；样本集不同时（§4.7）MissingRate 与 QC 都逐表型、
用该 trait 自己的 code 计数算 ⇒ **这个 bug 在我们的实现里不存在**。

### 9.7 GPU SPA

R 包侧有一个 `SAIGE_MT_GPUSPA=1` 的可选层挂在回落段上，
`spa_gpu_solve` 返回 rc≠0 时按「SPA 未收敛」退化到 normal 近似——
**这条兜底会静默改变数值**。
本设计**不移植这个行为**。如果将来加 GPU SPA 后端：失败时要么真回退 CPU SPA
（数值不变），要么直接报错退出。不许保留改变数值的静默兜底。

---

## 10. 明确不做的

| 不做 | 为什么 | 留什么接口 |
|---|---|---|
| **region / gene-based 多表型** | `mainRegionInCPP`（`main.cpp:4873` 起）的每个 region 要建 `P1Mat/P2Mat`、跑 SKAT-O/Burden/ACAT，per-trait 的量进得更深；而且 region 的并行维已经是 region。这是最大的空白，但它是**另一个设计**，不是本设计的扩展。 | `MTContext` 与 region 代码解耦：region 路径继续用 `ptr_gSAIGEobj = objs[0]`，P>1 且配了 `groupFile` 时**直接报错**（"multi-trait region testing is not supported"），不要静默只跑第一个 trait。 |
| **survival / Cox** | score/var 与 SPA 都是另一套（`SPA_survival`，主干里也没完整移植）。 | `TraitKind::Survival` 已经在 enum 里；`isBatchable()` 返回 false ⇒ 该 trait 自动走 per-pair 标量路径（Phase 1 的路径），**能出正确结果，只是不加速**。 |
| **稀疏 GRM 的批量** | var 要 per-marker 解一次 PCG（`getPCG1ofSigmaAndGtilde`），不是固定矩阵的收缩。 | 同上：`isBatchable()` 返回 false，走标量。`isFastTest=true` 的稀疏 GRM trait 仍然能批量（第一趟是稠密路径），只是显著位点走回落重算。 |
| **conditional analysis 的批量** | 要 per-marker 的 `G1tilde_P_G2tilde` 与 `m_VarInvMat_cond` 修正，表达不成「固定 per-trait 矩阵 × 基因型块」。而且 11 个 `m_*_cond` 成员是 per-trait 的，`assignConditionFactors`（`saige_test.cpp:1650`）要对每个 trait 各跑一次。 | `isBatchable()` 返回 false 走标量；`assign_conditionMarkers_factors`（`main.cpp:400`）改成对每个 `isCondition` 的 trait 各调一次。**这一条 Phase 1 就要做对**，否则条件分析的 MT 运行会用错模型。 |
| **`isnoadjCov` 的批量** | 另一套公式（不做 X 投影），批量段没写。收益也小（它本来就是快路径）。 | `isBatchable()` 返回 false 走标量。 |
| ~~不同样本集的 trait（子集 trait）~~ | **已做，见 §4.7**（PLINK）。BGEN/VCF/PGEN + 不同样本集、条件分析 + 不同样本集仍明确报错。 | `mtRequireSameSamples: true` 恢复「不同即报错」。 |
| **GxE** | 每个 marker 现搭交互项设计矩阵，X 随 marker 变 ⇒ 「X 与 marker 无关」这个前提直接失效，per-trait 缓存和 stack 都不成立。主干目前也没有 GxE。 | 无。将来若加 GxE，它与 MT 批量互斥。 |
| **降精度** | 全程 double，数值口径对 R。 | 无。这一条不许放宽。 |

---

## 11. 验收矩阵（不变量 I11 的具体化）

每一格的判据都是「MT 批量路径的输出 == 该 trait 单独跑的 golden」，
`nthreads: 1`，逐字节（回落对）/ p 值字符串零差异 + 1e-10 相对误差（批量对）。

**维度**

| 维度 | 取值 |
|---|---|
| LOCO | off / on / **混合**（部分 trait 静默回落全基因组） |
| trait 类型 | 全 binary / 全 quantitative / 混合 / 含一个 survival（验证被排除且结果正确） |
| p | 全相同 / **各 trait 不同**（触发 §9.1 的 `tl_X1` 路径） |
| varRatio | 单 VR / 分档 VR（各 trait 分档不同） |
| sparse GRM | off / on（该 trait 走标量） / on+isFastTest（该 trait 走批量+回落） |
| `isnoadjCov` | off / on（该 trait 走标量） |
| Firth | off / on（不同 `pCutoffforFirth`） |
| SPA | 不触发 / 触发（构造极端不平衡的 binary 表型） |
| ER | 不触发 / 触发（构造 MAC ≤ 4 的 marker，验证同 marker 的 quantitative trait 仍走批量） |
| conditional | off / on（该 trait 走标量） |
| `isMoreOutput` | off / on |
| 边界 marker | 全零列 / 高缺失（接近 `maxMissRate`）/ MAF 恰在 cutoff / p < 1e-300（§9.4）/ flip 触发（altFreq > 0.5） |
| 基因型格式 | plink / bgen / pgen / vcf（至少 plink + bgen 全铺，其余抽样） |
| P | 1 / 2 / 8 / 64 |
| 线程 | 1 / 8（输出必须相同） |
| `mtBlockSize` | 32 / 128 / 512（输出必须相同） |

**必须单独列出来的几格**（最容易出问题）：

1. **混合 LOCO**：P 个 trait 里有的用 `chr<N>/`、有的静默回落顶层。WARNING 打出来了吗？
   数值对吗？（依赖 `rda_to_arma.R` 的 LOCO 支持。）
2. **各 trait p 不同 + 多线程**：直接命中 §9.1。修复前应该能跑出 ASAN 报告，修复后干净。
3. **ER marker 上的 quantitative trait**：验证列分区（§3.2）没把它一起拖下水，
   且它的数值与 golden 一致。
4. **同一个 marker 上 SPA 触发和不触发的 trait 并存**：验证回落队列只带走该带的那些对。
5. **`pval == 0` 的 quantitative 对**：验证批量路径的 `"%.1fE%d"` 格式化（§9.4）。
6. **`fastRecomputeSameCtx` 为真的 trait**：验证第二趟被正确跳过、Firth 在第一趟内联跑了
   （`g_firthFitCalls` 的计数应等于候选数，无重复执行）。

**随机成分**：ER 的重采样。测试时固定 `erSeedStream`（它本来就是 marker 下标 + 1，
确定性），并验证 `nthreads` 变化不改输出。真实运行不引入额外随机源。

---

## 12. 一句话总结给实现者

结构上要做的事只有一件：**把 block 维从 B 扩成 B×P，
把 `SAIGEClass` 的标量成员（`m_res` / `m_mu2` / `m_X` / `m_XVX_inv_XV` / `m_XVX` / `m_S_a` / `m_tauvec[0]`）
换成一族 per-trait 的量，大的横向拼接、小的按 trait 下标存。**

正确性上要守的事也只有一件：**批量只算 normal 近似，
任何需要 SPA / Firth / fastTest 重算的对都退回去调那个 trait 自己的
`SAIGEClass::getMarkerPval`——同一段代码，不是等价重写。**
