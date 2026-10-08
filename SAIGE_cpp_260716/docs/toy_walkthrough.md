# Toy example：一步一步自己跑

目标：在你自己的文件夹里，从模拟数据开始，跑一次 step 1 和 step 2（4 个 binary trait，full GRM，GPU）。
每一步都写出原始命令、它读什么、写出什么、日志里该看到什么。数字是 2026-10-08 在 `saige-v100-dev` 上按同样命令跑出来的参考值。

跑完后文件夹长这样：

```
$ROOT/
├── input/
│   ├── geno.bed  geno.bim  geno.fam     基因型（PLINK 1），5,000 样本 × 5,000 marker
│   └── pheno.txt                        表型 + 协变量
├── step1/
│   ├── step1.yaml                       step 1 配置
│   ├── step1.log                        step 1 日志
│   ├── models/b1/ … models/b4/          每个表型一个零模型目录（step 2 读它）
│   ├── vr_b1.varianceRatio.txt …        每个表型一个 variance ratio（step 2 读它）
│   └── vr_b1.30markers.SAIGE.results.txt …   估 variance ratio 用的 30 个 marker 的结果（不用管）
└── step2/
    ├── step2.yaml                       step 2 配置
    ├── step2.log                        step 2 日志
    └── out/b1.txt … out/b4.txt          结果：每个表型一个文件，每个 marker 一行
```

---

## 0. 环境（每开一个新终端都要做一次）

```bash
# 你的工作文件夹（自己定）
export ROOT=/path/to/your/folder

# 程序（这台机器上已经编好的，基于 main 610ec5e9）
export SAIGE_HOME=/opt/saige/worktrees/usage-doc/SAIGE_cpp_260716
export S1=$SAIGE_HOME/step1_saige-null/saige-null       # step 1
export S2=$SAIGE_HOME/step2_saige-step2/saige-step2     # step 2
export PLINK2=/opt/saige/data/plink2                    # 只用来造模拟数据

# 依赖库
source ~/miniforge3/etc/profile.d/conda.sh
conda activate saige-build
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
export R_HOME=$CONDA_PREFIX/lib/R                       # step 1 运行时要用

mkdir -p $ROOT/input $ROOT/step1 $ROOT/step2/out        # step 2 不会自己建 out/
```

---

## 1. 造输入数据

### 1a. 基因型

```bash
cd $ROOT/input
$PLINK2 --dummy 5000 5000 0.01 acgt --seed 1 --make-bed --out geno0
# plink2 把所有 marker 都放在 1 号染色体、ID 是 dummy 名字；改成 snp1..snp5000、后一半放到 2 号染色体
awk 'BEGIN{OFS="\t"} {if (NR > 2500) $1 = 2; $2 = "snp" NR; $4 = 1000 * NR; print}' geno0.bim > geno.bim
mv geno0.bed geno.bed; mv geno0.fam geno.fam; rm -f geno0.*
```

读：无。写：`geno.bed`（2-bit 打包的基因型）、`geno.bim`（每个 marker 一行）、`geno.fam`（每个样本一行）。

```
$ head -2 geno.bim          染色体  marker ID  遗传距离  位置  allele1  allele2
1	snp1	0	1000	C	A
1	snp2	0	2000	C	G
$ head -2 geno.fam          FID  IID  父  母  性别  表型(不用)
0	per0	0	0	2	1
0	per1	0	0	2	1
```

### 1b. 表型 + 协变量

```bash
awk 'BEGIN{srand(7); OFS="\t"; print "IID","b1","b2","b3","b4","q1","q2","x1","x2"}
     function gauss(){ return sqrt(-2*log(1-rand()))*cos(6.283185307*rand()) }
     { x1 = gauss(); x2 = (rand() < 0.5) ? 1 : 0
       b1 = (rand() < 0.10 + 0.03*x1) ? 1 : 0
       b2 = (rand() < 0.30) ? 1 : 0
       b3 = (rand() < 0.05) ? 1 : 0
       b4 = (rand() < 0.10) ? "NA" : ((rand() < 0.20) ? 1 : 0)
       q1 = 0.5*x1 + gauss(); q2 = 0.3*x2 + gauss()
       printf "%s\t%d\t%d\t%d\t%s\t%.4f\t%.4f\t%.4f\t%d\n", $2, b1, b2, b3, b4, q1, q2, x1, x2 }' geno.fam > pheno.txt
```

格式要求：制表符分隔；第一行是列名；`IID` 列和 `.fam` 第 2 列对应；binary 用 0 = 对照、1 = 病例；缺失写 `NA`。
b1–b4 是 binary trait（b4 有约 10% 缺失），q1 q2 是 quantitative trait（这个例子不用），x1 x2 是协变量。

```
$ head -3 pheno.txt
IID	b1	b2	b3	b4	q1	q2	x1	x2
per0	0	1	0	0	-1.0620	0.3673	-0.4378	1
per1	0	1	0	NA	0.1620	1.3208	0.4618	0
```

---

## 2. Step 1：拟合零模型

### 2a. 写配置

```bash
cd $ROOT/step1
cat > step1.yaml <<YAML
paths:
  plinkFile: $ROOT/input/geno           # .bed/.bim/.fam 的前缀（不带扩展名）
  out_prefix: $ROOT/step1/models        # 零模型写到 models/<表型名>/
  out_prefix_vr: $ROOT/step1/vr         # variance ratio 写到 vr_<表型名>.varianceRatio.txt
  overwrite_varratio: true              # 允许重跑覆盖
design:
  csv: $ROOT/input/pheno.txt            # 表型文件
  iid_col: IID                          # 样本 ID 列
  y_cols: [b1, b2, b3, b4]              # 要拟合的表型（一次跑完）
  covar_cols: [x1, x2]                  # 协变量
fit:
  trait: binary
  loco: false
  nthreads: 8
  use_gpu: true
  firth_beta: true                      # 存进模型：step 2 对 p < 0.01 的对做 Firth
  p_cutoff_for_firth: 0.01
  spa_cutoff: 2.0
YAML
```

`<<YAML` 不加引号，`$ROOT` 会在写文件时展开成真实路径。写完 `cat step1.yaml` 看一眼。

### 2b. 运行

```bash
$S1 -c step1.yaml > step1.log 2>&1
```

命令结构：`saige-null -c <配置文件>`。参考耗时约 3 秒。

### 2c. 检查

```bash
grep -E "GPU tier|^Converged|^Variance ratio" step1.log
ls models/b1
cat vr_b1.varianceRatio.txt
```

应该看到：

```
[parallelCrossProd] GPU tier=4 enabled (...)       ← 用上了 GPU；没有这行就是在 CPU 上跑的（不会报错）
Converged: true                                    ← 每个表型都应该收敛
Variance ratio: $ROOT/step1/vr_b1.varianceRatio.txt
```

`models/b1/` 里是十几个 `.arma` 矩阵和 `nullmodel.json`，整个目录就是 step 2 要的「模型」。
`vr_b1.varianceRatio.txt` 参考值：

```
0.980001	null	1
0.979976	null_noXadj	1
```

（多线程跑 step 1，最后几位数字可能每次略有不同。）

---

## 3. Step 2：检验每个 marker

### 3a. 写配置

```bash
cd $ROOT/step2
{
cat <<YAML
genoType: plink
plinkFile: $ROOT/input/geno            # 和 step 1 同一份基因型
AlleleOrder: alt-first                 # 输出的 Allele2 = .bim 第 5 列
minMAF: 0
minMAC: 1
maxMissRate: 0.15
LOCO: false
isFirth: true                          # 只控制日志里的 Firth 汇总行；做不做 Firth 由 step 1 的 firth_beta 决定
MACCutoffforER: 4                      # MAC <= 4 且 |z| > 2 的对做精确检验（ER）
nThreads: 8
useGPU: true                           # GPU 的各个子开关默认全开
outputFormat: text
models:
YAML
for t in b1 b2 b3 b4; do
cat <<YAML
  - traitName: $t
    modelFile: $ROOT/step1/models/$t
    varianceRatioFile: $ROOT/step1/vr_$t.varianceRatio.txt
    outputFile: $ROOT/step2/out/$t.txt
YAML
done
} > step2.yaml
```

结构：上半部分是全局设置；`models:` 下每个表型一项，给出 step 1 的模型目录、variance ratio 文件和这个表型的输出文件。
写完 `cat step2.yaml` 看一眼。

### 3b. 运行

```bash
$S2 step2.yaml > step2.log 2>&1
```

命令结构：`saige-step2 <配置文件>`（step 2 没有 `-c`）。参考耗时约 3 秒。

### 3c. 检查

```bash
grep -E "useGPU|GPU coverage|device SPA|device Firth|device ER|markers were tested" step2.log
wc -l out/*.txt
head -3 out/b1.txt | column -t
```

日志参考：

```
useGPU: Tesla V100-SXM2-16GB, 16144 MiB, sm_70; fp64, ... 9984 markers per device batch ...
[b1] 5000 markers were tested (4769 on the GPU + 231 with the device SPA, 0 via the scalar CPU path; 49 Firth fits on the device).
...
[b4] 5000 markers were tested (4761 on the GPU + 238 with the device SPA, 1 via the scalar CPU path; 49 Firth fits on the device).
GPU coverage: 19999 / 20000 pairs (99.995%)
device SPA: 910 pairs ...
device Firth: 188 pairs ...
device ER: 0 pairs ...
```

没用上 GPU 时会打印 `useGPU: refused, running on the CPU (<原因>)`，照样出结果。

每个结果文件 5,001 行（表头 + 5,000 个 marker）。列：

| 列 | 意思 |
|---|---|
| CHR POS MarkerID | 染色体、位置、marker ID（来自 .bim） |
| Allele1 Allele2 | 两个等位基因；效应是 Allele2 的 |
| AC_Allele2 AF_Allele2 | Allele2 的计数和频率 |
| MissingRate | 这个 marker 的缺失率 |
| BETA SE | 效应和标准误（p < 0.01 时是 Firth 校正后的） |
| Tstat var | score 统计量和它的方差 |
| p.value | **最终 p 值** |
| p.value.NA | SPA 校正前的 p 值 |
| Is.SPA | 这一行有没有走 SPA |
| AF_case AF_ctrl | 病例 / 对照里的 Allele2 频率 |
| N_case N_ctrl | 病例 / 对照人数 |

参考（`out/b1.txt` 里 p < 0.01 的第一行）：

```
CHR POS   MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 MissingRate BETA      SE        Tstat    var     p.value      p.value.NA   Is.SPA AF_case  AF_ctrl  N_case N_ctrl
1   39000 snp39    C       T       7872.58    0.787258   0.0096      -0.238621 0.0796287 -36.4386 147.629 2.729456E-03 2.708778E-03 true   0.749802 0.791503 509    4491
```

---

## 4. 常见问题

| 现象 | 原因 |
|---|---|
| `plink2: command not found` | 没设 `PLINK2`，或者用了 `plink2` 而不是 `$PLINK2` |
| step 1 报找不到 R / libR.so | 没 `conda activate saige-build` 或没设 `R_HOME` |
| step 2 报写不了输出文件 | `$ROOT/step2/out` 没建 |
| step 2 报找不到模型 | `modelFile` 要指向 `models/<表型名>` 这个**目录** |
| step 1 日志里没有 `GPU tier=` | 在 CPU 上跑了（GPU 不可用时不报错） |
