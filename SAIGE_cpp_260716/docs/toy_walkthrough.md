# Toy example：一步一步自己跑

> 这一页是给开发机 `saige-v100-dev` 写的（路径都是这台机器上的）。在别的机器上请按 [Home](index.md) 的安装和教程走。

目标：在你自己的文件夹里，从模拟数据开始，跑一次 step 1 和 step 2（4 个 binary trait，full GRM，GPU）。
每一步都写出命令、它读什么、写出什么、日志里该看到什么。数字是 2026-10-08 在 `saige-v100-dev` 上按同样命令跑出来的参考值。
下面所有路径都相对于你的文件夹，所有命令都在这个文件夹里敲。

跑完后文件夹长这样：

```
你的文件夹/
├── saige-gpu-cpp                        程序（链接，见第 0 步）
├── input/
│   ├── geno.bed  geno.bim  geno.fam     基因型（PLINK 1），5,000 样本 × 5,000 marker
│   └── pheno.txt                        表型 + 协变量
├── step1.log                            step 1 日志
├── step1/
│   ├── step1.yaml                       这次 step 1 的全部设置（程序自动写的）
│   ├── models/b1/ … models/b4/          每个表型一个零模型目录（step 2 读它）
│   ├── vr_b1.varianceRatio.txt …        每个表型一个 variance ratio（step 2 读它）
│   └── vr_b1.30markers.SAIGE.results.txt …   估 variance ratio 用的 30 个 marker 的结果（不用管）
├── step2.log                            step 2 日志
└── step2/
    ├── step2.yaml                       这次 step 2 的全部设置（程序自动写的）
    └── b1.txt … b4.txt                  结果：每个表型一个文件，每个 marker 一行
```

---

## 0. 准备（一次）

```bash
mkdir -p 你的文件夹 && cd 你的文件夹
ln -s /opt/saige/worktrees/cli/SAIGE_cpp_260716/bin/saige-gpu-cpp .
./saige-gpu-cpp --version
```

`saige-gpu-cpp` 是 cli 分支 `make USE_CUDA=1 SM=70` 编出来的程序，它会调用同一个 `bin/` 目录里的
`saige-null`（step 1）和 `saige-step2`（step 2）。不需要 `conda activate`，也不需要设环境变量。

---

## 1. 造输入数据

### 1a. 基因型

```bash
mkdir -p input && cd input
/opt/saige/data/plink2 --dummy 5000 5000 0.01 acgt --seed 1 --make-bed --out geno0
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

回到你的文件夹（后面的命令都在这里敲）：

```bash
cd ..
```

---

## 2. Step 1：拟合零模型

### 2a. 运行

```bash
./saige-gpu-cpp step1 \
  --plinkFile input/geno \
  --phenoFile input/pheno.txt \
  --phenoCol b1,b2,b3,b4 \
  --covarColList x1,x2 \
  --traitType binary \
  --LOCO=FALSE \
  --nThreads 8 \
  --useGPU \
  --outDir step1 > step1.log 2>&1
```

| 参数 | 意思 |
|---|---|
| `--plinkFile input/geno` | 基因型 `.bed/.bim/.fam` 的前缀（不带扩展名） |
| `--phenoFile input/pheno.txt` | 表型文件；样本 ID 列默认叫 `IID`（`--sampleIDColinphenoFile` 可改） |
| `--phenoCol b1,b2,b3,b4` | 要拟合的表型，逗号隔开，一次跑完 |
| `--covarColList x1,x2` | 协变量 |
| `--traitType binary` | binary trait（quantitative trait 写 `quantitative`） |
| `--LOCO=FALSE` | 不做 LOCO（默认 TRUE，和 R 一样） |
| `--nThreads 8` | 线程数 |
| `--useGPU` | 用 GPU；没有 GPU 时自动在 CPU 上跑 |
| `--outDir step1` | 输出目录（不存在会自己建） |

参数名和 R SAIGE 的 `step1_fitNULLGLMM.R` 一样。全部参数和默认值：`./saige-gpu-cpp step1 --help`。
参考耗时约 2 秒。

程序先把所有设置写进 `step1/step1.yaml`，再在它上面运行 step 1 的引擎 `saige-null`。
所以 `step1.yaml` 就是这次运行的完整记录；文件开头几行注释写着原始命令，以及不经过 `saige-gpu-cpp`
直接重跑的命令（`saige-null -c step1/step1.yaml`）。同一个文件夹重跑 step 1 要加
`--IsOverwriteVarianceRatioFile=TRUE`，否则会停在已存在的 variance ratio 文件上（R 也是这样）。

### 2b. 检查

```bash
grep -E "GPU tier|^Converged|^Variance ratio" step1.log
ls step1/models/b1
cat step1/vr_b1.varianceRatio.txt
```

应该看到：

```
[parallelCrossProd] GPU tier=4 enabled (...)       ← 用上了 GPU；没有这行就是在 CPU 上跑的（不会报错）
Converged: true                                    ← 每个表型都应该收敛
Variance ratio: …/step1/vr_b1.varianceRatio.txt
```

`step1/models/b1/` 里是十几个 `.arma` 矩阵和 `nullmodel.json`，整个目录就是 step 2 要的「模型」。
`vr_b1.varianceRatio.txt` 参考值：

```
0.980001	null	1
0.979976	null_noXadj	1
```

（多线程跑 step 1，最后几位数字可能每次略有不同。）

---

## 3. Step 2：检验每个 marker

### 3a. 运行

```bash
./saige-gpu-cpp step2 \
  --step1Dir step1 \
  --plinkFile input/geno \
  --minMAF 0 \
  --minMAC 1 \
  --is_Firth_beta=TRUE \
  --pCutoffforFirth 0.01 \
  --nThreads 8 \
  --useGPU \
  --outDir step2 > step2.log 2>&1
```

| 参数 | 意思 |
|---|---|
| `--step1Dir step1` | step 1 的输出目录：里面每个表型都检验，模型和 variance ratio 自动找 |
| `--plinkFile input/geno` | 要检验的基因型（这里和 step 1 同一份） |
| `--minMAF 0 --minMAC 1` | marker 过滤 |
| `--is_Firth_beta=TRUE --pCutoffforFirth 0.01` | binary trait：p < 0.01 的做 Firth 校正 |
| `--nThreads 8 --useGPU` | 线程、GPU |
| `--outDir step2` | 输出目录（不存在会自己建）；每个表型一个 `step2/<表型>.txt` |

只想检验其中几个表型：加 `--phenoCol b1,b3`。参数名和 R SAIGE 的 `step2_SPAtests.R` 一样，
默认值也是 R 的（`--is_noadjCov=TRUE`、`--impute_method=best_guess`、`--is_fastTest=FALSE`、
`--SPAcutoff=2`、`--is_Firth_beta=FALSE`）。全部参数：`./saige-gpu-cpp step2 --help`。参考耗时约 2 秒。

和 step 1 一样，设置写在 `step2/step2.yaml`，引擎是 `saige-step2`。

### 3b. 检查

```bash
grep -E "useGPU|GPU coverage|device SPA|device Firth|device ER|markers were tested" step2.log
wc -l step2/*.txt
head -3 step2/b1.txt | column -t
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

参考（`step2/b1.txt` 里 p < 0.01 的第一行）：

```
CHR POS   MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 MissingRate BETA      SE        Tstat    var     p.value      p.value.NA   Is.SPA AF_case  AF_ctrl  N_case N_ctrl
1   39000 snp39    C       T       7872.58    0.787258   0.0096      -0.238621 0.0796287 -36.4386 147.629 2.729456E-03 2.708778E-03 true   0.749802 0.791503 509    4491
```

---

## 4. 常见问题

| 现象 | 原因 |
|---|---|
| `saige-gpu-cpp: error: unknown flag --phenoColl ... (did you mean --phenoCol?)` | 参数名拼错（区分大小写，和 R 一样） |
| `saige-gpu-cpp: error: --memoryChunk is an R SAIGE flag that ... does not support: ...` | R 有、这里没实现的参数；冒号后面是原因 |
| `saige-gpu-cpp: error: missing --outDir ...` | 缺必填参数 |
| `saige-gpu-cpp: error: --step1Dir step1: no step-1 models found ...` | `--step1Dir` 不是 step 1 的 `--outDir`，或者 step 1 失败了 |
| `Refusing to overwrite existing variance-ratio file` | 同一个 `--outDir` 重跑 step 1：加 `--IsOverwriteVarianceRatioFile=TRUE` |
| `saige-gpu-cpp: error: engine not found: ...` | `saige-gpu-cpp` 要和 `saige-null`、`saige-step2` 在同一个 `bin/` 里；用链接（`ln -s`），不要只拷贝它一个 |
| step 1 日志里没有 `GPU tier=` | 在 CPU 上跑了（GPU 不可用时不报错） |
