# 本地构建（saige-v100-dev，conda env RSAIGE_GPU）

上游此分支的构建假设 pixi/conda 打包环境，本地补两步：

1. plink2 静态库（Makevars 的 `-l:plink2_includes.a`）：
   ```
   git clone --depth 1 https://github.com/chrchang/plink-ng
   cd plink-ng/2.0/include && g++ -O2 -fPIC -std=c++14 -c *.cc
   ar rcs ../../../plink2_includes.a *.o
   ```
2. `src/Makevars` 的 `.pixi` include 路径已改为 `${CONDA_PREFIX}/include`
   （savvy/superlu/zstd 均来自 conda env）。

安装：`MAKEFLAGS=-j8 R CMD INSTALL --library=/opt/saige/Rlib-mt --no-docs .`
依赖补装：lintools（CRAN）。
