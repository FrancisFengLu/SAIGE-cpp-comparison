# SAIGE Step 2 -- C++ Standalone

Standalone C++ implementation of SAIGE-GENE Step 2 for gene/region-based and single-variant association testing. No R dependency required.

## Build Dependencies

### macOS (Homebrew)

```bash
brew install armadillo yaml-cpp boost superlu htslib zstd
```

### Ubuntu / Debian (apt)

```bash
sudo apt-get install -y \
  build-essential \
  libarmadillo-dev \
  libopenblas-dev \
  liblapack-dev \
  libyaml-cpp-dev \
  libboost-math-dev \
  libsuperlu-dev \
  libhts-dev \
  libzstd-dev \
  zlib1g-dev
```

A C++17 compiler is required (clang++ on macOS, g++ on Linux). The Makefile auto-detects the platform.

## Build

```bash
make -j4
```

This produces the `saige-step2` binary in the current directory.

## Usage

```bash
./saige-step2 config.yaml
```

## Input

- **Null model directory**: Step 1 output containing `nullmodel.json` and `.arma` binary files (mu, res, y, V, X, XVX, XVX_inv, XXVX_inv, XV, XVX_inv_XV, S_a).
- **Variance ratio file**: Tab-separated file from Step 1.
- **Genotype file**: one of the supported formats below.
- **Group file** (for region tests): 3-line-per-gene format (variants, annotations, weights).

## Supported Genotype Formats

| Format | Library | Config Key |
|--------|---------|------------|
| PLINK (.bed/.bim/.fam) | Built-in | `plinkFile` |
| VCF/BCF/VCF.GZ | htslib | `vcfFile`, `vcfField` |
| BGEN v1.2 | zstd + zlib | `bgenFile`, `bgenSampleFile` |
| PGEN (.pgen/.pvar/.psam) | Built-in (modes 0x01/0x02) | `pgenFile`, `pvarFile`, `psamFile` |
