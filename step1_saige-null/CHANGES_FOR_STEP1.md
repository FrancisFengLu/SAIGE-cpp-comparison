# Changes to SAIGE C++ Step 1 (March 25, 2026)

These changes make the Step 1 output compatible with the C++ Step 2 input format, and make the Makefile portable across macOS and Linux.

## 3 Files Modified

### 1. Makefile

**What changed:** Added automatic architecture detection so it builds on both macOS (Darwin) and Linux without manual edits.

- Added `UNAME_S := $(shell uname -s)` at the top
- Compiler: `clang++` on Darwin, `g++` on Linux
- `INCLUDES`: platform-conditional. Darwin uses Homebrew paths (`/opt/homebrew/...`), Linux uses standard system paths (`/usr/local/include`, `/usr/include/eigen3`)
- `LIBS`: platform-conditional. Darwin uses explicit `-L` Homebrew paths for openblas, superlu, tbb. Linux uses just library names.
- Removed the hardcoded temp cellar Eigen path (`/opt/homebrew/var/homebrew/tmp/.cellar/eigen/5.0.1/include/eigen3`) which was machine-specific
- All existing functionality preserved: R/Rcpp detection, pkg-config for armadillo/yaml-cpp, same target name (`saige-null`), same source files

### 2. null_model_engine.cpp

**What changed:** `.arma` file output paths changed from **dot-prefix** to **directory-based** convention, matching what Step 2 expects.

Before (dot-prefix):
```
<out_prefix>.mu.arma
<out_prefix>.res.arma
<out_prefix>.nullmodel.json
```

After (directory-based):
```
<out_prefix>/mu.arma
<out_prefix>/res.arma
<out_prefix>/nullmodel.json
```

Specific changes:
- `default_model_path()` (line ~320): `prefix + ".nullmodel.json"` -> `prefix + "/nullmodel.json"`
- Lines ~715-742: `prefix` renamed to `model_dir`, `ensure_parent_dir()` replaced with `fs::create_directories(model_dir)`, all 11 `.save()` calls changed from `prefix + ".NAME.arma"` to `model_dir + "/NAME.arma"`

The 11 .arma files: mu, res, y, V, S_a, X, XV, XVX, XVX_inv, XXVX_inv, XVX_inv_XV

### 3. glmm.cpp

**What changed:** `obj_noK.json` output path also changed to directory-based, for consistency.

- Added `#include <filesystem>` at the top
- `export_score_null_json()` (line ~122): `paths.out_prefix + ".obj_noK.json"` -> `paths.out_prefix + "/obj_noK.json"`, plus `std::filesystem::create_directories(paths.out_prefix)` before writing

## Why

Step 2 (`null_model_loader.cpp`) reads model files from a directory:
```cpp
data.mu = loadArmaVec(model_dir + "/mu.arma");
data.res = loadArmaVec(model_dir + "/res.arma");
// etc.
```

But Step 1 was writing them with a dot-prefix:
```cpp
mu_d.save(prefix + ".mu.arma", arma::arma_binary);
```

This meant Step 2 could not find Step 1's output. Now `out_prefix` in the Step 1 config is treated as a directory that contains all model output files.

## Config Impact

The `out_prefix` field in Step 1's YAML config (e.g., `config_test.yaml`) is now treated as a **directory path**, not a file prefix. For example:
```yaml
paths:
  out_prefix: /path/to/output/my_model
```
Will create directory `/path/to/output/my_model/` containing `nullmodel.json`, `obj_noK.json`, and all `.arma` files.

Step 2's `modelFile` config key should point to this same directory:
```yaml
modelFile: /path/to/output/my_model
```

## Tested

- Build: compiles cleanly on macOS (Darwin/arm64)
- Run: Step 1 converges, writes all 13 files to the output directory
- Integration: Step 2 successfully loads the Step 1 output and runs 128,868 single-variant tests to completion
