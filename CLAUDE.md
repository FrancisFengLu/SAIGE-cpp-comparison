# SAIGE R-to-C++ Comparison Project

## Project Goal

Make the C++ standalone SAIGE implementation produce **identical** results to the R SAIGE package for GLMM (Generalized Linear Mixed Model) null model fitting.

## Directory Structure

```
Jan_30_comparison/
├── AI_INSTRUCTIONS.txt    # Detailed AI agent instructions
├── SESSION_NOTES.txt      # Session-by-session changelog
├── PLAN.txt               # Work plan with checkboxes
├── test_sparse_grm.R      # R test script
├── reference/             # Documentation (HTML call graphs, analysis reports)
├── code_copy/
│   ├── SAIGE_isolated/    # Isolated R package (has own pixi environment)
│   │   ├── R/             # R source files
│   │   ├── src/           # Rcpp C++ source files
│   │   ├── extdata/input/ # Small test data (1000 samples, 128k markers)
│   │   └── pixi.toml      # R environment config
│   └── cpp_standalone/    # C++ standalone code
│       ├── glmm.cpp              # GLMM solvers (main fitting loop)
│       ├── saige_ai.cpp          # getCoefficients_cpp, getAIScore_cpp
│       ├── SAIGE_step1_fast.cpp  # PCG solver, GRM operations
│       ├── main.cpp              # Entry point
│       ├── config_test.yaml      # Test configuration
│       └── Makefile              # Build config
└── output/
    ├── checkpoints/       # R_CP*.rds and CPP_CP*.csv comparison files
    └── bypass/            # R-to-C++ value transfer files
```

## How to Build and Run

### Run Order: R first, then C++

**R must run first** because it generates bypass files (random vectors, QR transform) that C++ reads.

```bash
# Step 1: Run R (generates bypass files in output/bypass/)
cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison
~/.pixi/bin/pixi run --manifest-path=code_copy/SAIGE_isolated/pixi.toml Rscript test_sparse_grm.R 2>&1 | tee /tmp/r_output.txt

# Step 2: Build and run C++ (reads bypass files)
cd code_copy/cpp_standalone
make clean && make
./saige-null -c config_test.yaml 2>&1 | tee /tmp/cpp_output.txt
```

### Rebuild R Package
```bash
cd code_copy/SAIGE_isolated
~/.pixi/bin/pixi run --manifest-path=pixi.toml R CMD INSTALL --preclean .
```

## Bypasses

### Active: Random Vectors
- **File:** `output/bypass/random_vectors_seed10.csv` (30 vectors, N x 30 matrix)
- **Reason:** R and C++ RNG algorithms are fundamentally incompatible. Different random vectors change the Monte Carlo trace estimation (`tr(P*K)`), which directly affects Score and tau at every iteration.
- **Impact if disabled:** ~6.4% tau difference (0.2901 vs 0.3087). This is expected -- even different seeds within R would produce similar variation, since trace estimation is stochastic with only 30 vectors.
- R generates this file when run first; C++ reads it at startup in `GetTrace()`.

### Active: QR Transform
- **Files:** `output/bypass/R_qr_X1_transformed.csv`, `R_qr_qrr.csv`
- **Reason:** R and C++ QR decompositions differ slightly (sign conventions, column ordering).
- **Toggle:** Local `bool use_r_qr_bypass = true;` in `null_model_engine.cpp` (~line 361). Set to `false` to use C++ Eigen QR instead.
- **Impact when disabled:** ~2e-5 tau difference (negligible).
- R generates these files when run first; C++ loads them if `use_r_qr_bypass = true` and files exist.

## Critical Configuration

C++ now supports both solver paths matching R's `getPCG1ofSigmaAndVector`:
- `use_pcg_with_sparse_grm: false` (default) → direct sparse solve via `gen_spsolve_v4()` (matches R default `usePCGwithSparseGRM=FALSE`)
- `use_pcg_with_sparse_grm: true` → PCG iterative solver

Both R and C++ now default to the direct sparse solve. The R test script uses `usePCGwithSparseGRM = FALSE`.

## Current Status: VALUES MATCH (random vector bypass only)

| Metric | R | C++ (rand bypass) | C++ (both bypass) | C++ (no bypass) |
|--------|---|-------------------|-------------------|-----------------|
| tau[1] | 0.2901 | 0.2902 | 0.2901 | 0.3087 |
| Iterations | 4 | 4 | 4 | 4 |
| Diff vs R | — | ~2e-5 | ~4e-6 | ~6.4% |

Default configuration: random vector bypass ON, QR bypass OFF. The ~2e-5 difference is within float32 precision -- both implementations use `arma::fvec`/`arma::fmat` (32-bit) for core computations.

The ~6.4% difference without random vector bypass is expected and acceptable: Monte Carlo trace estimation with 30 vectors is inherently stochastic. Different seeds in R alone would produce similar variation.

## Bugs Found and Fixed (Feb 2, 2025)

### Bug 1: Inner IRLS Alpha Initialization (`glmm.cpp` ~line 455)
- **Problem:** C++ compared alpha against zeros; R compares against previous outer iteration's alpha
- **Fix:** Initialize `alpha_outer_prev` from `beta_init` instead of zeros

### Bug 2: Eta Overwrite (`glmm.cpp` line 715)
- **Problem:** Line `eta = (p > 0) ? (X * alpha + offset) : offset;` overwrote correct eta from getCoefficients
- **Why wrong:** In mixed models, `eta = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / W`, not simply `X * alpha`
- **Fix:** Removed the line

## Quantitative Trait Support (Feb 4, 2025)

Ported binary trait C++ code to support quantitative traits. See `code_copy/cpp_standalone/quantitative_conversion.txt` for full details.

**Status:** COMPLETE — Values match R.

| Metric | R | C++ (with bypass) | Difference |
|--------|---|-------------------|-----------|
| tau[0] | 0.210876 | 0.210876 | ~0 (float32 precision) |
| tau[1] | 0.474525 | 0.474526 | ~1e-6 |
| Iterations | 5 | 5 | Match |

**Key difference from binary:** Quantitative estimates 2 variance components (tau[0], tau[1] both free) vs binary (tau[0]=1 fixed). This means 2x2 AI matrix, 2 traces, 2 scores — but same algorithm. Initial tau is [1, 0] (not [0.5, 0.5]).

**Changes made:**
1. `GetTrace_q` — replaced Rcpp RNG with rademacher_vec + bypass loading (seed0)
2. `quant_glmm_solver` — full rewrite from binary_glmm_solver with 2D adaptations (inner IRLS, conservative first iteration, step halving, convergence after iter 1)
3. `getAIScore_q_cpp` — already structurally correct, added debug output
4. R's `GetTrace_q` (in SAIGE_isolated/src/) — added bypass file saving (random_vectors_seed0.csv)

### Bypass: Random Vectors (Quantitative)
- **File:** `output/bypass/random_vectors_seed0.csv` (30 vectors, N x 30 matrix)
- **Reason:** Same as binary — R and C++ RNG algorithms are incompatible
- R generates this file when run first; C++ reads it at startup in `GetTrace_q()`

**Test data:** `pheno_1000samples.txt_withdosages_withBothTraitTypes.txt` has `y_quantitative` column.

### hasCovariate QR Gating (Feb 4, 2025)

Added `hasCovariate` logic to `null_model_engine.cpp` to match R's QR transform gating. Binary trait with exactly 1 non-intercept covariate now correctly skips QR (matching R). Verified across 3 covariate configs: x1+x2, x1 only, none — all tau values match R.

## Rules for AI Agents

1. **Read AI_INSTRUCTIONS.txt and SESSION_NOTES.txt** for full context before making changes
2. **Use subagents** for code exploration and running code; main agent does reasoning only
3. **Write only to** `Jan_30_comparison/` and subdirectories; never modify `/SAIGE/` (original package)
4. **Do not change the algorithm** - fix bugs only, match R behavior exactly
5. **Acceptable differences:** <0.01% is floating-point precision; >1% is a bug
6. **Checkpoint strategy:** Both R and C++ output values independently to `output/checkpoints/`; compare outputs to verify match
7. **Bypass protocol:** If values cannot match independently, discuss with user first, then implement R-output-to-C++-read pattern and document in AI_INSTRUCTIONS.txt

## Key Functions

| Function | Purpose |
|----------|---------|
| `getCoefficients()` | Solve for fixed effects (alpha) using PCG |
| `getAIScore()` | Compute Score statistic and AI (Average Information) |
| `GetTrace()` | Monte Carlo estimation of tr(P*K) using 30 random vectors |
| `fitglmmaiRPCG()` | One AI-REML iteration (combines above) |

## Key Variables

- `tau[0]`, `tau[1]`: Variance components (tau[0]=1 fixed, tau[1] estimated)
- `alpha`: Fixed effect coefficients (in QR-transformed space during fitting)
- `eta`: Linear predictor (includes random effect contribution in mixed model)
- `Sigma_iY`, `Sigma_iX`: Sigma-inverse times Y and X
- `PY`: Projection of Y onto orthogonal complement of X
- `YPAPY`, `Trace`, `AI`: Components of AI-REML variance component update
- `qrr`: QR decomposition R matrix (for back-transforming alpha to original space)
