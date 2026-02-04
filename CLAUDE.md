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

### C++ Standalone
```bash
cd code_copy/cpp_standalone
make clean && make
./saige-null -c config_test.yaml 2>&1 | tee /tmp/cpp_output.txt
```

### R SAIGE
```bash
cd /Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison
~/.pixi/bin/pixi run --manifest-path=code_copy/SAIGE_isolated/pixi.toml Rscript test_sparse_grm.R 2>&1 | tee /tmp/r_output.txt
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

### Disabled: QR Transform
- **Files:** `output/bypass/R_qr_X1_transformed.csv.disabled`, `R_qr_qrr.csv.disabled`
- **Reason:** R and C++ QR decompositions differ slightly (sign conventions, column ordering), but the back-transformed coefficients and final tau converge to nearly identical values.
- **Impact when disabled:** ~2e-5 tau difference (negligible). C++ uses its own Eigen QR decomposition.
- **Status:** Disabled as of Feb 3, 2025. Can be re-enabled by removing `.disabled` suffix.

## Critical Configuration

The R test script **must** use `usePCGwithSparseGRM = TRUE` so R uses the PCG solver (matching C++). Without this, R uses a direct sparse solve that produces slightly different intermediate values.

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
