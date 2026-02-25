// null_model_engine.cpp
// ------------------------------------------------------------
// NullModelEngine: fits the null model (binary/quant/survival),
// computes offsets, and optionally runs LOCO in a single place.
// Glue points are provided to hook your existing C++ kernels
// (PCG/GLMM/LOCO) without R.
//
// Depends: Eigen3, (optional) TBB if you parallelize LOCO.
// ------------------------------------------------------------

#include "null_model_engine.hpp"
#include "saige_null.hpp"
#include "score.hpp"

#include <Eigen/Dense>
#include <Eigen/QR>
#include <armadillo>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

namespace fs = std::filesystem;

namespace saige {

// ======= Utility: map raw vectors to Eigen without copying =======
// NOTE: design.X is stored in ROW-MAJOR order (d.X[i*p + j] = element (i,j))
// We must use Eigen::RowMajor to correctly interpret the data layout.
using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
static inline Eigen::Map<const RowMajorMatrixXd>
map_mat(const std::vector<double>& buf, int n, int p) {
  if (buf.size() != static_cast<size_t>(n) * static_cast<size_t>(p)) {
    throw std::invalid_argument("Design.X buffer size does not match n*p");
  }
  return Eigen::Map<const RowMajorMatrixXd>(buf.data(), n, p);
}

static inline Eigen::Map<const Eigen::VectorXd>
map_vec(const std::vector<double>& buf, int n) {
  if (buf.size() != static_cast<size_t>(n)) {
    throw std::invalid_argument("Vector buffer size does not match n");
  }
  return Eigen::Map<const Eigen::VectorXd>(buf.data(), n);
}

static inline Eigen::VectorXd zeros(int n) { return Eigen::VectorXd::Zero(n); }

// ======= Minimal logging helper =======
static inline void log(const std::string& s) {
  // Replace with spdlog or your logger of choice
  std::fprintf(stderr, "[NullModelEngine] %s\n", s.c_str());
}

// ======= CSV Matrix Loading for R Bypass =======
static const std::string BYPASS_DIR = "/Users/francis/Desktop/Zhou_lab/SAIGE_gene_pixi/Jan_30_comparison/output/bypass/";

static bool bypass_files_exist() {
  std::string x1_path = BYPASS_DIR + "R_qr_X1_transformed.csv";
  std::string qrr_path = BYPASS_DIR + "R_qr_qrr.csv";
  return fs::exists(x1_path) && fs::exists(qrr_path);
}

// Parse a single CSV value (handles quoted values)
static double parse_csv_value(const std::string& s) {
  std::string cleaned = s;
  // Remove surrounding quotes if present
  if (!cleaned.empty() && cleaned.front() == '"') cleaned.erase(0, 1);
  if (!cleaned.empty() && cleaned.back() == '"') cleaned.pop_back();
  return std::stod(cleaned);
}

// Load a CSV matrix (skips header row, handles quoted column names)
static Eigen::MatrixXd load_csv_matrix(const std::string& filepath, int& out_rows, int& out_cols) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open CSV file: " + filepath);
  }

  std::vector<std::vector<double>> data;
  std::string line;
  bool first_line = true;

  while (std::getline(file, line)) {
    if (line.empty()) continue;
    
    // Skip header row (contains column names in quotes)
    if (first_line) {
      first_line = false;
      continue;
    }

    std::vector<double> row;
    std::stringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      row.push_back(parse_csv_value(cell));
    }
    if (!row.empty()) {
      data.push_back(row);
    }
  }

  if (data.empty()) {
    throw std::runtime_error("CSV file is empty or has only header: " + filepath);
  }

  out_rows = static_cast<int>(data.size());
  out_cols = static_cast<int>(data[0].size());

  Eigen::MatrixXd mat(out_rows, out_cols);
  for (int i = 0; i < out_rows; ++i) {
    if (static_cast<int>(data[i].size()) != out_cols) {
      throw std::runtime_error("Inconsistent column count in CSV at row " + std::to_string(i));
    }
    for (int j = 0; j < out_cols; ++j) {
      mat(i, j) = data[i][j];
    }
  }
  return mat;
}

// Load R's QR-transformed X matrix from bypass file
static Eigen::MatrixXd load_r_qr_x1_transformed(int expected_n) {
  std::string filepath = BYPASS_DIR + "R_qr_X1_transformed.csv";
  int rows, cols;
  Eigen::MatrixXd X = load_csv_matrix(filepath, rows, cols);
  
  if (rows != expected_n) {
    throw std::runtime_error("R bypass X1 has " + std::to_string(rows) + 
                            " rows but expected " + std::to_string(expected_n));
  }
  
  std::cout << "[R BYPASS] Loaded R_qr_X1_transformed.csv: " << rows << " x " << cols << std::endl;
  return X;
}

// Load R's QRR matrix from bypass file
static Eigen::MatrixXd load_r_qrr() {
  std::string filepath = BYPASS_DIR + "R_qr_qrr.csv";
  int rows, cols;
  Eigen::MatrixXd qrr = load_csv_matrix(filepath, rows, cols);
  std::cout << "[R BYPASS] Loaded R_qr_qrr.csv: " << rows << " x " << cols << std::endl;
  return qrr;
}

// ======= Baseline GLM (IRLS for logistic; OLS for gaussian) =======
struct BaselineGLMOut {
  Eigen::VectorXd beta;   // p
  Eigen::VectorXd eta;    // n
  Eigen::VectorXd mu;     // n
};

static BaselineGLMOut glm_gaussian(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
  // OLS via normal equations with LDLT (robust to collinearity)
  Eigen::VectorXd beta = X.colPivHouseholderQr().solve(y);
  Eigen::VectorXd eta  = X * beta;
  return {beta, eta, eta};
}

static BaselineGLMOut glm_logistic(const Eigen::MatrixXd& X,
                                   const Eigen::VectorXd& y,
                                   int maxit, double tol) {
  const int n = static_cast<int>(X.rows());
  const int p = static_cast<int>(X.cols());
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd eta  = X * beta;
  Eigen::VectorXd mu   = (1.0 / (1.0 + (-eta.array()).exp())).matrix();

  double dev_prev = std::numeric_limits<double>::infinity();
  for (int it = 0; it < maxit; ++it) {
    // Weights and working response
    Eigen::VectorXd w  = (mu.array() * (1.0 - mu.array())).matrix(); // n
    // Guard against zeros in w
    for (int i = 0; i < n; ++i) { if (w[i] < 1e-12) w[i] = 1e-12; }
    Eigen::VectorXd z  = eta + (y - mu).cwiseQuotient(w);

    // Weighted least squares step
    // Solve (X' W X) beta = X' W z
    Eigen::MatrixXd Xw = X.array().colwise() * w.array(); // n x p
    Eigen::MatrixXd XtWX = X.transpose() * Xw;            // p x p
    Eigen::VectorXd XtWz = X.transpose() * (w.asDiagonal() * z);

    Eigen::VectorXd delta = XtWX.ldlt().solve(XtWz) - beta;
    beta += delta;
    eta   = X * beta;
    mu    = (1.0 / (1.0 + (-eta.array()).exp())).matrix();

    // Deviance for convergence check
    double dev = 0.0;
    for (int i = 0; i < n; ++i) {
      double yi = y[i], mui = std::min(std::max(mu[i], 1e-12), 1.0 - 1e-12);
      dev += -2.0 * (yi * std::log(mui) + (1.0 - yi) * std::log(1.0 - mui));
    }
    if (std::abs(dev - dev_prev) < tol) break;
    dev_prev = dev;
  }
  return {beta, eta, mu};
}

// ======= QR transform with back-transform support =======
struct QRMap {
  Eigen::MatrixXd R;  // rank x rank (upper tri)
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P;
  int rank{0};
  bool scaled_sqrt_n{false};
  bool valid{false};
};

// Build X_qr with optional √n scaling; extract Q from the SAME QR object.
static Eigen::MatrixXd qr_transform(const Eigen::MatrixXd& X,
                                               QRMap& map,
                                               bool scale_sqrt_n)
{
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
  const int rank = qr.rank();
  map.rank = rank;
  map.R = qr.matrixR().topLeftCorner(rank, rank).template triangularView<Eigen::Upper>();
  map.P = qr.colsPermutation();
  map.scaled_sqrt_n = scale_sqrt_n;
  map.valid = true;

  // Thin Q from the same factorization
  Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(X.rows(), rank);

  if (scale_sqrt_n) {
    Q *= std::sqrt(static_cast<double>(X.rows()));
  }
  return Q;  // n x rank
}

// Backtransform coefficients (α_qr are coefficients in the Q-basis you returned above)
static Eigen::VectorXd qr_backtransform_beta(const Eigen::VectorXd& alpha_in_basis,
                                             const QRMap& map,
                                             int p_orig, int n_rows)
{
  if (!map.valid) return alpha_in_basis;

  // If Q was scaled by √n, match that here: effective RHS to R^{-1} is √n·γ (γ = coeffs for Q√n)
  Eigen::VectorXd alpha = alpha_in_basis;
  if (map.scaled_sqrt_n) {
    alpha *= std::sqrt(static_cast<double>(n_rows));
  }

  Eigen::VectorXd x = map.R.topLeftCorner(map.rank, map.rank)
                        .triangularView<Eigen::Upper>()
                        .solve(alpha);                  // x = R^{-1} alpha

  Eigen::VectorXd alpha_full = Eigen::VectorXd::Zero(p_orig);
  alpha_full.head(map.rank) = x;

  // Undo permutation back to original column order
  return map.P * alpha_full;
}

// Optional: recover keep_cols (original indices) like your Rcpp function exposed
static std::vector<int> qr_keep_cols_original(const QRMap& map) {
  std::vector<int> keep;
  keep.reserve(map.rank);
  const auto& idx = map.P.indices();   // length = p
  for (int k = 0; k < map.rank; ++k) keep.push_back(static_cast<int>(idx[k]));
  return keep; // these are original column indices of the pivoted-first rank columns
}
// ======= External solver hooks (plug in your native kernels) =======
// You can implement these elsewhere and link, or wrap your existing functions.
// The NullModelEngine will throw if a required hook is missing.

using BinarySolverFn = FitNullResult (*)(const Paths&,
                                         const FitNullConfig&,
                                         const Design& /*design used for GLMM*/,
                                         const std::vector<double>& /*offset*/,
                                         const std::vector<double>& /*beta_init (optional)*/);

using QuantSolverFn = FitNullResult (*)(const Paths&,
                                        const FitNullConfig&,
                                        const Design&,
                                        const std::vector<double>& /*offset*/,
                                        const std::vector<double>& /*beta_init (optional)*/);

static BinarySolverFn g_binary_solver = nullptr;
static QuantSolverFn  g_quant_solver  = nullptr;

void register_binary_solver(BinarySolverFn fn) { g_binary_solver = fn; }
void register_quant_solver (QuantSolverFn  fn) { g_quant_solver  = fn; }

// ======= LOCO runner hook (optional batch implementation) =======
// Provide a single entry point that executes LOCO across chromosomes;
// If you already have a native LOCO function, wrap it and assign here.

using LocoBatchFn = void (*)(const Paths&,
                             const FitNullConfig&,
                             const LocoRanges&,
                             const Design&,
                             const std::vector<double>& /*theta*/,
                             const std::vector<double>& /*alpha*/,
                             const std::vector<double>& /*offset*/);

static LocoBatchFn g_loco_batch = nullptr;
void register_loco_batch(LocoBatchFn fn) { g_loco_batch = fn; }

// ======= NullModelEngine impl =======

NullModelEngine::NullModelEngine(const Paths& paths,
                                 const FitNullConfig& cfg,
                                 const LocoRanges& chr)
  : paths_(paths), cfg_(cfg), chr_(chr) {}

static void ensure_parent_dir(const std::string& out_path) {
  fs::path p(out_path);
  if (!p.parent_path().empty()) fs::create_directories(p.parent_path());
}

static std::string default_model_path(const std::string& prefix) {
  return prefix + ".nullmodel.json"; // JSON stub; replace with your serializer
}

FitNullResult NullModelEngine::run(const Design& design_in_const) {
  // Make a mutable copy so we can update X if QR transform is used
  Design design_in = design_in_const;

  if (design_in.n <= 0) throw std::invalid_argument("Design.n must be > 0");
  if (design_in.p < 0)  throw std::invalid_argument("Design.p must be >= 0");
  if (design_in.y.size() != static_cast<size_t>(design_in.n)) {
    throw std::invalid_argument("Design.y length must equal n");
  }

  // --- Map inputs ---
  Eigen::Map<const Eigen::VectorXd> y  = map_vec(design_in.y, design_in.n);

  Eigen::VectorXd offset_in = design_in.offset.empty()
    ? zeros(design_in.n)
    : map_vec(design_in.offset, design_in.n);

  // --- Copy X from row-major buffer to column-major Eigen matrix ---
  // design.X is stored row-major (d.X[i*p + j] = element (i,j))
  Eigen::MatrixXd X;
  Eigen::MatrixXd X_orig;  // Keep original X before QR transform (for offset computation)
  if (design_in.p > 0) {
    X = map_mat(design_in.X, design_in.n, design_in.p);  // converts to col-major
    X_orig = X;  // save copy before potential QR transformation
  } else {
    X = Eigen::MatrixXd(design_in.n, 0);  // empty matrix for no covariates
    X_orig = X;
  }
  QRMap qrmap;
  const int p_orig = design_in.p;  // Store original p before QR might change it
  // Store R's QRR matrix for potential back-transformation
  Eigen::MatrixXd r_qrr_matrix;
  bool using_r_bypass = false;

  // --- Determine hasCovariate (match R's logic in SAIGE_fitGLMM_fast.R lines 1485-1617) ---
  // design_in.p includes the intercept column, so non-intercept covariate count = p - 1
  int n_covariates = design_in.p - 1;
  bool hasCovariate;
  if (n_covariates <= 0) {
    hasCovariate = false;
  } else if (cfg_.trait == "binary" && n_covariates == 1) {
    // R special case: binary trait with exactly 1 covariate -> skip QR
    hasCovariate = false;
  } else {
    hasCovariate = true;
  }

  if (!hasCovariate) {
    std::cout << "[QR] hasCovariate=false (trait=" << cfg_.trait
              << ", n_covariates=" << n_covariates
              << ") — skipping QR transform (matches R behavior)\n";
  }

  if (cfg_.covariate_qr && hasCovariate) {

    // === QR bypass: load R's QR-transformed X instead of using Eigen QR ===
    // Toggle this to true/false to enable/disable the R QR bypass
    bool use_r_qr_bypass = true;

    if (use_r_qr_bypass && bypass_files_exist()) {
      // Tentatively load to check column count matches current design
      auto X_bypass = load_r_qr_x1_transformed(design_in.n);
      auto qrr_bypass = load_r_qrr();

      if (static_cast<int>(X_bypass.cols()) != p_orig) {
        // Stale bypass files from a different covariate config — skip
        std::cout << "[R BYPASS] Column mismatch: bypass has " << X_bypass.cols()
                  << " cols but design has " << p_orig << " — falling back to C++ QR\n";
        use_r_qr_bypass = false;
      }

      if (use_r_qr_bypass) {
        std::cout << "\n========================================" << std::endl;
        std::cout << "[R BYPASS] Using R's QR transformation values!" << std::endl;
        std::cout << "========================================\n" << std::endl;

        X = X_bypass;
        r_qrr_matrix = qrr_bypass;

        // Set up minimal QRMap for bypass mode
        qrmap.rank = static_cast<int>(X.cols());
        qrmap.R = r_qrr_matrix;  // Store R's qrr for back-transform
        qrmap.valid = true;
        // R's transformed X is Q*sqrt(N), so back-transform needs sqrt(N) multiplication
        qrmap.scaled_sqrt_n = true;
        // R doesn't use column pivoting, so P is identity
        qrmap.P.setIdentity(p_orig);
        using_r_bypass = true;

        std::cout << "[R BYPASS] X dimensions: " << X.rows() << " x " << X.cols() << std::endl;
        std::cout << "[R BYPASS] QRR dimensions: " << r_qrr_matrix.rows() << " x " << r_qrr_matrix.cols() << std::endl;
      }
    }

    if (!use_r_qr_bypass) {
      // Use C++ Eigen's QR decomposition
      std::cout << "[C++ QR] R bypass files not found, using Eigen QR decomposition" << std::endl;
      X = qr_transform(X, qrmap, design_in.n); // X becomes orthonormal columns (rank cols)
    }


    // --- Debug print: dimensions + first few entries ---
    std::cout << std::setprecision(10);
    std::cout << "[QR fingerprint] X(0,0)=" << X(0,0) << " X(0,1)=" << X(0,1);
    if (X.cols() > 2) std::cout << " X(0,2)=" << X(0,2);
    std::cout << " X(1,0)=" << X(1,0) << " X(1,1)=" << X(1,1);
    if (X.cols() > 2) std::cout << " X(1,2)=" << X(1,2);
    std::cout << std::endl;
    std::cout << std::setprecision(6);  // restore default
    std::cout << "[design] after covariate QR transform: "
              << X.rows() << " x " << X.cols() << std::endl;

    int max_rows = std::min<int>(5, X.rows());
    int max_cols = std::min<int>(8, X.cols());

    // column headers
    std::cout << "      ";
    for (int j = 0; j < max_cols; ++j)
        std::cout << "col" << j << "\t";
    if (X.cols() > max_cols) std::cout << "...";
    std::cout << "\n";

    // preview first few rows
    for (int i = 0; i < max_rows; ++i) {
        std::cout << "row" << i << " ";
        for (int j = 0; j < max_cols; ++j)
            std::cout << std::fixed << std::setprecision(4) << X(i, j) << "\t";
        if (X.cols() > max_cols) std::cout << "...";
        std::cout << "\n";
    }
    if (X.rows() > max_rows)
        std::cout << "... (" << X.rows() << " rows total)\n";
      //

    // === UPDATE design_in to use QR-transformed X in GLMM loop (like R does) ===
    int new_p = static_cast<int>(X.cols());  // may be reduced if rank-deficient
    design_in.X.resize(design_in.n * new_p);
    // Copy X (col-major Eigen) to design_in.X (row-major)
    for (int i = 0; i < design_in.n; ++i) {
      for (int j = 0; j < new_p; ++j) {
        design_in.X[i * new_p + j] = X(i, j);
      }
    }
    design_in.p = new_p;
    std::cout << "[QR] Updated design_in.X to use transformed X ("
              << design_in.n << " x " << design_in.p << ")\n";
  }

  // --- Baseline GLM (on transformed or original X) ---
  BaselineGLMOut glm;
  const bool is_binary    = (cfg_.trait == "binary" || cfg_.trait == "survival");
  const bool is_quant     = (cfg_.trait == "quantitative");
  const int  maxit_glm    = std::max(5, cfg_.maxiter);
  const double tol_glm    = std::max(1e-10, cfg_.tol);

  if (design_in.p == 0) {
    // Intercept-only baseline (eta = offset, mu = link^{-1}(eta))
    Eigen::VectorXd eta  = offset_in;
    Eigen::VectorXd mu   = is_binary ? (1.0 / (1.0 + (-eta.array()).exp())).matrix() : eta;
    glm = {Eigen::VectorXd::Zero(0), eta, mu};
  } else if (is_binary) {
    glm = glm_logistic(X, y - offset_in, maxit_glm, tol_glm); // treat offset as prior eta, subtract here
    glm.eta.array() += offset_in.array();                     // restore total eta
    glm.mu  = (1.0 / (1.0 + (-glm.eta.array()).exp())).matrix();
  } else if (is_quant) {
    glm = glm_gaussian(X, y - offset_in);                    // linear link; offset handled by subtracting
    glm.eta.array() += offset_in.array();
    glm.mu = glm.eta;
  } else {
    throw std::invalid_argument("Unsupported trait: " + cfg_.trait);
  }


  // Back-transform beta if QR was used
  Eigen::VectorXd beta_cov;
  if (design_in.p == 0) {
    beta_cov = Eigen::VectorXd::Zero(0);
  } else if (cfg_.covariate_qr) {
    beta_cov = qr_backtransform_beta(glm.beta, qrmap, design_in.p, design_in.n);
  } else {
    beta_cov = glm.beta;
  }

  // ===== DEBUG: Print initial GLM coefficients =====
  std::cout << "\n===== C++ Initial GLM Coefficients =====" << std::endl;
  std::cout << "covariate_qr: " << (cfg_.covariate_qr ? "true" : "false") << std::endl;
  std::cout << "glm.beta (QR space): [";
  for (int i = 0; i < glm.beta.size(); ++i) {
    std::cout << glm.beta[i];
    if (i < glm.beta.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "beta_cov (original space): [";
  for (int i = 0; i < beta_cov.size(); ++i) {
    std::cout << beta_cov[i];
    if (i < beta_cov.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
  std::cout << "glm.eta[0:5]: " << glm.eta[0] << " " << glm.eta[1] << " " << glm.eta[2] << " " << glm.eta[3] << " " << glm.eta[4] << std::endl;
  std::cout << "glm.mu[0:5]: " << glm.mu[0] << " " << glm.mu[1] << " " << glm.mu[2] << " " << glm.mu[3] << " " << glm.mu[4] << std::endl;
  std::cout << "==========================================\n" << std::endl;


  // --- Compute final offset for GLMM ---
  // If covariate_offset=true, we fold fixed effects into offset and pass intercept-only to GLMM.
  // Otherwise, offset stays as provided (if any), and GLMM can refit fixed effects.
  // ensure offset_glmm has size n
  std::vector<double> offset_glmm;
  offset_glmm.reserve(design_in.n);

  if (!design_in.offset.empty()) {
      offset_glmm = design_in.offset;  // already aligned by preprocessing
  } else {
      offset_glmm.assign(design_in.n, 0.0);  // no offset given -> zeros
  }

  if (cfg_.covariate_offset) {
    // Recompute X * beta_cov on ORIGINAL design scale (X0) and add input offset
    if (design_in.p > 0) {
      Eigen::VectorXd xb = X_orig * beta_cov;
      for (int i = 0; i < design_in.n; ++i) offset_glmm[i] = offset_in[i] + xb[i];
    } else {
      for (int i = 0; i < design_in.n; ++i) offset_glmm[i] = offset_in[i];
    }
  } else {
    // Keep offset as input; GLMM will estimate fixed effects again.
    for (int i = 0; i < design_in.n; ++i) offset_glmm[i] = offset_in[i];
  }

  // (optional) sanity check
  if (offset_glmm.size() != static_cast<size_t>(design_in.n)) {
      throw std::runtime_error("offset_glmm must have length n; got " +
                              std::to_string(offset_glmm.size()) + " vs n=" +
                              std::to_string(design_in.n));
  }

  // Warm start is fine to keep as-is
  std::vector<double> beta_init;
  if (design_in.p > 0) {
      beta_init.assign(beta_cov.data(), beta_cov.data() + beta_cov.size());
  }

  // --- Call GLMM solver via hooks (you plug in your existing C++ kernels) ---
  FitNullResult out;
  if (is_binary) {
    if (!g_binary_solver) {
      throw std::runtime_error("Binary/survival GLMM solver not registered. Call register_binary_solver().");
    }
    // Optional warm start: pass coefficients matching current X space
    std::vector<double> beta_init;
    if (design_in.p > 0) {
      // If QR is enabled, design_in.X is in QR space, so pass glm.beta (QR space)
      // If QR is disabled, design_in.X is original, so pass beta_cov (original space)
      if (cfg_.covariate_qr && qrmap.valid) {
        beta_init.assign(glm.beta.data(), glm.beta.data() + glm.beta.size());
        std::cout << "[QR] Using glm.beta (QR space) as beta_init for GLMM solver\n";
      } else {
        beta_init.assign(beta_cov.data(), beta_cov.data() + beta_cov.size());
      }
    }

    out = g_binary_solver(paths_, cfg_, design_in, offset_glmm, beta_init);
    log("glmm solver Called-b") ;
  } else {
    if (!g_quant_solver) {
      throw std::runtime_error("Quantitative GLMM solver not registered. Call register_quant_solver().");
    }
    std::vector<double> beta_init;
    if (design_in.p > 0) {
      if (cfg_.covariate_qr && qrmap.valid) {
        beta_init.assign(glm.beta.data(), glm.beta.data() + glm.beta.size());
        std::cout << "[QR] Using glm.beta (QR space) as beta_init for GLMM solver\n";
      } else {
        beta_init.assign(beta_cov.data(), beta_cov.data() + beta_cov.size());
      }
    }
    out = g_quant_solver(paths_, cfg_, design_in, offset_glmm, beta_init);
  }

  // debug
  log("glmm solver Called") ;
  //

  // === Back-transform alpha from QR space to original space ===
  if (cfg_.covariate_qr && qrmap.valid && !out.alpha.empty()) {
    std::cout << "[QR] Back-transforming alpha from QR space to original space\n";
    std::cout << "[QR] alpha in QR space: [";
    for (size_t i = 0; i < out.alpha.size(); ++i) {
      std::cout << out.alpha[i];
      if (i < out.alpha.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";

    Eigen::VectorXd alpha_qr = Eigen::Map<Eigen::VectorXd>(out.alpha.data(), out.alpha.size());
    Eigen::VectorXd alpha_orig = qr_backtransform_beta(alpha_qr, qrmap, p_orig, design_in.n);

    out.alpha.resize(p_orig);
    for (int i = 0; i < p_orig; ++i) {
      out.alpha[i] = alpha_orig[i];
    }

    std::cout << "[QR] alpha in original space: [";
    for (size_t i = 0; i < out.alpha.size(); ++i) {
      std::cout << out.alpha[i];
      if (i < out.alpha.size() - 1) std::cout << ", ";
    }
    std::cout << "]\n";
  }

  // --- Attach offsets (ensure returned result has them) ---
  if (out.offset.empty()) {
    out.offset = std::move(offset_glmm);
  }
  out.loco        = cfg_.loco && chr_.enabled;
  out.lowmem_loco = cfg_.lowmem_loco;

  // --- LOCO batch (optional) ---
  if (out.loco && g_loco_batch) {
    log("Running LOCO batch");
    g_loco_batch(paths_, cfg_, chr_, design_in, out.theta, out.alpha, out.offset);
  }

  // --- Persist a compact JSON stub (optional; replace with your serializer) ---
  out.model_rda_path = default_model_path(paths_.out_prefix);
  try {
    ensure_parent_dir(out.model_rda_path);
    std::ofstream js(out.model_rda_path);
    js << "{\n";
    js << "  \"trait\": \"" << cfg_.trait << "\",\n";
    js << "  \"n\": " << design_in.n << ", \"p\": " << design_in.p << ",\n";
    js << "  \"theta\": [";
    for (size_t i = 0; i < out.theta.size(); ++i) js << (i ? "," : "") << out.theta[i];
    js << "],\n  \"alpha\": [";
    for (size_t i = 0; i < out.alpha.size(); ++i) js << (i ? "," : "") << out.alpha[i];
    js << "],\n  \"loco\": " << (out.loco ? "true" : "false")
       << ", \"lowmem_loco\": " << (out.lowmem_loco ? "true" : "false") << "\n";
    js << "}\n";
  } catch (...) {
    // Non-fatal: leave path set; upstream can decide how to handle
    log("Warning: failed to write model JSON; continuing.");
  }

  // --- Build obj_noK (ScoreNullPack) for Step 2 and save .arma binary files ---
  try {
    const int n_s = design_in.n;
    const int p_s = design_in.p;

    // Build arma matrices from design
    arma::fmat X_arma(n_s, p_s);
    for (int i = 0; i < n_s; ++i)
      for (int j = 0; j < p_s; ++j)
        X_arma(i, j) = static_cast<float>(design_in.X[i * p_s + j]);

    arma::fvec y_arma(n_s);
    for (int i = 0; i < n_s; ++i)
      y_arma(i) = static_cast<float>(design_in.y[i]);

    // Reconstruct mu from X*alpha + offset via link function
    arma::fvec alpha_arma(p_s);
    for (int i = 0; i < p_s; ++i)
      alpha_arma(i) = static_cast<float>(out.alpha[i]);

    arma::fvec eta_arma = X_arma * alpha_arma;
    if (!out.offset.empty()) {
      for (int i = 0; i < n_s; ++i)
        eta_arma(i) += static_cast<float>(out.offset[i]);
    }

    arma::fvec mu_arma(n_s);
    if (is_binary) {
      for (int i = 0; i < n_s; ++i) {
        float e = std::exp(eta_arma(i));
        mu_arma(i) = e / (1.0f + e);
      }
    } else {
      mu_arma = eta_arma;
    }

    // Build ScoreNull
    ScoreNull sn;
    if (is_binary) {
      sn = build_score_null_binary(X_arma, y_arma, mu_arma);
    } else {
      float tau0 = static_cast<float>(out.theta[0]);
      sn = build_score_null_quant(X_arma, y_arma, mu_arma, 1.0f / tau0);
    }

    // Convert to pack
    out.obj_noK = to_pack(sn, X_arma, y_arma, mu_arma, cfg_.trait);
    log("obj_noK (ScoreNullPack) populated");

    // --- Save .arma binary files for Step 2 ---
    std::string prefix = paths_.out_prefix;
    ensure_parent_dir(prefix + ".mu.arma");

    // Convert to double-precision arma vectors/matrices for saving
    arma::vec mu_d = arma::conv_to<arma::vec>::from(mu_arma);
    arma::vec res_d = arma::conv_to<arma::vec>::from(sn.res);
    arma::vec y_d = arma::conv_to<arma::vec>::from(y_arma);
    arma::vec V_d = arma::conv_to<arma::vec>::from(sn.V);
    arma::vec S_a_d = arma::conv_to<arma::vec>::from(sn.S_a);
    arma::mat X_d = arma::conv_to<arma::mat>::from(X_arma);
    arma::mat XV_d = arma::conv_to<arma::mat>::from(sn.XV);
    arma::mat XVX_d = arma::conv_to<arma::mat>::from(sn.XVX);
    arma::mat XVX_inv_d = arma::conv_to<arma::mat>::from(sn.XVX_inv);
    arma::mat XXVX_inv_d = arma::conv_to<arma::mat>::from(sn.XXVX_inv);
    arma::mat XVX_inv_XV_d = arma::conv_to<arma::mat>::from(sn.XVX_inv_XV);

    mu_d.save(prefix + ".mu.arma", arma::arma_binary);
    res_d.save(prefix + ".res.arma", arma::arma_binary);
    y_d.save(prefix + ".y.arma", arma::arma_binary);
    V_d.save(prefix + ".V.arma", arma::arma_binary);
    S_a_d.save(prefix + ".S_a.arma", arma::arma_binary);
    X_d.save(prefix + ".X.arma", arma::arma_binary);
    XV_d.save(prefix + ".XV.arma", arma::arma_binary);
    XVX_d.save(prefix + ".XVX.arma", arma::arma_binary);
    XVX_inv_d.save(prefix + ".XVX_inv.arma", arma::arma_binary);
    XXVX_inv_d.save(prefix + ".XXVX_inv.arma", arma::arma_binary);
    XVX_inv_XV_d.save(prefix + ".XVX_inv_XV.arma", arma::arma_binary);

    log("Saved .arma binary files for Step 2: mu, res, y, V, S_a, X, XV, XVX, XVX_inv, XXVX_inv, XVX_inv_XV");
  } catch (const std::exception& e) {
    log(std::string("Warning: failed to build obj_noK / save .arma files: ") + e.what());
  }

  return out;
}

} // namespace saige
