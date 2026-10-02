// Verbatim copy of src/firth_moment.hpp from the collaborator package
// saige-cpp-step2-source-demo-20260927 (GPLv3). Kept unmodified so the CPU
// reference for the max-step=1 estimator is exactly what that package runs.
#ifndef SAIGE_FIRTH_MOMENT_HPP
#define SAIGE_FIRTH_MOMENT_HPP

#include <armadillo>

#include <algorithm>
#include <cmath>
#include <limits>

namespace SAIGE {
namespace firth {

enum class FitFailure {
    none,
    invalid_input,
    singular_information,
    non_finite,
    line_search_failed
};

struct FitOptions {
    int max_iterations = 50;
    double max_step = 15.0;
    double gradient_tolerance = 1e-5;
    double parameter_tolerance = 1e-5;
    double determinant_tolerance =
        64.0 * std::numeric_limits<double>::epsilon();
    // Only large Fisher steps pay for an objective evaluation. Near the
    // solution, the original exact-equivalent moment iteration is retained.
    double line_search_trigger = 1.0;
    int max_halvings = 15;
    double objective_tolerance = 1e-10;
};

struct FitResult {
    double intercept = std::numeric_limits<double>::quiet_NaN();
    double beta = std::numeric_limits<double>::quiet_NaN();
    double se_beta = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    int line_search_evaluations = 0;
    int step_halvings = 0;
    bool converged = false;
    bool hit_iteration_limit = false;
    FitFailure failure = FitFailure::none;
};

inline double stable_sigmoid(double eta) noexcept {
    if (eta >= 0.0) {
        const double e = std::exp(-eta);
        return 1.0 / (1.0 + e);
    }
    const double e = std::exp(eta);
    return e / (1.0 + e);
}

inline double log_one_plus_exp(double eta) noexcept {
    if (eta > 0.0) return eta + std::log1p(std::exp(-eta));
    return std::log1p(std::exp(eta));
}

inline bool evaluate_penalized_objective(
    const arma::vec& z,
    const arma::vec& y,
    const arma::vec& offset,
    double alpha,
    double beta,
    bool use_firth,
    double determinant_tolerance,
    double& objective) {

    double log_likelihood = 0.0;
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    for (arma::uword i = 0; i < z.n_elem; ++i) {
        const double zi = z(i);
        const double eta = offset(i) + alpha + beta * zi;
        const double mu = stable_sigmoid(eta);
        const double w = mu * (1.0 - mu);
        log_likelihood += y(i) * eta - log_one_plus_exp(eta);
        s0 += w;
        s1 += w * zi;
        s2 += w * zi * zi;
    }
    const double determinant = s0 * s2 - s1 * s1;
    const double information_scale =
        std::max({1.0, std::abs(s0 * s2), std::abs(s1 * s1)});
    if (!std::isfinite(log_likelihood) || !std::isfinite(determinant) ||
        determinant <= determinant_tolerance * information_scale) {
        return false;
    }
    objective = log_likelihood;
    if (use_firth) objective += 0.5 * std::log(determinant);
    return std::isfinite(objective);
}

// Exact algebraic replacement for the QR/hat-vector Fisher-scoring iteration
// when X = [1, z]. The score, information, and leverage adjustment are exact
// moment identities. Safeguarded step-halving is activated only for large
// proposed steps and therefore changes the iteration path, not the objective.
inline FitResult fit_two_parameter_moment(
    const arma::vec& z,
    const arma::vec& y,
    const arma::vec& offset,
    bool use_firth,
    const arma::vec& init,
    const FitOptions& options = FitOptions()) {

    FitResult result;
    const arma::uword n = z.n_elem;
    if (n == 0 || y.n_elem != n || offset.n_elem != n || init.n_elem < 2 ||
        options.max_iterations <= 0 || options.max_step <= 0.0 ||
        options.max_halvings < 0 || options.line_search_trigger < 0.0) {
        result.failure = FitFailure::invalid_input;
        return result;
    }

    double alpha = init(0);
    double beta = init(1);
    double cov_beta_beta = std::numeric_limits<double>::quiet_NaN();

    for (int iteration = 0; iteration < options.max_iterations; ++iteration) {
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;
        double t0 = 0.0;
        double t1 = 0.0;
        double t2 = 0.0;
        double t3 = 0.0;
        double r0 = 0.0;
        double r1 = 0.0;
        double log_likelihood = 0.0;

        for (arma::uword i = 0; i < n; ++i) {
            const double zi = z(i);
            const double eta = offset(i) + alpha + beta * zi;
            const double mu = stable_sigmoid(eta);
            const double w = mu * (1.0 - mu);
            const double z2 = zi * zi;

            log_likelihood += y(i) * eta - log_one_plus_exp(eta);
            s0 += w;
            s1 += w * zi;
            s2 += w * z2;
            r0 += y(i) - mu;
            r1 += zi * (y(i) - mu);

            if (use_firth) {
                const double a = w * (0.5 - mu);
                t0 += a;
                t1 += a * zi;
                t2 += a * z2;
                t3 += a * z2 * zi;
            }
        }

        const double determinant = s0 * s2 - s1 * s1;
        const double information_scale =
            std::max({1.0, std::abs(s0 * s2), std::abs(s1 * s1)});
        if (!std::isfinite(determinant) ||
            determinant <= options.determinant_tolerance * information_scale) {
            result.failure = FitFailure::singular_information;
            return result;
        }

        double score_alpha = r0;
        double score_beta = r1;
        if (use_firth) {
            score_alpha +=
                (s2 * t0 - 2.0 * s1 * t1 + s0 * t2) / determinant;
            score_beta +=
                (s2 * t1 - 2.0 * s1 * t2 + s0 * t3) / determinant;
        }

        double delta_alpha =
            (s2 * score_alpha - s1 * score_beta) / determinant;
        double delta_beta =
            (-s1 * score_alpha + s0 * score_beta) / determinant;
        if (!std::isfinite(delta_alpha) || !std::isfinite(delta_beta) ||
            !std::isfinite(log_likelihood)) {
            result.failure = FitFailure::non_finite;
            return result;
        }

        double largest_step =
            std::max(std::abs(delta_alpha), std::abs(delta_beta));
        if (largest_step > options.max_step) {
            const double scale = options.max_step / largest_step;
            delta_alpha *= scale;
            delta_beta *= scale;
            largest_step = options.max_step;
        }

        if (largest_step > options.line_search_trigger) {
            double current_objective = log_likelihood;
            if (use_firth) current_objective += 0.5 * std::log(determinant);
            const double allowed_drop = options.objective_tolerance *
                (1.0 + std::abs(current_objective));
            bool accepted = false;
            double scale = 1.0;
            for (int hs = 0; hs <= options.max_halvings; ++hs) {
                double candidate_objective = 0.0;
                ++result.line_search_evaluations;
                if (evaluate_penalized_objective(
                        z, y, offset,
                        alpha + scale * delta_alpha,
                        beta + scale * delta_beta,
                        use_firth, options.determinant_tolerance,
                        candidate_objective) &&
                    candidate_objective >= current_objective - allowed_drop) {
                    delta_alpha *= scale;
                    delta_beta *= scale;
                    result.step_halvings += hs;
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!accepted) {
                result.failure = FitFailure::line_search_failed;
                return result;
            }
        }

        cov_beta_beta = s0 / determinant;
        alpha += delta_alpha;
        beta += delta_beta;
        result.iterations = iteration + 1;

        const bool small_step =
            std::max(std::abs(delta_alpha), std::abs(delta_beta)) <=
            options.parameter_tolerance;
        const bool small_score =
            std::max(std::abs(score_alpha), std::abs(score_beta)) <=
            options.gradient_tolerance;
        if (small_step && small_score) {
            result.converged = true;
            break;
        }
        if (result.iterations == options.max_iterations) {
            result.hit_iteration_limit = true;
        }
    }

    if (!std::isfinite(alpha) || !std::isfinite(beta) ||
        !std::isfinite(cov_beta_beta) || cov_beta_beta < 0.0) {
        result.failure = FitFailure::non_finite;
        return result;
    }

    result.intercept = alpha;
    result.beta = beta;
    result.se_beta = std::sqrt(cov_beta_beta);
    return result;
}

}  // namespace firth
}  // namespace SAIGE

#endif  // SAIGE_FIRTH_MOMENT_HPP
