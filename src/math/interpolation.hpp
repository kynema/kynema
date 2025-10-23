#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace kynema::math {

/**
 * @brief Computes weights for linear interpolation
 *
 * @param x Evaluation point
 * @param xs Interpolation nodes (sorted)
 * @param weights Output: weights for linear interpolation
 */
inline void LinearInterpWeights(double x, std::span<const double> xs, std::vector<double>& weights) {
    const auto n = xs.size();
    weights.assign(n, 0.);

    const auto lower = std::ranges::lower_bound(xs, x);
    // If x is less than the first node, first weight -> 1 and our work is done
    if (lower == xs.begin()) {
        weights.front() = 1.;
        return;
    }
    // If x is greater than the last node, last weight -> 1 and our work is done
    if (lower == xs.end()) {
        weights.back() = 1.;
        return;
    }
    // x is between two nodes, so compute weights for closest nodes
    const auto index = static_cast<unsigned>(std::distance(xs.begin(), lower));
    const auto lower_loc = xs[index - 1];
    const auto upper_loc = xs[index];
    const auto weight_upper = (x - lower_loc) / (upper_loc - lower_loc);
    weights[index - 1] = 1. - weight_upper;
    weights[index] = weight_upper;
}

/**
 * @brief Computes linear interpolation
 *
 * @param x Evaluation point
 * @param xs Value locations
 * @param values Values at given locations
 * @return Interpolated value at evaluation point
 */
inline double LinearInterp(double x, std::span<const double> xs, std::span<const double> values) {
    std::vector<double> weights(xs.size());
    LinearInterpWeights(x, xs, weights);
    return std::transform_reduce(
        std::begin(weights), std::end(weights), std::begin(values), 0., std::plus<>(),
        std::multiplies<>()
    );
}

/**
 * @brief Computes weights for Lagrange polynomial interpolation
 *
 * @param x Evaluation point
 * @param xs Interpolation nodes (sorted)
 * @param weights Output: weights for Lagrange polynomial interpolation
 */
inline void LagrangePolynomialInterpWeights(
    double x, std::span<const double> xs, std::vector<double>& weights
) {
    const auto n = xs.size();
    weights.assign(n, 1.);

    // Pre-compute (x - xs[m]) terms to avoid repeated calculations
    std::vector<double> x_minus_xm(n);
    for (auto m : std::views::iota(0U, n)) {
        x_minus_xm[m] = x - xs[m];
    }

    for (auto j : std::views::iota(0U, n)) {
        for (auto m : std::views::iota(0U, n)) {
            if (j != m) {
                weights[j] *= x_minus_xm[m] / (xs[j] - xs[m]);
            }
        }
    }
}

/**
 * @brief Computes weights for Lagrange polynomial derivative interpolation
 *
 * @param x Evaluation point
 * @param xs Interpolation nodes (sorted)
 * @param weights Output: weights for Lagrange polynomial derivative interpolation
 */
inline void LagrangePolynomialDerivWeights(
    double x, std::span<const double> xs, std::vector<double>& weights
) {
    const auto n = xs.size();
    weights.assign(n, 0.);

    for (auto i : std::views::iota(0U, n)) {
        auto xi = xs[i];
        auto weight = 0.;
        for (auto j : std::views::iota(0U, n) | std::views::filter([i](unsigned j) {
                          return j != i;
                      })) {
            auto range = std::views::iota(0U, n) | std::views::filter([i, j](unsigned k) {
                             return k != i && k != j;
                         }) |
                         std::views::common;
            auto prod = std::transform_reduce(
                std::begin(range), std::end(range), 1., std::multiplies<>(),
                [&xs, x, xi](auto k) {
                    return (x - xs[k]) / (xi - xs[k]);
                }
            );
            weight += prod / (xi - xs[j]);
        }
        weights[i] = weight;
    }
}

/**
 * @brief Calculates the value of Legendre polynomial of order n at point x
 * @details Uses the recurrence relation for Legendre polynomials:
 *          P_n(x) = ((2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x))/n
 *          Reference: Deville et al. (2002) "High-Order Methods for Incompressible Fluid Flow"
 *          DOI: https://doi.org/10.1017/CBO9780511546792, Eq. B.1.15, p.446
 *
 * @param n Order of the Legendre polynomial (n >= 0)
 * @param x Point at which to evaluate the polynomial, typically in [-1,1]
 * @return Value of the nth order Legendre polynomial at x
 */
inline double LegendrePolynomial(const size_t n, const double x) {
    // Base cases
    if (n == 0) {
        return 1.;
    }
    if (n == 1) {
        return x;
    }

    // Compute the nth Legendre polynomial iteratively
    double p_n_minus_2{1.};  // P_{n-2}(x)
    double p_n_minus_1{x};   // P_{n-1}(x)
    double p_n{0.};          // P_n(x)
    for (auto i : std::views::iota(2U, n + 1)) {
        const auto i_double = static_cast<double>(i);
        p_n = ((2. * i_double - 1.) * x * p_n_minus_1 - (i_double - 1.) * p_n_minus_2) / i_double;
        p_n_minus_2 = p_n_minus_1;
        p_n_minus_1 = p_n;
    }
    return p_n;
}

/**
 * @brief Generates Gauss-Lobatto-Legendre (GLL) points for spectral element discretization
 * @details Computes the GLL points, i.e. roots of the Legendre polynomial, using
 *          Newton-Raphson iteration. GLL points are optimal interpolation nodes for
 *          spectral methods.
 *
 * @param order Order of the polynomial interpolation (must be >= 1)
 * @return Vector of GLL points sorted in ascending order, size = order + 1
 * @throws std::invalid_argument if order < 1
 * @throws std::runtime_error if Newton-Raphson iteration fails to converge
 */
inline std::vector<double> GenerateGLLPoints(const size_t order) {
    constexpr auto max_iterations = 1000U;
    constexpr auto convergence_tolerance = std::numeric_limits<double>::epsilon();

    if (order < 1) {
        throw std::invalid_argument("Polynomial order must be >= 1");
    }

    const size_t n_nodes = order + 1;
    auto gll_points = std::vector<double>(n_nodes, 0.);
    gll_points.resize(n_nodes);

    // Set the endpoints fixed at [-1, 1]
    gll_points.front() = -1.;
    gll_points.back() = 1.;

    // Find interior GLL points (1, ..., order - 1) using Newton-Raphson iteration
    std::vector<double> legendre_poly(n_nodes, 0.);
    for (auto i : std::views::iota(1U, order)) {
        // Initial guess using Chebyshev-Gauss-Lobatto nodes
        auto x_it =
            -std::cos(static_cast<double>(i) * std::numbers::pi / static_cast<double>(order));

        bool converged{false};
        for ([[maybe_unused]] auto iteration : std::views::iota(0U, max_iterations)) {
            const auto x_old = x_it;

            // Compute Legendre polynomials up to order n
            for (auto k : std::views::iota(0U, n_nodes)) {
                legendre_poly[k] = LegendrePolynomial(k, x_it);
            }

            // Newton update: x_{n+1} = x_n - f(x_n)/f'(x_n)
            const auto numerator = (x_it * legendre_poly[n_nodes - 1]) - legendre_poly[n_nodes - 2];
            const auto denominator = static_cast<double>(n_nodes) * legendre_poly[n_nodes - 1];
            x_it -= numerator / denominator;

            // Check for convergence
            if (std::abs(x_it - x_old) <= convergence_tolerance) {
                converged = true;
                break;
            }
        }

        if (!converged) {
            throw std::runtime_error(
                "Newton-Raphson iteration failed to converge for GLL point index " +
                std::to_string(i)
            );
        }

        gll_points[i] = x_it;
    }

    return gll_points;
}

}  // namespace kynema::math
