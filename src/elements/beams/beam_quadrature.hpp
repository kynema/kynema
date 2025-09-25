#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include "math/gll_quadrature.hpp"
#include "math/interpolation.hpp"

namespace kynema::beams {

/**
 * @brief Creates a trapezoidal quadrature rule based on a given grid.
 *
 * This function generates a set of quadrature points and weights for trapezoidal integration
 * over a specified grid. The quadrature points are mapped to the interval [-1, 1].
 *
 * @param grid A span of grid points defining the integration domain.
 * @return A vector of arrays, where each array contains a quadrature point and its corresponding
 * weight.
 */
inline std::vector<std::array<double, 2>> CreateTrapezoidalQuadrature(std::span<const double> grid) {
    const auto n{grid.size()};
    const auto [grid_min, grid_max] = std::ranges::minmax(grid);
    const auto grid_range{grid_max - grid_min};
    auto quadrature = std::vector<std::array<double, 2>>{
        {-1., (grid[1] - grid[0]) / grid_range},
    };
    std::ranges::transform(
        std::views::iota(1U, n - 1), std::back_inserter(quadrature),
        [grid, gm = grid_min, grid_range](auto i) {
            return std::array{
                2. * (grid[i] - gm) / grid_range - 1., (grid[i + 1] - grid[i - 1]) / grid_range
            };
        }
    );
    quadrature.push_back({1., (grid[n - 1] - grid[n - 2]) / grid_range});
    return quadrature;
}

inline std::vector<std::array<double, 2>> CreateGaussLegendreLobattoQuadrature(
    std::span<const double> grid, size_t order
) {
    const auto n{grid.size()};
    const auto [grid_min, grid_max] = std::ranges::minmax(grid);
    const auto grid_range{grid_max - grid_min};
    const auto sectional_weights = math::GetGllWeights(order);
    const auto sectional_num_nodes = sectional_weights.size();
    const auto num_sections = (n - 1) / (sectional_num_nodes - 1);

    auto quadrature = std::vector<std::array<double, 2>>{};
    std::ranges::transform(
        grid, std::back_inserter(quadrature),
        [gm = grid_min, grid_range](auto grid_location) {
            return std::array{2. * (grid_location - gm) / grid_range - 1., 0.};
        }
    );

    auto section_index = 0UL;
    for ([[maybe_unused]] auto section : std::views::iota(0U, num_sections)) {
        const auto section_range =
            grid[section_index + sectional_num_nodes - 1] - grid[section_index];
        const auto weight_scaling = section_range / grid_range;
        for (auto node : std::views::iota(0U, sectional_num_nodes)) {
            quadrature[section_index + node][1] += sectional_weights[node] * weight_scaling;
        }
        section_index += sectional_num_nodes - 1UL;
    }
    return quadrature;
}

/**
 * @brief Creates Gauss-Legendre (GL) quadrature points and weights on [-1, 1]
 *
 * @details GL quadrature provides optimal accuracy for polynomial integration.
 *          The points are roots of P_n(x) where P_n is the nth Legendre polynomial.
 *          GL quadrature does NOT include the endpoints (-1, 1).
 *
 * @param order Number of quadrature points (n >= 1). Returns n quadrature points.
 * @return Vector of {point, weight} pairs for GL quadrature, sorted by point value
 */
inline std::vector<std::array<double, 2>> CreateGaussLegendreQuadrature(const size_t order) {
    constexpr auto max_iterations = 1000U;
    constexpr auto convergence_tolerance = std::numeric_limits<double>::epsilon();

    // If order is 0, return the base case
    if (order < 1) {
        return {{0., 2.}};  // point -> 0, weight->2
    }

    const size_t n_points{order};               // GL has n points for order n
    const size_t n_unique{(n_points + 1) / 2};  // Unique points due to symmetry
    std::vector<double> points(n_points, 0.);   // Quadrature points
    std::vector<double> weights(n_points, 0.);  // Quadrature weights

    // Find Gauss-Legendre (GL) points using Newton-Raphson method
    for (auto i_GL_point : std::views::iota(0U, n_unique)) {
        // Initial guess
        auto x_it = std::cos(
            std::numbers::pi * (static_cast<double>(i_GL_point) + 0.75) /
            (static_cast<double>(n_points) + 0.5)
        );

        // Newton-Raphson iteration to find root of P_n(x)
        bool converged = false;
        for ([[maybe_unused]] auto iteration : std::views::iota(0U, max_iterations)) {
            const auto x_old = x_it;  // Store old value for convergence check

            // Evaluate Legendre polynomials P_{n-1}(x) and P_n(x)
            const auto p_n_minus_1 = kynema::math::LegendrePolynomial(n_points - 1, x_it);
            const auto p_n = kynema::math::LegendrePolynomial(n_points, x_it);

            // Compute derivative P'_n(x) using recurrence relation
            const auto p_n_prime =
                static_cast<double>(n_points) * (x_it * p_n - p_n_minus_1) / (x_it * x_it - 1.);

            // Newton update: x_{n+1} = x_n - f(x_n)/f'(x_n)
            x_it -= p_n / p_n_prime;

            // Check convergence
            if (std::abs(x_it - x_old) <= convergence_tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error(
                "Newton-Raphson failed to converge for GL point " + std::to_string(i_GL_point) +
                " of order " + std::to_string(n_points)
            );
        }

        // Store symmetric points
        points[i_GL_point] = -x_it;                // left side
        points[n_points - 1 - i_GL_point] = x_it;  // right side

        // Compute GL weights: w = 2 / ((1 - x²) * [P'_n(x)]²)
        const auto p_n_minus_1 = kynema::math::LegendrePolynomial(n_points - 1, x_it);
        const auto p_n = kynema::math::LegendrePolynomial(n_points, x_it);
        const auto p_n_prime =
            static_cast<double>(n_points) * (x_it * p_n - p_n_minus_1) / (x_it * x_it - 1.);
        const auto weight = 2. / ((1. - x_it * x_it) * p_n_prime * p_n_prime);

        // Store symmetric weights
        weights[i_GL_point] = weight;
        weights[n_points - 1 - i_GL_point] = weight;
    }

    // Build and sort quadrature pairs
    auto quadrature = std::vector<std::array<double, 2>>(n_points);
    for (auto i_GL_point : std::views::iota(0U, n_points)) {
        quadrature[i_GL_point] = {points[i_GL_point], weights[i_GL_point]};
    }
    std::ranges::sort(quadrature, [](const auto& a, const auto& b) {
        return a[0] < b[0];
    });
    return quadrature;
}

}  // namespace kynema::beams
