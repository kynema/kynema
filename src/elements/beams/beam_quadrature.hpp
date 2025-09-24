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

/**
 * @brief Creates Gauss-Legendre quadrature points and weights on [-1, 1]
 *
 * @param order Number of quadrature points (n >= 1). Uses P_n for n points.
 * @return Vector of {point, weight} pairs for Gauss-Legendre quadrature
 */
inline std::vector<std::array<double, 2>> CreateGaussLegendreQuadrature(const size_t order) {
    constexpr auto max_iterations = 1000U;
    constexpr auto convergence_tolerance = std::numeric_limits<double>::epsilon();

    // If order is 0, return the base case
    if (order < 1) {
        return {{0., 2.}};  // point -> 0, weight->2
    }

    const size_t n{order};        // Total number of quadrature points to be calculated
    const size_t m{(n + 1) / 2};  // Number of unique Gauss-Legendre points needed due to symmetry
    std::vector<double> points(n, 0.);   // Quadrature points
    std::vector<double> weights(n, 0.);  // Quadrature weights

    // Find the Gauss-Legendre points using Newton-Raphson method
    for (auto i_GL_point : std::views::iota(0U, m)) {
        // Initial guess
        auto x_it = std::cos(
            std::numbers::pi * (static_cast<double>(i_GL_point) + 0.75) /
            (static_cast<double>(n) + 0.5)
        );

        // Newton-Raphson iterations using P_n, P_{n-1}, and P'_n relation where
        // P_n      ->  nth order Legendre polynomial
        // P_{n-1}  ->  (n-1)th order Legendre polynomial
        // P'_n     ->  derivative of the nth order Legendre polynomial
        bool converged{false};
        for ([[maybe_unused]] auto it : std::views::iota(0U, max_iterations)) {
            const auto p_n_minus_1 = kynema::math::LegendrePolynomial(n - 1, x_it);
            const auto p_n = kynema::math::LegendrePolynomial(n, x_it);
            const auto p_n_prime =
                static_cast<double>(n) * (x_it * p_n - p_n_minus_1) / (x_it * x_it - 1.);

            // Newton update: x_{n+1} = x_n - f(x_n)/f'(x_n)
            const auto x_old = x_it;  // Store old value for convergence check
            x_it -= p_n / p_n_prime;  // Newton-Raphson update

            // Check for convergence
            if (std::abs(x_it - x_old) <= convergence_tolerance) {
                converged = true;
                break;
            }
        }
        if (!converged) {
            throw std::runtime_error(
                "Newton-Raphson iteration failed to converge for GLL point index " +
                std::to_string(i_GL_point)
            );
        }

        // Symmetric nodes i.e. points_i = -x_it and points_{n-1-i} = x_it
        points[i_GL_point] = -x_it;
        points[n - 1 - i_GL_point] = x_it;

        // Weight: w = 2 / ((1 - x_it^2) [P'_n(x_it)]^2)
        const auto p_n_minus_1 = kynema::math::LegendrePolynomial(n - 1, x_it);
        const auto p_n = kynema::math::LegendrePolynomial(n, x_it);
        const auto p_n_prime =
            static_cast<double>(n) * (x_it * p_n - p_n_minus_1) / (x_it * x_it - 1.);
        const auto weight = 2. / ((1. - x_it * x_it) * p_n_prime * p_n_prime);

        // Symmetric nodes i.e. weights_i = weights_{n-1-i}
        weights[i_GL_point] = weight;
        weights[n - 1 - i_GL_point] = weight;
    }

    // Return the quadrature points and weights by sorting the points in ascending order
    auto quadrature = std::vector<std::array<double, 2>>(n);
    for (auto i_GL_point : std::views::iota(0U, n)) {
        quadrature[i_GL_point] = {points[i_GL_point], weights[i_GL_point]};
    }
    std::ranges::sort(quadrature, [](const auto& a, const auto& b) {
        return a[0] < b[0];
    });
    return quadrature;
}

}  // namespace kynema::beams
