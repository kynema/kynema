#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <ranges>
#include <span>
#include <vector>

#include "math/gl_quadrature.hpp"
#include "math/gll_quadrature.hpp"

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
                (2. * (grid[i] - gm) / grid_range) - 1., (grid[i + 1] - grid[i - 1]) / grid_range
            };
        }
    );
    quadrature.push_back({1., (grid[n - 1] - grid[n - 2]) / grid_range});
    return quadrature;
}

inline std::vector<std::array<double, 2>> CreateGaussLegendreLobattoQuadrature(
    std::span<const double> grid, std::span<const double> original_grid, size_t order
) {
    const auto n{grid.size()};
    const auto [grid_min, grid_max] = std::ranges::minmax(grid);
    const auto grid_range{grid_max - grid_min};
    const auto sectional_weights = math::GetGllWeights(order);
    const auto sectional_num_nodes = sectional_weights.size();
    const auto num_sections = (n - 1) / (sectional_num_nodes - 1);

    auto quadrature = std::vector<std::array<double, 2>>{};
    std::ranges::transform(
        grid, std::back_inserter(quadrature), [gm = grid_min, grid_range](auto grid_location) {
            return std::array{(2. * (grid_location - gm) / grid_range) - 1., 0.};
        }
    );

    auto section_index = 0UL;
    for ([[maybe_unused]] auto section : std::views::iota(0U, num_sections)) {
        const auto section_range = original_grid[section + 1] - original_grid[section];
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
 * @param grid grid locations of each quadrature point
 * @param original_grid grid locations of the origionally specified sections
 * @param order Number of quadrature points (n >= 1). Returns n quadrature points.
 * @return Vector of {point, weight} pairs for GL quadrature, sorted by point value
 */
inline std::vector<std::array<double, 2>> CreateGaussLegendreQuadrature(
    std::span<const double> grid, std::span<const double> original_grid, size_t order
) {
    const auto n{grid.size()};
    const auto [grid_min, grid_max] = std::ranges::minmax(original_grid);
    const auto grid_range{grid_max - grid_min};
    const auto sectional_weights = math::GetGlWeights(order);
    const auto sectional_num_nodes = sectional_weights.size();
    const auto num_sections = n / sectional_num_nodes;

    auto quadrature = std::vector<std::array<double, 2>>{};
    std::ranges::transform(
        grid, std::back_inserter(quadrature), [gm = grid_min, grid_range](auto grid_location) {
            return std::array{(2. * (grid_location - gm) / grid_range) - 1., 0.};
        }
    );

    auto section_index = 0UL;
    for ([[maybe_unused]] auto section : std::views::iota(0U, num_sections)) {
        const auto section_range = original_grid[section + 1] - original_grid[section];
        const auto weight_scaling = section_range / grid_range;
        for (auto node : std::views::iota(0U, sectional_num_nodes)) {
            quadrature[section_index + node][1] += sectional_weights[node] * weight_scaling;
        }
        section_index += sectional_num_nodes;
    }
    return quadrature;
}

}  // namespace kynema::beams
