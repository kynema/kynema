#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <span>
#include <vector>

#include "math/gll_quadrature.hpp"

namespace kynema::beams {

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

}  // namespace kynema::beams
