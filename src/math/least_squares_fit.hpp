#pragma once

#include <array>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

#include "interpolation.hpp"

namespace kynema::math {

/**
 * @brief Maps input geometric locations -> normalized domain using linear mapping
 *
 * @param geom_locations Input geometric locations (typically in domain [0, 1]),
 *                       sorted in ascending order
 * @return std::vector<double> Mapped/normalized evaluation points in domain [-1, 1]
 */
inline std::vector<double> MapGeometricLocations(std::span<const double> geom_locations) {
    // Get first and last points of the input domain (assumed to be sorted)
    const double domain_start = geom_locations.front();
    const double domain_end = geom_locations.back();
    if (domain_end == domain_start) {
        throw std::invalid_argument(
            "Invalid geometric locations: domain start and end points are equal."
        );
    }

    // Map each point from domain -> [-1, 1]
    std::vector<double> mapped_locations(geom_locations.size());
    const auto domain_span = domain_end - domain_start;
    std::ranges::transform(
        geom_locations, std::begin(mapped_locations),
        [domain_start, domain_span](auto geom_location) {
            return (2. * (geom_location - domain_start) / domain_span) - 1.;
        }
    );
    return mapped_locations;
}

/**
 * @brief Computes shape function matrices ϕg relating points ξb to ξg
 * At least two input points are required and it is assumed that there are more
 * output points than input points.
 *
 * @param input_points Input points, ξb, in [-1, 1]
 * @param output_points Output points, ξg, in [-1, 1]
 * @return Shape function matrix
 */
inline std::vector<std::vector<double>> ComputeShapeFunctionValues(
    std::span<const double> input_points, std::span<const double> output_points
) {
    // Number of points in input and output arrays
    const auto num_input_points = input_points.size();
    const auto num_output_points = output_points.size();

    // Compute weights for the shape functions and its derivatives at the evaluation points
    std::vector<double> weights(num_output_points, 0.);

    // Create shape function interpolation matrix
    auto shape_functions = std::vector<std::vector<double>>(
        num_output_points, std::vector<double>(num_input_points, 0.)
    );
    for (auto input_point : std::views::iota(0U, num_input_points)) {
        math::LagrangePolynomialInterpWeights(input_points[input_point], output_points, weights);
        for (auto output_point : std::views::iota(0U, num_output_points)) {
            shape_functions[output_point][input_point] = weights[output_point];
        }
    }

    return shape_functions;
}

/**
 * @brief Computes shape function derivatives dϕg relating points ξb to ξg
 * At least two input points are required and it is assumed that there are more
 * output points than input points.
 *
 * @param input_points Input points, ξb, in [-1, 1]
 * @param output_points Output points, ξg, in [-1, 1]
 * @return Shape function derivative matrix
 */
inline std::vector<std::vector<double>> ComputeShapeFunctionDerivatives(
    std::span<const double> input_points, std::span<const double> output_points
) {
    // Number of points in input and output arrays
    const auto num_input_points = input_points.size();
    const auto num_output_points = output_points.size();

    // Compute weights for the shape functions and its derivatives at the evaluation points
    std::vector<double> weights(num_output_points, 0.);

    // Create shape function derivative matrix
    auto derivative_functions = std::vector<std::vector<double>>(
        num_output_points, std::vector<double>(num_input_points, 0.)
    );
    for (auto input_point : std::views::iota(0U, num_input_points)) {
        math::LagrangePolynomialDerivWeights(input_points[input_point], output_points, weights);
        for (auto output_point : std::views::iota(0U, num_output_points)) {
            derivative_functions[output_point][input_point] = weights[output_point];
        }
    }

    return derivative_functions;
}

/**
 * @brief Performs least squares fitting to determine interpolation coefficients
 * @details Performs least squares fitting to determine interpolation coefficients
 *          by solving a dense linear system [A][X] = [B], where [A] is the shape
 *          function matrix (p x n), [B] is the input points (n x 3), and [X] is the
 *          interpolation coefficients (p x 3)
 *
 * @param shape_functions Shape function matrix (p x n)
 * @param points_to_fit x,y,z coordinates of the points to fit (n x 3)
 * @return Interpolation coefficients (p x 3)
 */
std::vector<std::array<double, 3>> PerformLeastSquaresFitting(
    std::span<const std::vector<double>> shape_functions,
    std::span<const std::array<double, 3>> points_to_fit
);
}  // namespace kynema::math
