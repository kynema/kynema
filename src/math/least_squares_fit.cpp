#if defined __NVCC__
#pragma nv_diag_once 20014, 20015
#pragma nv_diag_suppress 20014, 20015
#endif
#include "least_squares_fit.hpp"

#include <Eigen/Dense>

namespace kynema::math {
std::vector<std::array<double, 3>> PerformLeastSquaresFitting(
    std::span<const std::vector<double>> shape_functions,
    std::span<const std::array<double, 3>> points_to_fit
) {
    const auto p = static_cast<unsigned>(shape_functions.size());
    const auto n = shape_functions.front().size();
    if (std::ranges::any_of(shape_functions, [n](const auto& row) {
            return row.size() != n;
        })) {
        throw std::invalid_argument("Inconsistent number of columns in shape_functions.");
    }

    auto S = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(p, n);
    for (auto j : std::views::iota(0U, shape_functions.front().size())) {
        for (auto i : std::views::iota(0U, p)) {
            S(i, j) = shape_functions[i][j];
        }
    }

    auto P = Eigen::Matrix<double, Eigen::Dynamic, 3>(n, 3);
    for (auto j : std::views::iota(0U, 3U)) {
        for (auto i : std::views::iota(0U, points_to_fit.size())) {
            P(i, j) = points_to_fit[i][j];
        }
    }

    auto A = (S * S.transpose()).eval();
    A(0, 0) = 1.;
    A(p - 1U, p - 1U) = 1.;
    for (auto i : std::views::iota(0U, p - 1U)) {
        A(0, i + 1U) = 0.;
        A(p - 1U, i) = 0.;
    }

    auto B = (S * P).eval();
    for (auto dim : std::views::iota(0U, 3U)) {
        B(0, dim) = points_to_fit[0][dim];
        B(p - 1U, dim) = points_to_fit[n - 1][dim];
    }

    auto lu =
        Eigen::PartialPivLU<Eigen::Ref<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>(A);
    auto x = lu.solve(B).eval();
    auto result = std::vector<std::array<double, 3>>(p);
    for (auto i : std::views::iota(0U, p)) {
        for (auto j : std::views::iota(0U, 3U)) {
            result[i][j] = x(i, j);
        }
    }

    return result;
}
}  // namespace kynema::math
