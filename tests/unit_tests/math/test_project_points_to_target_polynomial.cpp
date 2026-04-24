#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "math/project_points_to_target_polynomial.hpp"

namespace kynema_fmb::tests {

TEST(ProjectPointsToTargetPolynomialTest, Project2ndOrderTo4thOrderPolynomial) {
    const size_t num_input_pts = 3;
    const size_t num_output_pts = 5;
    const auto input_points = std::vector<std::array<double, 3>>{
        {0., 0., 0.},      // pt 1
        {2.5, -0.25, 0.},  // pt 2
        {5., 1., -1.}      // pt 3
    };

    const auto output_points =
        math::ProjectPointsToTargetPolynomial(num_input_pts, num_output_pts, input_points);

    // Expected projected points from Mathematica
    const std::vector<std::array<double, 3>> expected_projected_points = {
        {0., 0., 0.},  // First point - same as input pt 1
        {0.8633658232300573, -0.2558982639254172, 0.11304112106827431},  // Interior point 1
        {2.5, -0.25, 0.},  // Second point - same as input pt 2
        {4.1366341767699435, 0.39875540678255994, -0.5416125496397028},  // Interior point 2
        {5., 1., -1.}  // Last point - same as input pt 3
    };

    ASSERT_EQ(output_points.size(), expected_projected_points.size());
    for (auto i : std::views::iota(0U, output_points.size())) {
        for (auto j : std::views::iota(0U, 3U)) {
            EXPECT_NEAR(output_points[i][j], expected_projected_points[i][j], 1e-6)
                << "Mismatch at projected point [" << i << "][" << j << "]";
        }
    }
}

}  // namespace kynema_fmb::tests
