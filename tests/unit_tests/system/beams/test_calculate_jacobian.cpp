#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "elements/beams/calculate_jacobian.hpp"
#include "test_calculate.hpp"

namespace kynema::beams::tests {

TEST(CalculateJacobian, LinearElement) {
    constexpr size_t num_elems{1};
    constexpr size_t num_nodes{2};
    constexpr size_t num_qps{1};

    const auto num_nodes_per_elem = CreateView<size_t[1]>("num_nodes", std::array{num_nodes});
    const auto num_qps_per_elem = CreateView<size_t[1]>("num_qps", std::array{num_qps});
    const auto shape_derivative =
        CreateView<double[1][2][1]>("shape_derivative", std::array{-0.5, 0.5});
    const auto node_position_rotation = CreateView<double[1][2][7]>(
        "node_position_rotation", std::array{-1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.}
    );
    const auto qp_position_derivative = Kokkos::View<double[1][1][3]>("position_derivative");
    const auto qp_jacobian = Kokkos::View<double[1][1]>("qp_jacobian");

    const auto calculate_jacobian = CalculateJacobian<Kokkos::DefaultExecutionSpace>{
        .num_nodes_per_element = num_nodes_per_elem,
        .num_qps_per_element = num_qps_per_elem,
        .shape_derivative = shape_derivative,
        .node_position_rotation = node_position_rotation,
        .qp_position_derivative = qp_position_derivative,
        .qp_jacobian = qp_jacobian
    };
    Kokkos::parallel_for("calculate_jacobian", 1, calculate_jacobian);

    auto host_jacobian = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_jacobian);

    // For a linear element of length 2., the Jacobian should be 1.
    ASSERT_EQ(host_jacobian.extent(0), num_elems);  // 1 element
    ASSERT_EQ(host_jacobian.extent(1), num_qps);    // 1 quadrature point
    EXPECT_DOUBLE_EQ(host_jacobian(0, 0), 1.);

    auto host_position_derivative =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_position_derivative);

    // Should be a unit vector in x direction
    ASSERT_EQ(host_position_derivative.extent(0), num_elems);  // 1 element
    ASSERT_EQ(host_position_derivative.extent(1), num_qps);    // 1 quadrature point
    ASSERT_EQ(host_position_derivative.extent(2), 3);          // 3 dimensions
    EXPECT_DOUBLE_EQ(host_position_derivative(0, 0, 0), 1.);
    EXPECT_DOUBLE_EQ(host_position_derivative(0, 0, 1), 0.);
    EXPECT_DOUBLE_EQ(host_position_derivative(0, 0, 2), 0.);
}

TEST(CalculateJacobian, FourthOrderElement) {
    constexpr size_t num_elems{1};
    constexpr size_t num_nodes{5};  // 4th order = 5 nodes
    constexpr size_t num_qps{5};    // 5 quadrature points

    const auto num_nodes_per_elem =
        CreateView<size_t[num_elems]>("num_nodes", std::array{num_nodes});
    const auto num_qps_per_elem = CreateView<size_t[num_elems]>("num_qps", std::array{num_qps});
    const auto node_position_rotation = CreateView<double[num_elems][num_nodes][7]>(
        "node_position_rotation",
        std::array{
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.16237631096713473,
            0.17578464768961147,
            0.1481911137890286,
            0.,
            0.,
            0.,
            0.,
            0.25,
            1.,
            1.1875,
            0.,
            0.,
            0.,
            0.,
            -0.30523345382427747,
            2.4670724951675314,
            2.953849702537502,
            0.,
            0.,
            0.,
            0.,
            -1.,
            3.5,
            4.,
            0.,
            0.,
            0.,
            0.
        }
    );
    const auto shape_derivative = CreateView<double[num_elems][num_nodes][num_qps]>(
        "shape_derivative",
        std::array{
            -3.70533645359145147,     -0.528715267980273795, 0.375000000000000111,
            -0.243518021129600221,    0.144236409367993451,  4.33282116876392998,
            -1.0976579678283287,      -1.33658457769545347,  0.74973857001328692,
            -0.420676230427679154,    -0.903924536232162623, 2.13259378469228889,
            -2.22044604925031308e-16, -2.13259378469228977,  0.90392453623216229,
            0.420676230427679099,     -0.749738570013287031, 1.33658457769545347,
            1.09765796782832847,      -4.33282116876392909,  -0.144236409367993451,
            0.243518021129600221,     -0.375000000000000111, 0.528715267980273684,
            3.70533645359145192,
        }
    );
    const auto qp_position_derivative =
        Kokkos::View<double[num_elems][num_nodes][3]>("position_derivative");
    const auto qp_jacobian = Kokkos::View<double[num_elems][num_nodes]>("qp_jacobian");

    const auto calculate_jacobian = CalculateJacobian<Kokkos::DefaultExecutionSpace>{
        .num_nodes_per_element = num_nodes_per_elem,
        .num_qps_per_element = num_qps_per_elem,
        .shape_derivative = shape_derivative,
        .node_position_rotation = node_position_rotation,
        .qp_position_derivative = qp_position_derivative,
        .qp_jacobian = qp_jacobian
    };
    Kokkos::parallel_for("calculate_jacobian", 1, calculate_jacobian);

    const auto host_qp_jacobian =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_jacobian);

    // Expected jacobians at each quadrature point (from BeamDyn)
    ASSERT_EQ(host_qp_jacobian.extent(0), num_elems);  // 1 element
    ASSERT_EQ(host_qp_jacobian.extent(1), num_qps);    // 5 quadrature points

    constexpr std::array<double, 5> expected_jacobians = {
        0.671587005850145, 1.509599209717606, 2.861380785564898, 4.097191592895187, 4.880926263217582
    };
    for (auto qp : std::views::iota(0U, num_qps)) {
        EXPECT_NEAR(host_qp_jacobian(0, qp), expected_jacobians[qp], 1e-12);
    }

    const auto host_qp_position_derivative =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_position_derivative);

    // Verify that position derivatives are unit vectors
    ASSERT_EQ(host_qp_position_derivative.extent(0), num_elems);  // 1 element
    ASSERT_EQ(host_qp_position_derivative.extent(1), num_qps);    // 5 quadrature points
    ASSERT_EQ(host_qp_position_derivative.extent(2), 3);          // 3 dimensions

    for (auto qp : std::views::iota(0U, num_qps)) {
        const auto magnitude = std::sqrt(
            (host_qp_position_derivative(0, qp, 0) * host_qp_position_derivative(0, qp, 0)) +
            (host_qp_position_derivative(0, qp, 1) * host_qp_position_derivative(0, qp, 1)) +
            (host_qp_position_derivative(0, qp, 2) * host_qp_position_derivative(0, qp, 2))
        );
        EXPECT_NEAR(magnitude, 1., 1e-12);  // unit vector
    }
}

}  // namespace kynema::beams::tests
