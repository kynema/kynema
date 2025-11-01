#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "elements/beams/interpolate_QP_vector.hpp"
#include "test_interpolate_QP.hpp"

namespace {

inline auto create_node_u_dot_OneNode() {
    return kynema::beams::tests::CreateView<double[1][1][6]>(
        "node_u_dot", std::array{1., 2., 3., 4., 5., 6.}
    );
}

inline auto create_node_u_dot_TwoNode() {
    return kynema::beams::tests::CreateView<double[1][2][6]>(
        "node_u_dot", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.}
    );
}

}  // namespace

namespace kynema::beams::tests {

TEST(InterpolateQPVectorTests, OneNodeOneQP) {
    constexpr auto num_qp = size_t{1U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeOneQP();
    auto node_u_dot = create_node_u_dot_OneNode();
    auto qp_u_dot = Kokkos::View<double[1][num_qp][3]>("qp_u_dot");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPVector<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_vector = Kokkos::subview(node_u_dot, Kokkos::ALL, Kokkos::ALL, Kokkos::pair(0, 3)),
            .qp_vector = qp_u_dot
        }
    );
    auto qp_u_dot_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u_dot);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 2), 6., tolerance);
}

TEST(InterpolateQPVectorTests, OneNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeTwoQP();
    auto node_u_dot = create_node_u_dot_OneNode();
    auto qp_u_dot = Kokkos::View<double[1][num_qp][3]>("qp_u_dot");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPVector<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_vector = Kokkos::subview(node_u_dot, Kokkos::ALL, Kokkos::ALL, Kokkos::pair(0, 3)),
            .qp_vector = qp_u_dot
        }
    );
    auto qp_u_dot_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u_dot);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 2), 6., tolerance);

    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 0), 4., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 1), 8., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 2), 12., tolerance);
}

TEST(InterpolateQPVectorTests, TwoNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{2U};
    auto shape_interp = create_shape_interp_TwoNodeTwoQP();
    auto node_u_dot = create_node_u_dot_TwoNode();
    auto qp_u_dot = Kokkos::View<double[1][num_qp][3]>("qp_u_dot");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPVector<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_vector = Kokkos::subview(node_u_dot, Kokkos::ALL, Kokkos::ALL, Kokkos::pair(0, 3)),
            .qp_vector = qp_u_dot
        }
    );
    auto qp_u_dot_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u_dot);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 0), 30., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 1), 36., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 0, 2), 42., tolerance);

    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 0), 38., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 1), 46., tolerance);
    EXPECT_NEAR(qp_u_dot_mirror(0, 1, 2), 54., tolerance);
}

}  // namespace kynema::beams::tests
