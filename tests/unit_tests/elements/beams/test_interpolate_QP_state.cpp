#include <array>
#include <cmath>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "elements/beams/interpolate_QP_state.hpp"
#include "test_interpolate_QP.hpp"

namespace {
inline auto create_node_u_OneNode() {
    return kynema::beams::tests::CreateView<double[1][1][7]>(
        "node_u", std::array{1., 2., 3., 4., 5., 6., 7.}
    );
}

inline auto create_node_u_TwoNode() {
    return kynema::beams::tests::CreateView<double[1][2][7]>(
        "node_u", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.}
    );
}

}  // namespace

namespace kynema::beams::tests {

TEST(InterpolateQPStateTests, u_OneNodeOneQP) {
    constexpr auto num_qp = size_t{1U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeOneQP();
    auto node_u = create_node_u_OneNode();
    auto qp_u = Kokkos::View<double[1][num_qp][3]>("qp_u");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_u<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_u = qp_u
        }
    );
    auto qp_u_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 2), 6., tolerance);
}

TEST(InterpolateQPStateTests, uprime_OneNodeOneQP) {
    constexpr auto num_qp = size_t{1U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_deriv = create_shape_deriv_OneNodeOneQP();
    auto jacobian = create_jacobian_OneQP();
    auto node_u = create_node_u_OneNode();
    auto qp_uprime = Kokkos::View<double[1][num_qp][3]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_uprime<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_uprime = qp_uprime
        }
    );
    auto qp_uprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_uprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 2), 6., tolerance);
}

TEST(InterpolateQPStateTests, r_OneNodeOneQP) {
    constexpr auto num_qp = size_t{1U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeOneQP();
    auto node_u = create_node_u_OneNode();
    auto qp_r = Kokkos::View<double[1][num_qp][4]>("qp_r");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_r<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_r = qp_r
        }
    );
    auto qp_r_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_r);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_r_mirror(0, 0, 0), 8. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 1), 10. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 2), 12. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 3), 14. / (6. * std::sqrt(14.)), tolerance);
}

TEST(InterpolateQPStateTests, rprime_OneNodeOneQP) {
    constexpr auto num_qp = size_t{1U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_deriv = create_shape_deriv_OneNodeOneQP();
    auto jacobian = create_jacobian_OneQP();
    auto node_u = create_node_u_OneNode();
    auto qp_rprime = Kokkos::View<double[1][num_qp][4]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_rprime<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_rprime = qp_rprime
        }
    );
    auto qp_rprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_rprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 0), 8., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 1), 10., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 2), 12., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 3), 14., tolerance);
}

TEST(InterpolateQPStateTests, u_OneNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeTwoQP();
    auto node_u = create_node_u_OneNode();
    auto qp_u = Kokkos::View<double[1][num_qp][3]>("qp_u");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_u<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_u = qp_u
        }
    );
    auto qp_u_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 2), 6., tolerance);

    EXPECT_NEAR(qp_u_mirror(0, 1, 0), 4., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 1, 1), 8., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 1, 2), 12., tolerance);
}

TEST(InterpolateQPStateTests, uprime_OneNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_deriv = create_shape_deriv_OneNodeTwoQP();
    auto jacobian = create_jacobian_TwoQP();
    auto node_u = create_node_u_OneNode();
    auto qp_uprime = Kokkos::View<double[1][num_qp][3]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_uprime<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_uprime = qp_uprime
        }
    );
    auto qp_uprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_uprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 0), 2., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 1), 4., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 2), 6., tolerance);

    EXPECT_NEAR(qp_uprime_mirror(0, 1, 0), 3., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 1, 1), 6., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 1, 2), 9., tolerance);
}

TEST(InterpolateQPStateTests, r_OneNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_interp = create_shape_interp_OneNodeTwoQP();
    auto node_u = create_node_u_OneNode();
    auto qp_r = Kokkos::View<double[1][num_qp][4]>("qp_r");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_r<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_r = qp_r
        }
    );
    auto qp_r_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_r);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_r_mirror(0, 0, 0), 8. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 1), 10. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 2), 12. / (6. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 3), 14. / (6. * std::sqrt(14.)), tolerance);

    EXPECT_NEAR(qp_r_mirror(0, 1, 0), 16. / (12. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 1), 20. / (12. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 2), 24. / (12. * std::sqrt(14.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 3), 28. / (12. * std::sqrt(14.)), tolerance);
}

TEST(InterpolateQPStateTests, rprime_OneNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{1U};
    auto shape_deriv = create_shape_deriv_OneNodeTwoQP();
    auto jacobian = create_jacobian_TwoQP();
    auto node_u = create_node_u_OneNode();
    auto qp_rprime = Kokkos::View<double[1][num_qp][4]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_rprime<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_rprime = qp_rprime
        }
    );
    auto qp_rprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_rprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 0), 8., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 1), 10., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 2), 12., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 3), 14., tolerance);

    EXPECT_NEAR(qp_rprime_mirror(0, 1, 0), 12., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 1), 15., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 2), 18., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 3), 21., tolerance);
}

TEST(InterpolateQPStateTests, u_TwoNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{2U};
    auto shape_interp = create_shape_interp_TwoNodeTwoQP();
    auto node_u = create_node_u_TwoNode();
    auto qp_u = Kokkos::View<double[1][num_qp][3]>("qp_u");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_u<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_u = qp_u
        }
    );
    auto qp_u_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_u);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_u_mirror(0, 0, 0), 34., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 1), 40., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 0, 2), 46., tolerance);

    EXPECT_NEAR(qp_u_mirror(0, 1, 0), 43., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 1, 1), 51., tolerance);
    EXPECT_NEAR(qp_u_mirror(0, 1, 2), 59., tolerance);
}

TEST(InterpolateQPStateTests, uprime_TwoNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{2U};
    auto shape_deriv = create_shape_deriv_TwoNodeTwoQP();
    auto jacobian = create_jacobian_TwoQP();
    auto node_u = create_node_u_TwoNode();
    auto qp_uprime = Kokkos::View<double[1][num_qp][3]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_uprime<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_uprime = qp_uprime
        }
    );
    auto qp_uprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_uprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 0), 34., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 1), 40., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 0, 2), 46., tolerance);

    EXPECT_NEAR(qp_uprime_mirror(0, 1, 0), 43., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 1, 1), 51., tolerance);
    EXPECT_NEAR(qp_uprime_mirror(0, 1, 2), 59., tolerance);
}

TEST(InterpolateQPStateTests, r_TwoNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{2U};
    auto shape_interp = create_shape_interp_TwoNodeTwoQP();
    auto node_u = create_node_u_TwoNode();
    auto qp_r = Kokkos::View<double[1][num_qp][4]>("qp_r");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_r<Kokkos::DefaultExecutionSpace>{
            .element = 0,
            .num_nodes = num_nodes,
            .shape_interp = shape_interp,
            .node_u = node_u,
            .qp_r = qp_r
        }
    );
    auto qp_r_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_r);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_r_mirror(0, 0, 0), 52. / (2. * std::sqrt(3766.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 1), 58. / (2. * std::sqrt(3766.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 2), 64. / (2. * std::sqrt(3766.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 0, 3), 70. / (2. * std::sqrt(3766.)), tolerance);

    EXPECT_NEAR(qp_r_mirror(0, 1, 0), 67. / (14. * std::sqrt(129.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 1), 75. / (14. * std::sqrt(129.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 2), 83. / (14. * std::sqrt(129.)), tolerance);
    EXPECT_NEAR(qp_r_mirror(0, 1, 3), 91. / (14. * std::sqrt(129.)), tolerance);
}

TEST(InterpolateQPStateTests, rprime_TwoNodeTwoQP) {
    constexpr auto num_qp = size_t{2U};
    constexpr auto num_nodes = size_t{2U};
    auto shape_deriv = create_shape_deriv_TwoNodeTwoQP();
    auto jacobian = create_jacobian_TwoQP();
    auto node_u = create_node_u_TwoNode();
    auto qp_rprime = Kokkos::View<double[1][num_qp][4]>("qp_uprime");
    Kokkos::parallel_for(
        num_qp,
        InterpolateQPState_rprime<Kokkos::DefaultExecutionSpace>{
            .element = 0U,
            .num_nodes = num_nodes,
            .shape_deriv = shape_deriv,
            .qp_jacobian = jacobian,
            .node_u = node_u,
            .qp_rprime = qp_rprime
        }
    );
    auto qp_rprime_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), qp_rprime);
    constexpr auto tolerance = 1.e-16;
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 0), 52., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 1), 58., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 2), 64., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 0, 3), 70., tolerance);

    EXPECT_NEAR(qp_rprime_mirror(0, 1, 0), 67., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 1), 75., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 2), 83., tolerance);
    EXPECT_NEAR(qp_rprime_mirror(0, 1, 3), 91., tolerance);
}

}  // namespace kynema::beams::tests
