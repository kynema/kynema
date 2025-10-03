#include <array>
#include <cstddef>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "constraints/calculate_revolute_joint_output.hpp"
#include "create_view.hpp"

namespace kynema::tests {

struct ExecuteCalculateRevoluteJointOutput {
    int i_constraint;
    Kokkos::View<size_t*>::const_type target_node_index;
    Kokkos::View<double* [3][3]>::const_type axes;
    Kokkos::View<double* [7]>::const_type node_x0;
    Kokkos::View<double* [7]>::const_type node_u;
    Kokkos::View<double* [6]>::const_type node_udot;
    Kokkos::View<double* [6]>::const_type node_uddot;
    Kokkos::View<double* [3]> outputs;

    KOKKOS_FUNCTION
    void operator()(int) const {
        constraints::CalculateRevoluteJointOutput<Kokkos::DefaultExecutionSpace>{
            i_constraint, target_node_index, axes, node_x0, node_u, node_udot, node_uddot, outputs
        }();
    }
};

TEST(CalculateRevoluteJointOutputTests, OneConstraint) {
    const auto target_node_index =
        CreateView<size_t[1]>("target_node_index", std::array<size_t, 1>{1UL});
    const auto axes =
        CreateView<double[1][3][3]>("axes", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9.});
    const auto node_x0 = CreateView<double[2][7]>(
        "node_x0", std::array{0., 0., 0., 0., 0., 0., 0., 11., 12., 13., 14., 15., 16., 17.}
    );
    const auto node_u = CreateView<double[2][7]>(
        "node_u", std::array{0., 0., 0., 0., 0., 0., 0., 18., 19., 20., 21., 22., 23., 24.}
    );
    const auto node_udot = CreateView<double[2][6]>(
        "node_udot", std::array{0., 0., 0., 0., 0., 0., 25., 26., 27., 28., 29., 30.}
    );
    const auto node_uddot = CreateView<double[2][6]>(
        "node_uddot", std::array{0., 0., 0., 0., 0., 0., 31., 32., 33., 34., 35., 36.}
    );

    const auto outputs = Kokkos::View<double[1][3]>("outputs");

    Kokkos::parallel_for(
        "CalculateRevoluteJointOutput", 1,
        ExecuteCalculateRevoluteJointOutput{
            0, target_node_index, axes, node_x0, node_u, node_udot, node_uddot, outputs
        }
    );

    const auto outputs_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), outputs);

    constexpr auto outputs_exact_data = std::array{0., 358792., 433384.};
    const auto outputs_exact =
        Kokkos::View<double[1][3], Kokkos::HostSpace>::const_type(outputs_exact_data.data());

    for (auto i : std::views::iota(0, 3)) {
        EXPECT_NEAR(outputs_mirror(0, i), outputs_exact(0, i), 1.e-15);
    }
}

}  // namespace kynema::tests
