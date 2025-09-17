#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "constraints/calculate_revolute_joint_force.hpp"
#include "create_view.hpp"

namespace kynema::tests {

struct ExecuteCalculateRevoluteJointForce {
    Kokkos::View<double[3][3]>::const_type axes;
    Kokkos::View<double[7]>::const_type constraint_inputs;
    Kokkos::View<double[7]>::const_type node_u;
    Kokkos::View<double[6]> residual_terms;

    KOKKOS_FUNCTION
    void operator()(int) const {
        constraints::CalculateRevoluteJointForce<Kokkos::DefaultExecutionSpace>::invoke(
            axes, constraint_inputs, node_u, residual_terms
        );
    }
};

TEST(CalculateRevoluteJointForceTests, OneConstraint) {
    const auto axes =
        CreateView<double[3][3]>("axes", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9.});
    const auto constraint_inputs =
        CreateView<double[7]>("constraint_inputs", std::array{4., 5., 6., 7., 8., 9., 10.});
    const auto node_u =
        CreateView<double[7]>("node_u", std::array{11., 12., 13., 14., 15., 16., 17.});

    const auto residual_terms = Kokkos::View<double[6]>("residual_terms");

    Kokkos::parallel_for(
        "CalculateRevoluteJointConstraint", 1,
        ExecuteCalculateRevoluteJointForce{axes, constraint_inputs, node_u, residual_terms}
    );

    const auto residual_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_terms);

    constexpr auto residual_terms_exact_data = std::array{0., 0., 0., 11032., 4816., 8008.};
    const auto residual_terms_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(residual_terms_exact_data.data());

    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(residual_terms_mirror(i), residual_terms_exact(i), 1.e-15);
    }
}

}  // namespace kynema::tests
