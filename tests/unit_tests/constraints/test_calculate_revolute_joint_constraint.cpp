#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "constraints/calculate_revolute_joint_constraint.hpp"
#include "create_view.hpp"

namespace kynema::tests {

struct ExecuteCalculateRevoluteJointConstraint {
    Kokkos::View<double[3]>::const_type X0;
    Kokkos::View<double[3][3]>::const_type axes;
    Kokkos::View<double[7]>::const_type base_node_u;
    Kokkos::View<double[7]>::const_type target_node_u;
    Kokkos::View<double[6]> residual_terms;
    Kokkos::View<double[6][6]> base_gradient_terms;
    Kokkos::View<double[6][6]> target_gradient_terms;

    KOKKOS_FUNCTION
    void operator()(int) const {
        constraints::CalculateRevoluteJointConstraint<Kokkos::DefaultExecutionSpace>::invoke(
            X0, axes, base_node_u, target_node_u, residual_terms, base_gradient_terms,
            target_gradient_terms
        );
    }
};

TEST(CalculateRevoluteJointConstraintTests, OneConstraint) {
    const auto X0 = CreateView<double[3]>("X0", std::array{1., 2., 3.});
    const auto axes =
        CreateView<double[3][3]>("axes", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9.});
    const auto base_node_u =
        CreateView<double[7]>("base_node_u", std::array{18., 19., 20., 21., 22., 23., 24.});
    const auto target_node_u =
        CreateView<double[7]>("target_node_u", std::array{11., 12., 13., 14., 15., 16., 17.});

    const auto residual_terms = Kokkos::View<double[6]>("residual_terms");
    const auto base_gradient_terms = Kokkos::View<double[6][6]>("base_gradient_terms");
    const auto target_gradient_terms = Kokkos::View<double[6][6]>("target_gradient_terms");

    Kokkos::parallel_for(
        "CalculateRevoluteJointConstraint", 1,
        ExecuteCalculateRevoluteJointConstraint{
            X0, axes, base_node_u, target_node_u, residual_terms, base_gradient_terms,
            target_gradient_terms
        }
    );

    const auto residual_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_terms);

    constexpr auto residual_terms_exact_data =
        std::array{-5900., -2385., -4162., 97314000., 62379744., 0.};
    const auto residual_terms_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(residual_terms_exact_data.data());

    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(residual_terms_mirror(i), residual_terms_exact(i), 1.e-15);
    }

    const auto base_gradient_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), base_gradient_terms);

    constexpr auto base_gradient_terms_exact_data = std::array{
        -1., 0.,  0.,  0.,         0.,         0.,         // Row 1
        0.,  -1., 0.,  0.,         0.,         0.,         // Row 2
        0.,  0.,  -1., 0.,         0.,         0.,         // Row 3
        0.,  0.,  0.,  -10930136., -15850520., 24566248.,  // Row 4
        0.,  0.,  0.,  -5585804.,  -8091272.,  12549292.,  // Row 5
        0.,  0.,  0.,  0.,         0.,         0.          // Row 6
    };
    const auto base_gradient_terms_exact = Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(
        base_gradient_terms_exact_data.data()
    );

    for (auto i : std::views::iota(0, 6)) {
        for (auto j : std::views::iota(0, 6)) {
            EXPECT_NEAR(base_gradient_terms_mirror(i, j), base_gradient_terms_exact(i, j), 1.e-15);
        }
    }

    const auto target_gradient_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), target_gradient_terms);

    for (auto i : std::views::iota(0, 6)) {
        for (auto j : std::views::iota(0, 6)) {
            EXPECT_NEAR(
                target_gradient_terms_mirror(i, j), -base_gradient_terms_exact(i, j), 1.e-15
            );
        }
    }
}

}  // namespace kynema::tests
