#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "constraints/calculate_prescribed_bc_constraint.hpp"
#include "create_view.hpp"

namespace kynema::tests {

struct ExecuteCalculatePrescribedBCConstraint {
    Kokkos::View<double[3]>::const_type X0;
    Kokkos::View<double[7]>::const_type inputs;
    Kokkos::View<double[7]>::const_type node_u;
    Kokkos::View<double[6]> residual_terms;
    Kokkos::View<double[6][6]> target_gradient_terms;

    KOKKOS_FUNCTION
    void operator()(int) const {
        constraints::CalculatePrescribedBCConstraint<Kokkos::DefaultExecutionSpace>::invoke(
            X0, inputs, node_u, residual_terms, target_gradient_terms
        );
    }
};

TEST(CalculatePrescribedBCConstraintTests, OneConstraint) {
    const auto X0 = CreateView<double[3]>("X0", std::array{1., 2., 3.});
    const auto constraint_inputs =
        CreateView<double[7]>("constraint_inputs", std::array{4., 5., 6., 7., 8., 9., 10.});
    const auto node_u =
        CreateView<double[7]>("node_u", std::array{11., 12., 13., 14., 15., 16., 17.});

    const auto residual_terms = Kokkos::View<double[6]>("residual_terms");
    const auto target_gradient_terms = Kokkos::View<double[6][6]>("target_gradient_terms");

    Kokkos::parallel_for(
        "CalculatePrescribedBCConstraint", 1,
        ExecuteCalculatePrescribedBCConstraint{
            X0, constraint_inputs, node_u, residual_terms, target_gradient_terms
        }
    );

    const auto residual_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_terms);

    constexpr auto residual_terms_exact_data = std::array{
        -790., -411., -620., -50.666666666666657, -7.1054273576010019e-15, -101.33333333333343
    };
    const auto residual_terms_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(residual_terms_exact_data.data());

    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(residual_terms_mirror(i), residual_terms_exact(i), 1.e-12);
    }

    const auto target_gradient_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), target_gradient_terms);

    // clang-format off
    constexpr auto target_gradient_terms_exact_data = std::array{
        1., 0., 0., 0., 0., 0.,  // Row 1
        0., 1., 0., 0., 0., 0.,  // Row 2
        0., 0., 1., 0., 0., 0.,  // Row 3
        0., 0., 0., 961.99999999999977, 50.666666666666714, -1.3333333333333379,  // Row 4
        0., 0., 0., -50.666666666666714, 962.6666666666664, 25.333333333333329,   // Row 5
        0., 0., 0., -1.3333333333333308, -25.333333333333329, 959.99999999999977  // Row 6
    };
    // clang-format on
    const auto target_gradient_terms_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(
            target_gradient_terms_exact_data.data()
        );

    for (auto i : std::views::iota(0, 6)) {
        for (auto j : std::views::iota(0, 6)) {
            EXPECT_NEAR(
                target_gradient_terms_mirror(i, j), target_gradient_terms_exact(i, j), 1.e-12
            );
        }
    }
}

}  // namespace kynema::tests
