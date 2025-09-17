#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "constraints/calculate_fixed_bc_constraint.hpp"
#include "create_view.hpp"

namespace kynema::tests {

struct ExecuteCalculateFixedBCConstraint {
    Kokkos::View<double[3]>::const_type X0;
    Kokkos::View<double[7]>::const_type node_u;
    Kokkos::View<double[6]> residual_terms;
    Kokkos::View<double[6][6]> target_gradient_terms;

    KOKKOS_FUNCTION
    void operator()(int) const {
        constraints::CalculateFixedBCConstraint<Kokkos::DefaultExecutionSpace>::invoke(
            X0, node_u, residual_terms, target_gradient_terms
        );
    }
};

TEST(CalculateFixedBCConstraintTests, OneConstraint) {
    const auto X0 = CreateView<double[3]>("X0", std::array{1., 2., 3.});
    const auto node_u =
        CreateView<double[7]>("node_u", std::array{11., 12., 13., 14., 15., 16., 17.});

    const auto residual_terms = Kokkos::View<double[6]>("residual_terms");
    const auto target_gradient_terms = Kokkos::View<double[6][6]>("target_gradient_terms");

    Kokkos::parallel_for(
        "CalculateFixedBCConstraint", 1,
        ExecuteCalculateFixedBCConstraint{X0, node_u, residual_terms, target_gradient_terms}
    );

    const auto residual_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), residual_terms);

    constexpr auto residual_terms_exact_data = std::array{11., 12., 13., 420., 448., 476.};
    const auto residual_terms_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(residual_terms_exact_data.data());

    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(residual_terms_mirror(i), residual_terms_exact(i), 1.e-15);
    }

    const auto target_gradient_terms_mirror =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), target_gradient_terms);

    constexpr auto target_gradient_terms_exact_data = std::array{
        1., 0., 0., 0.,    0.,    0.,     // Row 1
        0., 1., 0., 0.,    0.,    0.,     // Row 2
        0., 0., 1., 0.,    0.,    0.,     // Row 3
        0., 0., 0., -29.,  -478., -31.,   // Row 4
        0., 0., 0., -2.,   -60.,  -482.,  // Row 5
        0., 0., 0., -479., -62.,  -93.    // Row 6
    };
    const auto target_gradient_terms_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(
            target_gradient_terms_exact_data.data()
        );

    for (auto i : std::views::iota(0, 6)) {
        for (auto j : std::views::iota(0, 6)) {
            EXPECT_NEAR(
                target_gradient_terms_mirror(i, j), target_gradient_terms_exact(i, j), 1.e-15
            );
        }
    }
}

}  // namespace kynema::tests
