#include <array>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "create_view.hpp"
#include "state/calculate_displacement.hpp"

namespace kynema::tests {
TEST(CalculateDisplacement, OneNode) {
    constexpr auto h = 2.;
    const auto q_delta = CreateView<double[1][6]>("q_delta", std::array{1., 2., 3., 4., 5., 6.});
    const auto q_prev =
        CreateView<double[1][7]>("q_prev", std::array{7., 8., 9., 10., 11., 12., 13.});

    const auto q = Kokkos::View<double[1][7]>("q");

    Kokkos::parallel_for(
        "CalculateDisplacement", 1,
        state::CalculateDisplacement<Kokkos::DefaultExecutionSpace>{
            .h = h, .q_delta = q_delta, .q_prev = q_prev, .q = q
        }
    );

    constexpr auto q_exact_data = std::array{9.,
                                             12.,
                                             15.,
                                             -20.510952969318009,
                                             -6.4827969278626361,
                                             -5.141528734877415,
                                             -6.696180594260384};
    const auto q_exact =
        Kokkos::View<double[1][7], Kokkos::HostSpace>::const_type(q_exact_data.data());
    const auto q_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), q);
    for (auto i : std::views::iota(0, 7)) {
        EXPECT_NEAR(q_mirror(0, i), q_exact(0, i), 1.e-14);
    }
}
}  // namespace kynema::tests
