#include <array>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "create_view.hpp"
#include "state/update_global_position.hpp"

namespace kynema::tests {
TEST(UpdateGlobalPosition, OneNode) {
    const auto q = CreateView<double[1][7]>("q", std::array{1., 2., 3., 4., 5., 6., 7.});
    const auto x0 = CreateView<double[1][7]>("x0", std::array{8., 9., 10., 11., 12., 13., 14.});

    const auto x = Kokkos::View<double[1][7]>("x");

    Kokkos::parallel_for(
        "UpdateGlobalPosition", 1,
        state::UpdateGlobalPosition<Kokkos::DefaultExecutionSpace>{.q = q, .x0 = x0, .x = x}
    );

    constexpr auto x_exact_data = std::array{9., 11., 13., -192., 96., 132., 126.};
    const auto x_exact =
        Kokkos::View<double[1][7], Kokkos::HostSpace>::const_type(x_exact_data.data());
    const auto x_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);
    for (auto i : std::views::iota(0, 7)) {
        EXPECT_NEAR(x_mirror(0, i), x_exact(0, i), 1.e-14);
    }
}

}  // namespace kynema::tests
