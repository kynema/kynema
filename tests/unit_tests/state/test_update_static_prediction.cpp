#include <array>
#include <cstddef>
#include <ranges>
#include <string>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "create_view.hpp"
#include "dof_management/freedom_signature.hpp"
#include "state/update_static_prediction.hpp"

namespace kynema::tests {
TEST(UpdateStaticPrediction, TwoNodes) {
    constexpr auto h = 2.;

    const auto node_freedom_allocation_table = CreateView<dof::FreedomSignature[2]>(
        "nfat",
        std::array{dof::FreedomSignature::AllComponents, dof::FreedomSignature::AllComponents}
    );
    const auto node_freedom_map_table = CreateView<size_t[2]>("nfmt", std::array{0UL, 6UL});
    const auto x_delta = CreateLeftView<double[12][1]>(
        "x_delta", std::array{2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 22., 24.}
    );

    const auto q_delta = Kokkos::View<double[2][6]>("q_delta");

    Kokkos::parallel_for(
        "UpdateStaticPrediction", 2,
        state::UpdateStaticPrediction<Kokkos::DefaultExecutionSpace>{
            .h = h,
            .node_freedom_allocation_table = node_freedom_allocation_table,
            .node_freedom_map_table = node_freedom_map_table,
            .x_delta = x_delta,
            .q_delta = q_delta
        }
    );

    constexpr auto q_delta_exact_data =
        std::array{1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.};
    const auto q_delta_exact =
        Kokkos::View<double[2][6], Kokkos::HostSpace>::const_type(q_delta_exact_data.data());
    const auto q_delta_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), q_delta);
    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(q_delta_mirror(0, i), q_delta_exact(0, i), 1.e-14);
    }
    for (auto i : std::views::iota(0, 6)) {
        EXPECT_NEAR(q_delta_mirror(1, i), q_delta_exact(1, i), 1.e-14);
    }
}

}  // namespace kynema::tests
