#include <array>
#include <ranges>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

#include "create_view.hpp"
#include "state/calculate_next_state.hpp"

namespace {

inline void CompareWithExpected(
    const Kokkos::View<const double**>::host_mirror_type& result,
    const Kokkos::View<const double**, Kokkos::HostSpace>& expected
) {
    for (auto i : std::views::iota(0U, result.extent(0))) {
        for (auto j : std::views::iota(0U, result.extent(1))) {
            EXPECT_NEAR(result(i, j), expected(i, j), 1.e-14);
        }
    }
}

}  // namespace
namespace kynema::tests {

TEST(CalculateNextState, OneNode) {
    constexpr auto h = 2.;
    constexpr auto alpha_f = 3.;
    constexpr auto alpha_m = 4.;
    constexpr auto beta = 5.;
    constexpr auto gamma = 6.;

    const auto q_delta =
        CreateMutableView<double[1][6]>("q_delta", std::array{1., 2., 3., 4., 5., 6.});
    const auto v = CreateMutableView<double[1][6]>("v", std::array{1., 2., 3., 4., 5., 6.});
    const auto vd = CreateMutableView<double[1][6]>("vd", std::array{7., 8., 9., 10., 11., 12.});
    const auto a = CreateMutableView<double[1][6]>("a", std::array{13., 14., 15., 16., 17., 18.});

    Kokkos::parallel_for(
        "CalculateNextState", 1,
        state::CalculateNextState<Kokkos::DefaultExecutionSpace>{
            .h = h,
            .alpha_f = alpha_f,
            .alpha_m = alpha_m,
            .beta = beta,
            .gamma = gamma,
            .q_delta = q_delta,
            .v = v,
            .vd = vd,
            .a = a
        }
    );

    constexpr auto q_delta_exact_data = std::array{-12.666666666666661, -17.333333333333339, -22.,
                                                   -26.666666666666661, -31.333333333333339, -36.};
    const auto q_delta_exact =
        Kokkos::View<double[1][6], Kokkos::HostSpace>::const_type(q_delta_exact_data.data());
    const auto q_delta_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), q_delta);
    CompareWithExpected(q_delta_mirror, q_delta_exact);

    constexpr auto v_exact_data = std::array{-4.9999999999999929, -10.000000000000007, -15.,
                                             -19.999999999999993, -25.000000000000007, -30.};
    const auto v_exact =
        Kokkos::View<double[1][6], Kokkos::HostSpace>::const_type(v_exact_data.data());
    const auto v_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
    CompareWithExpected(v_mirror, v_exact);

    constexpr auto vd_exact_data = std::array{0., 0., 0., 0., 0., 0.};
    const auto vd_exact =
        Kokkos::View<double[1][6], Kokkos::HostSpace>::const_type(vd_exact_data.data());
    const auto vd_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vd);
    CompareWithExpected(vd_mirror, vd_exact);

    constexpr auto a_exact_data = std::array{10.333333333333334, 10.666666666666666, 11.,
                                             11.333333333333334, 11.666666666666666, 12.};
    const auto a_exact =
        Kokkos::View<double[1][6], Kokkos::HostSpace>::const_type(a_exact_data.data());
    const auto a_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
    CompareWithExpected(a_mirror, a_exact);
}
}  // namespace kynema::tests
