#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_P_D2.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateP_D2() {
    const auto Duu = kynema::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto xr_prime =
        kynema::beams::tests::CreateView<double[3]>("xr_prime", std::array{62., 63., 64.});
    const auto u_prime =
        kynema::beams::tests::CreateView<double[3]>("u_prime", std::array{65., 66., 67.});
    const auto omega =
        kynema::beams::tests::CreateView<double[3]>("omega", std::array{47., 48., 49.});
    const auto eps_dot =
        kynema::beams::tests::CreateView<double[3]>("eps_dot", std::array{56., 57., 58.});
    const auto kappa_dot =
        kynema::beams::tests::CreateView<double[3]>("kappa_dot", std::array{59., 60., 61.});

    const auto P_D2 = Kokkos::View<double[6][6]>("P_D2");

    Kokkos::parallel_for(
        "CalculateP_D2", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateP_D2<Kokkos::DefaultExecutionSpace>::invoke(
                xr_prime, u_prime, omega, eps_dot, kappa_dot, Duu, P_D2
            );
        }
    );

    constexpr auto P_D2_exact_data = std::array{
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //
        -0.0, 5458.0, -3352.0, -848.0, 1696.0, -848.0, //
        -5458.0, -0.0, 1246.0, 1696.0, -3392.0, 1696.0, //
        3352.0, -1246.0, -0.0, -848.0, 1696.0, -848.0, //
    };
    const auto P_D2_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(P_D2_exact_data.data());

    const auto P_D2_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), P_D2);
    kynema::beams::tests::CompareWithExpected(P_D2_mirror, P_D2_exact);
}

}  // namespace

namespace kynema::tests {

TEST(CalculateP_D2Tests, OneNode) {
    TestCalculateP_D2();
}

}  // namespace kynema::tests
