#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_G_D2.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateG_D2() {
    const auto Duu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto r = kynema_fmb::beams::tests::CreateView<double[4]>("r", std::array{40., 41., 42., 43.});
    const auto xr_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("xr_prime", std::array{62., 63., 64.});
    const auto u_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("u_prime", std::array{65., 66., 67.});
    const auto omega =
        kynema_fmb::beams::tests::CreateView<double[3]>("omega", std::array{47., 48., 49.});
    const auto kappa =
        kynema_fmb::beams::tests::CreateView<double[3]>("kappa", std::array{53., 54., 55.});
    const auto eps_dot =
        kynema_fmb::beams::tests::CreateView<double[3]>("eps_dot", std::array{56., 57., 58.});
    const auto kappa_dot =
        kynema_fmb::beams::tests::CreateView<double[3]>("kappa_dot", std::array{59., 60., 61.});

    const auto G_D2 = Kokkos::View<double[6][6]>("G_D2");

    Kokkos::parallel_for(
        "CalculateG_D2", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateG_D2<Kokkos::DefaultExecutionSpace>::invoke(
                r, xr_prime, u_prime, kappa, Duu, G_D2
            );
        }
    );

    constexpr auto G_D2_exact_data = std::array{
        0.0, 0.0, 0.0, 0.0,        0.0,         0.0,          //
        0.0, 0.0, 0.0, 0.0,        0.0,         0.0,          //
        0.0, 0.0, 0.0, 0.0,        0.0,         0.0,          //
        0.0, 0.0, 0.0, -6291004.0, 20202488.0,  -13931644.0,  //
        0.0, 0.0, 0.0, 12582008.0, -40404976.0, 27863288.0,   //
        0.0, 0.0, 0.0, -6291004.0, 20202488.0,  -13931644.0,  //
    };
    const auto G_D2_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(G_D2_exact_data.data());

    const auto G_D2_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), G_D2);
    kynema_fmb::beams::tests::CompareWithExpected(G_D2_mirror, G_D2_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateG_D2Tests, OneNode) {
    TestCalculateG_D2();
}

}  // namespace kynema_fmb::tests
