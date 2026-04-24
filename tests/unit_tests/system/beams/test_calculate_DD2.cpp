#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_D_D2.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateD_D2() {
    const auto Duu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto xr_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("xr_prime", std::array{62., 63., 64.});
    const auto u_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("u_prime", std::array{65., 66., 67.});

    const auto D_D2 = Kokkos::View<double[6][6]>("D_D2");

    Kokkos::parallel_for(
        "CalculateD_D2", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateD_D2<Kokkos::DefaultExecutionSpace>::invoke(
                xr_prime, u_prime, Duu, D_D2
            );
        }
    );

    constexpr auto D_D2_exact_data = std::array{
        0.0,    0.0,    0.0,    0.0,    0.0,    0.0,     //
        0.0,    0.0,    0.0,    0.0,    0.0,    0.0,     //
        0.0,    0.0,    0.0,    0.0,    0.0,    0.0,     //
        -760.0, -758.0, -756.0, -754.0, -752.0, -750.0,  //
        1520.0, 1516.0, 1512.0, 1508.0, 1504.0, 1500.0,  //
        -760.0, -758.0, -756.0, -754.0, -752.0, -750.0,  //
    };
    const auto D_D2_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(D_D2_exact_data.data());

    const auto D_D2_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D_D2);
    kynema_fmb::beams::tests::CompareWithExpected(D_D2_mirror, D_D2_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateD_D2Tests, OneNode) {
    TestCalculateD_D2();
}

}  // namespace kynema_fmb::tests
