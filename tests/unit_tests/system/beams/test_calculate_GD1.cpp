#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_G_D1.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateG_D1() {
    const auto Duu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto r = kynema_fmb::beams::tests::CreateView<double[4]>("r", std::array{40., 41., 42., 43.});
    const auto xr_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("xr_prime", std::array{62., 63., 64.});
    const auto kappa =
        kynema_fmb::beams::tests::CreateView<double[3]>("kappa", std::array{53., 54., 55.});
    const auto G_D1 = Kokkos::View<double[6][6]>("G_D1");

    Kokkos::parallel_for(
        "CalculateG_D1", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateG_D1<Kokkos::DefaultExecutionSpace>::invoke(
                r, xr_prime, kappa, Duu, G_D1
            );
        }
    );

    constexpr auto G_D1_exact_data = std::array{
        0.0, 0.0, 0.0, -423257.0, 816274.0, -403097.0,  //
        0.0, 0.0, 0.0, -380387.0, 670054.0, -299747.0,  //
        0.0, 0.0, 0.0, -337517.0, 523834.0, -196397.0,  //
        0.0, 0.0, 0.0, -294647.0, 377614.0, -93047.0,   //
        0.0, 0.0, 0.0, -251777.0, 231394.0, 10303.0,    //
        0.0, 0.0, 0.0, -208907.0, 85174.0,  113653.0,   //
    };
    const auto G_D1_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(G_D1_exact_data.data());

    const auto G_D1_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), G_D1);
    kynema_fmb::beams::tests::CompareWithExpected(G_D1_mirror, G_D1_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateG_D1Tests, OneNode) {
    TestCalculateG_D1();
}

}  // namespace kynema_fmb::tests
