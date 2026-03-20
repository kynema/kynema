#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_K_D1.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateK_D1() {
    const auto Duu = kynema::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto r = kynema::beams::tests::CreateView<double[4]>("r", std::array{40., 41., 42., 43.});
    const auto xr_prime =
        kynema::beams::tests::CreateView<double[3]>("xr_prime", std::array{62., 63., 64.});
    const auto omega =
        kynema::beams::tests::CreateView<double[3]>("omega", std::array{47., 48., 49.});
    const auto kappa =
        kynema::beams::tests::CreateView<double[3]>("kappa", std::array{53., 54., 55.});
    const auto eps_dot =
        kynema::beams::tests::CreateView<double[3]>("eps_dot", std::array{56., 57., 58.});
    const auto kappa_dot =
        kynema::beams::tests::CreateView<double[3]>("kappa_dot", std::array{59., 60., 61.});

    const auto K_D1 = Kokkos::View<double[6][6]>("K_D1");

    Kokkos::parallel_for(
        "CalculateK_D1", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateK_D1<Kokkos::DefaultExecutionSpace>::invoke(
                r, xr_prime, omega, kappa, eps_dot, kappa_dot, Duu, K_D1
            );
        }
    );

    constexpr auto K_D1_exact_data = std::array{
        0.0, 0.0, 0.0, 61038794.0, 1126704.0, -58800314.0,  //
        0.0, 0.0, 0.0, 53071816.0, 975002.0,  -51126612.0,  //
        0.0, 0.0, 0.0, 45119106.0, 827512.0,  -43458754.0,  //
        0.0, 0.0, 0.0, 37154234.0, 694290.0,  -35799320.0,  //
        0.0, 0.0, 0.0, 29180938.0, 536270.0,  -28112982.0,  //
        0.0, 0.0, 0.0, 21240864.0, 382462.0,  -20451442.0,  //
    };
    const auto K_D1_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(K_D1_exact_data.data());

    const auto K_D1_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), K_D1);
    kynema::beams::tests::CompareWithExpected(K_D1_mirror, K_D1_exact);
}

}  // namespace

namespace kynema::tests {

TEST(CalculateK_D1Tests, OneNode) {
    TestCalculateK_D1();
}

}  // namespace kynema::tests
