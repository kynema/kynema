#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_force_FD2.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateForceFD2() {
    const auto Duu = kynema::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );

    const auto xr_prime =
        kynema::beams::tests::CreateView<double[3]>("xr_prime", std::array{37., 38., 39.});

    const auto u_prime =
        kynema::beams::tests::CreateView<double[3]>("u_prime", std::array{40., 41., 42.});

    const auto strain_dot = kynema::beams::tests::CreateView<double[6]>(
        "strain_dot", std::array{43., 44., 45., 46., 47., 48.}
    );

    const auto FD2 = Kokkos::View<double[6]>("FD2");

    Kokkos::parallel_for(
        "CalculateForceFD2", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateForceFD2<Kokkos::DefaultExecutionSpace>::invoke(
                Duu, xr_prime, u_prime, strain_dot, FD2
            );
        }
    );

    constexpr auto FD2_exact_data = std::array{0., 0., 0., -124180., 248360., -124180.};
    const auto FD2_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(FD2_exact_data.data());

    const auto FD2_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FD2);
    kynema::beams::tests::CompareWithExpected(FD2_mirror, FD2_exact);
}

}  // namespace

namespace kynema::tests {

TEST(CalculateForceFD2Tests, OneNode) {
    TestCalculateForceFD2();
}

}  // namespace kynema::tests
