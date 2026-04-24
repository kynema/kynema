#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_force_FD1.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateForceFD1() {
    const auto Duu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );

    const auto strain_dot = kynema_fmb::beams::tests::CreateView<double[6]>(
        "strain_dot", std::array{37., 38., 39., 40., 41., 42.}
    );

    const auto FD1 = Kokkos::View<double[6]>("FD1");

    Kokkos::parallel_for(
        "CalculateForceFD1", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateForceFD1<Kokkos::DefaultExecutionSpace>::invoke(
                Duu, strain_dot, FD1
            );
        }
    );

    constexpr auto FD1_exact_data = std::array{847., 2269., 3691., 5113., 6535., 7957.};
    const auto FD1_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(FD1_exact_data.data());

    const auto FD1_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FD1);
    kynema_fmb::beams::tests::CompareWithExpected(FD1_mirror, FD1_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateForceFD1Tests, OneNode) {
    TestCalculateForceFD1();
}

}  // namespace kynema_fmb::tests
