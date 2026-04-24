#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_force_FC.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateForceFC() {
    const auto Cuu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Cuu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );

    const auto strain = kynema_fmb::beams::tests::CreateView<double[6]>(
        "strain", std::array{37., 38., 39., 40., 41., 42.}
    );

    const auto FC = Kokkos::View<double[6]>("FC");

    Kokkos::parallel_for(
        "CalculateForceFC", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateForceFC<Kokkos::DefaultExecutionSpace>::invoke(Cuu, strain, FC);
        }
    );

    constexpr auto FC_exact_data = std::array{847., 2269., 3691., 5113., 6535., 7957.};
    const auto FC_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(FC_exact_data.data());

    const auto FC_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FC);
    kynema_fmb::beams::tests::CompareWithExpected(FC_mirror, FC_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateForceFCTests, OneNode) {
    TestCalculateForceFC();
}

}  // namespace kynema_fmb::tests
