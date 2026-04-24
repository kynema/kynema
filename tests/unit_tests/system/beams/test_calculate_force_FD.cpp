#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_force_FD.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateForceFD() {
    const auto x0pupSS = kynema_fmb::beams::tests::CreateView<double[3][3]>(
        "x0pupSS", std::array{1., 2., 3., 4., 5., 6., 7., 8., 9.}
    );
    const auto FC = kynema_fmb::beams::tests::CreateView<double[6]>(
        "FC", std::array{10., 11., 12., 13., 14., 15.}
    );

    const auto FD = Kokkos::View<double[6]>("FD");

    Kokkos::parallel_for(
        "CalculateForceFD", 1, KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateForceFD<Kokkos::DefaultExecutionSpace>::invoke(
                x0pupSS, FC, FD
            );
        }
    );

    constexpr auto FD_exact_data = std::array{0., 0., 0., 138., 171., 204.};
    const auto FD_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(FD_exact_data.data());

    const auto FD_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FD);
    kynema_fmb::beams::tests::CompareWithExpected(FD_mirror, FD_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateForceFDTests, OneNode) {
    TestCalculateForceFD();
}

}  // namespace kynema_fmb::tests
