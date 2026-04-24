#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_strain.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateStrain() {
    const auto x0_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("x0_prime", std::array{1., 2., 3.});
    const auto u_prime =
        kynema_fmb::beams::tests::CreateView<double[3]>("u_prime", std::array{4., 5., 6.});
    const auto r = kynema_fmb::beams::tests::CreateView<double[4]>("r", std::array{7., 8., 9., 10.});
    const auto r_prime =
        kynema_fmb::beams::tests::CreateView<double[4]>("r_prime", std::array{11., 12., 13., 14.});

    const auto strain = Kokkos::View<double[6]>("strain");

    Kokkos::parallel_for(
        "CalculateStrain", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateStrain<Kokkos::DefaultExecutionSpace>::invoke(
                x0_prime, u_prime, r, r_prime, strain
            );
        }
    );

    constexpr auto strain_exact_data = std::array{-793., -413., -621., -16., 0., -32.};
    const auto strain_exact =
        Kokkos::View<const double[6], Kokkos::HostSpace>(strain_exact_data.data());

    const auto strain_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), strain);
    kynema_fmb::beams::tests::CompareWithExpected(strain_mirror, strain_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateStrainTests, OneNode) {
    TestCalculateStrain();
}

}  // namespace kynema_fmb::tests
