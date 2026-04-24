#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/test_calculate.hpp"
#include "system/springs/calculate_force_coefficients.hpp"

namespace {

void TestCalculateForceCoefficient1_ThreeElements() {
    const auto k = kynema_fmb::beams::tests::CreateView<double[3]>("k", std::array{100., 200., 300.});
    const auto l_ref = kynema_fmb::beams::tests::CreateView<double[3]>("l_ref", std::array{1., 5., 3.});
    const auto l = kynema_fmb::beams::tests::CreateView<double[3]>("l", std::array{2., 4., 6.});

    const auto c1 = Kokkos::View<double[3]>("c1");

    Kokkos::parallel_for(
        "CalculateForceCoefficients", 3,
        KOKKOS_LAMBDA(const size_t i_elem) {
            c1(i_elem) = kynema_fmb::springs::CalculateForceCoefficient1<Kokkos::DefaultExecutionSpace>(
                k(i_elem), l_ref(i_elem), l(i_elem)
            );
        }
    );

    constexpr auto c1_exact_data = std::array{-50., 50., -150.};
    const auto c1_exact =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(c1_exact_data.data());

    const auto c1_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), c1);
    kynema_fmb::beams::tests::CompareWithExpected(c1_result, c1_exact);
}

void TestCalculateForceCoefficient2_ThreeElements() {
    const auto k = kynema_fmb::beams::tests::CreateView<double[3]>("k", std::array{100., 200., 300.});
    const auto l_ref = kynema_fmb::beams::tests::CreateView<double[3]>("l_ref", std::array{1., 5., 3.});
    const auto l = kynema_fmb::beams::tests::CreateView<double[3]>("l", std::array{2., 4., 6.});

    const auto c2 = Kokkos::View<double[3]>("c2");

    Kokkos::parallel_for(
        "CalculateForceCoefficients", 3,
        KOKKOS_LAMBDA(const size_t i_elem) {
            c2(i_elem) = kynema_fmb::springs::CalculateForceCoefficient2<Kokkos::DefaultExecutionSpace>(
                k(i_elem), l_ref(i_elem), l(i_elem)
            );
        }
    );

    constexpr auto c2_exact_data = std::array{12.5, 15.625, 25. / 6.};
    const auto c2_exact =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(c2_exact_data.data());

    const auto c2_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), c2);
    kynema_fmb::beams::tests::CompareWithExpected(c2_result, c2_exact);
}
}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateForceCoefficients1Test, ThreeElements) {
    TestCalculateForceCoefficient1_ThreeElements();
}

TEST(CalculateForceCoefficients2Test, ThreeElements) {
    TestCalculateForceCoefficient2_ThreeElements();
}

}  // namespace kynema_fmb::tests
