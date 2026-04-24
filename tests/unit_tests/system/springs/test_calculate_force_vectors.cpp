#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/test_calculate.hpp"
#include "system/springs/calculate_force_vectors.hpp"

namespace {

void TestCalculateForceVectors_ThreeElements_1() {
    const auto f0 = Kokkos::View<double[3]>("f0");

    Kokkos::parallel_for(
        "CalculateForceVectors", 1,
        KOKKOS_LAMBDA(const size_t) {
            constexpr auto r0_data = Kokkos::Array<double, 3>{1., 2., 3.};
            const auto r0 = Kokkos::View<double[3]>::const_type(r0_data.data());
            constexpr auto c10 = 2.;

            kynema_fmb::springs::CalculateForceVectors<Kokkos::DefaultExecutionSpace>::invoke(
                r0, c10, f0
            );
        }
    );

    constexpr auto f0_exact_data = std::array{2., 4., 6.};
    const auto f0_exact =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(f0_exact_data.data());

    const auto f0_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f0);
    kynema_fmb::beams::tests::CompareWithExpected(f0_result, f0_exact);
}

void TestCalculateForceVectors_ThreeElements_2() {
    const auto f1 = Kokkos::View<double[3]>("f1");

    Kokkos::parallel_for(
        "CalculateForceVectors", 1,
        KOKKOS_LAMBDA(const size_t) {
            constexpr auto r1_data = Kokkos::Array<double, 3>{4., 5., 6.};
            const auto r1 = Kokkos::View<double[3]>::const_type(r1_data.data());
            constexpr auto c11 = -1.;

            kynema_fmb::springs::CalculateForceVectors<Kokkos::DefaultExecutionSpace>::invoke(
                r1, c11, f1
            );
        }
    );

    constexpr auto f1_exact_data = std::array{-4., -5., -6.};
    const auto f1_exact =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(f1_exact_data.data());

    const auto f1_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f1);
    kynema_fmb::beams::tests::CompareWithExpected(f1_result, f1_exact);
}

void TestCalculateForceVectors_ThreeElements_3() {
    const auto f2 = Kokkos::View<double[3]>("f2");

    Kokkos::parallel_for(
        "CalculateForceVectors", 1,
        KOKKOS_LAMBDA(const size_t) {
            constexpr auto r2_data = Kokkos::Array<double, 3>{7., 8., 9.};
            const auto r2 = Kokkos::View<double[3]>::const_type(r2_data.data());
            constexpr auto c12 = 0.5;

            kynema_fmb::springs::CalculateForceVectors<Kokkos::DefaultExecutionSpace>::invoke(
                r2, c12, f2
            );
        }
    );

    constexpr auto f2_exact_data = std::array{3.5, 4., 4.5};
    const auto f2_exact =
        Kokkos::View<double[3], Kokkos::HostSpace>::const_type(f2_exact_data.data());

    const auto f2_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f2);
    kynema_fmb::beams::tests::CompareWithExpected(f2_result, f2_exact);
}
}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateForceVectorsTests, ThreeElements_1) {
    TestCalculateForceVectors_ThreeElements_1();
}

TEST(CalculateForceVectorsTests, ThreeElements_2) {
    TestCalculateForceVectors_ThreeElements_2();
}

TEST(CalculateForceVectorsTests, ThreeElements_3) {
    TestCalculateForceVectors_ThreeElements_3();
}
}  // namespace kynema_fmb::tests
