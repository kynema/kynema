#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/springs/calculate_length.hpp"

namespace {
void TestCalclateLength_1() {
    auto l = 0.;
    Kokkos::parallel_reduce(
        "CalculateLength", 1,
        KOKKOS_LAMBDA(size_t, double& length) {
            constexpr auto r0_data = Kokkos::Array<double, 3>{1., 0., 0.};
            const auto r0 = Kokkos::View<double[3]>::const_type(r0_data.data());

            length = kynema_fmb::springs::CalculateLength<Kokkos::DefaultExecutionSpace>(r0);
        },
        l
    );

    EXPECT_NEAR(l, 1., 1.e-16);
}

void TestCalclateLength_2() {
    auto l = 0.;
    Kokkos::parallel_reduce(
        "CalculateLength", 1,
        KOKKOS_LAMBDA(size_t, double& length) {
            constexpr auto r0_data = Kokkos::Array<double, 3>{3., 4., 0.};
            const auto r0 = Kokkos::View<double[3]>::const_type(r0_data.data());

            length = kynema_fmb::springs::CalculateLength<Kokkos::DefaultExecutionSpace>(r0);
        },
        l
    );

    EXPECT_NEAR(l, 5., 1.e-16);
}

void TestCalclateLength_3() {
    auto l = 0.;
    Kokkos::parallel_reduce(
        "CalculateLength", 1,
        KOKKOS_LAMBDA(size_t, double& length) {
            constexpr auto r0_data = Kokkos::Array<double, 3>{1., 2., 2.};
            const auto r0 = Kokkos::View<double[3]>::const_type(r0_data.data());

            length = kynema_fmb::springs::CalculateLength<Kokkos::DefaultExecutionSpace>(r0);
        },
        l
    );

    EXPECT_NEAR(l, 3., 1.e-16);
}
}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateLengthTests, ThreeElements_1) {
    TestCalclateLength_1();
}

TEST(CalculateLengthTests, ThreeElements_2) {
    TestCalclateLength_2();
}

TEST(CalculateLengthTests, ThreeElements_3) {
    TestCalclateLength_3();
}

}  // namespace kynema_fmb::tests
