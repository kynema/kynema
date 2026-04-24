#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/test_calculate.hpp"
#include "system/springs/calculate_stiffness_matrix.hpp"

namespace {

void TestCalculateStiffnessMatrixTests_OneElement() {
    auto a = Kokkos::View<double[3][3]>("a");

    Kokkos::parallel_for(
        "CalculateStiffnessMatrix", 1, KOKKOS_LAMBDA(size_t) {
            constexpr auto c1 = 2.;
            constexpr auto c2 = 1.;
            constexpr auto l = 1.;

            constexpr auto r_data = Kokkos::Array<double, 3>{1., 2., 3.};
            const auto r = Kokkos::View<double[3]>::const_type(r_data.data());

            kynema_fmb::springs::CalculateStiffnessMatrix<Kokkos::DefaultExecutionSpace>::invoke(
                c1, c2, r, l, a
            );
        }
    );

    // Expected a matrix: diag_term * I - c2 * r_tilde * r_tilde
    // where diag_term = c1 - c2 * l^2 = 2. - 1. * 1. = 1.
    constexpr auto a_exact_data = std::array{
        14., -2., -3.,  // row 1
        -2., 11., -6.,  // row 2
        -3., -6., 6.    // row 3
    };

    const auto a_exact = Kokkos::View<const double[3][3], Kokkos::HostSpace>(a_exact_data.data());

    const auto a_result = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
    Kokkos::deep_copy(a_result, a);

    kynema_fmb::beams::tests::CompareWithExpected(a_result, a_exact);
}

}  // namespace
namespace kynema_fmb::tests {

TEST(CalculateStiffnessMatrixTests, OneElement) {
    TestCalculateStiffnessMatrixTests_OneElement();
}

}  // namespace kynema_fmb::tests
