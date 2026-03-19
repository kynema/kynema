#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_D_D1.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateD_D1() {
    const auto Duu = kynema::beams::tests::CreateView<double[6][6]>(
        "Duu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto omega =
        kynema::beams::tests::CreateView<double[3]>("omega", std::array{47., 48., 49.});

    const auto D_D1 = Kokkos::View<double[6][6]>("D_D1");

    Kokkos::parallel_for(
        "CalculateD_D1", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema::beams::CalculateD_D1<Kokkos::DefaultExecutionSpace>::invoke(omega, Duu, D_D1);
        }
    );

    constexpr auto D_D1_exact_data = std::array{
        0.0, 0.0, 0.0, -43.0, 86.0, -43.0,  //
        0.0, 0.0, 0.0, -37.0, 74.0, -37.0,  //
        0.0, 0.0, 0.0, -31.0, 62.0, -31.0,  //
        0.0, 0.0, 0.0, -25.0, 50.0, -25.0,  //
        0.0, 0.0, 0.0, -19.0, 38.0, -19.0,  //
        0.0, 0.0, 0.0, -13.0, 26.0, -13.0,  //
    };
    const auto D_D1_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(D_D1_exact_data.data());

    const auto D_D1_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D_D1);
    kynema::beams::tests::CompareWithExpected(D_D1_mirror, D_D1_exact);
}

}  // namespace

namespace kynema::tests {

TEST(CalculateD_D1Tests, OneNode) {
    TestCalculateD_D1();
}

}  // namespace kynema::tests
