#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_Quu.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateQuu() {
    const auto Cuu = kynema_fmb::beams::tests::CreateView<double[6][6]>(
        "Cuu", std::array{1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12.,
                          13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24.,
                          25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36.}
    );
    const auto x0pupSS = kynema_fmb::beams::tests::CreateView<double[3][3]>(
        "x0pupSS", std::array{37., 38., 39., 40., 41., 42., 43., 44., 45.}
    );
    const auto M_tilde = kynema_fmb::beams::tests::CreateView<double[3][3]>(
        "M_tilde", std::array{46., 47., 48., 49., 50., 51., 52., 53., 54.}
    );
    const auto N_tilde = kynema_fmb::beams::tests::CreateView<double[3][3]>(
        "N_tilde", std::array{55., 56., 57., 58., 59., 60., 61., 62., 63.}
    );

    const auto Quu = Kokkos::View<double[6][6]>("Quu");

    Kokkos::parallel_for(
        "CalculateQuu", 1, KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateQuu<Kokkos::DefaultExecutionSpace>::invoke(
                Cuu, x0pupSS, N_tilde, Quu
            );
        }
    );

    constexpr auto Quu_exact_data =
        std::array{0., 0., 0., 0.,      0.,      0.,      0., 0., 0., 0.,      0.,      0.,
                   0., 0., 0., 0.,      0.,      0.,      0., 0., 0., 113262., 116130., 118998.,
                   0., 0., 0., 115986., 118923., 121860., 0., 0., 0., 118710., 121716., 124722.};
    const auto Quu_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(Quu_exact_data.data());

    const auto Quu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Quu);
    kynema_fmb::beams::tests::CompareWithExpected(Quu_mirror, Quu_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateQuuTests, OneNode) {
    TestCalculateQuu();
}

}  // namespace kynema_fmb::tests
