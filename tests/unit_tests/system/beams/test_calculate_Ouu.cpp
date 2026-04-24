#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/beams/calculate_Ouu.hpp"
#include "test_calculate.hpp"

namespace {

void TestCalculateOuu() {
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

    const auto Ouu = Kokkos::View<double[6][6]>("Ouu");

    Kokkos::parallel_for(
        "CalculateOuu", 1,
        KOKKOS_LAMBDA(size_t) {
            kynema_fmb::beams::CalculateOuu<Kokkos::DefaultExecutionSpace>::invoke(
                Cuu, x0pupSS, M_tilde, N_tilde, Ouu
            );
        }
    );

    constexpr auto Ouu_exact_data =
        std::array{0., 0., 0., 191.,  196.,  201.,  0., 0., 0., 908.,  931.,  954.,
                   0., 0., 0., 1625., 1666., 1707., 0., 0., 0., 2360., 2419., 2478.,
                   0., 0., 0., 3077., 3154., 3231., 0., 0., 0., 3794., 3889., 3984.};
    const auto Ouu_exact =
        Kokkos::View<double[6][6], Kokkos::HostSpace>::const_type(Ouu_exact_data.data());

    const auto Ouu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Ouu);
    kynema_fmb::beams::tests::CompareWithExpected(Ouu_mirror, Ouu_exact);
}

}  // namespace

namespace kynema_fmb::tests {

TEST(CalculateOuuTests, OneNode) {
    TestCalculateOuu();
}

}  // namespace kynema_fmb::tests
