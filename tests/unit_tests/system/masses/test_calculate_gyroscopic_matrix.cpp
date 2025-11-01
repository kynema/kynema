#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/masses/calculate_gyroscopic_matrix.hpp"
#include "test_calculate.hpp"

namespace kynema::tests {

struct ExecuteCalculateGyroscopicMatrix {
    double mass;
    Kokkos::View<double[3]>::const_type omega;
    Kokkos::View<double[3][3]>::const_type omega_tilde;
    Kokkos::View<double[3][3]>::const_type rho;
    Kokkos::View<double[3]>::const_type eta;
    Kokkos::View<double[6][6]> Guu;

    KOKKOS_FUNCTION
    void operator()(size_t) const {
        masses::CalculateGyroscopicMatrix<Kokkos::DefaultExecutionSpace>::invoke(
            mass, omega, eta, rho, omega_tilde, Guu
        );
    }
};

TEST(CalculateGyroscopicMatrixMassesTests, OneNode) {
    const double mass = 1.;

    const auto omega = CreateView<double[3]>("omega", std::array{40., 41., 42.});
    const auto omega_tilde = CreateView<double[3][3]>(
        "omega_tilde", std::array{46., 47., 48., 49., 50., 51., 52., 53., 54.}
    );
    const auto rho =
        CreateView<double[3][3]>("rho", std::array{55., 56., 57., 58., 59., 60., 61., 62., 63.});
    const auto eta = CreateView<double[3]>("eta", std::array{64., 65., 66.});

    const auto Guu = Kokkos::View<double[6][6]>("Guu");

    Kokkos::parallel_for(
        "CalculateGyroscopicMatrix", 1,
        ExecuteCalculateGyroscopicMatrix{
            .mass = mass,
            .omega = omega,
            .omega_tilde = omega_tilde,
            .rho = rho,
            .eta = eta,
            .Guu = Guu
        }
    );

    constexpr auto Guu_exact_data =
        std::array{0., 0., 0., 18.,   10301., -9734., 0., 0., 0., -10322., -30.,   9182.,
                   0., 0., 0., 9764., -9191., 12.,    0., 0., 0., 8184.,   15953., 1207.,
                   0., 0., 0., 1078., 8856.,  15896., 0., 0., 0., 16487.,  2497.,  9546.};
    const auto Guu_exact =
        Kokkos::View<const double[6][6], Kokkos::HostSpace>(Guu_exact_data.data());

    const auto Guu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Guu);
    CompareWithExpected(Guu_mirror, Guu_exact);
}

}  // namespace kynema::tests
