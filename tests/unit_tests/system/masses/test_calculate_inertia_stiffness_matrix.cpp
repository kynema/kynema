#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/masses/calculate_inertia_stiffness_matrix.hpp"
#include "test_calculate.hpp"

namespace kynema::tests {

struct ExecuteCalculateInertiaStiffnessMatrix {
    double mass;
    Kokkos::View<double[3]>::const_type u_ddot;
    Kokkos::View<double[3]>::const_type omega;
    Kokkos::View<double[3]>::const_type omega_dot;
    Kokkos::View<double[3][3]>::const_type omega_tilde;
    Kokkos::View<double[3][3]>::const_type omega_dot_tilde;
    Kokkos::View<double[3][3]>::const_type rho;
    Kokkos::View<double[3]>::const_type eta;
    Kokkos::View<double[6][6]> Kuu;

    KOKKOS_FUNCTION
    void operator()(size_t) const {
        masses::CalculateInertiaStiffnessMatrix<Kokkos::DefaultExecutionSpace>::invoke(
            mass, u_ddot, omega, omega_dot, eta, rho, omega_tilde, omega_dot_tilde, Kuu
        );
    }
};

TEST(CalculateInertiaStiffnessMatrixMassesTests, OneNode) {
    const double mass = 1.;

    const auto u_ddot = CreateView<double[3]>("u_ddot", std::array{37., 38., 39.});
    const auto omega = CreateView<double[3]>("omega", std::array{40., 41., 42.});
    const auto omega_dot = CreateView<double[3]>("omega_dot", std::array{43., 44., 45.});
    const auto eta_tilde = CreateView<double[3][3]>(
        "eta_tilde", std::array{46., 47., 48., 49., 50., 51., 52., 53., 54.}
    );
    const auto omega_tilde = CreateView<double[3][3]>(
        "omega_tilde", std::array{46., 47., 48., 49., 50., 51., 52., 53., 54.}
    );
    const auto omega_dot_tilde = CreateView<double[3][3]>(
        "omega_dot_tilde", std::array{55., 56., 57., 58., 59., 60., 61., 62., 63.}
    );
    const auto rho =
        CreateView<double[3][3]>("rho", std::array{64., 65., 66., 67., 68., 69., 70., 71., 72.});
    const auto eta = CreateView<double[3]>("eta", std::array{73., 74., 75.});

    const auto Kuu = Kokkos::View<double[6][6]>("Kuu");

    Kokkos::parallel_for(
        "CalculateInertiaStiffnessMatrix", 1,
        ExecuteCalculateInertiaStiffnessMatrix{
            .mass = mass,
            .u_ddot = u_ddot,
            .omega = omega,
            .omega_dot = omega_dot,
            .omega_tilde = omega_tilde,
            .omega_dot_tilde = omega_dot_tilde,
            .rho = rho,
            .eta = eta,
            .Kuu = Kuu
        }
    );

    constexpr auto Kuu_exact_data = std::array<double, 36>{
        0., 0., 0., 3396.,    -6792.,   3396.,    0., 0., 0., 3609.,    -7218.,   3609.,
        0., 0., 0., 3822.,    -7644.,   3822.,    0., 0., 0., 1407766., 1481559., 1465326.,
        0., 0., 0., 1496300., 1558384., 1576048., 0., 0., 0., 1604122., 1652877., 1652190.
    };
    const auto Kuu_exact =
        Kokkos::View<const double[6][6], Kokkos::HostSpace>(Kuu_exact_data.data());

    const auto Kuu_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kuu);
    CompareWithExpected(Kuu_mirror, Kuu_exact);
}

}  // namespace kynema::tests
