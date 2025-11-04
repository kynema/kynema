#include <array>
#include <cstddef>
#include <string>

#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>
#include <gtest/gtest.h>

#include "system/masses/calculate_inertial_force.hpp"
#include "test_calculate.hpp"

namespace kynema::tests {

struct ExecuteCalculateInertialForces {
    double mass;
    Kokkos::View<double[3]>::const_type u_ddot;
    Kokkos::View<double[3]>::const_type omega;
    Kokkos::View<double[3]>::const_type omega_dot;
    Kokkos::View<double[3][3]>::const_type eta_tilde;
    Kokkos::View<double[3][3]>::const_type omega_tilde;
    Kokkos::View<double[3][3]>::const_type omega_dot_tilde;
    Kokkos::View<double[3][3]>::const_type rho;
    Kokkos::View<double[3]>::const_type eta;
    Kokkos::View<double[6]> FI;

    KOKKOS_FUNCTION
    void operator()(size_t) const {
        masses::CalculateInertialForce<Kokkos::DefaultExecutionSpace>::invoke(
            mass, u_ddot, omega, omega_dot, eta, eta_tilde, rho, omega_tilde, omega_dot_tilde, FI
        );
    }
};

TEST(CalculateInertialForcesTestsMasses, OneNode) {
    constexpr auto mass = 1.;

    const auto u_ddot = CreateView<double[3]>("u_ddot", std::array{37., 38., 39.});
    const auto omega = CreateView<double[3]>("omega", std::array{40., 41., 42.});
    const auto omega_dot = CreateView<double[3]>("omega_dot", std::array{43., 44., 45.});
    const auto eta_tilde = CreateView<double[3][3]>(
        "eta_tilde", std::array{46., 47., 48., 49., 50., 51., 52., 53., 54.}
    );
    const auto rho =
        CreateView<double[3][3]>("rho", std::array{55., 56., 57., 58., 59., 60., 61., 62., 63.});
    const auto eta = CreateView<double[3]>("eta", std::array{64., 65., 66.});
    const auto omega_tilde = CreateView<double[3][3]>(
        "omega_tilde", std::array{0., -42., 41., 42., 0., -40., -41., 40., 0.}
    );
    const auto omega_dot_tilde = CreateView<double[3][3]>(
        "omega_dot_tilde", std::array{0., -45., 44., 45., 0., -43., -44., 43., 0.}
    );

    const auto FI = Kokkos::View<double[6]>("FI");

    Kokkos::parallel_for(
        "CalculateInertialForces", 1,
        ExecuteCalculateInertialForces{
            .mass = mass,
            .u_ddot = u_ddot,
            .omega = omega,
            .omega_dot = omega_dot,
            .eta_tilde = eta_tilde,
            .omega_tilde = omega_tilde,
            .omega_dot_tilde = omega_dot_tilde,
            .rho = rho,
            .eta = eta,
            .FI = FI
        }
    );

    constexpr auto FI_exact_data = std::array{-2984., 32., 2922., 20624., -2248., 22100.};
    const auto FI_exact =
        Kokkos::View<double[6], Kokkos::HostSpace>::const_type(FI_exact_data.data());

    const auto FI_mirror = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), FI);
    CompareWithExpected(FI_mirror, FI_exact);
}

}  // namespace kynema::tests
