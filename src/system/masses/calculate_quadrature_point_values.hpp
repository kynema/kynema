#pragma once

#include <Kokkos_Core.hpp>

#include "system/masses/calculate_gravity_force.hpp"
#include "system/masses/calculate_gyroscopic_matrix.hpp"
#include "system/masses/calculate_inertia_stiffness_matrix.hpp"
#include "system/masses/calculate_inertial_force.hpp"
#include "system/masses/calculate_mass_matrix_components.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema_fmb::masses {

template <typename DeviceType>
struct CalculateQuadraturePointValues {
    using TeamPolicy = typename Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    double beta_prime;
    double gamma_prime;

    View<double* [7]> Q;
    View<double* [6]> V;
    View<double* [6]> A;

    ConstView<double* [6][6]> tangent;
    ConstView<size_t*> node_state_indices;
    ConstView<double[3]> gravity;
    ConstView<double* [6][6]> qp_Mstar;
    ConstView<double* [7]> node_x0;

    View<double* [6]> residual_vector_terms;
    View<double* [6][6]> system_matrix_terms;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;

        const auto index = node_state_indices(element);

        // Allocate scratch views
        const auto u_ddot_data = Array<double, 3>{A(index, 0), A(index, 1), A(index, 2)};
        const auto omega_dot_data = Array<double, 3>{A(index, 3), A(index, 4), A(index, 5)};
        const auto omega_data = Array<double, 3>{V(index, 3), V(index, 4), V(index, 5)};
        auto Muu_data = Array<double, 36>{};
        auto eta_data = Array<double, 3>{};
        auto rho_data = Array<double, 9>{};
        auto omega_tilde_data = Array<double, 9>{};
        auto omega_dot_tilde_data = Array<double, 9>{};

        // Set up Views
        const auto u_ddot = ConstView<double[3]>(u_ddot_data.data());
        const auto omega_dot = ConstView<double[3]>(omega_dot_data.data());
        const auto omega = ConstView<double[3]>(omega_data.data());
        auto Muu = View<double[6][6]>(Muu_data.data());
        auto eta = View<double[3]>(eta_data.data());
        auto rho = View<double[3][3]>(rho_data.data());
        auto omega_tilde = View<double[3][3]>(omega_tilde_data.data());
        auto omega_dot_tilde = View<double[3][3]>(omega_dot_tilde_data.data());

        // Do the math
        {
            auto Mstar_data = Array<double, 36>{};
            const auto r0_data = Array<double, 4>{
                node_x0(element, 3), node_x0(element, 4), node_x0(element, 5), node_x0(element, 6)
            };
            const auto r_data = Array<double, 4>{Q(index, 3), Q(index, 4), Q(index, 5), Q(index, 6)};
            auto xr_data = Array<double, 4>{};

            const auto Mstar = View<double[6][6]>(Mstar_data.data());
            const auto r0 = ConstView<double[4]>(r0_data.data());
            const auto r = ConstView<double[4]>(r_data.data());
            const auto xr = View<double[4]>(xr_data.data());

            CopyMatrix::invoke(subview(qp_Mstar, element, ALL, ALL), Mstar);
            math::QuaternionCompose(r, r0, xr);
            math::VecTilde(omega, omega_tilde);
            math::VecTilde(omega_dot, omega_dot_tilde);

            RotateSectionMatrix<DeviceType>::invoke(xr, Mstar, Muu);

            CalculateEta<DeviceType>(Muu, eta);
            CalculateRho<DeviceType>(Muu, rho);
        }

        {
            auto eta_tilde_data = Array<double, 9>{};
            auto Fg_data = Array<double, 6>{};
            auto Fi_data = Array<double, 6>{};

            auto eta_tilde = View<double[3][3]>(eta_tilde_data.data());
            auto Fg = View<double[6]>(Fg_data.data());
            auto Fi = View<double[6]>(Fi_data.data());

            const auto mass = Muu(0, 0);

            math::VecTilde(eta, eta_tilde);

            CalculateGravityForce<DeviceType>::invoke(mass, gravity, eta_tilde, Fg);
            CalculateInertialForce<DeviceType>::invoke(
                mass, u_ddot, omega, omega_dot, eta, eta_tilde, rho, omega_tilde, omega_dot_tilde, Fi
            );

            for (auto component = 0U; component < 6U; ++component) {
                residual_vector_terms(element, component) = Fi(component) - Fg(component);
            }
        }

        {
            auto Guu_data = Array<double, 36>{};
            auto Kuu_data = Array<double, 36>{};
            auto T_data = Array<double, 36>{};
            auto STpI_data = Array<double, 36>{};

            auto Guu = View<double[6][6]>(Guu_data.data());
            auto Kuu = View<double[6][6]>(Kuu_data.data());
            auto T = View<double[6][6]>(T_data.data());
            auto STpI = View<double[6][6]>(STpI_data.data());

            const auto mass = Muu(0, 0);

            CalculateGyroscopicMatrix<DeviceType>::invoke(mass, omega, eta, rho, omega_tilde, Guu);
            CalculateInertiaStiffnessMatrix<DeviceType>::invoke(
                mass, u_ddot, omega, omega_dot, eta, rho, omega_tilde, omega_dot_tilde, Kuu
            );

            CopyMatrix::invoke(Kokkos::subview(tangent, index, Kokkos::ALL, Kokkos::ALL), T);

            for (auto component_1 = 0U; component_1 < 6U; ++component_1) {
                for (auto component_2 = 0U; component_2 < 6U; ++component_2) {
                    STpI(component_1, component_2) = beta_prime * Muu(component_1, component_2) +
                                                     gamma_prime * Guu(component_1, component_2);
                }
            }

            Gemm::invoke(1., Kuu, T, 1., STpI);

            CopyMatrix::invoke(
                STpI, Kokkos::subview(system_matrix_terms, element, Kokkos::ALL, Kokkos::ALL)
            );
        }
    }
};

}  // namespace kynema_fmb::masses
