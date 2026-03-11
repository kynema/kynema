#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "interpolate_to_quadrature_point_for_inertia.hpp"
#include "system/masses/calculate_gravity_force.hpp"
#include "system/masses/calculate_gyroscopic_matrix.hpp"
#include "system/masses/calculate_inertia_stiffness_matrix.hpp"
#include "system/masses/calculate_inertial_force.hpp"
#include "system/masses/calculate_mass_matrix_components.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateQuadraturePointInertialValues {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    size_t element;

    ConstLeftView<double**> shape_interp;
    ConstView<double[3]> gravity;
    ConstView<double** [4]> qp_r0;
    ConstView<double** [6][6]> qp_Mstar;
    ConstView<double* [7]> node_u;
    ConstView<double* [6]> node_u_dot;
    ConstView<double* [6]> node_u_ddot;

    View<double* [6]> qp_Fi;
    View<double* [6]> qp_Fg;
    View<double* [6][6]> qp_Muu;
    View<double* [6][6]> qp_Guu;
    View<double* [6][6]> qp_Kuu;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using CopyVector = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 1>;

        const auto r0_data = Array<double, 4>{
            qp_r0(element, qp, 0), qp_r0(element, qp, 1), qp_r0(element, qp, 2),
            qp_r0(element, qp, 3)
        };
        auto r_data = Array<double, 4>{};
        auto xr_data = Array<double, 4>{};
        auto u_ddot_data = Array<double, 3>{};
        auto omega_data = Array<double, 3>{};
        auto omega_dot_data = Array<double, 3>{};
        auto Mstar_data = Array<double, 36>{};

        auto eta_data = Array<double, 3>{};
        auto eta_tilde_data = Array<double, 9>{};
        auto rho_data = Array<double, 9>{};
        auto omega_tilde_data = Array<double, 9>{};
        auto omega_dot_tilde_data = Array<double, 9>{};
        auto FI_data = Array<double, 6>{};
        auto FG_data = Array<double, 6>{};
        auto Muu_data = Array<double, 36>{};
        auto Guu_data = Array<double, 36>{};
        auto Kuu_data = Array<double, 36>{};

        const auto r0 = ConstView<double[4]>(r0_data.data());
        const auto r = View<double[4]>(r_data.data());
        const auto xr = View<double[4]>(xr_data.data());
        const auto u_ddot = View<double[3]>(u_ddot_data.data());
        const auto omega = View<double[3]>(omega_data.data());
        const auto omega_dot = View<double[3]>(omega_dot_data.data());

        const auto eta = View<double[3]>(eta_data.data());
        const auto eta_tilde = View<double[3][3]>(eta_tilde_data.data());
        const auto rho = View<double[3][3]>(rho_data.data());
        const auto omega_tilde = View<double[3][3]>(omega_tilde_data.data());
        const auto omega_dot_tilde = View<double[3][3]>(omega_dot_tilde_data.data());
        const auto FI = View<double[6]>(FI_data.data());
        const auto FG = View<double[6]>(FG_data.data());
        const auto Mstar = View<double[6][6]>(Mstar_data.data());
        const auto Muu = View<double[6][6]>(Muu_data.data());
        const auto Guu = View<double[6][6]>(Guu_data.data());
        const auto Kuu = View<double[6][6]>(Kuu_data.data());

        CopyMatrix::invoke(subview(qp_Mstar, element, qp, ALL, ALL), Mstar);
        beams::InterpolateToQuadraturePointForInertia<DeviceType>::invoke(
            subview(shape_interp, ALL, qp), node_u, node_u_dot, node_u_ddot, r, u_ddot, omega,
            omega_dot
        );

        math::QuaternionCompose(r, r0, xr);
        masses::RotateSectionMatrix<DeviceType>::invoke(xr, Mstar, Muu);

        const auto mass = Muu(0, 0);
        masses::CalculateEta<DeviceType>(Muu, eta);
        math::VecTilde(eta, eta_tilde);
        masses::CalculateRho<DeviceType>(Muu, rho);

        math::VecTilde(omega, omega_tilde);
        math::VecTilde(omega_dot, omega_dot_tilde);

        masses::CalculateInertialForce<DeviceType>::invoke(
            mass, u_ddot, omega, omega_dot, eta, eta_tilde, rho, omega_tilde, omega_dot_tilde, FI
        );
        masses::CalculateGravityForce<DeviceType>::invoke(mass, gravity, eta_tilde, FG);

        masses::CalculateGyroscopicMatrix<DeviceType>::invoke(
            mass, omega, eta, rho, omega_tilde, Guu
        );
        masses::CalculateInertiaStiffnessMatrix<DeviceType>::invoke(
            mass, u_ddot, omega, omega_dot, eta, rho, omega_tilde, omega_dot_tilde, Kuu
        );

        CopyVector::invoke(FI, subview(qp_Fi, qp, ALL));
        CopyVector::invoke(FG, subview(qp_Fg, qp, ALL));
        CopyMatrix::invoke(Muu, subview(qp_Muu, qp, ALL, ALL));
        CopyMatrix::invoke(Guu, subview(qp_Guu, qp, ALL, ALL));
        CopyMatrix::invoke(Kuu, subview(qp_Kuu, qp, ALL, ALL));
    }
};

}  // namespace kynema::beams
