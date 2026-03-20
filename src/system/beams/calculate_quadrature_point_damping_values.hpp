#pragma once

#include <ranges>

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "calculate_D_D1.hpp"
#include "calculate_D_D2.hpp"
#include "calculate_G_D1.hpp"
#include "calculate_G_D2.hpp"
#include "calculate_K_D1.hpp"
#include "calculate_K_D2.hpp"
#include "calculate_P_D2.hpp"
#include "calculate_force_FD1.hpp"
#include "calculate_force_FD2.hpp"
#include "calculate_strain.hpp"
#include "calculate_strain_dot.hpp"
#include "interpolate_to_quadrature_point_for_damping.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema::beams {
template <typename DeviceType>
struct CalculateQuadraturePointDampingValues {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;
    using Gemm = KokkosBatched::SerialGemm<
        KokkosBatched::Trans::NoTranspose, KokkosBatched::Trans::NoTranspose,
        KokkosBatched::Algo::Gemm::Default>;

    size_t element;
    ConstView<double[6]> mu;
    ConstView<double*> qp_jacobian;
    ConstLeftView<double**> shape_interp;
    ConstLeftView<double**> shape_deriv;
    ConstView<double** [4]> qp_r0;
    ConstView<double** [3]> qp_x0_prime;
    ConstView<double** [6][6]> qp_Cstar;
    ConstView<double* [7]> node_u;
    ConstView<double* [6]> node_u_dot;

    View<double* [6]> qp_FD1;
    View<double* [6]> qp_FD2;
    View<double* [6][6]> qp_Duu;
    View<double* [6][6]> qp_DD1;
    View<double* [6][6]> qp_DD2;
    View<double* [6][6]> qp_GD1;
    View<double* [6][6]> qp_GD2;
    View<double* [6][6]> qp_PD2;
    View<double* [6][6]> qp_KD1;
    View<double* [6][6]> qp_KD2;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        using Kokkos::ALL;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using CopyVector = KokkosBatched::SerialCopy<KokkosBatched::Trans::NoTranspose, 1>;

        // If mu is all zeros, zero the output and skip the calculation
        if (mu(0) == 0.0 && mu(1) == 0.0 && mu(2) == 0.0 && mu(3) == 0.0 && mu(4) == 0.0 &&
            mu(5) == 0.0) {
            KokkosBlas::SerialSet::invoke(0., subview(qp_FD1, qp, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_FD2, qp, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_Duu, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_DD1, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_DD2, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_GD1, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_GD2, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_PD2, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_KD1, qp, ALL, ALL));
            KokkosBlas::SerialSet::invoke(0., subview(qp_KD2, qp, ALL, ALL));
            return;
        }

        const auto r0_data = Array<double, 4>{
            qp_r0(element, qp, 0), qp_r0(element, qp, 1), qp_r0(element, qp, 2),
            qp_r0(element, qp, 3)
        };
        const auto r0 = ConstView<double[4]>(r0_data.data());

        const auto xr_prime_data = Array<double, 3>{
            qp_x0_prime(element, qp, 0), qp_x0_prime(element, qp, 1), qp_x0_prime(element, qp, 2)
        };
        const auto xr_prime = ConstView<double[3]>(xr_prime_data.data());

        auto xr_data = Array<double, 4>{};
        const auto xr = View<double[4]>(xr_data.data());

        auto u_data = Array<double, 3>{};
        const auto u = View<double[3]>(u_data.data());
        auto u_prime_data = Array<double, 3>{};
        const auto u_prime = View<double[3]>(u_prime_data.data());

        auto r_data = Array<double, 4>{};
        const auto r = View<double[4]>(r_data.data());
        auto r_prime_data = Array<double, 4>{};
        const auto r_prime = View<double[4]>(r_prime_data.data());

        auto u_dot_data = Array<double, 3>{};
        const auto u_dot = View<double[3]>(u_dot_data.data());
        auto u_dot_prime_data = Array<double, 3>{};
        const auto u_dot_prime = View<double[3]>(u_dot_prime_data.data());

        auto omega_data = Array<double, 3>{};
        const auto omega = View<double[3]>(omega_data.data());
        auto omega_prime_data = Array<double, 3>{};
        const auto omega_prime = View<double[3]>(omega_prime_data.data());

        auto strain_data = Array<double, 6>{};
        const auto strain = View<double[6]>(strain_data.data());
        auto strain_dot_data = Array<double, 6>{};
        const auto strain_dot = View<double[6]>(strain_dot_data.data());

        auto Cstar_data = Array<double, 36>{};
        const auto Cstar = View<double[6][6]>(Cstar_data.data());
        auto Dstar_data = Array<double, 36>{};
        const auto Dstar = View<double[6][6]>(Dstar_data.data());

        auto mu_diag_data = Array<double, 36>{};
        const auto mu_diag = View<double[6][6]>(mu_diag_data.data());

        auto Cuu_data = Array<double, 36>{};
        const auto Cuu = View<double[6][6]>(Cuu_data.data());
        auto Duu_data = Array<double, 36>{};
        const auto Duu = View<double[6][6]>(Duu_data.data());

        auto FD1_data = Array<double, 6>{};
        const auto FD1 = View<double[6]>(FD1_data.data());
        auto FD2_data = Array<double, 6>{};
        const auto FD2 = View<double[6]>(FD2_data.data());

        auto DD1_data = Array<double, 36>{};
        const auto DD1 = View<double[6][6]>(DD1_data.data());
        auto DD2_data = Array<double, 36>{};
        const auto DD2 = View<double[6][6]>(DD2_data.data());
        auto GD1_data = Array<double, 36>{};
        const auto GD1 = View<double[6][6]>(GD1_data.data());
        auto GD2_data = Array<double, 36>{};
        const auto GD2 = View<double[6][6]>(GD2_data.data());
        auto PD2_data = Array<double, 36>{};
        const auto PD2 = View<double[6][6]>(PD2_data.data());
        auto KD1_data = Array<double, 36>{};
        const auto KD1 = View<double[6][6]>(KD1_data.data());
        auto KD2_data = Array<double, 36>{};
        const auto KD2 = View<double[6][6]>(KD2_data.data());

        CopyMatrix::invoke(subview(qp_Cstar, element, qp, ALL, ALL), Cstar);

        // Build damping matrix Dstar from Cstar and mu_diag
        KokkosBlas::SerialSet::invoke(0., mu_diag);
        for (auto i = 0U; i < 6U; ++i) {
            mu_diag(i, i) = mu(i);
        }
        Gemm::invoke(1., mu_diag, Cstar, 0., Dstar);

        beams::InterpolateToQuadraturePointForDamping<DeviceType>::invoke(
            qp_jacobian(qp), subview(shape_interp, ALL, qp), subview(shape_deriv, ALL, qp), node_u,
            node_u_dot, u, r, u_prime, r_prime, u_dot, omega, u_dot_prime, omega_prime
        );
        math::QuaternionCompose(r, r0, xr);

        // Rotate stiffness and damping matrices into the global coordinate system
        masses::RotateSectionMatrix<DeviceType>::invoke(xr, Cstar, Cuu);
        masses::RotateSectionMatrix<DeviceType>::invoke(xr, Dstar, Duu);

        // Calculate strain and extract curvature (kappa)
        beams::CalculateStrain<DeviceType>::invoke(xr_prime, u_prime, r, r_prime, strain);
        auto kappa_data = Array<double, 3>{};
        const auto kappa = View<double[3]>(kappa_data.data());
        CopyVector::invoke(subview(strain, make_pair(3UL, 6UL)), kappa);

        // Calculate strain rate and extract the strain rate components
        beams::CalculateStrainDot<DeviceType>::invoke(
            kappa, r, omega, u_dot_prime, omega_prime, xr_prime, strain_dot
        );
        auto eps_dot_data = Array<double, 3>{};
        const auto eps_dot = View<double[3]>(eps_dot_data.data());
        CopyVector::invoke(subview(strain_dot, make_pair(0UL, 3UL)), eps_dot);
        auto kappa_dot_data = Array<double, 3>{};
        const auto kappa_dot = View<double[3]>(kappa_dot_data.data());
        CopyVector::invoke(subview(strain_dot, make_pair(3UL, 6UL)), kappa_dot);

        // Calculate damping forces and copy to quadrature point values
        CalculateForceFD1<DeviceType>::invoke(Duu, strain_dot, FD1);
        CalculateForceFD2<DeviceType>::invoke(Duu, xr_prime, u_prime, strain_dot, FD2);
        CopyVector::invoke(FD1, subview(qp_FD1, qp, ALL));
        CopyVector::invoke(FD2, subview(qp_FD2, qp, ALL));

        // Calculate damping tangent matrices and copy to quadrature point values
        CalculateD_D1<DeviceType>::invoke(omega, Duu, DD1);
        CalculateD_D2<DeviceType>::invoke(xr_prime, u_prime, Duu, DD2);
        CalculateG_D1<DeviceType>::invoke(r, xr_prime, kappa, Duu, GD1);
        CalculateG_D2<DeviceType>::invoke(r, xr_prime, u_prime, kappa, Duu, GD2);
        CalculateP_D2<DeviceType>::invoke(xr_prime, u_prime, omega, eps_dot, kappa_dot, Duu, PD2);
        CalculateK_D1<DeviceType>::invoke(r, xr_prime, omega, kappa, eps_dot, kappa_dot, Duu, KD1);
        CalculateK_D2<DeviceType>::invoke(r, xr_prime, omega, kappa, eps_dot, kappa_dot, Duu, KD2);
        CopyMatrix::invoke(Duu, subview(qp_Duu, qp, ALL, ALL));
        CopyMatrix::invoke(DD1, subview(qp_DD1, qp, ALL, ALL));
        CopyMatrix::invoke(DD2, subview(qp_DD2, qp, ALL, ALL));
        CopyMatrix::invoke(GD1, subview(qp_GD1, qp, ALL, ALL));
        CopyMatrix::invoke(GD2, subview(qp_GD2, qp, ALL, ALL));
        CopyMatrix::invoke(PD2, subview(qp_PD2, qp, ALL, ALL));
        CopyMatrix::invoke(KD1, subview(qp_KD1, qp, ALL, ALL));
        CopyMatrix::invoke(KD2, subview(qp_KD2, qp, ALL, ALL));
    }
};

}  // namespace kynema::beams
