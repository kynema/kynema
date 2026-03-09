#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "calculate_strain_dot.hpp"
#include "calculate_D_D1.hpp"
#include "calculate_D_D2.hpp"
#include "calculate_G_D1.hpp"
#include "calculate_G_D2.hpp"
#include "calculate_K_D1.hpp"
#include "calculate_K_D2.hpp"
#include "calculate_P_D2.hpp"
#include "calculate_force_FD1.hpp"
#include "calculate_force_FD2.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema::beams {
template <typename DeviceType>
struct CalculateDampingQuadraturePointValues {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;
    template <typename ValueType>
    using ConstLeftView = typename LeftView<ValueType>::const_type;

    size_t element;

    ConstView<double*> qp_jacobian;
    ConstLeftView<double**> shape_interp;
    ConstLeftView<double**> shape_deriv;
    ConstView<double** [4]> qp_r0;
    ConstView<double** [3]> qp_x0_prime;
    ConstView<double** [6][6]> qp_Cstar;
    ConstView<double* [7]> node_u;
    ConstView<double* [7]> node_v;

    View<double* [6]> qp_FD1;
    View<double* [6]> qp_FD2;
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

        const auto r0_data = Array<double, 4>{
            qp_r0(element, qp, 0), qp_r0(element, qp, 1), qp_r0(element, qp, 2),
            qp_r0(element, qp, 3)
        };
        const auto x0_prime_data = Array<double, 3>{
            qp_x0_prime(element, qp, 0), qp_x0_prime(element, qp, 1), qp_x0_prime(element, qp, 2)
        };
        auto xr_data = Array<double, 4>{};
        auto u_data = Array<double, 3>{};
        auto u_prime_data = Array<double, 3>{};
        auto r_data = Array<double, 4>{};
        auto r_prime_data = Array<double, 4>{};
        auto Cstar_data = Array<double, 36>{};

        auto strain_data = Array<double, 6>{};
        auto x0pupSS_data = Array<double, 9>{};
        auto M_tilde_data = Array<double, 9>{};
        auto N_tilde_data = Array<double, 9>{};
        auto FD1_data = Array<double, 6>{};
        auto FD2_data = Array<double, 6>{};
        auto Cuu_data = Array<double, 36>{};
        // auto Ouu_data = Array<double, 36>{};
        // auto Puu_data = Array<double, 36>{};
        // auto Quu_data = Array<double, 36>{};

        const auto r0 = ConstView<double[4]>(r0_data.data());
        const auto x0_prime = ConstView<double[3]>(x0_prime_data.data());
        const auto xr = View<double[4]>(xr_data.data());
        const auto u = View<double[3]>(u_data.data());
        const auto u_prime = View<double[3]>(u_prime_data.data());
        const auto r = View<double[4]>(r_data.data());
        const auto r_prime = View<double[4]>(r_prime_data.data());
        const auto strain = View<double[6]>(strain_data.data());
        const auto x0pupSS = View<double[3][3]>(x0pupSS_data.data());
        const auto M_tilde = View<double[3][3]>(M_tilde_data.data());
        const auto N_tilde = View<double[3][3]>(N_tilde_data.data());
        const auto FD1 = View<double[6]>(FD1_data.data());
        const auto FD2 = View<double[6]>(FD2_data.data());
        const auto Cstar = View<double[6][6]>(Cstar_data.data());
        const auto Cuu = View<double[6][6]>(Cuu_data.data());
        // const auto Ouu = View<double[6][6]>(Ouu_data.data());
        // const auto Puu = View<double[6][6]>(Puu_data.data());
        // const auto Quu = View<double[6][6]>(Quu_data.data());

        CopyMatrix::invoke(subview(qp_Cstar, element, qp, ALL, ALL), Cstar);

        // beams::InterpolateToQuadraturePointForStiffness<DeviceType>::invoke(
        //     qp_jacobian(qp), subview(shape_interp, ALL, qp), subview(shape_deriv, ALL, qp), node_u,
        //     u, r, u_prime, r_prime
        // );
        math::QuaternionCompose(r, r0, xr);

        masses::RotateSectionMatrix<DeviceType>::invoke(xr, Cstar, Cuu);

        CalculateForceFD1<DeviceType>::invoke();
        CalculateForceFD2<DeviceType>::invoke();
        // CalculateD_D1<DeviceType>::invoke(omega, D, D_D1);
        CalculateD_D2<DeviceType>::invoke();
        CalculateG_D1<DeviceType>::invoke();
        CalculateG_D2<DeviceType>::invoke();
    }
};

}  // namespace kynema::beams
