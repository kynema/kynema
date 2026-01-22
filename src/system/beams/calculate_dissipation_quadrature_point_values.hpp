#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "calculate_force_FD1_viscoelastic.hpp"
#include "calculate_force_FD2.hpp"
#include "calculate_D_D1.hpp"
#include "calculate_D_D2.hpp"
#include "calculate_G_D1.hpp"
#include "calculate_G_D2.hpp"
#include "calculate_P_D2.hpp"
#include "calculate_K_D1.hpp"
#include "calculate_K_D2.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema::beams {
template <typename DeviceType>
struct CalculateDissipationQuadraturePointValues {
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
    ConstView<double**[4]> qp_r0;
    ConstView<double** [3]> qp_x0_prime;
    ConstView<double* [7]> node_u;

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

        const auto r0 = ConstView<double[4]>(r0_data.data());
        const auto x0_prime = ConstView<double[3]>(x0_prime_data.data());

        CalculateForceFD1Viscoelastic<DeviceType>::invoke(h, tau_i, kv_i, rr0, strain_dot_n, strain_dot_n1, visco_hist, FD1);
        CalculateForceFD2<DeviceType>::invoke();
        CalculateD_D1<DeviceType>::invoke(omega, D, D_D1);
        CalculateD_D2<DeviceType>::invoke();
        CalculateG_D1<DeviceType>::invoke();
        CalculateG_D2<DeviceType>::invoke();
    }
};

}
