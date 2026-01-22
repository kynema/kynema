#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <Kokkos_Core.hpp>

#include "calculate_Ouu.hpp"
#include "calculate_Puu.hpp"
#include "calculate_Quu.hpp"
#include "calculate_force_FC.hpp"
#include "calculate_force_FD.hpp"
#include "calculate_strain.hpp"
#include "calculate_temporary_variables.hpp"
#include "interpolate_to_quadrature_point_for_stiffness.hpp"
#include "system/masses/rotate_section_matrix.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateStiffnessQuadraturePointValues {
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

    View<double* [6]> qp_Fc;
    View<double* [6]> qp_Fd;
    View<double* [6][6]> qp_Cuu;
    View<double* [6][6]> qp_Ouu;
    View<double* [6][6]> qp_Puu;
    View<double* [6][6]> qp_Quu;

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
        auto FC_data = Array<double, 6>{};
        auto FD_data = Array<double, 6>{};
        auto Cuu_data = Array<double, 36>{};
        auto Ouu_data = Array<double, 36>{};
        auto Puu_data = Array<double, 36>{};
        auto Quu_data = Array<double, 36>{};

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
        const auto FC = View<double[6]>(FC_data.data());
        const auto FD = View<double[6]>(FD_data.data());
        const auto Cstar = View<double[6][6]>(Cstar_data.data());
        const auto Cuu = View<double[6][6]>(Cuu_data.data());
        const auto Ouu = View<double[6][6]>(Ouu_data.data());
        const auto Puu = View<double[6][6]>(Puu_data.data());
        const auto Quu = View<double[6][6]>(Quu_data.data());

        CopyMatrix::invoke(subview(qp_Cstar, element, qp, ALL, ALL), Cstar);

        beams::InterpolateToQuadraturePointForStiffness<DeviceType>::invoke(
            qp_jacobian(qp), subview(shape_interp, ALL, qp), subview(shape_deriv, ALL, qp), node_u,
            u, r, u_prime, r_prime
        );
        math::QuaternionCompose(r, r0, xr);

        masses::RotateSectionMatrix<DeviceType>::invoke(xr, Cstar, Cuu);

        beams::CalculateStrain<DeviceType>::invoke(x0_prime, u_prime, r, r_prime, strain);
        beams::CalculateTemporaryVariables<DeviceType>::invoke(x0_prime, u_prime, x0pupSS);
        beams::CalculateForceFC<DeviceType>::invoke(Cuu, strain, FC);
        beams::CalculateForceFD<DeviceType>::invoke(x0pupSS, FC, FD);

        math::VecTilde(subview(FC, make_pair(0, 3)), N_tilde);
        math::VecTilde(subview(FC, make_pair(3, 6)), M_tilde);

        beams::CalculateOuu<DeviceType>::invoke(Cuu, x0pupSS, M_tilde, N_tilde, Ouu);
        beams::CalculatePuu<DeviceType>::invoke(Cuu, x0pupSS, N_tilde, Puu);
        beams::CalculateQuu<DeviceType>::invoke(Cuu, x0pupSS, N_tilde, Quu);

        CopyVector::invoke(FC, subview(qp_Fc, qp, ALL));
        CopyVector::invoke(FD, subview(qp_Fd, qp, ALL));

        CopyMatrix::invoke(Cuu, subview(qp_Cuu, qp, ALL, ALL));
        CopyMatrix::invoke(Ouu, subview(qp_Ouu, qp, ALL, ALL));
        CopyMatrix::invoke(Puu, subview(qp_Puu, qp, ALL, ALL));
        CopyMatrix::invoke(Quu, subview(qp_Quu, qp, ALL, ALL));
    }
};

}  // namespace kynema::beams
