#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema_fmb::masses {

template <typename DeviceType>
struct CalculateInertiaStiffnessMatrix {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double mass, const ConstView<double[3]>& u_ddot, const ConstView<double[3]>& omega,
        const ConstView<double[3]>& omega_dot, const ConstView<double[3]>& eta,
        const ConstView<double[3][3]>& rho, const ConstView<double[3][3]>& omega_tilde,
        const ConstView<double[3][3]>& omega_dot_tilde, const View<double[6][6]>& Kuu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using GemmDefault = KokkosBatched::Algo::Gemm::Default;
        using GemvDefault = KokkosBatched::Algo::Gemv::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, GemmDefault>;
        using GemmNT = KokkosBatched::SerialGemm<NoTranspose, Transpose, GemmDefault>;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, GemvDefault>;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto v1 = Array<double, 3>{};
        auto V1 = View<double[3]>(v1.data());
        auto m1 = Array<double, 9>{};
        auto M1 = View<double[3][3]>(m1.data());
        auto m2 = Array<double, 9>{};
        auto M2 = View<double[3][3]>(m2.data());

        KokkosBlas::serial_axpy(1., omega_dot_tilde, M1);
        GemmNN::invoke(1., omega_tilde, omega_tilde, 1., M1);
        KokkosBlas::serial_axpy(mass, eta, V1);
        math::VecTilde(V1, M2);
        auto Kuu_12 = subview(Kuu, make_pair(0, 3), make_pair(3, 6));
        GemmNT::invoke(1., M1, M2, 0., Kuu_12);
        math::VecTilde(u_ddot, M1);
        auto Kuu_22 = subview(Kuu, make_pair(3, 6), make_pair(3, 6));
        GemmNN::invoke(1., rho, omega_dot_tilde, 0., Kuu_22);
        GemmNN::invoke(1., M1, M2, 1., Kuu_22);
        Gemv::invoke(1., rho, omega_dot, 0., V1);
        math::VecTilde(V1, M2);
        KokkosBlas::serial_axpy(-1., M2, Kuu_22);
        Gemv::invoke(1., rho, omega, 0., V1);
        math::VecTilde(V1, M1);
        GemmNN::invoke(1., rho, omega_tilde, -1., M1);
        GemmNN::invoke(1., omega_tilde, M1, 1., Kuu_22);
    }
};
}  // namespace kynema_fmb::masses
