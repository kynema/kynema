#pragma once

#include <KokkosBatched_Copy_Decl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema_fmb::masses {

template <typename DeviceType>
struct CalculateGyroscopicMatrix {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double mass, const ConstView<double[3]>& omega, const ConstView<double[3]>& eta,
        const ConstView<double[3][3]>& rho, const ConstView<double[3][3]>& omega_tilde,
        const Kokkos::View<double[6][6], DeviceType>& Guu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using GemmDefault = KokkosBatched::Algo::Gemm::Default;
        using GemvDefault = KokkosBatched::Algo::Gemv::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, GemmDefault>;
        using GemmNT = KokkosBatched::SerialGemm<NoTranspose, Transpose, GemmDefault>;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, GemvDefault>;
        using CopyMatrix = KokkosBatched::SerialCopy<>;
        using CopyMatrixTranspose = KokkosBatched::SerialCopy<Transpose>;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto v1 = Array<double, 3>{};
        auto V1 = View<double[3]>(v1.data());
        auto v2 = Array<double, 3>{};
        auto V2 = View<double[3]>(v2.data());
        auto m1 = Array<double, 9>{};
        auto M1 = View<double[3][3]>(m1.data());

        // omega.tilde() * m * eta.tilde().t() + (omega.tilde() * m * eta).tilde().t()
        auto Guu_12 = subview(Guu, make_pair(0, 3), make_pair(3, 6));
        KokkosBlas::serial_axpy(mass, eta, V1);
        Gemv::invoke(1., omega_tilde, V1, 0., V2);
        math::VecTilde(V2, M1);
        CopyMatrixTranspose::invoke(M1, Guu_12);

        math::VecTilde(V1, M1);
        GemmNT::invoke(1., omega_tilde, M1, 1., Guu_12);
        // Guu_22 = omega.tilde() * rho - (rho * omega).tilde()
        Gemv::invoke(1., rho, omega, 0., V1);
        math::VecTilde(V1, M1);
        auto Guu_22 = subview(Guu, make_pair(3, 6), make_pair(3, 6));
        CopyMatrix::invoke(M1, Guu_22);
        GemmNN::invoke(1., omega_tilde, rho, -1., Guu_22);
    }
};
}  // namespace kynema_fmb::masses
