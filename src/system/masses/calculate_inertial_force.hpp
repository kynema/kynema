#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::masses {

template <typename DeviceType>
struct CalculateInertialForce {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double mass, const ConstView<double[3]>& u_ddot, const ConstView<double[3]>& omega,
        const ConstView<double[3]>& omega_dot, const ConstView<double[3]>& eta,
        const ConstView<double[3][3]>& eta_tilde, const ConstView<double[3][3]>& rho,
        const ConstView<double[3][3]>& omega_tilde, const ConstView<double[3][3]>& omega_dot_tilde,
        const View<double[6]>& FI
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using GemmDefault = KokkosBatched::Algo::Gemm::Default;
        using GemvDefault = KokkosBlas::Algo::Gemv::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, GemmDefault>;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, GemvDefault>;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto v1 = Array<double, 3>{};
        auto V1 = View<double[3]>(v1.data());
        auto m1 = Array<double, 9>{};
        auto M1 = View<double[3][3]>(m1.data());

        // Compute first 3 components of FI
        // FI_1 = m * u_ddot + (omega_dot_tilde + omega_tilde * omega_tilde) * m * eta
        auto FI_1 = subview(FI, make_pair(0, 3));
        Gemm::invoke(mass, omega_tilde, omega_tilde, 0., M1);
        KokkosBlas::serial_axpy(mass, omega_dot_tilde, M1);
        Gemv::invoke(1., M1, eta, 0., FI_1);
        KokkosBlas::serial_axpy(mass, u_ddot, FI_1);

        // Compute last 3 components of FI
        // FI_2 = m * eta_tilde * u_ddot + rho * omega_dot + omega_tilde * rho * omega
        auto FI_2 = subview(FI, make_pair(3, 6));
        KokkosBlas::serial_axpy(mass, u_ddot, V1);
        Gemv::invoke(1., eta_tilde, V1, 0., FI_2);
        Gemv::invoke(1., rho, omega_dot, 1., FI_2);
        Gemm::invoke(1., omega_tilde, rho, 0., M1);
        Gemv::invoke(1., M1, omega, 1., FI_2);
    }
};
}  // namespace kynema_fmb::masses
