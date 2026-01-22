#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateP_D2 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& u_prime, const ConstView<double[3]>& omega,const ConstView<double[3]>& eps_dot, const ConstView<double[3]>& kappa_dot, const ConstView<double[3][3]>& D, const View<double[6][6]>& P_D2) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Gemv = KokkosBlas::SerialGemv<KokkosBlas::Trans::NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., P_D2);

        auto omega_tilde_data = Kokkos::Array<double, 9>{};
        auto omega_tilde = View<double[3][3]>(omega_tilde_data.data());
        math::VecTilde(omega, omega_tilde);

        auto D11 = subview(D, make_pair(0, 3), make_pair(0, 3));
        auto D11_eps_dot_data = Kokkos::Array<double, 3>{};
        auto D11_eps_dot = View<double[3]>(D11_eps_dot_data.data());
        Gemv::invoke(1., D11, eps_dot, 0., D11_eps_dot);
        auto D11_eps_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D11_eps_dot_tilde = View<double[3][3]>(D11_eps_dot_tilde_data.data());
        math::VecTilde(D11_eps_dot, D11_eps_dot_tilde);

        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        auto D12_kappa_dot_data = Kokkos::Array<double, 3>{};
        auto D12_kappa_dot = View<double[3]>(D12_kappa_dot_data.data());
        Gemv::invoke(1., D12, kappa_dot, 0., D12_kappa_dot);
        auto D12_kappa_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D12_kappa_dot_tilde = View<double[3][3]>(D12_kappa_dot_tilde_data.data());
        math::VecTilde(D12_kappa_dot, D12_kappa_dot_tilde);

        auto P_D2_21 = subview(P_D2, make_pair(3, 6), make_pair(0, 3));
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                P_D2_21(i, j) = -D11_eps_dot_tilde(i, j) - D12_kappa_dot_tilde(i, j);
            }
        }

        auto D12_omega_tilde_data = Kokkos::Array<double, 9>{};
        auto D12_omega_tilde = View<double[3][3]>(D12_omega_tilde_data.data());
        GemmNN::invoke(1., D12, omega_tilde, 0., D12_omega_tilde);

        auto xr_u_prime_data = Kokkos::Array<double, 3>{};
        auto xr_u_prime = View<double[3]>(xr_u_prime_data.data());
        for (auto component = 0; component < 3; ++component) {
            xr_u_prime(component) = xr_prime(component) + u_prime(component);
        }
        auto tilde_xp_up_data = Kokkos::Array<double, 9>{};
        auto tilde_xp_up = View<double[3][3]>(tilde_xp_up_data.data());
        math::VecTilde(xr_u_prime, tilde_xp_up);

        auto P_D2_22 = subview(P_D2, make_pair(3, 6), make_pair(3, 6));
        GemmTN::invoke(1., tilde_xp_up, D12_omega_tilde, 0., P_D2_22);
    }
};
}
