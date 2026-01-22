#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateK_D2 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[4]>& r, const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& omega, const ConstView<double[3]>& kappa, const ConstView<double[3]>& eps_dot, const ConstView<double[3]>& kappa_dot, const ConstView<double[3][3]>& D, const View<double[6][6]>& K_D2) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Gemv = KokkosBlas::SerialGemv<KokkosBlas::Trans::NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., K_D2);

        auto omega_tilde_data = Kokkos::Array<double, 9>{};
        auto omega_tilde = View<double[3][3]>(omega_tilde_data.data());
        math::VecTilde(omega, omega_tilde);

        auto kappa_tilde_data = Kokkos::Array<double, 0>{};
        auto kappa_tilde = View<double[3][3]>(kappa_tilde_data.data());
        math::VecTilde(kappa, kappa_tilde);

        auto eps_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto eps_dot_tilde = View<double[3][3]>(eps_dot_tilde_data.data());
        math::VecTilde(eps_dot, eps_dot_tilde);

        auto kappa_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto kappa_dot_tilde = View<double[3][3]>(kappa_dot_tilde_data.data());
        math::VecTilde(kappa_dot, kappa_dot_tilde);

        auto D11 = subview(D, make_pair(0, 3), make_pair(0, 3));
        auto D11_eps_dot_data = Kokkos::Array<double, 3>{};
        auto D11_eps_dot = View<double[3]>(D11_eps_dot_data.data());
        Gemv::invoke(1., D11, eps_dot, 0., D11_eps_dot);
        auto D11_eps_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D11_eps_dot_tilde = View<double[3][3]>(D11_eps_dot_tilde_data.data());
        math::VecTilde(D11_eps_dot, D11_eps_dot_tilde);

        auto D11_times_eps_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D11_times_eps_dot_tilde = View<double[3][3]>(D11_times_eps_dot_tilde_data.data());
        GemmNN::invoke(1., D11, eps_dot_tilde, 0., D11_times_eps_dot_tilde);

        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        auto D12_kappa_dot_data = Kokkos::Array<double, 3>{};
        auto D12_kappa_dot = View<double[3]>(D12_kappa_dot_data.data());
        Gemv::invoke(1., D12, kappa_dot, 0., D12_kappa_dot);
        auto D12_kappa_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D12_kappa_dot_tilde = View<double[3][3]>(D12_kappa_dot_tilde_data.data());
        math::VecTilde(D12_kappa_dot, D12_kappa_dot_tilde);

        auto D12_times_kappa_dot_tilde_data = Kokkos::Array<double, 9>{};
        auto D12_times_kappa_dot_tilde = View<double[3][3]>(D12_times_kappa_dot_tilde_data.data());
        GemmNN::invoke(1., D12, kappa_dot_tilde, 0., D12_times_kappa_dot_tilde);

        auto D11_omega_tilde_data = Kokkos::Array<double, 9>{};
        auto D11_omega_tilde = View<double[3][3]>(D11_omega_tilde_data.data());
        GemmNN::invoke(1., D11, omega_tilde, 0., D11_omega_tilde);

        auto D12_omega_tilde_data = Kokkos::Array<double, 9>{};
        auto D12_omega_tilde = View<double[3][3]>(D12_omega_tilde_data.data());
        GemmNN::invoke(1., D12, omega_tilde, 0., D12_omega_tilde);

        auto r_xr_prime_data = Kokkos::Array<double, 3>{};
        auto r_xr_prime = View<double[3]>(r_xr_prime_data.data());
        math::RotateVectorByQuaternion(r, xr_prime, r_xr_prime);
        auto r_xr_prime_tilde_data = Kokkos::Array<double, 9>{};
        auto r_xr_prime_tilde = View<double[3][3]>(r_xr_prime_tilde_data.data());
        math::VecTilde(r_xr_prime, r_xr_prime_tilde);

        auto K_D2_22 = subview(K_D2, make_pair(0, 3), make_pair(3, 6));
        GemmNN::invoke(1., D11_omega_tilde, r_xr_prime_tilde, 0., K_D2_22);
        GemmNN::invoke(-1., D12_omega_tilde, kappa_tilde, 1., K_D2_22);
        for (auto i =  0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                K_D2_22(i, j) += -D11_times_eps_dot_tilde(i, j) + D11_eps_dot_tilde(i, j) - D12_times_kappa_dot_tilde(i, j) + D12_kappa_dot_tilde(i, j);
            }
        }
    }
};
}
