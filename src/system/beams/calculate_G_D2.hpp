#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateG_D2 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[4]>& r, const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& u_prime, const ConstView<double[3]>& kappa, const ConstView<double[3][3]>& D, const View<double[6][6]>& G_D2) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., G_D2);

        auto r_xr_prime_data = Kokkos::Array<double, 3>{};
        auto r_xr_prime = View<double[3]>(r_xr_prime_data.data());
        math::RotateVectorByQuaternion(r, xr_prime, r_xr_prime);
        auto r_xr_prime_tilde_data = Kokkos::Array<double, 9>{};
        auto r_xr_prime_tilde = View<double[3][3]>(r_xr_prime_tilde_data.data());
        math::VecTilde(r_xr_prime, r_xr_prime_tilde);

        auto kappa_tilde_data = Kokkos::Array<double, 9>{};
        auto kappa_tilde = View<double[3][3]>(kappa_tilde_data.data());
        math::VecTilde(kappa, kappa_tilde);

        auto xr_u_prime_data = Kokkos::Array<double, 3>{};
        auto xr_u_prime = View<double[3]>(xr_u_prime_data.data());
        for (auto component = 0; component < 3; ++component) {
            xr_u_prime(component) = xr_prime(component) + u_prime(component);
        }
        auto tilde_xp_up_data = Kokkos::Array<double, 9>{};
        auto tilde_xp_up = View<double[3][3]>(tilde_xp_up_data.data());
        math::VecTilde(xr_u_prime, tilde_xp_up);

        auto d_tmp_data = Kokkos::Array<double, 9>{};
        auto d_tmp = View<double[3][3]>(d_tmp_data.data());
        auto D11 = subview(D, make_pair(0, 3), make_pair(0, 3));
        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        GemmNN::invoke(1., D11, r_xr_prime_tilde, 0., d_tmp);
        GemmNN::invoke(-1., D12, kappa_tilde, 1., d_tmp);

        auto G_D2_22 = subview(G_D2, make_pair(3, 6), make_pair(3, 6));
        GemmTN::invoke(1., tilde_xp_up, d_tmp, 0., G_D2_22);
    }
};
}
