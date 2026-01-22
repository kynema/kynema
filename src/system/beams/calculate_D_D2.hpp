#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateD_D2 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& u_prime, const ConstView<double[6][6]>& D, const View<double[6][6]>& D_D2) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., D_D2);

        auto xr_u_prime_data = Kokkos::Array<double, 3>{};
        auto xr_u_prime = View<double[3]>(xr_u_prime_data.data());
        for (auto component = 0; component < 3; ++component) {
            xr_u_prime(component) = xr_prime(component) + u_prime(component);
        }
        auto tilde_xp_up_data = Kokkos::Array<double, 9>{};
        auto tilde_xp_up = View<double[3][3]>(tilde_xp_up_data.data());
        math::VecTilde(xr_u_prime, tilde_xp_up);

        auto D11 = subview(D, make_pair(0, 3), make_pair(0, 3));
        auto D_D2_21 = subview(D, make_pair(3, 6), make_pair(0, 3));
        GemmTN::invoke(1., tilde_xp_up, D11, 0., D_D2_21);

        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        auto D_D2_22 = subview(D, make_pair(3, 6), make_pair(3, 6));
        GemmTN::invoke(1., tilde_xp_up, D12, 0., D_D2_22);
    }
};
}
