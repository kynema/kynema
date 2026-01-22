#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateG_D1 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(const ConstView<double[4]>& r, const ConstView<double[3]>& xr_prime, const ConstView<double[3]>& kappa, const ConstView<double[6][6]>& D, const View<double[6][6]>& G_D1) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., G_D1);

        auto r_xr_prime_data = Kokkos::Array<double, 3>{};
        auto r_xr_prime = View<double[3]>(r_xr_prime_data.data());
        math::RotateVectorByQuaternion(r, xr_prime, r_xr_prime);
        auto r_xr_prime_tilde_data = Kokkos::Array<double, 9>{};
        auto r_xr_prime_tilde = View<double[3][3]>(r_xr_prime_tilde_data.data());
        math::VecTilde(r_xr_prime, r_xr_prime_tilde);

        auto kappa_tilde_data = Kokkos::Array<double, 9>{};
        auto kappa_tilde = View<double[3][3]>(kappa_tilde_data.data());
        math::VecTilde(kappa, kappa_tilde);

        auto D11 = subview(D, make_pair(0, 3), make_pair(0, 3));
        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        auto G_D1_12 = subview(G_D1, make_pair(0, 3), make_pair(3, 6));
        GemmNN::invoke(1., D11, r_xr_prime_tilde, 0., G_D1_12);
        GemmNN::invoke(-1., D12, kappa_tilde, 1., G_D1_12);

        auto D21 = subview(D, make_pair(3, 6), make_pair(0, 3));
        auto D22 = subview(D, make_pair(3, 6), make_pair(3, 6));
        auto G_D1_22 = subview(G_D1, make_pair(3, 6), make_pair(3,6));
        GemmNN::invoke(1., D21, r_xr_prime_tilde, 0., G_D1_22);
        GemmNN::invoke(-1., D22, kappa_tilde, 1., G_D1_22);
    }
};
}
