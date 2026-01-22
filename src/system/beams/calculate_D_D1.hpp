#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

#include "math/vector_operations.hpp"

namespace kynema::beams {

template <typename DeviceType>
struct CalculateD_D1 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& omega, const ConstView<double[6][6]>& D,
        const View<double[6][6]>& D_D1
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::SerialSet::invoke(0., D_D1);

        auto tilde_omega_data = Kokkos::Array<double, 9>{};
        auto tilde_omega = View<double[3][3]>(tilde_omega_data.data());
        math::VecTilde(omega, tilde_omega);

        auto D12 = subview(D, make_pair(0, 3), make_pair(3, 6));
        auto D_D1_12 = subview(D_D1, make_pair(0, 3), make_pair(3, 6));
        GemmNN::invoke(1., D12, tilde_omega, 0., D_D1_12);

        auto D22 = subview(D, make_pair(3, 6), make_pair(3, 6));
        auto D_D1_22 = subview(D_D1, make_pair(3, 6), make_pair(3, 6));
        GemmNN::invoke(1., D22, tilde_omega, 0., D_D1_22);
    }
};
}  // namespace kynema::beams
