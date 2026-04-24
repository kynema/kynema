#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateOuu {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[6][6]>& Cuu, const ConstView<double[3][3]>& x0pupSS,
        const ConstView<double[3][3]>& M_tilde, const ConstView<double[3][3]>& N_tilde,
        const View<double[6][6]>& Ouu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using Gemm = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto C11 = subview(Cuu, make_pair(0, 3), make_pair(0, 3));
        auto C21 = subview(Cuu, make_pair(3, 6), make_pair(0, 3));
        KokkosBlas::SerialSet::invoke(0., Ouu);
        auto Ouu_12 = subview(Ouu, make_pair(0, 3), make_pair(3, 6));
        KokkosBlas::serial_axpy(1., N_tilde, Ouu_12);
        Gemm::invoke(1., C11, x0pupSS, -1., Ouu_12);
        auto Ouu_22 = subview(Ouu, make_pair(3, 6), make_pair(3, 6));
        KokkosBlas::serial_axpy(1., M_tilde, Ouu_22);
        Gemm::invoke(1., C21, x0pupSS, -1., Ouu_22);
    }
};
}  // namespace kynema_fmb::beams
