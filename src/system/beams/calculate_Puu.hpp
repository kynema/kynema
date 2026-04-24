#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculatePuu {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[6][6]>& Cuu, const ConstView<double[3][3]>& x0pupSS,
        const ConstView<double[3][3]>& N_tilde, const View<double[6][6]>& Puu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto C11 = subview(Cuu, make_pair(0, 3), make_pair(0, 3));
        auto C12 = subview(Cuu, make_pair(0, 3), make_pair(3, 6));
        KokkosBlas::SerialSet::invoke(0., Puu);
        auto Puu_21 = subview(Puu, make_pair(3, 6), make_pair(0, 3));
        KokkosBlas::serial_axpy(1., N_tilde, Puu_21);
        GemmTN::invoke(1., x0pupSS, C11, 1., Puu_21);
        auto Puu_22 = subview(Puu, make_pair(3, 6), make_pair(3, 6));
        GemmTN::invoke(1., x0pupSS, C12, 0., Puu_22);
    }
};
}  // namespace kynema_fmb::beams
