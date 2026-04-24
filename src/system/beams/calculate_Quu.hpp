#pragma once

#include <KokkosBatched_Gemm_Decl.hpp>
#include <KokkosBlas.hpp>
#include <KokkosBlas1_set.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateQuu {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[6][6]>& Cuu, const ConstView<double[3][3]>& x0pupSS,
        const ConstView<double[3][3]>& N_tilde, const View<double[6][6]>& Quu
    ) {
        using NoTranspose = KokkosBatched::Trans::NoTranspose;
        using Transpose = KokkosBatched::Trans::Transpose;
        using Default = KokkosBatched::Algo::Gemm::Default;
        using GemmNN = KokkosBatched::SerialGemm<NoTranspose, NoTranspose, Default>;
        using GemmTN = KokkosBatched::SerialGemm<Transpose, NoTranspose, Default>;
        using Kokkos::Array;
        using Kokkos::make_pair;
        using Kokkos::subview;

        auto m1 = Array<double, 9>{};
        auto M1 = View<double[3][3]>(m1.data());
        auto C11 = subview(Cuu, make_pair(0, 3), make_pair(0, 3));
        KokkosBlas::SerialSet::invoke(0., Quu);
        KokkosBlas::serial_axpy(1., N_tilde, M1);
        GemmNN::invoke(1., C11, x0pupSS, -1., M1);
        auto Quu_22 = subview(Quu, make_pair(3, 6), make_pair(3, 6));
        GemmTN::invoke(1., x0pupSS, M1, 0., Quu_22);
    }
};
}  // namespace kynema_fmb::beams
