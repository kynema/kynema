#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::masses {

template <typename DeviceType>
struct CalculateGravityForce {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        double mass, const ConstView<double[3]>& gravity, const ConstView<double[3][3]>& eta_tilde,
        const View<double[6]>& FG
    ) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        KokkosBlas::serial_axpy(mass, gravity, subview(FG, make_pair(0, 3)));
        Gemv::invoke(1., eta_tilde, subview(FG, make_pair(0, 3)), 0., subview(FG, make_pair(3, 6)));
    }
};
}  // namespace kynema_fmb::masses
