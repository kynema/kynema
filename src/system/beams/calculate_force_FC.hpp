#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateForceFC {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[6][6]>& Cuu, const ConstView<double[6]>& strain,
        const View<double[6]>& FC
    ) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        Gemv::invoke(1., Cuu, strain, 0., FC);
    }
};
}  // namespace kynema_fmb::beams
