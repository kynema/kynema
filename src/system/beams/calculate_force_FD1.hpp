#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema::beams {

template <typename DeviceType>
struct CalculateForceFD1 {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[6][6]>& D, const ConstView<double[6]>& strain_dot,
        const View<double[6]>& FD1
    ) {
        using NoTranspose = KokkosBlas::Trans::NoTranspose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<NoTranspose, Default>;
        Gemv::invoke(1., D, strain_dot, 0., FD1);
    }
};

}  // namespace kynema::beams
