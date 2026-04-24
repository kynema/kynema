#pragma once

#include <KokkosBlas.hpp>
#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

template <typename DeviceType>
struct CalculateForceFD {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3][3]>& x0pupSS, const ConstView<double[6]>& FC,
        const View<double[6]>& FD
    ) {
        using Transpose = KokkosBlas::Trans::Transpose;
        using Default = KokkosBlas::Algo::Gemv::Default;
        using Gemv = KokkosBlas::SerialGemv<Transpose, Default>;
        using Kokkos::make_pair;
        using Kokkos::subview;

        for (auto component = 0U; component < 6U; ++component) {
            FD(component) = 0.;
        }
        Gemv::invoke(1., x0pupSS, subview(FC, make_pair(0, 3)), 0., subview(FD, make_pair(3, 6)));
    }
};
}  // namespace kynema_fmb::beams
