#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::springs {

template <typename DeviceType>
struct CalculateDistanceComponents {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& x0, const ConstView<double[3]>& u1,
        const ConstView<double[3]>& u2, const View<double[3]>& r
    ) {
        for (auto component = 0; component < 3; ++component) {
            r(component) = x0(component) + u2(component) - u1(component);
        }
    }
};
}  // namespace kynema_fmb::springs
