#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::springs {

template <typename DeviceType>
struct CalculateForceVectors {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    KOKKOS_FUNCTION static void invoke(
        const ConstView<double[3]>& r, double c1, const View<double[3]>& f
    ) {
        for (auto component = 0; component < 3; ++component) {
            f(component) = c1 * r(component);
        }
    }
};
}  // namespace kynema_fmb::springs
