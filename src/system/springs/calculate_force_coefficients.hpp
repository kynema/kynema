#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::springs {

template <typename DeviceType>
KOKKOS_INLINE_FUNCTION double CalculateForceCoefficient1(double k, double l_ref, double l) {
    return k * (l_ref / l - 1.);
}

template <typename DeviceType>
KOKKOS_INLINE_FUNCTION double CalculateForceCoefficient2(double k, double l_ref, double l) {
    return k * l_ref / (l * l * l);
}

}  // namespace kynema_fmb::springs
