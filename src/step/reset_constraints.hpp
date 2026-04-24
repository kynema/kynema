#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/constraints.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void ResetConstraints(Constraints<DeviceType>& constraints) {
    auto region = Kokkos::Profiling::ScopedRegion("Reset Constraints");
    Kokkos::deep_copy(constraints.lambda, 0.);
    Kokkos::deep_copy(constraints.system_residual_terms, 0.);
}

}  // namespace kynema_fmb::step
