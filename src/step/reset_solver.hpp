#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "solver/solver.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void ResetSolver(Solver<DeviceType>& solver) {
    auto region = Kokkos::Profiling::ScopedRegion("Reset Solver");
    Kokkos::deep_copy(solver.A.values, 0.);
    Kokkos::deep_copy(solver.b, 0.);
}

}  // namespace kynema_fmb::step
