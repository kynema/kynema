#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "elements/elements.hpp"
#include "state/state.hpp"
#include "update_system_variables_beams.hpp"
#include "update_system_variables_masses.hpp"
#include "update_system_variables_springs.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void UpdateSystemVariables(
    StepParameters& parameters, const Elements<DeviceType>& elements, State<DeviceType>& state
) {
    auto region = Kokkos::Profiling::ScopedRegion("Update System Variables");

    if (elements.beams.num_elems > 0) {
        UpdateSystemVariablesBeams(parameters, elements.beams, state);
    }
    if (elements.masses.num_elems > 0) {
        UpdateSystemVariablesMasses(parameters, elements.masses, state);
    }
    if (elements.springs.num_elems > 0) {
        UpdateSystemVariablesSprings(elements.springs, state);
    }
}

}  // namespace kynema_fmb::step
