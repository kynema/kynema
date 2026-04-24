#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "elements/springs/springs.hpp"
#include "state/state.hpp"
#include "system/springs/calculate_quadrature_point_values.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void UpdateSystemVariablesSprings(
    const Springs<DeviceType>& springs, State<DeviceType>& state
) {
    auto region = Kokkos::Profiling::ScopedRegion("Update System Variables Springs");

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto range_policy = RangePolicy(0, springs.num_elems);

    Kokkos::parallel_for(
        "springs::CalculateQuadraturePointValues", range_policy,
        springs::CalculateQuadraturePointValues<DeviceType>{
            state.q, springs.node_state_indices, springs.x0, springs.l_ref, springs.k,
            springs.residual_vector_terms, springs.stiffness_matrix_terms
        }
    );
}

}  // namespace kynema_fmb::step
