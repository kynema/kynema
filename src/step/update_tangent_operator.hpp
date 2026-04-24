#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "state/state.hpp"
#include "step_parameters.hpp"
#include "system/calculate_tangent_operator.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void UpdateTangentOperator(StepParameters& parameters, State<DeviceType>& state) {
    auto region = Kokkos::Profiling::ScopedRegion("Update Tangent Operator");

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;

    Kokkos::parallel_for(
        "CalculateTangentOperator", RangePolicy(0, state.num_system_nodes),
        system::CalculateTangentOperator<DeviceType>{
            parameters.h,
            state.q_delta,
            state.tangent,
        }
    );
}

}  // namespace kynema_fmb::step
