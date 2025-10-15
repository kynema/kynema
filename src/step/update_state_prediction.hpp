#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "solver/solver.hpp"
#include "state/calculate_displacement.hpp"
#include "state/state.hpp"
#include "state/update_dynamic_prediction.hpp"
#include "state/update_static_prediction.hpp"
#include "step_parameters.hpp"

namespace kynema::step {

/**
 * @brief Updates the predicted next state values, based on computed solver solution, solver.x
 *
 * @tparam DeviceType The Kokkos Device where the solver and state data structures reside
 *
 * @param parameters A struct containing the control parameters for time stepping
 * @param solver A struct containing the linear system, already solved for this iteration
 * @param state A struct containing the state information at each node
 */
template <typename DeviceType>
inline void UpdateStatePrediction(
    StepParameters& parameters, const Solver<DeviceType>& solver, State<DeviceType>& state
) {
    auto region = Kokkos::Profiling::ScopedRegion("Update State Prediction");
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto range_policy = RangePolicy(0, solver.num_system_nodes);

    // Update the state prediction for a DYNAMIC solve
    if (parameters.is_dynamic_solve) {
        Kokkos::parallel_for(
            "UpdateDynamicPrediction", range_policy,
            state::UpdateDynamicPrediction<DeviceType>{
                parameters.h,
                parameters.beta_prime,
                parameters.gamma_prime,
                state.node_freedom_allocation_table,
                state.node_freedom_map_table,
                solver.x,
                state.q_delta,
                state.v,
                state.vd,
            }
        );
    } else {
        // Update the state prediction for a STATIC solve
        Kokkos::parallel_for(
            "UpdateStaticPrediction", range_policy,
            state::UpdateStaticPrediction<DeviceType>{
                parameters.h,
                state.node_freedom_allocation_table,
                state.node_freedom_map_table,
                solver.x,
                state.q_delta,
            }
        );
    }

    Kokkos::parallel_for(
        "CalculateDisplacement", range_policy,
        state::CalculateDisplacement<DeviceType>{
            parameters.h,
            state.q_delta,
            state.q_prev,
            state.q,
        }
    );
}

}  // namespace kynema::step
