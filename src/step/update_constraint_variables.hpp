#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/calculate_constraint_residual_gradient.hpp"
#include "constraints/constraints.hpp"
#include "state/state.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void UpdateConstraintVariables(
    State<DeviceType>& state, Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("Update Constraint Variables");

    if (constraints.num_constraints == 0) {
        return;
    }

    constraints.UpdateViews();

    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto range_policy = RangePolicy(0, constraints.num_constraints);
    Kokkos::parallel_for(
        "CalculateConstraintResidualGradient", range_policy,
        constraints::CalculateConstraintResidualGradient<DeviceType>{
            constraints.type, constraints.base_node_index, constraints.target_node_index,
            constraints.X0, constraints.axes, constraints.input, constraints.lambda, state.tangent,
            state.q, constraints.residual_terms, constraints.base_lambda_residual_terms,
            constraints.target_lambda_residual_terms, constraints.system_residual_terms,
            constraints.base_gradient_terms, constraints.target_gradient_terms,
            constraints.base_gradient_transpose_terms, constraints.target_gradient_transpose_terms
        }
    );
}

}  // namespace kynema_fmb::step
