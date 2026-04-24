#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/constraints.hpp"
#include "constraints/update_lambda_prediction.hpp"
#include "solver/solver.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void UpdateConstraintPrediction(
    Solver<DeviceType>& solver, Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("Update Constraint Prediction");
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    auto range_policy = RangePolicy(0, constraints.num_constraints);
    Kokkos::parallel_for(
        "UpdateLambdaPrediction", range_policy,
        constraints::UpdateLambdaPrediction<DeviceType>{
            solver.num_system_dofs, constraints.row_range, solver.x, constraints.lambda
        }
    );
}

}  // namespace kynema_fmb::step
