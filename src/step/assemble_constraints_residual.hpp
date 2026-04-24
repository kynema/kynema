#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/constraints.hpp"
#include "solver/contribute_constraints_system_residual_to_vector.hpp"
#include "solver/contribute_lambda_to_vector.hpp"
#include "solver/copy_constraints_residual_to_vector.hpp"
#include "solver/solver.hpp"

namespace kynema_fmb::step {

template <typename DeviceType>
inline void AssembleConstraintsResidual(
    Solver<DeviceType>& solver, Constraints<DeviceType>& constraints
) {
    auto resid_region = Kokkos::Profiling::ScopedRegion("Assemble Constraints Residual");

    if (constraints.num_constraints == 0) {
        return;
    }

    auto range_policy =
        Kokkos::RangePolicy<typename DeviceType::execution_space>(0, constraints.num_constraints);

    Kokkos::parallel_for(
        "ContributeConstraintsSystemResidualToVector", range_policy,
        solver::ContributeConstraintsSystemResidualToVector<DeviceType>{
            constraints.target_node_freedom_table, constraints.target_active_dofs,
            constraints.system_residual_terms, solver.b
        }
    );

    Kokkos::parallel_for(
        "ContributeLambdaToVector", range_policy,
        solver::ContributeLambdaToVector<DeviceType>{
            constraints.base_node_freedom_signature, constraints.target_node_freedom_signature,
            constraints.base_node_freedom_table, constraints.target_node_freedom_table,
            constraints.base_lambda_residual_terms, constraints.target_lambda_residual_terms,
            solver.b
        }
    );

    Kokkos::parallel_for(
        "CopyConstraintsResidualToVector", range_policy,
        solver::CopyConstraintsResidualToVector<DeviceType>{
            solver.num_system_dofs, constraints.row_range, constraints.residual_terms, solver.b
        }
    );
}

}  // namespace kynema_fmb::step
