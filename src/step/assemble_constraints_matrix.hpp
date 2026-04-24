#pragma once

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "constraints/constraints.hpp"
#include "solver/copy_constraints_to_sparse_matrix.hpp"
#include "solver/copy_constraints_transpose_to_sparse_matrix.hpp"
#include "solver/solver.hpp"

namespace kynema_fmb::step {
template <typename DeviceType>
inline void AssembleConstraintsMatrix(
    Solver<DeviceType>& solver, Constraints<DeviceType>& constraints
) {
    auto region = Kokkos::Profiling::ScopedRegion("Assemble Constraints Matrix");

    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;

    auto constraint_policy =
        TeamPolicy(static_cast<int>(constraints.num_constraints), Kokkos::AUTO());

    Kokkos::parallel_for(
        "CopyConstraintsToSparseMatrix", constraint_policy,
        solver::CopyConstraintsToSparseMatrix<typename Solver<DeviceType>::CrsMatrixType>{
            solver.num_system_dofs, constraints.row_range, constraints.base_node_freedom_signature,
            constraints.target_node_freedom_signature, constraints.base_node_freedom_table,
            constraints.target_node_freedom_table, constraints.base_gradient_terms,
            constraints.target_gradient_terms, solver.A
        }
    );

    Kokkos::parallel_for(
        "CopyConstraintsTransposeToSparseMatrix", constraint_policy,
        solver::CopyConstraintsTransposeToSparseMatrix<typename Solver<DeviceType>::CrsMatrixType>{
            solver.num_system_dofs, constraints.row_range, constraints.base_node_freedom_signature,
            constraints.target_node_freedom_signature, constraints.base_node_freedom_table,
            constraints.target_node_freedom_table, constraints.base_gradient_transpose_terms,
            constraints.target_gradient_transpose_terms, solver.A
        }
    );
}
}  // namespace kynema_fmb::step
