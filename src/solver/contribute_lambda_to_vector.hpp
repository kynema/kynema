#pragma once

#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb::solver {

/**
 * @brief A kernel which contributes the constraint Lagrange multiplier terms to the correct
 * locations in the global RHS vector.
 */
template <typename DeviceType>
struct ContributeLambdaToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<dof::FreedomSignature*> base_node_freedom_signature;
    ConstView<dof::FreedomSignature*> target_node_freedom_signature;

    ConstView<size_t* [6]> base_node_freedom_table;
    ConstView<size_t* [6]> target_node_freedom_table;

    ConstView<double* [6]> base_lambda_residual_terms;
    ConstView<double* [6]> target_lambda_residual_terms;

    LeftView<double* [1]> R;

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto base_num_dofs = count_active_dofs(base_node_freedom_signature(constraint));
        const auto first_base_dof = base_node_freedom_table(constraint, 0);
        for (auto component = 0U; component < base_num_dofs; ++component) {
            Kokkos::atomic_add(
                &R(first_base_dof + component, 0), base_lambda_residual_terms(constraint, component)
            );
        }

        const auto target_num_dofs = count_active_dofs(target_node_freedom_signature(constraint));
        const auto first_target_dof = target_node_freedom_table(constraint, 0);
        for (auto component = 0U; component < target_num_dofs; ++component) {
            Kokkos::atomic_add(
                &R(first_target_dof + component, 0),
                target_lambda_residual_terms(constraint, component)
            );
        }
    }
};

}  // namespace kynema_fmb::solver
