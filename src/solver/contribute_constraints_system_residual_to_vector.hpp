#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the system residual contributions for a constraint's target
 * node into the correct location of the global RHS vector.
 */
template <typename DeviceType>
struct ContributeConstraintsSystemResidualToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<size_t* [6]> target_node_freedom_table;
    ConstView<size_t*> target_active_dofs;
    ConstView<double* [6]> system_residual_terms;
    LeftView<double* [1]> residual;

    KOKKOS_FUNCTION
    void operator()(size_t i_constraint) const {
        const auto num_dofs = target_active_dofs(i_constraint);
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        for (auto component = 0U; component < num_dofs; ++component) {
            if constexpr (force_atomic) {
                Kokkos::atomic_add(
                    &residual(target_node_freedom_table(i_constraint, component), 0),
                    system_residual_terms(i_constraint, component)
                );
            } else {
                residual(target_node_freedom_table(i_constraint, component), 0) +=
                    system_residual_terms(i_constraint, component);
            }
        }
    }
};

}  // namespace kynema_fmb::solver
