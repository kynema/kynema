#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which copies the residual contributions for each constraint to their correct
 * location in the global RHS vector.
 */
template <typename DeviceType>
struct CopyConstraintsResidualToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    size_t start_row;
    ConstView<Kokkos::pair<size_t, size_t>*> row_range;
    ConstView<double* [6]> constraint_residual_terms;
    LeftView<double* [1]> residual;

    KOKKOS_FUNCTION
    void operator()(size_t constraint) const {
        const auto first_row = row_range(constraint).first + start_row;
        const auto num_rows = row_range(constraint).second - row_range(constraint).first;
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        for (auto component = 0U; component < num_rows; ++component) {
            if constexpr (force_atomic) {
                Kokkos::atomic_add(
                    &residual(first_row + component, 0),
                    constraint_residual_terms(constraint, component)
                );
            } else {
                residual(first_row + component, 0) =
                    constraint_residual_terms(constraint, component);
            }
        }
    }
};
}  // namespace kynema_fmb::solver
