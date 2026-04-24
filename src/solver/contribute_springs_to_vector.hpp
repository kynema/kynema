#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the residual contributions computed at each of the nodes in a spring
 * element into the correct location of the global RHS vector.
 */
template <typename DeviceType>
struct ContributeSpringsToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<size_t* [2][3]> element_freedom_table;
    ConstView<double* [2][3]> elements;
    LeftView<double* [1]> vector;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        for (auto component = 0U; component < element_freedom_table.extent(2); ++component) {
            if constexpr (force_atomic) {
                Kokkos::atomic_add(
                    &vector(element_freedom_table(element, 0, component), 0),
                    elements(element, 0, component)
                );
                Kokkos::atomic_add(
                    &vector(element_freedom_table(element, 1, component), 0),
                    elements(element, 1, component)
                );
            } else {
                vector(element_freedom_table(element, 0, component), 0) +=
                    elements(element, 0, component);
                vector(element_freedom_table(element, 1, component), 0) +=
                    elements(element, 1, component);
            }
        }
    };
};

}  // namespace kynema_fmb::solver
