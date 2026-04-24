#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::solver {

/**
 * @brief A Kernel which sums the residual contributions computed at a mass element's node
 * into the correct location of the global RHS vector.
 */
template <typename DeviceType>
struct ContributeMassesToVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using LeftView = Kokkos::View<ValueType, Kokkos::LayoutLeft, DeviceType>;

    ConstView<size_t* [6]> element_freedom_table;
    ConstView<double* [6]> elements;
    LeftView<double* [1]> vector;

    KOKKOS_FUNCTION
    void operator()(size_t element) const {
        constexpr auto force_atomic =
            !std::is_same_v<typename DeviceType::execution_space, Kokkos::Serial>;
        for (auto component = 0U; component < element_freedom_table.extent(1); ++component) {
            if constexpr (force_atomic) {
                Kokkos::atomic_add(
                    &vector(element_freedom_table(element, component), 0),
                    elements(element, component)
                );
            } else {
                vector(element_freedom_table(element, component), 0) += elements(element, component);
            }
        }
    }
};

}  // namespace kynema_fmb::solver
