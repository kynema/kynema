#pragma once

#include <Kokkos_Core.hpp>

namespace kynema_fmb::beams {

/**
 * @brief A  Kernel which interpolates a vector quantity from nodes on a given element
 * to a quadrature point given a basis function
 */
template <typename DeviceType>
struct InterpolateQPVector {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;
    template <typename ValueType>
    using StrideView = Kokkos::View<ValueType, Kokkos::LayoutStride, DeviceType>;
    template <typename ValueType>
    using ConstStrideView = typename StrideView<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double***> shape_interp;
    ConstStrideView<double** [3]> node_vector;
    View<double** [3]> qp_vector;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        // const auto j = first_qp + qp;
        auto local_total = Kokkos::Array<double, 3>{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto phi = shape_interp(element, node, qp);
            for (auto component = 0U; component < 3U; ++component) {
                local_total[component] += node_vector(element, node, component) * phi;
            }
        }
        for (auto component = 0U; component < 3U; ++component) {
            qp_vector(element, qp, component) = local_total[component];
        }
    }
};

}  // namespace kynema_fmb::beams
