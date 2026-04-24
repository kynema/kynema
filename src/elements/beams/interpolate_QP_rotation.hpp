#pragma once

#include <Kokkos_Core.hpp>

#include "interpolation_operations.hpp"

namespace kynema_fmb::beams {

/**
 * @brief A Kernel which interpolates a rotation quaternion on a given element from its nodes
 * to all of it quadrature points.
 */
template <typename DeviceType>
struct InterpolateQPRotation {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t*> num_qps_per_element;
    ConstView<double***> shape_interpolation;        // Num Nodes x Num Quadrature points
    ConstView<double** [7]> node_position_rotation;  // Node global position vector
    View<double** [4]> qp_rotation;                  // quadrature point rotation

    KOKKOS_FUNCTION
    void operator()(int element) const {
        using Kokkos::ALL;
        using Kokkos::make_pair;
        using Kokkos::subview;

        const auto num_nodes = num_nodes_per_element(element);
        const auto num_qps = num_qps_per_element(element);
        const auto shape_interp = subview(
            shape_interpolation, element, make_pair(size_t{0U}, num_nodes),
            make_pair(size_t{0U}, num_qps)
        );
        const auto node_rot = subview(
            node_position_rotation, element, make_pair(size_t{0U}, num_nodes), make_pair(3, 7)
        );
        const auto qp_rot = subview(qp_rotation, element, make_pair(size_t{0U}, num_qps), ALL);

        InterpQuaternion(shape_interp, node_rot, qp_rot);
    }
};

}  // namespace kynema_fmb::beams
