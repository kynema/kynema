#pragma once

#include <Kokkos_Core.hpp>

#include "interpolation_operations.hpp"

namespace kynema_fmb::beams {

/**
 * @brief Functor to calculate Jacobians and unit tangent vectors at quadrature points for beam
 * elements
 *
 * This functor performs two main operations for a provided element:
 * - Calculates the Jacobian (J = |dx/dξ|) at each quadrature point, which represents
 *   the mapping between parametric (ξ) and physical (x) space
 * - Computes normalized position derivatives to obtain unit tangent vectors
 *
 * @note The functor modifies qp_position_derivative in place, normalizing it to create unit tangent
 * vectors
 */
template <typename DeviceType>
struct CalculateJacobian {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;
    ConstView<size_t*> num_qps_per_element;
    ConstView<double***> shape_derivative;           //< num_elems x num_nodes x num_qps
    ConstView<double** [7]> node_position_rotation;  //< num_elems x num_nodes x 7
    View<double** [3]> qp_position_derivative;       //< output: num_elems x num_qps x 3
    View<double**> qp_jacobian;                      //< output: num_elems x num_qps

    KOKKOS_FUNCTION
    void operator()(int element) const {
        using Kokkos::ALL;
        using Kokkos::make_pair;
        using Kokkos::pow;
        using Kokkos::sqrt;
        using Kokkos::subview;

        const auto num_nodes = num_nodes_per_element(element);
        const auto num_qps = num_qps_per_element(element);
        const auto shape_deriv = subview(
            shape_derivative, element, make_pair(size_t{0U}, num_nodes),
            make_pair(size_t{0U}, num_qps)
        );
        const auto qp_pos_deriv =
            subview(qp_position_derivative, element, make_pair(size_t{0U}, num_qps), ALL);
        const auto node_pos = subview(
            node_position_rotation, element, make_pair(size_t{0U}, num_nodes), make_pair(0U, 3U)
        );
        const auto qp_jacob = subview(qp_jacobian, element, make_pair(size_t{0U}, num_qps));

        // Interpolate position derivatives at quadrature points using shape functions
        // qp_pos_deriv = Σ(dN/dξ * node_pos)
        InterpVector3(shape_deriv, node_pos, qp_pos_deriv);

        for (auto qp = 0U; qp < num_qps; ++qp) {
            // Calculate Jacobian - this is a scalar that represents the "stretching" factor between
            // parametric (ξ) and physical (x) space
            // J = |dx/dξ| = sqrt((dx/dξ)² + (dy/dξ)² + (dz/dξ)²)
            const auto jacobian = sqrt(
                pow(qp_pos_deriv(qp, 0), 2.) + pow(qp_pos_deriv(qp, 1), 2.) +
                pow(qp_pos_deriv(qp, 2), 2.)
            );
            qp_jacob(qp) = jacobian;

            // Normalize position derivatives by Jacobian to get unit tangent vector that points
            // in the direction of curve/beam
            for (auto component = 0U; component < 3U; ++component) {
                qp_pos_deriv(qp, component) /= jacobian;
            }
        }
    }
};

}  // namespace kynema_fmb::beams
