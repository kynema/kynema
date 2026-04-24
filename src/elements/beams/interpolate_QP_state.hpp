#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::beams {

/**
 * @brief Interpolates the displacement (u) part of the state at a quadrature point
 */
template <typename DeviceType>
struct InterpolateQPState_u {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double***> shape_interp;
    ConstView<double** [7]> node_u;
    View<double** [3]> qp_u;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        auto local_total = Kokkos::Array<double, 3>{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto phi = shape_interp(element, node, qp);
            for (auto component = 0U; component < 3U; ++component) {
                local_total[component] += node_u(element, node, component) * phi;
            }
        }
        for (auto component = 0U; component < 3U; ++component) {
            qp_u(element, qp, component) = local_total[component];
        }
    }
};

/**
 * @brief Interpolates the displacement derivative (u') part of the state at a quadrature point
 */
template <typename DeviceType>
struct InterpolateQPState_uprime {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double***> shape_deriv;
    ConstView<double**> qp_jacobian;
    ConstView<double** [7]> node_u;
    View<double** [3]> qp_uprime;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        const auto jacobian = qp_jacobian(element, qp);
        auto local_total = Kokkos::Array<double, 3>{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto dphi = shape_deriv(element, node, qp);
            for (auto component = 0U; component < 3U; ++component) {
                local_total[component] += node_u(element, node, component) * dphi / jacobian;
            }
        }
        for (auto component = 0U; component < 3U; ++component) {
            qp_uprime(element, qp, component) = local_total[component];
        }
    }
};

/**
 * @brief Interpolates the rotation (r) part of the state at a quadrature point
 */
template <typename DeviceType>
struct InterpolateQPState_r {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double***> shape_interp;
    ConstView<double** [7]> node_u;
    View<double** [4]> qp_r;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        auto local_total = Kokkos::Array<double, 4>{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto phi = shape_interp(element, node, qp);
            for (auto component = 0U; component < 4U; ++component) {
                local_total[component] += node_u(element, node, component + 3) * phi;
            }
        }

        for (auto component = 0U; component < 4U; ++component) {
            qp_r(element, qp, component) = math::NormalizeQuaternion(local_total)[component];
        }
    }
};

/**
 * @brief Interpolates the rotation derivative (r') part of the state at a quadrature point
 */
template <typename DeviceType>
struct InterpolateQPState_rprime {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    size_t num_nodes;
    ConstView<double***> shape_deriv;
    ConstView<double**> qp_jacobian;
    ConstView<double** [7]> node_u;
    View<double** [4]> qp_rprime;

    KOKKOS_FUNCTION
    void operator()(size_t qp) const {
        const auto jacobian = qp_jacobian(element, qp);
        auto local_total = Kokkos::Array<double, 4>{};
        for (auto node = 0U; node < num_nodes; ++node) {
            const auto dphi = shape_deriv(element, node, qp);
            for (auto component = 0U; component < 4U; ++component) {
                local_total[component] += node_u(element, node, component + 3) * dphi / jacobian;
            }
        }
        for (auto component = 0U; component < 4U; ++component) {
            qp_rprime(element, qp, component) = local_total[component];
        }
    }
};

}  // namespace kynema_fmb::beams
