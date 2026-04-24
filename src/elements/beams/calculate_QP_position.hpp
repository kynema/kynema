#pragma once

#include <Kokkos_Core.hpp>

#include "math/quaternion_operations.hpp"

namespace kynema_fmb::beams {

/**
 * @brief Functor to calculate current position and orientation at quadrature points
 *
 * This functor updates the current position and orientation (qp_x_) of quadrature points
 * by combining:
 * - Initial position (qp_x0_) with displacement (qp_u_) for translational components
 * - Initial orientation (qp_r0_) with rotation (qp_r_) for rotational components using quaternions
 *
 * @param element Element index
 * @param qp_x0_ Initial quadrature point positions (num_elems x num_qps x 3)
 * @param qp_u_ Displacements at quadrature points (num_elems x num_qps x 3)
 * @param qp_r0_ Initial orientations as quaternions (num_elems x num_qps x 4)
 * @param qp_r_ Rotations as quaternions (num_elems x num_qps x 4)
 * @param qp_x_ Output current positions and orientations (num_elems x num_qps x 7)
 *              where [0:3] = position, [3:7] = orientation quaternion
 */
template <typename DeviceType>
struct CalculateQPPosition {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    size_t element;
    ConstView<double** [3]> qp_x0_;
    ConstView<double** [3]> qp_u_;
    ConstView<double** [4]> qp_r0_;
    ConstView<double** [4]> qp_r_;
    View<double** [7]> qp_x_;

    KOKKOS_FUNCTION void operator()(int qp) const {
        using Kokkos::ALL;
        using Kokkos::subview;

        // Calculate current position
        qp_x_(element, qp, 0) = qp_x0_(element, qp, 0) + qp_u_(element, qp, 0);
        qp_x_(element, qp, 1) = qp_x0_(element, qp, 1) + qp_u_(element, qp, 1);
        qp_x_(element, qp, 2) = qp_x0_(element, qp, 2) + qp_u_(element, qp, 2);

        // Calculate current orientation
        auto RR0_data = Kokkos::Array<double, 4>{};
        auto RR0 = View<double[4]>(RR0_data.data());
        math::QuaternionCompose(
            subview(qp_r_, element, qp, ALL), subview(qp_r0_, element, qp, ALL), RR0
        );
        qp_x_(element, qp, 3) = RR0(0);
        qp_x_(element, qp, 4) = RR0(1);
        qp_x_(element, qp, 5) = RR0(2);
        qp_x_(element, qp, 6) = RR0(3);
    }
};

}  // namespace kynema_fmb::beams
