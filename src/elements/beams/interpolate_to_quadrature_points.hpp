#pragma once

#include <Kokkos_Core.hpp>

#include "calculate_QP_position.hpp"
#include "interpolate_QP_state.hpp"
#include "interpolate_QP_vector.hpp"

namespace kynema_fmb::beams {

/**
 * @brief Interpolates various quantities from nodes to quadrature points for beam elements
 *
 * This functor handles the interpolation of multiple quantities from nodal points to
 * quadrature points (QPs) for beam elements, including:
 * - Displacements (u) and their derivatives
 * - Rotations (r) and their derivatives
 * - Velocities (u_dot) and angular velocities (omega)
 * - Accelerations (u_ddot) and angular accelerations (omega_dot)
 * - Final positions (x) of quadrature points
 */
template <typename DeviceType>
struct InterpolateToQuadraturePoints {
    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;
    using member_type = typename TeamPolicy::member_type;
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;
    template <typename ValueType>
    using ConstView = typename View<ValueType>::const_type;

    ConstView<size_t*> num_nodes_per_element;  //< Number of Nodes
    ConstView<size_t*> num_qps_per_element;    //< Number of QPs per element
    ConstView<double***> shape_interp;         //< shape functions at QPs
    ConstView<double***> shape_deriv;          //< shape function derivatives at QPs
    ConstView<double**> qp_jacobian;           //< Jacobian at QPs
    ConstView<double** [7]> node_u;            //< Nodal displacements
    ConstView<double** [6]> node_u_dot;        //< Nodal velocities
    ConstView<double** [6]> node_u_ddot;       //< Nodal accelerations
    ConstView<double** [3]> qp_x0;             //< Initial positions at QPs
    ConstView<double** [4]> qp_r0;             //< Initial rotations at QPs

    // Output quantities at quadrature points
    View<double** [3]> qp_u;          //< Interpolated displacements at QPs
    View<double** [3]> qp_uprime;     //< Displacement derivatives at QPs
    View<double** [4]> qp_r;          //< Interpolated rotations at QPs
    View<double** [4]> qp_rprime;     //< Rotation derivatives at QPs
    View<double** [3]> qp_u_dot;      //< Interpolated velocities at QPs
    View<double** [3]> qp_omega;      //< Interpolated angular velocities at QPs
    View<double** [3]> qp_u_ddot;     //< Interpolated accelerations at QPs
    View<double** [3]> qp_omega_dot;  //< Interpolated angular accelerations at QPs
    View<double** [7]> qp_x;          //< Final positions of quadrature points

    KOKKOS_FUNCTION
    void operator()(member_type member) const {
        using Kokkos::ALL;
        using Kokkos::make_pair;
        using Kokkos::parallel_for;

        const auto element = static_cast<size_t>(member.league_rank());
        const auto num_nodes = num_nodes_per_element(element);
        const auto num_qps = num_qps_per_element(element);

        auto qp_range = Kokkos::TeamThreadRange(member, num_qps);

        parallel_for(
            qp_range,
            InterpolateQPState_u<DeviceType>{element, num_nodes, shape_interp, node_u, qp_u}
        );
        parallel_for(
            qp_range, InterpolateQPState_uprime<DeviceType>{
                          element, num_nodes, shape_deriv, qp_jacobian, node_u, qp_uprime
                      }
        );
        parallel_for(
            qp_range,
            InterpolateQPState_r<DeviceType>{element, num_nodes, shape_interp, node_u, qp_r}
        );
        parallel_for(
            qp_range, InterpolateQPState_rprime<DeviceType>{
                          element, num_nodes, shape_deriv, qp_jacobian, node_u, qp_rprime
                      }
        );
        parallel_for(
            qp_range, InterpolateQPVector<DeviceType>{
                          element, num_nodes, shape_interp,
                          subview(node_u_dot, ALL, ALL, make_pair(0, 3)), qp_u_dot
                      }
        );
        parallel_for(
            qp_range, InterpolateQPVector<DeviceType>{
                          element, num_nodes, shape_interp,
                          subview(node_u_dot, ALL, ALL, make_pair(3, 6)), qp_omega
                      }
        );
        parallel_for(
            qp_range, InterpolateQPVector<DeviceType>{
                          element, num_nodes, shape_interp,
                          subview(node_u_ddot, ALL, ALL, make_pair(0, 3)), qp_u_ddot
                      }
        );
        parallel_for(
            qp_range, InterpolateQPVector<DeviceType>{
                          element, num_nodes, shape_interp,
                          subview(node_u_ddot, ALL, ALL, make_pair(3, 6)), qp_omega_dot
                      }
        );
        parallel_for(
            qp_range, CalculateQPPosition<DeviceType>{element, qp_x0, qp_u, qp_r0, qp_r, qp_x}
        );
    }
};

}  // namespace kynema_fmb::beams
