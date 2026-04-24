#pragma once
#include <span>

#include <Kokkos_Core.hpp>

#include "beams.hpp"
#include "beams_input.hpp"
#include "calculate_jacobian.hpp"
#include "interpolate_QP_position.hpp"
#include "interpolate_QP_rotation.hpp"
#include "interpolate_to_quadrature_points.hpp"
#include "populate_element_views.hpp"

namespace kynema_fmb {

/**
 * @brief Creates a beams data structure and initializes its data
 *
 * @tparam DeviceData The Kokkos Device where the newly created Beams structure will reside
 *
 * @param beams_input A BeamsInput object defining the beam elements
 * @param nodes A vector defining all of the nodes in the problem
 * @return A fully initialized Beams data structure
 *
 */
template <typename DeviceType>
inline Beams<DeviceType> CreateBeams(const BeamsInput& beams_input, std::span<const Node> nodes) {
    using Kokkos::ALL;
    using Kokkos::create_mirror_view;
    using Kokkos::deep_copy;
    using Kokkos::make_pair;
    using Kokkos::parallel_for;
    using Kokkos::subview;
    using Kokkos::WithoutInitializing;
    using RangePolicy = Kokkos::RangePolicy<typename DeviceType::execution_space>;
    using TeamPolicy = Kokkos::TeamPolicy<typename DeviceType::execution_space>;

    Beams<DeviceType> beams(
        beams_input.NumElements(), beams_input.MaxElemNodes(), beams_input.MaxElemQuadraturePoints()
    );

    auto host_gravity = create_mirror_view(WithoutInitializing, beams.gravity);

    auto host_num_nodes_per_element =
        create_mirror_view(WithoutInitializing, beams.num_nodes_per_element);
    auto host_num_qps_per_element =
        create_mirror_view(WithoutInitializing, beams.num_qps_per_element);
    auto host_element_mu = create_mirror_view(WithoutInitializing, beams.element_mu);
    auto host_node_state_indices = create_mirror_view(WithoutInitializing, beams.node_state_indices);
    auto host_node_x0 = create_mirror_view(WithoutInitializing, beams.node_x0);
    auto host_node_u = create_mirror_view(WithoutInitializing, beams.node_u);
    auto host_node_u_dot = create_mirror_view(WithoutInitializing, beams.node_u_dot);
    auto host_node_u_ddot = create_mirror_view(WithoutInitializing, beams.node_u_ddot);

    auto host_qp_weight = create_mirror_view(WithoutInitializing, beams.qp_weight);
    auto host_qp_Mstar = create_mirror_view(WithoutInitializing, beams.qp_Mstar);
    auto host_qp_Cstar = create_mirror_view(WithoutInitializing, beams.qp_Cstar);

    auto host_shape_interp = create_mirror_view(WithoutInitializing, beams.shape_interp);
    auto host_shape_deriv = create_mirror_view(WithoutInitializing, beams.shape_deriv);

    host_gravity(0) = beams_input.gravity[0];
    host_gravity(1) = beams_input.gravity[1];
    host_gravity(2) = beams_input.gravity[2];

    for (auto element : std::views::iota(0U, beams_input.NumElements())) {
        // Get number of nodes and quadrature points in element
        const auto num_nodes = beams_input.elements[element].node_ids.size();
        const auto num_qps = beams_input.elements[element].quadrature.size();

        // Create element indices and set in host mirror
        host_num_nodes_per_element(element) = num_nodes;
        host_num_qps_per_element(element) = num_qps;

        // Set element damping coefficients
        for (auto i : std::views::iota(0U, 6U)) {
            host_element_mu(element, i) = beams_input.elements[element].mu[i];
        }

        // Populate beam node->state indices
        for (auto node : std::views::iota(0U, num_nodes)) {
            host_node_state_indices(element, node) = beams_input.elements[element].node_ids[node];
        }

        // Populate views for this element
        beams::PopulateNodeX0(
            beams_input.elements[element], nodes,
            subview(host_node_x0, element, make_pair(0UL, num_nodes), ALL)
        );
        beams::PopulateQPWeight(
            beams_input.elements[element], subview(host_qp_weight, element, make_pair(0UL, num_qps))
        );
        beams::PopulateShapeFunctionValues(
            beams_input.elements[element], nodes,
            subview(host_shape_interp, element, make_pair(0UL, num_nodes), make_pair(0UL, num_qps))
        );
        beams::PopulateShapeFunctionDerivatives(
            beams_input.elements[element], nodes,
            subview(host_shape_deriv, element, make_pair(0UL, num_nodes), make_pair(0UL, num_qps))
        );
        beams::PopulateQPMstar(
            beams_input.elements[element],
            subview(host_qp_Mstar, element, make_pair(0UL, num_qps), ALL, ALL)
        );
        beams::PopulateQPCstar(
            beams_input.elements[element],
            subview(host_qp_Cstar, element, make_pair(0UL, num_qps), ALL, ALL)
        );
    }

    deep_copy(beams.num_nodes_per_element, host_num_nodes_per_element);
    deep_copy(beams.num_qps_per_element, host_num_qps_per_element);
    deep_copy(beams.element_mu, host_element_mu);
    deep_copy(beams.gravity, host_gravity);
    deep_copy(beams.node_state_indices, host_node_state_indices);
    deep_copy(beams.node_x0, host_node_x0);
    deep_copy(beams.node_u, host_node_u);
    deep_copy(beams.node_u_dot, host_node_u_dot);
    deep_copy(beams.node_u_ddot, host_node_u_ddot);
    deep_copy(beams.qp_weight, host_qp_weight);
    deep_copy(beams.qp_Mstar, host_qp_Mstar);
    deep_copy(beams.qp_Cstar, host_qp_Cstar);
    deep_copy(beams.shape_interp, host_shape_interp);
    deep_copy(beams.shape_deriv, host_shape_deriv);

    deep_copy(beams.node_FX, 0.);
    deep_copy(beams.qp_Fe, 0.);

    auto range_policy = RangePolicy(0, beams.num_elems);

    parallel_for(
        "InterpolateQPPosition", range_policy,
        beams::InterpolateQPPosition<DeviceType>{
            beams.num_nodes_per_element, beams.num_qps_per_element, beams.shape_interp,
            beams.node_x0, beams.qp_x0
        }
    );

    parallel_for(
        "InterpolateQPRotation", range_policy,
        beams::InterpolateQPRotation<DeviceType>{
            beams.num_nodes_per_element, beams.num_qps_per_element, beams.shape_interp,
            beams.node_x0, beams.qp_r0
        }
    );

    parallel_for(
        "CalculateJacobian", range_policy,
        beams::CalculateJacobian<DeviceType>{
            beams.num_nodes_per_element,
            beams.num_qps_per_element,
            beams.shape_deriv,
            beams.node_x0,
            beams.qp_x0_prime,
            beams.qp_jacobian,
        }
    );

    auto team_policy = TeamPolicy(static_cast<int>(beams.num_elems), Kokkos::AUTO());

    parallel_for(
        "InterpolateToQuadraturePoints", team_policy,
        beams::InterpolateToQuadraturePoints<DeviceType>{
            beams.num_nodes_per_element, beams.num_qps_per_element, beams.shape_interp,
            beams.shape_deriv, beams.qp_jacobian, beams.node_u, beams.node_u_dot, beams.node_u_ddot,
            beams.qp_x0, beams.qp_r0, beams.qp_u, beams.qp_u_prime, beams.qp_r, beams.qp_r_prime,
            beams.qp_u_dot, beams.qp_omega, beams.qp_u_ddot, beams.qp_omega_dot, beams.qp_x
        }
    );

    return beams;
}

}  // namespace kynema_fmb
