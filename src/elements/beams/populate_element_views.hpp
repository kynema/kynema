#pragma once

#include <ranges>
#include <span>
#include <vector>

#include "beam_element.hpp"
#include "math/interpolation.hpp"
#include "model/node.hpp"

namespace kynema::beams {

/// @brief Populate the node initial position and orientation
inline void PopulateNodeX0(
    const BeamElement& elem, std::span<const Node> nodes,
    const Kokkos::View<double* [7], Kokkos::LayoutStride, Kokkos::HostSpace>& node_x0
) {
    for (auto node : std::views::iota(0U, elem.node_ids.size())) {
        for (auto component : std::views::iota(0U, 7U)) {
            node_x0(node, component) = nodes[elem.node_ids[node]].x0[component];
        }
    }
}

/// @brief Populate the integration weights at each quadrature point
inline void PopulateQPWeight(
    const BeamElement& elem,
    const Kokkos::View<double*, Kokkos::LayoutStride, Kokkos::HostSpace>& qp_weight
) {
    for (auto qp : std::views::iota(0U, elem.quadrature.size())) {
        qp_weight(qp) = elem.quadrature[qp][1];
    }
}

/// @brief Map node positions from [0,1] to [-1,1]
inline std::vector<double> MapNodePositions(const BeamElement& elem, std::span<const Node> nodes) {
    std::vector<double> node_xi(elem.node_ids.size());
    std::ranges::transform(elem.node_ids, std::begin(node_xi), [&](auto node_id) {
        return (2. * nodes[node_id].s) - 1.;
    });
    return node_xi;
}

/// @brief Populate shape function values at each quadrature point
inline void PopulateShapeFunctionValues(
    const BeamElement& elem, std::span<const Node> nodes,
    const Kokkos::View<double**, Kokkos::LayoutStride, Kokkos::HostSpace>& shape_interp
) {
    const auto node_xi = MapNodePositions(elem, nodes);

    std::vector<double> weights;
    for (auto qp : std::views::iota(0U, elem.quadrature.size())) {
        auto qp_xi = elem.quadrature[qp][0];

        math::LagrangePolynomialInterpWeights(qp_xi, node_xi, weights);
        for (auto node : std::views::iota(0U, node_xi.size())) {
            shape_interp(node, qp) = weights[node];
        }
    }
}

/// @brief Populate shape function derivatives at each quadrature point
inline void PopulateShapeFunctionDerivatives(
    const BeamElement& elem, std::span<const Node> nodes,
    const Kokkos::View<double**, Kokkos::LayoutStride, Kokkos::HostSpace>& shape_deriv
) {
    const auto node_xi = MapNodePositions(elem, nodes);

    std::vector<double> weights;
    for (auto qp : std::views::iota(0U, elem.quadrature.size())) {
        auto qp_xi = elem.quadrature[qp][0];

        math::LagrangePolynomialDerivWeights(qp_xi, node_xi, weights);
        for (auto node : std::views::iota(0U, node_xi.size())) {
            shape_deriv(node, qp) = weights[node];
        }
    }
}

/// @brief Map section positions from [0,1] to [-1,1]
inline std::vector<double> MapSectionPositions(const BeamElement& elem) {
    std::vector<double> section_xi(elem.sections.size());
    std::ranges::transform(elem.sections, std::begin(section_xi), [](const auto& section) {
        return (2. * section.position) - 1.;
    });
    return section_xi;
}

/// @brief Populate mass matrix values at each quadrature point
inline void PopulateQPMstar(
    const BeamElement& elem,
    const Kokkos::View<double* [6][6], Kokkos::LayoutStride, Kokkos::HostSpace>& qp_Mstar
) {
    const auto section_xi = MapSectionPositions(elem);
    std::vector<double> weights(elem.sections.size());

    for (auto qp : std::views::iota(0U, elem.quadrature.size())) {
        auto qp_xi = elem.quadrature[qp][0];
        math::LinearInterpWeights(qp_xi, section_xi, weights);
        for (auto component_1 : std::views::iota(0U, 6U)) {
            for (auto component_2 : std::views::iota(0U, 6U)) {
                qp_Mstar(qp, component_1, component_2) = 0.;
            }
        }
        for (auto section : std::views::iota(0U, section_xi.size())) {
            for (auto component_1 : std::views::iota(0U, 6U)) {
                for (auto component_2 : std::views::iota(0U, 6U)) {
                    qp_Mstar(qp, component_1, component_2) +=
                        elem.sections[section].M_star[component_1][component_2] * weights[section];
                }
            }
        }
    }
}

/// @brief Populate stiffness matrix values at each quadrature point
inline void PopulateQPCstar(
    const BeamElement& elem,
    const Kokkos::View<double* [6][6], Kokkos::LayoutStride, Kokkos::HostSpace>& qp_Cstar
) {
    const auto section_xi = MapSectionPositions(elem);
    std::vector<double> weights(elem.sections.size());

    for (auto qp : std::views::iota(0U, elem.quadrature.size())) {
        auto qp_xi = elem.quadrature[qp][0];
        math::LinearInterpWeights(qp_xi, section_xi, weights);
        for (auto component_1 : std::views::iota(0U, 6U)) {
            for (auto component_2 : std::views::iota(0U, 6U)) {
                qp_Cstar(qp, component_1, component_2) = 0.;
            }
        }
        for (auto section : std::views::iota(0U, section_xi.size())) {
            for (auto component_1 : std::views::iota(0U, 6U)) {
                for (auto component_2 : std::views::iota(0U, 6U)) {
                    qp_Cstar(qp, component_1, component_2) +=
                        elem.sections[section].C_star[component_1][component_2] * weights[section];
                }
            }
        }
    }
}
}  // namespace kynema::beams
