#pragma once

#include <algorithm>
#include <array>
#include <numeric>
#include <span>
#include <vector>

#include "beam_element.hpp"

namespace kynema_fmb {

/**
 * @brief Represents the input data for creating flexible beams
 *
 * This struct encapsulates all necessary input parameters for instantiating flex beams in
 * kynema i.e. the beam elements and environmental factors such as gravity.
 * It also provides some utilities for computing properties such as
 * - total number of nodes/quadrature points
 * - maximum number of nodes/quadrature points per element etc.
 */
struct BeamsInput {
    std::vector<BeamElement> elements;  //< Elements in the beam
    std::array<double, 3> gravity;      //< Gravity vector

    BeamsInput(std::span<const BeamElement> elems, std::span<const double, 3> g)
        : gravity({g[0], g[1], g[2]}) {
        elements.assign(std::begin(elems), std::end(elems));
    }

    /// Returns the number of elements in the beam
    [[nodiscard]] size_t NumElements() const { return elements.size(); }

    /// Computes the sum of a value across all elements
    template <typename Accessor>
    [[nodiscard]] size_t ComputeSum(Accessor accessor) const {
        return std::transform_reduce(
            elements.begin(), elements.end(), size_t{0U}, std::plus{},
            [&accessor](const BeamElement& e) {
                return accessor(e);
            }
        );
    }

    /// Computes the maximum of a value across all elements
    template <typename Accessor>
    [[nodiscard]] size_t ComputeMax(Accessor accessor) const {
        return std::transform_reduce(
            elements.begin(), elements.end(), size_t{0U},
            [](auto a, auto b) {
                return std::max(a, b);
            },
            [&accessor](const auto& e) {
                return accessor(e);
            }
        );
    }

    /// Returns the total number of nodes in the beam
    [[nodiscard]] size_t NumNodes() const {
        return ComputeSum([](const BeamElement& e) {
            return e.node_ids.size();
        });
    }

    /// Returns the total number of quadrature points in the beam
    [[nodiscard]] size_t NumQuadraturePoints() const {
        return ComputeSum([](const BeamElement& e) {
            return e.quadrature.size();
        });
    }

    /// Returns the maximum number of nodes in any element of the beam
    [[nodiscard]] size_t MaxElemNodes() const {
        return ComputeMax([](const BeamElement& e) {
            return e.node_ids.size();
        });
    }

    /// Returns the maximum number of quadrature points in any element of the beam
    [[nodiscard]] size_t MaxElemQuadraturePoints() const {
        return ComputeMax([](const BeamElement& e) {
            return e.quadrature.size();
        });
    }
};

}  // namespace kynema_fmb
