#pragma once

#include <array>
#include <cmath>

namespace kynema_fmb {

/**
 * @brief Spring element represents a constitutively linear spring connecting two nodes and defined
 * by its scalar stiffness and undeformed length
 */
struct SpringElement {
    size_t ID;                       // Element identifier
    std::array<size_t, 2> node_ids;  // 2 node IDs (start and end points of spring)
    double stiffness;                // Spring stiffness coefficient
    double undeformed_length;        // Reference/undeformed length of spring

    SpringElement(size_t id, std::array<size_t, 2> n_ids, double k, double l0)
        : ID(id), node_ids(n_ids), stiffness(k), undeformed_length(l0) {}
};

}  // namespace kynema_fmb
