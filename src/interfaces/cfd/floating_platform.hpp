#pragma once

#include <cstddef>
#include <vector>

#include "interfaces/cfd/mooring_line.hpp"
#include "interfaces/cfd/node_data.hpp"

namespace kynema_fmb::interfaces::cfd {

//------------------------------------------------------------------------------
// Data Structures
//------------------------------------------------------------------------------

/**
 * @brief A structure describing a point mass and its attached mooring lines
 */
struct FloatingPlatform {
    /// @brief Flag indicating platform is active in model
    bool active = false;

    /// @brief Platform node data (ID, motion, loads)
    NodeData node;

    /// @brief Platform mass element identifier
    size_t mass_element_id;

    /// @brief Mooring line data array
    std::vector<MooringLine> mooring_lines;
};

}  // namespace kynema_fmb::interfaces::cfd
