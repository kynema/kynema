#pragma once

#include <string>
#include <vector>

#include "utilities/netcdf/node_state_writer.hpp"

namespace kynema::interfaces::components {

/**
 * @brief A configuration object used to create the Outputs object
 */
struct OutputsConfig {
    /// @brief Output file path for NetCDF results (empty = no outputs will be written)
    std::string output_file_path;

    /// @brief Vector of state component prefixes to enable for writing
    std::vector<std::string> output_state_prefixes{"x", "u", "v", "a", "f"};

    /// @brief Number of timesteps to buffer before auto-flush
    size_t buffer_size{util::NodeStateWriter::kDefaultBufferSize};

    /// @brief Check if outputs are enabled
    [[nodiscard]] bool Enabled() const { return !this->output_file_path.empty(); }
};

}  // namespace kynema::interfaces::components
