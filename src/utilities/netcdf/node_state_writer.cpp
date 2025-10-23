#include "node_state_writer.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace kynema::util {
NodeStateWriter::NodeStateWriter(const std::string& file_path, bool create, size_t num_nodes)
    : file_(file_path, create), num_nodes_(num_nodes) {
    // Define dimensions for time and nodes
    // NOTE: We are assuming these dimensions were not added previously
    const int time_dim =
        file_.AddDimension("time", NC_UNLIMITED);  // Unlimited timesteps can be added
    const int node_dim = file_.AddDimension("nodes", num_nodes);
    const std::vector<int> dimensions = {time_dim, node_dim};

    // Define variables for each state component
    this->DefineStateVariables("x", dimensions, true);   // Position (x, y, z, w, i, j, k)
    this->DefineStateVariables("u", dimensions, true);   // Displacement (x, y, z, w, i, j, k)
    this->DefineStateVariables("v", dimensions, false);  // Velocity (x, y, z, i, j, k)
    this->DefineStateVariables("a", dimensions, false);  // Acceleration (x, y, z, i, j, k)
    this->DefineStateVariables("f", dimensions, false);  // Force (x, y, z, i, j, k)

    // Define variables for deformation
    this->DefineStateVariables("deformation", dimensions, true);  // Deformation (x, y, z)
}

void NodeStateWriter::WriteStateDataAtTimestep(
    size_t timestep, const std::string& component_prefix, const std::vector<double>& x,
    const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& i,
    const std::vector<double>& j, const std::vector<double>& k, const std::vector<double>& w
) const {
    // Validate the component prefix - must be one of the valid prefixes:
    // "x" -> position
    // "u" -> displacement
    // "v" -> velocity
    // "a" -> acceleration
    // "f" -> force
    static const std::array<std::string_view, 5> valid_prefixes = {"x", "u", "v", "a", "f"};
    if (std::ranges::none_of(valid_prefixes, [&](const auto& prefix) {
            return prefix == component_prefix;
        })) {
        throw std::invalid_argument("Invalid component prefix: " + component_prefix);
    }

    // Validate vector sizes - must be the same for all components
    const size_t size = x.size();
    if (y.size() != size || z.size() != size || i.size() != size || j.size() != size ||
        k.size() != size) {
        throw std::invalid_argument("All vectors must have the same size");
    }

    // Write data to variables
    const std::vector<size_t> start = {timestep, 0};  // start at the current timestep and node 0
    const std::vector<size_t> count = {
        1, x.size()
    };  // write one timestep worth of data for all nodes
    file_.WriteVariableAt(component_prefix + "_x", start, count, x);
    file_.WriteVariableAt(component_prefix + "_y", start, count, y);
    file_.WriteVariableAt(component_prefix + "_z", start, count, z);
    file_.WriteVariableAt(component_prefix + "_i", start, count, i);
    file_.WriteVariableAt(component_prefix + "_j", start, count, j);
    file_.WriteVariableAt(component_prefix + "_k", start, count, k);

    // Write w component only for position and displacement
    if (!w.empty() && (component_prefix == "x" || component_prefix == "u")) {
        file_.WriteVariableAt(component_prefix + "_w", start, count, w);
    }
}

void NodeStateWriter::WriteDeformationDataAtTimestep(
    size_t timestep, const std::vector<double>& x, const std::vector<double>& y,
    const std::vector<double>& z
) const {
    // Validate vector sizes - must be the same for all components
    const size_t size = x.size();
    if (y.size() != size || z.size() != size) {
        throw std::invalid_argument("All vectors must have the same size");
    }

    // Write data to variables
    const std::vector<size_t> start = {timestep, 0};  // start at the current timestep and node 0
    const std::vector<size_t> count = {
        1, x.size()
    };  // write one timestep worth of data for all nodes
    file_.WriteVariableAt("deformation_x", start, count, x);
    file_.WriteVariableAt("deformation_y", start, count, y);
    file_.WriteVariableAt("deformation_z", start, count, z);
}

const NetCDFFile& NodeStateWriter::GetFile() const {
    return file_;
}

size_t NodeStateWriter::GetNumNodes() const {
    return num_nodes_;
}

void NodeStateWriter::DefineStateVariables(
    const std::string& prefix, const std::vector<int>& dimensions, bool has_w
) {
    (void)file_.AddVariable<double>(prefix + "_x", dimensions);
    (void)file_.AddVariable<double>(prefix + "_y", dimensions);
    (void)file_.AddVariable<double>(prefix + "_z", dimensions);
    (void)file_.AddVariable<double>(prefix + "_i", dimensions);
    (void)file_.AddVariable<double>(prefix + "_j", dimensions);
    (void)file_.AddVariable<double>(prefix + "_k", dimensions);

    if (has_w) {
        (void)file_.AddVariable<double>(prefix + "_w", dimensions);
    }
}
}  // namespace kynema::util
