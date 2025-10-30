#include "node_state_writer.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <unordered_map>

namespace kynema::util {

NodeStateWriter::NodeStateWriter(
    const std::string& file_path, bool create, size_t num_nodes,
    const std::vector<std::string>& enabled_state_prefixes, size_t buffer_size
)
    : file_(file_path, create),
      num_nodes_(num_nodes),
      enabled_state_prefixes_(enabled_state_prefixes),
      buffer_size_(buffer_size) {
    //-----------------------------------
    // define variables
    //-----------------------------------

    // Define dimensions for time and nodes
    // NOTE: We are assuming these dimensions were not added previously
    const int time_dim =
        file_.AddDimension("time", NC_UNLIMITED);  // Unlimited timesteps can be added
    const int node_dim = file_.AddDimension("nodes", num_nodes);
    const std::vector<int> dimensions = {time_dim, node_dim};

    // Define variables for each state component
    for (const auto& prefix : enabled_state_prefixes_) {
        const bool has_w = (prefix == "x" || prefix == "u");  // position and displacement have w
        this->DefineStateVariables(prefix, dimensions, has_w);
    }

    //-----------------------------------
    // write in chunks
    //-----------------------------------

    // Set chunking for the state variables
    // position -> x, displacement -> u, velocity -> v, acceleration -> a, force -> f
    const size_t chunk_size =
        std::max(this->buffer_size_, size_t(1));  // minimum chunk size is 1 timestep
    const std::array<size_t, 2> chunking = {chunk_size, this->num_nodes_};
    const std::array<std::string_view, 7> components = {"x", "y", "z", "i", "j", "k", "w"};

    for (const auto& prefix : this->enabled_state_prefixes_) {
        for (const auto& comp : components) {
            const std::string var_name = std::string(prefix) + "_" + std::string(comp);
            // Skip w for velocity, acceleration, and force
            if (comp == "w" && (prefix == "v" || prefix == "a" || prefix == "f")) {
                continue;
            }
            this->file_.SetChunking(var_name, std::span<const size_t>(chunking));
        }
    }

    // Reserve buffer capacity to avoid reallocations during writing
    if (this->buffer_size_ > 0U) {
        for (auto& buffer : this->state_buffers_) {
            buffer.reserve(this->buffer_size_);
        }
    }
}

NodeStateWriter::~NodeStateWriter() {
    // flush all remaining buffered data before destructing the object
    if (this->buffer_size_ > 0U) {
        this->FlushAllBuffers();
    }
}

void NodeStateWriter::WriteStateDataAtTimestep(
    size_t timestep, const std::string& component_prefix, const std::vector<double>& x,
    const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& i,
    const std::vector<double>& j, const std::vector<double>& k, const std::vector<double>& w
) {
    // Validate the component prefix - must be one of the valid prefixes:
    // "x" -> position
    // "u" -> displacement
    // "v" -> velocity
    // "a" -> acceleration
    // "f" -> force
    static const std::array<std::string_view, 5> valid_prefixes = {"x", "u", "v", "a", "f"};
    if (std::none_of(valid_prefixes.begin(), valid_prefixes.end(), [&](const auto& prefix) {
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

    //-----------------------------------
    // No buffering
    //-----------------------------------
    if (this->buffer_size_ == 0U) {
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
        return;
    }

    //-----------------------------------
    // write data w/ buffering
    //-----------------------------------
    const size_t component_index = GetComponentIndex(component_prefix);
    auto& buffer = this->state_buffers_[component_index];
    buffer.push_back(StateTimestepData{timestep, component_prefix, x, y, z, i, j, k, w});
    if (buffer.size() >= this->buffer_size_) {
        this->FlushStateBuffer(component_index);
    }
}

const NetCdfFile& NodeStateWriter::GetFile() const {
    return file_;
}

size_t NodeStateWriter::GetNumNodes() const {
    return num_nodes_;
}

void NodeStateWriter::Close() {
    // flush all remaining buffered data before closing the file
    if (this->buffer_size_ > 0U) {
        this->FlushAllBuffers();
    }
    this->file_.Close();
}

void NodeStateWriter::Open() {
    this->file_.Open();
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

size_t NodeStateWriter::GetComponentIndex(const std::string& prefix) {
    // Map: "x" -> 0, "u" -> 1, "v" -> 2, "a" -> 3, "f" -> 4
    static const std::unordered_map<std::string, size_t> prefix_to_index{
        {"x", 0}, {"u", 1}, {"v", 2}, {"a", 3}, {"f", 4}
    };
    auto it = prefix_to_index.find(prefix);
    if (it != prefix_to_index.end()) {
        return it->second;
    }
    throw std::invalid_argument("Unknown component prefix: " + prefix);
}

void NodeStateWriter::FlushStateBuffer(size_t component_index) {
    // Get the buffer for the component
    auto& buffer = this->state_buffers_[component_index];
    if (buffer.empty()) {
        return;  // nothing to flush
    }

    // Write data to file
    const std::string& prefix = buffer.front().component_prefix;
    const std::array<std::string, 6> var_names{prefix + "_x", prefix + "_y", prefix + "_z",
                                               prefix + "_i", prefix + "_j", prefix + "_k"};
    bool write_w{false};  // Write _w component only for position and displacement (x, u)
    if (!buffer.front().w.empty() && (prefix == "x" || prefix == "u")) {
        write_w = true;
    }
    for (const auto& record : buffer) {
        const std::vector<size_t> start = {record.timestep, 0};
        const std::vector<size_t> count = {1, record.x.size()};
        file_.WriteVariableAt(var_names[0], start, count, record.x);
        file_.WriteVariableAt(var_names[1], start, count, record.y);
        file_.WriteVariableAt(var_names[2], start, count, record.z);
        file_.WriteVariableAt(var_names[3], start, count, record.i);
        file_.WriteVariableAt(var_names[4], start, count, record.j);
        file_.WriteVariableAt(var_names[5], start, count, record.k);
        if (write_w && !record.w.empty()) {
            file_.WriteVariableAt(prefix + "_w", start, count, record.w);
        }
    }
    buffer.clear();
}

void NodeStateWriter::FlushAllBuffers() {
    // There are 5 state component types: x, u, v, a, f -> flush all
    for (size_t component_index = 0U; component_index < 5U; ++component_index) {
        this->FlushStateBuffer(component_index);
    }
    this->file_.Sync();
}

}  // namespace kynema::util
