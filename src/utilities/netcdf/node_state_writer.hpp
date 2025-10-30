#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "netcdf_file.hpp"

namespace kynema::util {

/**
 * @brief Class for writing Kynema nodal state data to NetCDF-based output files
 *
 * This class handles the writing of nodal state data for Kynema simulations to NetCDF format.
 * It manages the output of:
 *   - Position (x, y, z, w, i, j, k)
 *   - Displacement (x, y, z, w, i, j, k)
 *   - Velocity (x, y, z, i, j, k)
 *   - Acceleration (x, y, z, i, j, k)
 *   - Force (x, y, z, i, j, k)
 *
 * @note Each item is stored as a separate variable in the NetCDF file, organized by
 * timestep and node index. The file structure uses an unlimited time dimension to allow
 * for continuous writing of timesteps during simulation.
 *
 * @note The class includes buffering to improve write performance by batching multiple
 * timesteps together before writing to disk.
 */
class NodeStateWriter {
public:
    /// Default buffer size (number of timesteps to accumulate before auto-flush, 0 = no buffering)
    static constexpr size_t kDefaultBufferSize{0};

    /**
     * @brief Constructor to create a NodeStateWriter object
     *
     * @param file_path Path to the output NetCDF file
     * @param create Whether to create a new file or open an existing one
     * @param num_nodes Number of nodes in the simulation
     * @param enabled_state_prefixes Vector of state component prefixes to enable for writing
     *                              (default: all state components i.e. {"x", "u", "v", "a", "f"})
     * @param buffer_size Number of timesteps to accumulate before auto-flush (0 = no buffering)
     */
    NodeStateWriter(
        const std::string& file_path, bool create, size_t num_nodes,
        const std::vector<std::string>& enabled_state_prefixes = {"x", "u", "v", "a", "f"},
        size_t buffer_size = kDefaultBufferSize
    );

    /// @brief Destructor to flush any remaining buffered data
    ~NodeStateWriter();

    // Rule of five -> explicitly delete the copy ctor, copy assignment operator,
    // move ctor, and move assignment operator
    NodeStateWriter(const NodeStateWriter&) = delete;
    NodeStateWriter& operator=(const NodeStateWriter&) = delete;
    NodeStateWriter(NodeStateWriter&&) = delete;
    NodeStateWriter& operator=(NodeStateWriter&&) = delete;

    /**
     * @brief Writes state data for a specific timestep
     *
     * @param timestep Current timestep index
     * @param component_prefix Prefix for the component
     * @param x Data for component 1
     * @param y Data for component 2
     * @param z Data for component 3
     * @param i Data for component 4
     * @param j Data for component 5
     * @param k Data for component 6
     * @param w Data for component 7 (optional, only used for position and displacement)
     */
    void WriteStateDataAtTimestep(
        size_t timestep, const std::string& component_prefix, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& i,
        const std::vector<double>& j, const std::vector<double>& k,
        const std::vector<double>& w = std::vector<double>()
    );

    /// @brief Get the NetCDF file object
    [[nodiscard]] const NetCdfFile& GetFile() const;

    /// @brief Get the number of nodes with state data in output file
    [[nodiscard]] size_t GetNumNodes() const;

    /// @brief Manually flush and close the underlying NetCDF file
    void Close();

    /// @brief Manually (re)open the underlying NetCDF file
    void Open();

private:
    NetCdfFile file_;                                  //< NetCDF file object for writing output data
    size_t num_nodes_;                                 //< number of nodes in the simulation
    std::vector<std::string> enabled_state_prefixes_;  //< prefixes for the state components to write
    size_t buffer_size_;  //< number of timesteps to accumulate before auto-flush

    /**
     * @brief Defines variables for a state component (position, velocity, etc.)
     *
     * @param prefix Prefix for the variable names
     * @param dimensions Vector of dimension IDs
     * @param has_w Whether the component has a w component
     */
    void DefineStateVariables(
        const std::string& prefix, const std::vector<int>& dimensions, bool has_w
    );

    //-----------------------------------
    // support for buffering
    //-----------------------------------

    /// Structure to hold one timestep's worth of state data
    struct StateTimestepData {
        size_t timestep;
        std::string component_prefix;             //< prefix for the component
        std::vector<double> x, y, z, i, j, k, w;  //< data for the component
    };

    std::vector<StateTimestepData>
        state_buffers_[5];  //< buffer for each state component type  -- x, u, v, a, f

    /**
     * @brief Gets the index of a state component by its prefix
     * @param prefix The string prefix ("x", "u", "v", "a", or "f")
     * @return The index corresponding to the prefix [0, 1, 2, 3, 4]
     */
    [[nodiscard]] static size_t GetComponentIndex(const std::string& prefix);

    /**
     * @brief Flushes buffered data for a state component to disk
     * @param component_index Index of the component (0=x, 1=u, 2=v, 3=a, 4=f)
     */
    void FlushStateBuffer(size_t component_index);

    /**
     * @brief Flushes all buffered data and syncs the file
     */
    void FlushAllBuffers();
};

}  // namespace kynema::util
