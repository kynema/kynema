#pragma once

#include <memory>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>

#include "utilities/netcdf/node_state_writer.hpp"
#include "utilities/netcdf/time_series_writer.hpp"

namespace kynema::interfaces {

template <typename DeviceType>
struct HostState;

/**
 * @brief Handles writing state data to disk as simulation outputs and provides a means for
 * post-processing e.g. visualization
 */
class Outputs {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /**
     * @brief Constructor taking an output file and location
     *
     * @param output_file name of the output file
     * @param num_nodes number of nodes to be written to the file
     * @param time_series_file optional name of the file with time-series data (empty string = no
     *                         time series will be written)
     * @param enabled_state_prefixes which state component prefixes to enable for writing
     *        (default: all states i.e. {"x", "u", "v", "a", "f"})
     * @param enable_deformation whether to write deformation data (default: false)
     * @param buffer_size number of timesteps to buffer before auto-flush (default: 0 = no buffering)
     */
    Outputs(
        const std::string& output_file, size_t num_nodes, const std::string& time_series_file = "",
        const std::vector<std::string>& enabled_state_prefixes = {"x", "u", "v", "a", "f"},
        bool enable_deformation = false,
        size_t buffer_size = util::NodeStateWriter::kDefaultBufferSize
    );

    /**
     * @brief Gets a reference to the NodeStateWriter for direct usage
     */
    [[nodiscard]] std::unique_ptr<util::NodeStateWriter>& GetOutputWriter();

    /**
     * @brief Gets a reference to the TimeSeriesWriter for direct usage
     */
    [[nodiscard]] std::unique_ptr<util::TimeSeriesWriter>& GetTimeSeriesWriter();

    /**
     * @brief Write node state outputs to NetCDF file at specified timestep
     *
     * @param host_state an updated host_state object with the current state loaded
     * @param timestep The timestep number to write data to
     */
    void WriteNodeOutputsAtTimestep(const HostState<DeviceType>& host_state, size_t timestep);

    /**
     * @brief Write rotor time-series data at specified timestep
     *
     * @param timestep The timestep number to write data to
     * @param name The name of the variable to write (e.g., "azimuth_angle", "rotor_speed")
     * @param value The current value of the variable
     */
    void WriteValueAtTimestep(size_t timestep, const std::string& name, double value);

    /// @brief Manually close the underlying NetCDF files
    void Close();

    /// @brief Manually (re)open the underlying NetCDF files
    void Open();

private:
    std::unique_ptr<util::NodeStateWriter> output_writer_;  ///< Output writer
    size_t num_nodes_;  ///< Number of nodes to be written in the output file
    std::unique_ptr<util::TimeSeriesWriter> time_series_writer_;  ///< Time series writer
    std::vector<std::string> enabled_state_prefixes_;  ///< Which state components to write

    // Pre-allocated vectors for storing output data
    std::vector<double> x_data_;
    std::vector<double> y_data_;
    std::vector<double> z_data_;
    std::vector<double> i_data_;
    std::vector<double> j_data_;
    std::vector<double> k_data_;
    std::vector<double> w_data_;

    /**
     * @brief Helper to check if a state component is enabled
     * @param prefix The state prefix to check ("x", "u", "v", "a", "f")
     * @return true if the state is enabled for writing
     */
    [[nodiscard]] bool IsStateComponentEnabled(const std::string& prefix) const;
};

}  // namespace kynema::interfaces
