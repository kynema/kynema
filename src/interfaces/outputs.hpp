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
 * @brief A class that handles writing state data to disk for analysis/visualization
 */
class Outputs {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /// @brief Enum for selecting where on elements to write the outputs
    enum class OutputLocation : std::uint8_t {
        kNodes = 0,  ///< Write outputs at node locations
        kQPs = 1     ///< Write outputs at quadrature points
    };

    /**
     * @brief Constructor taking an output file and location
     *
     * @param output_file The name of the output file
     * @param num_nodes the number of nodes to be written to the file
     * @param location where the output will be written (quadrature points vs nodes)
     */
    Outputs(
        const std::string& output_file, size_t num_nodes,
        OutputLocation location = OutputLocation::kNodes
    );

    /**
     * @brief Constructor taking an output file, time-series file, and location
     *
     * @param output_file The name of the output file
     * @param time_series_file The name of the file with time-series data
     * @param num_nodes the number of nodes to be written to the file
     * @param location where the output will be written (quadrature points vs nodes)
     */
    Outputs(
        const std::string& output_file, const std::string& time_series_file, size_t num_nodes,
        OutputLocation location = OutputLocation::kNodes
    );

    /**
     * @brief Gets a reference to the NodeStateWriter for direct usage
     *
     * @return A reference to the NodeStateWriter
     */
    [[nodiscard]] std::unique_ptr<util::NodeStateWriter>& GetOutputWriter();

    /**
     * @brief Gets a reference to the TimeSeriesWriter for direct usage
     *
     * @return A reference to the TimeSeriesWriter
     */
    [[nodiscard]] std::unique_ptr<util::TimeSeriesWriter>& GetTimeSeriesWriter();

    /**
     * @brief Returns the output location as set by the constructor
     *
     * @return the OutputLocation
     */
    [[nodiscard]] OutputLocation GetLocation() const;

    /**
     * @brief Write node state outputs to NetCDF file at specified timestep
     *
     * @param host_state an updated host_state object with the current state loaded
     * @param timestep The timestep number to write data to
     */
    void WriteNodeOutputsAtTimestep(const HostState<DeviceType>& host_state, size_t timestep);

    /// @brief
    /**
     * @brief Write rotor time-series data at specified timestep
     *
     * @param timestep The timestep number to write data to
     * @param name The name of the variable to write (e.g., "azimuth_angle", "rotor_speed")
     * @param value The current value of the variable
     */
    void WriteValueAtTimestep(const size_t timestep, const std::string& name, const double value);

private:
    std::unique_ptr<util::NodeStateWriter> output_writer_;        ///< Output writer
    std::unique_ptr<util::TimeSeriesWriter> time_series_writer_;  ///< Time series writer
    size_t num_nodes_;         ///< Number of nodes to be written in the output file
    OutputLocation location_;  ///< Output writing location in element

    // Pre-allocated vectors for storing output data
    std::vector<double> x_data_;
    std::vector<double> y_data_;
    std::vector<double> z_data_;
    std::vector<double> i_data_;
    std::vector<double> j_data_;
    std::vector<double> k_data_;
    std::vector<double> w_data_;
};

}  // namespace kynema::interfaces
