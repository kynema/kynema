#pragma once

#include <span>
#include <string>
#include <vector>

#include "netcdf_file.hpp"

namespace kynema::util {

/**
 * @brief Class for writing time-series data to NetCDF file
 */
class TimeSeriesWriter {
public:
    /**
     * @brief Constructor to create a TimeSeriesWriter object
     *
     * @param file_path Path to the output NetCDF file
     * @param create Whether to create a new file or open an existing one
     */
    explicit TimeSeriesWriter(const std::string& file_path, bool create = true);

    /**
     * @brief Writes multiple values for a time-series variable at a specific timestep
     *
     * @param variable_name Name of the variable to write
     * @param timestep Current timestep index
     * @param values Vector of values to write at the current timestep
     */
    void WriteValuesAtTimestep(
        const std::string& variable_name, size_t timestep, std::span<const double> values
    );

    /**
     * @brief Writes a single value for a time-series variable at a specific timestep
     *
     * @param variable_name Name of the variable to write
     * @param timestep Current timestep index
     * @param value Value to write at the current timestep
     */
    void WriteValueAtTimestep(
        const std::string& variable_name, size_t timestep, const double& value
    );

    /// @brief Get the NetCDF file object
    [[nodiscard]] const NetCdfFile& GetFile() const;

private:
    NetCdfFile file_;
    int time_dim_;
};

}  // namespace kynema::util
