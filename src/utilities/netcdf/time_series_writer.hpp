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
    /// Default buffer size (number of rows to accumulate before auto-flush, 0 = no buffering)
    static constexpr size_t kDefaultBufferSize{0};

    /**
     * @brief Constructor to create a TimeSeriesWriter object for writing time-series data with a
     * single variable
     *
     * @param file_path Path to the output NetCDF file
     * @param create Whether to create a new file or open an existing one
     */
    explicit TimeSeriesWriter(const std::string& file_path, bool create = true);

    /**
     * @brief Constructor to create a TimeSeriesWriter with predefined channels and buffering
     *
     * @param file_path Path to the output NetCDF file
     * @param create Whether to create a new file or open an existing one
     * @param channel_names Names of the channels (columns) stored in the time-series variable
     * @param channel_units Optional units for each channel (same size as channel_names if provided)
     * @param buffer_size Number of rows to buffer before auto-flush (0 = no buffering)
     */
    TimeSeriesWriter(
        const std::string& file_path, bool create, const std::vector<std::string>& channel_names,
        const std::vector<std::string>& channel_units = {}, size_t buffer_size = kDefaultBufferSize
    );

    /// @brief Destructor flushes any remaining buffered rows
    ~TimeSeriesWriter();

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

    /**
     * @brief Writes a full row (all channels) at a specific timestep
     *
     * @param timestep Current timestep index
     * @param row Vector of values to write at the current timestep
     */
    void WriteRowAtTimestep(size_t timestep, std::span<const double> row);

    /// @brief Flushes any remaining buffered rows
    void Flush();

    /// @brief Manually closes the underlying NetCDF file and flush any remaining buffered rows
    void Close();

    /// @brief Manually (re)opens the underlying NetCDF file
    void Open();

    /// @brief Gets the NetCDF file object
    [[nodiscard]] const NetCdfFile& GetFile() const;

    /// @brief Gets the channel names
    [[nodiscard]] const std::vector<std::string>& GetChannelNames() const;

    /// @brief Gets the channel units
    [[nodiscard]] const std::vector<std::string>& GetChannelUnits() const;

    /// @brief Gets the number of channels
    [[nodiscard]] size_t GetNumChannels() const;

    /// @brief Gets the buffer size
    [[nodiscard]] size_t GetBufferSize() const;

private:
    NetCdfFile file_;

    //-----------------------------------
    // channel output variables
    //-----------------------------------

    std::vector<std::string> channel_names_;
    std::vector<std::string> channel_units_;
    size_t num_channels_;
    size_t buffer_size_;

    //-----------------------------------
    // support for buffering
    //-----------------------------------

    /// Structure to hold buffered time-series data
    struct BufferedData {
        std::vector<double> data;  //< Flattened row-major buffer
        size_t start_timestep{0};  //< First timestep in the buffer
        size_t num_rows{0};        //< Number of complete rows currently in buffer
    };

    BufferedData buffer_;  //< Buffer for time-series data

    /**
     * @brief Flushes the buffered rows to the NetCDF file
     */
    void FlushBuffer();
};

}  // namespace kynema::util
