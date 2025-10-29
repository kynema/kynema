#pragma once

#include "outputs_config.hpp"

namespace kynema::interfaces::components {

/**
 * @brief A builder class for building a outputs configuration
 */
struct OutputsBuilder {
    /**
     * @brief Sets the output file path for NetCDF results
     *
     * @details If the file path is empty, no output will be written to disk. This must be set
     *          to enable automatic output operations at each timestep.
     *
     * @param output_file_path Path where the output NetCDF file will be written
     * @return Reference to this OutputsBuilder object, enabling method chaining
     */
    OutputsBuilder& SetOutputFilePath(const std::string& output_file_path) {
        config.output_file_path = output_file_path;
        return *this;
    }

    /**
     * @brief Sets the output state prefixes to be written to the output file
     *
     * @param prefixes Vector of variable name prefixes
     * @return A reference to this outputs builder object to allow chaining
     */
    OutputsBuilder& SetOutputStatePrefixes(const std::vector<std::string>& prefixes) {
        config.output_state_prefixes = prefixes;
        return *this;
    }

    /**
     * @brief Sets the number of output timesteps to buffer before flushing to disk
     *
     * @param buffer_size Number of timesteps to buffer
     * @return A reference to this outputs builder object to allow chaining
     */
    OutputsBuilder& SetOutputBufferSize(std::size_t buffer_size) {
        config.buffer_size = buffer_size;
        return *this;
    }

    /**
     * @brief Creates an OutputsConfig object based on the previously set parameters
     *
     * @return The completed OutputsConfig object
     */
    [[nodiscard]] OutputsConfig Config() const { return this->config; }

private:
    OutputsConfig config;
};

}  // namespace kynema::interfaces::components
