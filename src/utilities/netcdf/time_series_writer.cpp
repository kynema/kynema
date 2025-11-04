#include "time_series_writer.hpp"

#include <stdexcept>

namespace kynema::util {
TimeSeriesWriter::TimeSeriesWriter(const std::string& file_path, bool create)
    : file_(file_path, create) {
    // Check if the "time" dimension already exists in the file
    try {
        this->time_dim_ = file_.GetDimensionId("time");
    } catch (const std::runtime_error&) {
        this->time_dim_ =
            file_.AddDimension("time", NC_UNLIMITED);  // Unlimited timesteps can be added
    }
}

void TimeSeriesWriter::WriteValuesAtTimestep(
    const std::string& variable_name, size_t timestep, std::span<const double> values
) {
    // Check if the variable already exists in the file
    try {
        (void)file_.GetVariableId(variable_name);
    } catch (const std::runtime_error&) {
        int value_dim{-1};
        const std::string value_dim_name = variable_name + "_dimension";
        // Check if the value dimension already exists in the variable
        try {
            value_dim = file_.GetDimensionId(value_dim_name);
        } catch (const std::runtime_error&) {
            value_dim = file_.AddDimension(value_dim_name, values.size());
        }

        // Add the variable to the file
        const std::vector<int> dimensions = {time_dim_, value_dim};
        (void)file_.AddVariable<double>(variable_name, dimensions);
    }

    // Write the values to the time-series variable
    const std::vector<size_t> start = {timestep, 0};  // Start at the current timestep and value 0
    const std::vector<size_t> count = {
        1, values.size()
    };  // Write one timestep worth of data for all values
    file_.WriteVariableAt(variable_name, start, count, values);
}

void TimeSeriesWriter::WriteValueAtTimestep(
    const std::string& variable_name, size_t timestep, const double& value
) {
    WriteValuesAtTimestep(variable_name, timestep, std::array{value});
}

const NetCdfFile& TimeSeriesWriter::GetFile() const {
    return file_;
}
}  // namespace kynema::util
