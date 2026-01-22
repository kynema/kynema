#include "time_series_writer.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace kynema::util {

TimeSeriesWriter::TimeSeriesWriter(const std::string& file_path, bool create)
    : file_(file_path, create), num_channels_{0}, buffer_size_{kDefaultBufferSize} {
    (void)file_.AddDimension("time", NC_UNLIMITED);
}

TimeSeriesWriter::TimeSeriesWriter(
    const std::string& file_path, bool create, const std::vector<std::string>& channel_names,
    const std::vector<std::string>& channel_units, size_t buffer_size
)
    : file_(file_path, create),
      channel_names_(channel_names),
      channel_units_(channel_units),
      num_channels_(channel_names.size()),
      buffer_size_(buffer_size) {
    //-----------------------------------
    // validate input parameters
    //-----------------------------------

    if (this->channel_names_.empty()) {
        throw std::invalid_argument(
            "TimeSeriesWriter must be constructed with at least one output channel"
        );
    }

    if (!this->channel_units_.empty() && this->channel_units_.size() != this->num_channels_) {
        throw std::invalid_argument(
            "TimeSeriesWriter must be constructed with output channel_units size matching "
            "output channel_names size"
        );
    }

    //-----------------------------------
    // define variables
    //-----------------------------------

    // Define the time dimension
    // NOTE: We are assuming this dimension was not added previously
    const int time_dim =
        file_.AddDimension("time", NC_UNLIMITED);  // Unlimited timesteps can be added

    // Define output channel dimension and metadata
    const int channel_dim = file_.AddDimension("channel", this->num_channels_);

    // Add channel metadata variables
    (void)file_.AddVariable<std::string>("channel_names", std::array<int, 1>{channel_dim});
    file_.WriteVariable("channel_names", this->channel_names_);
    if (!this->channel_units_.empty()) {
        (void)file_.AddVariable<std::string>("channel_units", std::array<int, 1>{channel_dim});
        file_.WriteVariable("channel_units", this->channel_units_);
    }

    // Define the data variable [time, channel] - this is a 2D variable with the time and channel
    // dimensions where channel is a vector of outputs
    const std::vector<int> dimensions = {time_dim, channel_dim};
    (void)file_.AddVariable<double>("time_series_channel_outputs", dimensions);

    //-----------------------------------
    // write in chunks
    //-----------------------------------

    // Configure chunking: time-major, write full channel width at once
    const size_t chunk_size =
        std::max(this->buffer_size_, size_t(1));  // minimum chunk size is 1 timestep
    const std::array<size_t, 2> chunking{chunk_size, this->num_channels_};
    file_.SetChunking("time_series_channel_outputs", std::span<const size_t>(chunking));

    //-----------------------------------
    // initialize buffer
    //-----------------------------------

    // Reserve buffer capacity to avoid reallocations during writing
    if (this->buffer_size_ > 0U) {
        this->buffer_.data.reserve(this->buffer_size_ * this->num_channels_);
    }
}

TimeSeriesWriter::~TimeSeriesWriter() {
    // flush any remaining buffered rows before destructing the object
    this->Flush();
}

void TimeSeriesWriter::WriteValuesAtTimestep(
    const std::string& variable_name, size_t timestep, std::span<const double> values
) {
    // Check if the variable already exists in the file
    try {
        (void)file_.GetVariableId(variable_name);
    } catch (const std::runtime_error&) {
        // Get the time dimension id
        const int time_dim = file_.GetDimensionId("time");

        // Get or create the value dimension id for this variable
        const std::string value_dim_name = variable_name + "_dimension";
        int value_dim{-1};
        try {
            value_dim = file_.GetDimensionId(value_dim_name);
        } catch (const std::runtime_error&) {
            value_dim = file_.AddDimension(value_dim_name, values.size());
        }

        // Add the variable to the file [time, value_dim]
        const std::vector<int> dimensions = {time_dim, value_dim};
        (void)file_.AddVariable<double>(variable_name, dimensions);
    }

    // Write the values to the time-series variable
    const std::vector<size_t> start = {timestep, 0};
    const std::vector<size_t> count = {1, values.size()};
    file_.WriteVariableAt(variable_name, start, count, values);
}

void TimeSeriesWriter::WriteValueAtTimestep(
    const std::string& variable_name, size_t timestep, const double& value
) {
    WriteValuesAtTimestep(variable_name, timestep, std::array{value});
}

void TimeSeriesWriter::WriteRowAtTimestep(size_t timestep, std::span<const double> row) {
    // Validate inputs
    if (row.size() != this->num_channels_) {
        throw std::invalid_argument("Row size must match number of channels to be written");
    }

    //-----------------------------------
    // No buffering
    //-----------------------------------

    // No buffering requested -> write one row at a time
    if (this->buffer_size_ == 0U) {
        const std::array<size_t, 2> start{timestep, 0};
        const std::array<size_t, 2> count{1, this->num_channels_};
        this->file_.WriteVariableAt("time_series_channel_outputs", start, count, row);
        return;
    }

    //-----------------------------------
    // write data w/ buffering
    //-----------------------------------

    // Ensure contiguous timesteps in the buffer
    if (this->buffer_.num_rows == 0U) {
        // buffer is empty -> start new batch at this timestep
        this->buffer_.start_timestep = timestep;
    } else if (timestep != this->buffer_.start_timestep + this->buffer_.num_rows) {
        // timestep is not contiguous -> flush current buffer and restart
        this->FlushBuffer();
        this->buffer_.start_timestep = timestep;
    }

    // Append row to buffer
    this->buffer_.data.insert(this->buffer_.data.end(), row.begin(), row.end());
    ++this->buffer_.num_rows;

    // Auto flush when buffer reaches capacity
    if (this->buffer_.num_rows >= this->buffer_size_) {
        this->FlushBuffer();
    }
}

void TimeSeriesWriter::FlushBuffer() {
    if (this->buffer_.num_rows == 0U) {
        return;  // nothing to flush
    }

    // Write the buffered data to the file
    const std::array<size_t, 2> start{this->buffer_.start_timestep, 0};
    const std::array<size_t, 2> count{this->buffer_.num_rows, this->num_channels_};
    const size_t total_values = this->buffer_.num_rows * this->num_channels_;
    this->file_.WriteVariableAt(
        "time_series_channel_outputs", start, count,
        std::span<const double>(this->buffer_.data.data(), total_values)
    );

    // Reset the buffer
    this->buffer_.num_rows = 0U;
    // no need to clear the data, it will be overwritten next time
}

void TimeSeriesWriter::Flush() {
    this->FlushBuffer();
    file_.Sync();
}

void TimeSeriesWriter::Close() {
    // flush any remaining buffered rows before closing the file
    if (this->buffer_size_ > 0U) {
        this->Flush();
    }
    this->file_.Close();
}

void TimeSeriesWriter::Open() {
    this->file_.Open();
}

const NetCdfFile& TimeSeriesWriter::GetFile() const {
    return this->file_;
}

const std::vector<std::string>& TimeSeriesWriter::GetChannelNames() const {
    return this->channel_names_;
}

const std::vector<std::string>& TimeSeriesWriter::GetChannelUnits() const {
    return this->channel_units_;
}

size_t TimeSeriesWriter::GetNumChannels() const {
    return this->num_channels_;
}

size_t TimeSeriesWriter::GetBufferSize() const {
    return this->buffer_size_;
}

}  // namespace kynema::util
