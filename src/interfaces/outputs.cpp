#include "outputs.hpp"

#include <ranges>

#include "host_state.hpp"

namespace kynema::interfaces {

Outputs::Outputs(
    const std::string& output_file, size_t num_nodes, const std::string& time_series_file,
    const std::vector<std::string>& enabled_state_prefixes, size_t buffer_size
)
    : output_writer_(std::make_unique<util::NodeStateWriter>(
          output_file, true, num_nodes, enabled_state_prefixes, buffer_size
      )),
      num_nodes_(num_nodes),
      time_series_writer_(
          time_series_file.empty() ? nullptr
                                   : std::make_unique<util::TimeSeriesWriter>(time_series_file, true)
      ),
      enabled_state_prefixes_(enabled_state_prefixes) {
    this->x_data_.resize(num_nodes_);
    this->y_data_.resize(num_nodes_);
    this->z_data_.resize(num_nodes_);
    this->i_data_.resize(num_nodes_);
    this->j_data_.resize(num_nodes_);
    this->k_data_.resize(num_nodes_);
    this->w_data_.resize(num_nodes_);
}

std::unique_ptr<util::NodeStateWriter>& Outputs::GetOutputWriter() {
    return this->output_writer_;
}

std::unique_ptr<util::TimeSeriesWriter>& Outputs::GetTimeSeriesWriter() {
    return this->time_series_writer_;
}

bool Outputs::IsStateComponentEnabled(const std::string& prefix) const {
    return std::find(
               this->enabled_state_prefixes_.begin(), this->enabled_state_prefixes_.end(), prefix
           ) != this->enabled_state_prefixes_.end();
}

void Outputs::WriteNodeOutputsAtTimestep(const HostState<DeviceType>& host_state, size_t timestep) {
    if (!this->output_writer_) {
        return;
    }

    // Position data
    if (this->IsStateComponentEnabled("x")) {
        for (auto node : std::views::iota(0U, num_nodes_)) {
            this->x_data_[node] = host_state.x(node, 0);
            this->y_data_[node] = host_state.x(node, 1);
            this->z_data_[node] = host_state.x(node, 2);
            this->w_data_[node] = host_state.x(node, 3);
            this->i_data_[node] = host_state.x(node, 4);
            this->j_data_[node] = host_state.x(node, 5);
            this->k_data_[node] = host_state.x(node, 6);
        }
        this->output_writer_->WriteStateDataAtTimestep(
            timestep, "x", this->x_data_, this->y_data_, this->z_data_, this->i_data_, this->j_data_,
            this->k_data_, this->w_data_
        );
    }

    // Displacement data
    if (this->IsStateComponentEnabled("u")) {
        for (auto node : std::views::iota(0U, num_nodes_)) {
            this->x_data_[node] = host_state.q(node, 0);
            this->y_data_[node] = host_state.q(node, 1);
            this->z_data_[node] = host_state.q(node, 2);
            this->w_data_[node] = host_state.q(node, 3);
            this->i_data_[node] = host_state.q(node, 4);
            this->j_data_[node] = host_state.q(node, 5);
            this->k_data_[node] = host_state.q(node, 6);
        }
        this->output_writer_->WriteStateDataAtTimestep(
            timestep, "u", this->x_data_, this->y_data_, this->z_data_, this->i_data_, this->j_data_,
            this->k_data_, this->w_data_
        );
    }

    // Velocity data
    if (this->IsStateComponentEnabled("v")) {
        for (auto node : std::views::iota(0U, num_nodes_)) {
            this->x_data_[node] = host_state.v(node, 0);
            this->y_data_[node] = host_state.v(node, 1);
            this->z_data_[node] = host_state.v(node, 2);
            this->i_data_[node] = host_state.v(node, 3);
            this->j_data_[node] = host_state.v(node, 4);
            this->k_data_[node] = host_state.v(node, 5);
        }
        this->output_writer_->WriteStateDataAtTimestep(
            timestep, "v", this->x_data_, this->y_data_, this->z_data_, this->i_data_, this->j_data_,
            this->k_data_
        );
    }

    // Acceleration data
    if (this->IsStateComponentEnabled("a")) {
        for (auto node : std::views::iota(0U, num_nodes_)) {
            this->x_data_[node] = host_state.vd(node, 0);
            this->y_data_[node] = host_state.vd(node, 1);
            this->z_data_[node] = host_state.vd(node, 2);
            this->i_data_[node] = host_state.vd(node, 3);
            this->j_data_[node] = host_state.vd(node, 4);
            this->k_data_[node] = host_state.vd(node, 5);
        }
        this->output_writer_->WriteStateDataAtTimestep(
            timestep, "a", this->x_data_, this->y_data_, this->z_data_, this->i_data_, this->j_data_,
            this->k_data_
        );
    }

    // Force data
    if (this->IsStateComponentEnabled("f")) {
        for (auto node : std::views::iota(0U, num_nodes_)) {
            this->x_data_[node] = host_state.f(node, 0);
            this->y_data_[node] = host_state.f(node, 1);
            this->z_data_[node] = host_state.f(node, 2);
            this->i_data_[node] = host_state.f(node, 3);
            this->j_data_[node] = host_state.f(node, 4);
            this->k_data_[node] = host_state.f(node, 5);
        }
        this->output_writer_->WriteStateDataAtTimestep(
            timestep, "f", this->x_data_, this->y_data_, this->z_data_, this->i_data_, this->j_data_,
            this->k_data_
        );
    }
}

void Outputs::WriteValueAtTimestep(size_t timestep, const std::string& name, double value) {
    if (this->time_series_writer_) {
        this->time_series_writer_->WriteValueAtTimestep(name, timestep, value);
    }
}

void Outputs::Close() {
    if (this->output_writer_) {
        this->output_writer_->Close();
    }
    if (this->time_series_writer_) {
        this->time_series_writer_->Close();
    }
}

void Outputs::Open() {
    if (this->output_writer_) {
        this->output_writer_->Open();
    }
    if (this->time_series_writer_) {
        this->time_series_writer_->Open();
    }
}

}  // namespace kynema::interfaces
