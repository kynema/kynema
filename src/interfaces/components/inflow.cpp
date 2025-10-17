#include "inflow.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <stdexcept>

namespace kynema::interfaces::components {

std::array<double, 3> UniformFlowParameters::Velocity(const std::array<double, 3>& position) const {
    // Calculate horizontal velocity using power law wind shear
    // vh = v_ref * (z / z_ref)^alpha
    const double vh = velocity_horizontal * std::pow(position[2] / height_reference, shear_vertical);

    // Apply horizontal direction to the velocity vector
    const double sin_flow_angle = std::sin(flow_angle_horizontal);
    const double cos_flow_angle = std::cos(flow_angle_horizontal);
    return {vh * cos_flow_angle, -vh * sin_flow_angle, 0.};
}

std::array<double, 3> UniformFlow::Velocity(
    [[maybe_unused]] double t, const std::array<double, 3>& position
) const {
    // If no data is available, throw an error
    if (data.empty()) {
        throw std::runtime_error("No uniform flow data available");
    }

    // If there is only one data point or time is less than the first data point's time,
    // use the first uniform flow parameters
    if (data.size() == 1 || t <= data.front().time) {
        return data[0].Velocity(position);
    }

    // If time is greater than the last data point's time, use the last uniform flow parameters
    if (t >= data.back().time) {
        return data.back().Velocity(position);
    }

    const auto time_iterator = std::find_if(std::cbegin(data), std::cend(data), [t](auto d) {
        return d.time > t;
    });

    const auto t_index = static_cast<std::size_t>(std::distance(std::cbegin(data), time_iterator));

    const double alpha =
        (t - data[t_index - 1].time) / (data[t_index].time - data[t_index - 1].time);

    // Create interpolated UniformFlowParameters at given time
    const auto t_data = UniformFlowParameters{
        t,
        data[t_index - 1].velocity_horizontal * (1 - alpha) +
            data[t_index].velocity_horizontal * alpha,
        data[t_index - 1].height_reference * (1 - alpha) + data[t_index].height_reference * alpha,
        data[t_index - 1].shear_vertical * (1 - alpha) + data[t_index].shear_vertical * alpha,
        data[t_index - 1].flow_angle_horizontal * (1 - alpha) +
            data[t_index].flow_angle_horizontal * alpha
    };

    return t_data.Velocity(position);
}

Inflow Inflow::SteadyWind(double vh, double z_ref, double alpha, double angle_h) {
    return Inflow(
        InflowType::Uniform, UniformFlow{std::vector<UniformFlowParameters>{
                                 {0., vh, z_ref, alpha, angle_h}  // Use aggregate initialization
                             }}
    );
}

std::array<double, 3> Inflow::Velocity(double t, const std::array<double, 3>& position) const {
    if (type == InflowType::Uniform) {
        return uniform_flow.Velocity(t, position);
    }
    throw std::runtime_error("Unknown inflow type");
}

}  // namespace kynema::interfaces::components
