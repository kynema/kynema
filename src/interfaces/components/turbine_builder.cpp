#include "turbine_builder.hpp"

#include "model/model.hpp"
#include "turbine.hpp"

namespace kynema::interfaces::components {
[[nodiscard]] const TurbineInput& TurbineBuilder::Input() {
    // Add the blade inputs from the blade builders
    this->input.blades.clear();
    for (const auto& builder : this->blade_builders) {
        this->input.blades.push_back(builder.Input());
    }

    // Get the tower input from the tower builder
    this->input.tower = this->tower_builder.Input();

    return input;
}

[[nodiscard]] Turbine TurbineBuilder::Build(Model& model) {
    // build turbine and return
    return {this->Input(), model};
}

[[nodiscard]] components::BeamBuilder& TurbineBuilder::Blade(size_t blade_index) {
    // Ensure we have enough blade builders
    if (blade_index >= this->blade_builders.size()) {
        this->blade_builders.resize(blade_index + 1);
    }
    return this->blade_builders[blade_index];
}

[[nodiscard]] components::BeamBuilder& TurbineBuilder::Tower() {
    // return tower builder
    return this->tower_builder;
}

TurbineBuilder& TurbineBuilder::SetYawBearingInertiaMatrix(
    const std::array<std::array<double, 6>, 6>& matrix
) {
    this->input.yaw_bearing_inertia_matrix = matrix;
    this->input.yaw_bearing_inertia_matrix_set = true;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetHubInertiaMatrix(
    const std::array<std::array<double, 6>, 6>& matrix
) {
    this->input.hub_inertia_matrix = matrix;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetNacelleInertiaMatrix(
    const std::array<std::array<double, 6>, 6>& matrix, const std::array<double, 3>& cm_offset
) {
    this->input.nacelle_inertia_matrix = matrix;
    this->input.nacelle_cm_offset = cm_offset;
    this->input.nacelle_inertia_matrix_set = true;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetTowerBasePosition(const std::array<double, 7>& position) {
    this->input.tower_base_position = position;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetTowerAxisToRotorApex(double distance) {
    this->input.tower_axis_to_rotor_apex = distance;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetTowerTopToRotorApex(double height) {
    this->input.tower_top_to_rotor_apex = height;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetRotorApexToHub(double distance) {
    this->input.rotor_apex_to_hub = distance;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetHubDiameter(double diameter) {
    this->input.hub_diameter = diameter;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetShaftTiltAngle(double angle) {
    this->input.shaft_tilt_angle = angle;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetNacelleYawAngle(double angle) {
    this->input.nacelle_yaw_angle = angle;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetConeAngle(double angle) {
    this->input.cone_angle = angle;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetBladePitchAngle(double angle) {
    this->input.blade_pitch_angle = angle;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetAzimuthAngle(double angle) {
    this->input.azimuth_angle = angle;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetRotorSpeed(double speed) {
    this->input.rotor_speed = speed;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetGearBoxRatio(double ratio) {
    this->input.gear_box_ratio = ratio;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetGeneratorPower(double power) {
    this->input.generator_power = power;
    return *this;
}

TurbineBuilder& TurbineBuilder::SetHubWindSpeed(double speed) {
    this->input.hub_wind_speed = speed;
    return *this;
}
}  // namespace kynema::interfaces::components
