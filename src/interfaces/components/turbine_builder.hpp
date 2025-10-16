#pragma once

#include "beam_builder.hpp"
#include "turbine_input.hpp"

namespace kynema::interfaces::components {

class Turbine;

/// @brief Builder class for creating Turbine objects with a fluent interface pattern
class TurbineBuilder {
public:
    TurbineBuilder() = default;

    /**
     * @brief Get the current turbine input configuration
     * @return Reference to the current turbine input
     */
    [[nodiscard]] const TurbineInput& Input();

    /**
     * @brief Build a Turbine object from the current configuration
     * @param model The model to associate with this turbine
     * @return A new Turbine object
     */
    [[nodiscard]] Turbine Build(Model& model);

    //--------------------------------------------------------------------------
    // Build components
    //--------------------------------------------------------------------------

    /**
     * @brief Get reference to builder for a specific blade
     * @param blade_index The index of the blade
     * @return Reference to the blade builder
     */
    [[nodiscard]] components::BeamBuilder& Blade(size_t blade_index);

    /**
     * @brief Get reference to builder for the tower
     * @return Reference to the tower builder
     */
    [[nodiscard]] components::BeamBuilder& Tower();

    /**
     * @brief Set the yaw bearing inertia matrix (6x6)
     * @param matrix The inertia matrix to set, includes yaw bearing and nacelle mass with inertia
     * about yaw bearing i.e. system_inertia_tt in WindIO
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetYawBearingInertiaMatrix(const std::array<std::array<double, 6>, 6>& matrix);

    /**
     * @brief Set the hub inertia matrix (6x6)
     * @param matrix The inertia matrix to set, includes hub assembly mass and inertia
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetHubInertiaMatrix(const std::array<std::array<double, 6>, 6>& matrix);

    /**
     * @brief Set the nacelle inertia matrix (6x6) and center of mass offset (3x1)
     * @param matrix The inertia matrix to set, includes mass and inertia
     * @param cm_offset The nacelle center of mass offset from the yaw bearing (meters)
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetNacelleInertiaMatrix(
        const std::array<std::array<double, 6>, 6>& matrix, const std::array<double, 3>& cm_offset
    );

    //--------------------------------------------------------------------------
    // Build geometric configuration of the turbine
    //--------------------------------------------------------------------------

    /**
     * @brief Set the position of the tower base node
     * @param position Array containing position/orientation [x,y,z,qw,qx,qy,qz]
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetTowerBasePosition(const std::array<double, 7>& position);

    /**
     * @brief Set the distance from tower axis to hub i.e. distance from tower axis -> rotor apex
     * (meters)
     * @param distance The distance to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetTowerAxisToRotorApex(double distance);

    /**
     * @brief Set the hub height above the tower top i.e. distrance from tower top -> rotor apex
     * (meters)
     * @param height The hub height to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetTowerTopToRotorApex(double height);

    /**
     * @brief Distance from rotor apex to hub center of mass (meters)
     * @param distance The distance to set (meters)
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetRotorApexToHub(double distance);

    /**
     * @brief Set the hub diameter (meters)
     * @param diameter The hub diameter to set (meters)
     * @return Reference to this builder for method chaining
     */
    TurbineBuilder& SetHubDiameter(double diameter);

    //--------------------------------------------------------------------------
    // Build initial operating conditions of the turbine
    //--------------------------------------------------------------------------

    /**
     * @brief Set the shaft tilt angle (degrees)
     * @param angle The shaft tilt angle to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetShaftTiltAngle(double angle);

    /**
     * @brief Set the initial nacelle yaw angle (radians)
     * @param angle The nacelle yaw angle to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetNacelleYawAngle(double angle);

    /**
     * @brief Set the initial cone angle (radians)
     * @param angle The cone angle to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetConeAngle(double angle);

    /**
     * @brief Set the initial blade pitch angle (radians)
     * @param angle The blade pitch angle to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetBladePitchAngle(double angle);

    /**
     * @brief Set the azimuth angle (radiams)
     * @param angle The azimuth angle to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetAzimuthAngle(double angle);

    /**
     * @brief Set the initial rotor speed (rad/s)
     * @param speed The rotor speed to set
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetRotorSpeed(double speed);

    /**
     * @brief Set the gear box ratio
     * @param ratio ratio between generator speed and rotor speed
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetGearBoxRatio(double ratio);

    /**
     * @brief Set the generator power
     * @param power initial generator power (W)
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetGeneratorPower(double power);

    /**
     * @brief Set the controller wind speed
     * @param speed The wind speed at the hub
     * @return Reference to the builder for method chaining
     */
    TurbineBuilder& SetHubWindSpeed(double speed);

private:
    TurbineInput input;                       ///< turbine configuration being built
    std::vector<BeamBuilder> blade_builders;  ///< builders for the blade components
    BeamBuilder tower_builder;                ///< builder for the tower component
};

}  // namespace kynema::interfaces::components
