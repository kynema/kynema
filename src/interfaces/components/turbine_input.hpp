#pragma once

#include <array>
#include <vector>

#include "beam_input.hpp"

namespace kynema::interfaces::components {

/**
 * @brief Complete input specification for a turbine
 *
 * This structure defines all necessary parameters to configure a complete wind turbine model
 * including the structural components (blades, tower, hub, nacelle), geometric relationships,
 * inertial properties, and initial operating conditions.
 */
struct TurbineInput {
    //--------------------------------------------------------------------------
    // Structural component inputs
    //--------------------------------------------------------------------------

    /**
     * @brief Configuration for the turbine blades
     * @details Each blade is represented as a flexible beam. Vector contains BeamInput for each
     * blade of the turbine.
     */
    std::vector<BeamInput> blades;

    /**
     * @brief Configuration for the tower structure
     * @details Defines the tower as a beam element from base to top
     */
    BeamInput tower;

    //--------------------------------------------------------------------------
    // Inertial properties inputs
    //--------------------------------------------------------------------------

    /**
     * @brief Yaw bearing inertia matrix (6x6)
     * @details Includes yaw bearing and nacelle mass with inertia about yaw bearing
     */
    std::array<std::array<double, 6>, 6> yaw_bearing_inertia_matrix{};

    /**
     * @brief Flag indicating if yaw_bearing inertia matrix has been set
     */
    bool yaw_bearing_inertia_matrix_set{false};

    /**
     * @brief Hub inertia matrix (6x6)
     * @details Includes hub assembly mass and inertia properties
     */
    std::array<std::array<double, 6>, 6> hub_inertia_matrix{};

    /**
     * @brief Nacelle inertia matrix (6x6)
     * @details Includes nacelle mass and inertia properties
     */
    std::array<std::array<double, 6>, 6> nacelle_inertia_matrix{};

    /**
     * @brief Flag indicating if nacelle inertia matrix has been set
     */
    bool nacelle_inertia_matrix_set{false};

    /**
     * @brief Nacelle center of mass offset from the yaw bearing (meters)
     */
    std::array<double, 3> nacelle_cm_offset{};

    //--------------------------------------------------------------------------
    // Geometric configuration inputs
    //--------------------------------------------------------------------------

    /**
     * @brief Position of the tower base node in the global coordinate system
     * @details Defines the initial location/orientation of the entire turbine.
     */
    std::array<double, 7> tower_base_position{0., 0., 0., 1., 0., 0., 0.};

    /**
     * @brief Horizontal distance between the tower axis -> rotor apex (meters)
     */
    double tower_axis_to_rotor_apex{0.};

    /**
     * @brief Vertical distance between the tower top -> rotor apex (meters)
     */
    double tower_top_to_rotor_apex{0.};

    /**
     * @brief Distance from rotor apex -> hub center of mass (meters)
     * @details Along the rotor shaft axis, accounts for hub geometry. Also known as HubCM.
     */
    double rotor_apex_to_hub{0.};

    /**
     * @brief Hub diameter (meters)
     * @details Distance from rotor apex to blade root node.
     * @note Default value is near-zero to avoid singularities
     */
    double hub_diameter{1e-7};

    //--------------------------------------------------------------------------
    // Initial operating condition inputs
    //--------------------------------------------------------------------------

    /**
     * @brief Initial nacelle yaw angle (radians)
     * @details Rotation of nacelle about the tower/yaw axis.
     */
    double nacelle_yaw_angle{0.};

    /**
     * @brief Shaft tilt angle (radians)
     * @details Angle between the rotor shaft and the horizontal plane. Positive angle tilts shaft
     * up.
     */
    double shaft_tilt_angle{0.};

    /**
     * @brief Blade cone angle (radians)
     * @details Angle between the blade and the plane perpendicular to shaft
     */
    double cone_angle{0.};

    /**
     * @brief Initial blade pitch angle (radians)
     * @details Rotation of blades about their longitudinal/pitch axis. 0 means blades are
     * parallel to the tower.
     */
    double blade_pitch_angle{0.};

    /**
     * @brief Current rotor azimuth angle (radians)
     * @details Angular position of rotor about shaft axis. 0 means blade 1 is at 12 o'clock
     * position.
     */
    double azimuth_angle{0.};

    /**
     * @brief Initial rotor rotational speed (rad/s)
     * @details Rotational velocity about the shaft axis
     */
    double rotor_speed{0.};

    //--------------------------------------------------------------------------
    // Controller inputs
    //--------------------------------------------------------------------------

    /**
     * @brief Gearbox ratio of the turbine (-)
     * @details Ratio between generator speed and rotor speed
     */
    double gear_box_ratio{1.};

    /**
     * @brief Initial generator power (W)
     * @details Power output of the generator
     */
    double generator_power{0.};

    /**
     * @brief Hub height wind speed (m/s)
     * @details Wind speed at the hub height
     */
    double hub_wind_speed{0.};
};

}  // namespace kynema::interfaces::components
