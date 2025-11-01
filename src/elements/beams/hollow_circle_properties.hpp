#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>

#include "beam_section.hpp"
#include "generate_sectional_properties.hpp"

namespace kynema::beams {

/**
 * @brief Struct containing geometric properties for a hollow circular cross-section
 */
struct HollowCircleProperties {
    double area;  ///< Cross-sectional area [Length^2]
    double Ixx;   ///< Second moment of area about x-axis [Length^4]
    double Iyy;   ///< Second moment of area about y-axis [Length^4]
    double J;     ///< Polar moment of inertia [Length^4]
    double kx;    ///< Shear correction factor in x direction [dimensionless]
    double ky;    ///< Shear correction factor in y direction [dimensionless]
};

/**
 * @brief Calculates geometric properties for a hollow circular cross-section
 *
 * @param outer_diameter Outer diameter of the hollow circle [Length]
 * @param wall_thickness Wall thickness [Length]
 * @param nu Poisson's ratio for shear correction factor calculation [dimensionless]
 *
 * @return HollowCircleProperties struct containing geometric properties
 */
static HollowCircleProperties CalculateHollowCircleProperties(
    double outer_diameter, double wall_thickness, double nu = 0.33
) {
    //--------------------------------------------------------------------------
    // Check input arguments
    //--------------------------------------------------------------------------
    if (wall_thickness >= outer_diameter / 2.) {
        throw std::invalid_argument("Wall thickness must be less than outer radius");
    }

    //--------------------------------------------------------------------------
    // Calculate geometric properties
    //--------------------------------------------------------------------------
    const double outer_radius = outer_diameter / 2.;
    const double inner_radius = outer_radius - wall_thickness;

    // Calculate geometric properties using difference of circles
    const double area = std::numbers::pi * (std::pow(outer_radius, 2) - std::pow(inner_radius, 2));
    const double Ixx =
        std::numbers::pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 4.;
    const double Iyy = Ixx;  // Circular symmetry
    const double J = std::numbers::pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 2.;

    // Shear correction factors for hollow circular sections
    const double kx = (6. * (1. + nu)) / (7. + 6. * nu);  // Timoshenko-Ehrenfest beam theory
    const double ky = kx;                                 // Circular symmetry

    return HollowCircleProperties{.area = area, .Ixx = Ixx, .Iyy = Iyy, .J = J, .kx = kx, .ky = ky};
}

/**
 * @brief Generates a BeamSection with 6x6 mass and stiffness matrices for a hollow circular
 * cross-section
 *
 * This is a convenience function specifically for hollow circular sections commonly
 * used in wind turbine towers. It calculates the geometric properties and then
 * calls the general GenerateMassMatrix and GenerateStiffnessMatrix functions.
 *
 * @param s Normalized position along beam element [dimensionless]
 * @param E Young's modulus [Force/Length²]
 * @param G Shear modulus [Force/Length²]
 * @param rho Material density [Mass/Length³]
 * @param outer_diameter Outer diameter of the hollow circle [Length]
 * @param wall_thickness Wall thickness [Length]
 * @param nu Poisson's ratio [dimensionless]
 * @param x_C x-coordinate of elastic centroid relative to reference point [Length]
 * @param y_C y-coordinate of elastic centroid relative to reference point [Length]
 * @param theta_p Rotation angle from reference axes to principal bending axes [radians]
 * @param x_S x-coordinate of shear center relative to reference point [Length]
 * @param y_S y-coordinate of shear center relative to reference point [Length]
 * @param theta_s Rotation angle from reference axes to principal shear axes [radians]
 * @param x_G x-coordinate of center of gravity relative to reference point [Length]
 * @param y_G y-coordinate of center of gravity relative to reference point [Length]
 * @param theta_i Rotation angle from reference axes to principal inertia axes [radians]
 *
 * @return BeamSection object containing stiffness matrix, mass matrix, and position
 *
 * @note For hollow circular sections, the elastic centroid, shear center, and center
 *       of gravity all coincide at the geometric center, so offset parameters are typically zero.
 * @note Principal axes align with the reference axes due to circular symmetry.
 */
static BeamSection GenerateHollowCircleSection(
    double s, double E, double G, double rho, double outer_diameter, double wall_thickness,
    double nu, double x_C = 0., double y_C = 0., double theta_p = 0., double x_S = 0.,
    double y_S = 0., double theta_s = 0., double x_G = 0., double y_G = 0., double theta_i = 0
) {
    auto properties = CalculateHollowCircleProperties(outer_diameter, wall_thickness, nu);

    // Calculate mass properties
    const double m = rho * properties.area;   // Mass per unit length
    const double I_x = rho * properties.Ixx;  // Mass moment of inertia about x-axis
    const double I_y = rho * properties.Iyy;  // Mass moment of inertia about y-axis
    const double I_p = I_x + I_y;             // Polar mass moment of inertia

    // Calculate stiffness properties
    const double EA = E * properties.area;
    const double EI_x = E * properties.Ixx;
    const double EI_y = E * properties.Iyy;
    const double GKt = G * properties.J;
    const double GA = G * properties.area;

    // Generate mass and stiffness matrices and return BeamSection
    const auto mass = GenerateMassMatrix(m, I_x, I_y, I_p, x_G, y_G, theta_i);
    const auto stiffness = GenerateStiffnessMatrix(
        EA, EI_x, EI_y, GKt, GA, properties.kx, properties.ky, x_C, y_C, theta_p, x_S, y_S, theta_s
    );

    return {s, mass, stiffness};
}

}  // namespace kynema::beams
