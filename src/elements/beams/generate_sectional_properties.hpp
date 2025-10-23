#pragma once

#include <array>
#include <cmath>

namespace kynema::beams {

/**
 * @brief Generates a 6x6 cross-sectional stiffness matrix for use in beam elements
 *
 * This function constructs the cross-sectional stiffness matrix that relates the generalized
 * strains to the generalized forces/moments in a beam cross-section. The matrix accounts for
 * coupling between axial, bending, shear, and torsional deformations due to offset between
 * the elastic centroid and shear center locations. Returns a 6x6 stiffness matrix at the cross
 * section origin 'O', based on inputs at the centroid 'C' and shear center 'S' locations.
 *
 * @param EA Axial stiffness (E×A) [Force]
 * @param EI_x Bending stiffness around local x-axis in principal frame and at the
               centroid 'C' [Force × Length²]
 * @param EI_y Bending stiffness around local y-axis in principal frame and at the
               centroid 'C' [Force × Length²]
 * @param GKt Torsional stiffness around principal frame and at shear center 'S'
              (G×J) [Force × Length²]
 * @param GA Shear stiffness around principal frame and at shear center 'S'
              (G×A) [Force]
 * @param kxs Shear correction factor in x-direction [dimensionless]
 * @param kys Shear correction factor in y-direction [dimensionless]
 * @param x_C x-coordinate of elastic centroid relative to reference point i.e. distance
             FROM 'O' TO 'C' [Length]
 * @param y_C y-coordinate of elastic centroid relative to reference point i.e. distance
             FROM 'O' TO 'C' [Length]
 * @param theta_p Rotation angle (around z) FROM reference axes TO principal bending axes [radians]
 * @param x_S x-coordinate of shear center relative to reference point i.e. distance
             FROM 'O' TO 'S' [Length]
 * @param y_S y-coordinate of shear center relative to reference point i.e. distance
             FROM 'O' TO 'S' [Length]
 * @param theta_s Rotation angle (around z) FROM reference axes TO principal shear axes [radians]
 *
 * @return 6x6 cross-sectional stiffness matrix
 */
static std::array<std::array<double, 6>, 6> GenerateStiffnessMatrix(
    double EA, double EI_x, double EI_y, double GKt, double GA, double kxs, double kys, double x_C,
    double y_C, double theta_p, double x_S, double y_S, double theta_s
) {
    const double cos_theta_p = std::cos(theta_p);
    const double sin_theta_p = std::sin(theta_p);
    const double cos_theta_s = std::cos(theta_s);
    const double sin_theta_s = std::sin(theta_s);

    // Calculate bending stiffness in principal frame
    const double bending_xx =
        (EI_x * cos_theta_p * cos_theta_p) + (EI_y * sin_theta_p * sin_theta_p);
    const double bending_yy =
        (EI_x * sin_theta_p * sin_theta_p) + (EI_y * cos_theta_p * cos_theta_p);
    const double bending_xy = (EI_y - EI_x) * sin_theta_p * cos_theta_p;

    // Calculate shear stiffness in principal frame
    const double shear_xx = GA * (kxs * cos_theta_s * cos_theta_s + kys * sin_theta_s * sin_theta_s);
    const double shear_yy = GA * (kxs * sin_theta_s * sin_theta_s + kys * cos_theta_s * cos_theta_s);
    const double shear_xy = GA * (kys - kxs) * sin_theta_s * cos_theta_s;

    //--------------------------------------------------------------------------
    // Assemble stiffness matrix by blocks
    //--------------------------------------------------------------------------
    std::array<std::array<double, 6>, 6> stiffness_matrix = {};  // initialized to zeros

    // Shear-shear coupling (rows 0-1, cols 0-1)
    stiffness_matrix[0][0] = shear_xx;
    stiffness_matrix[0][1] = -shear_xy;
    stiffness_matrix[1][0] = -shear_xy;
    stiffness_matrix[1][1] = shear_yy;

    // Axial stiffness (row 2, col 2)
    stiffness_matrix[2][2] = EA;

    // Axial-bending coupling due to centroid offset (row 2, cols 3-4)
    stiffness_matrix[2][3] = EA * y_C;
    stiffness_matrix[2][4] = -EA * x_C;

    // Bending-axial coupling due to centroid offset (rows 3-4, col 2)
    stiffness_matrix[3][2] = EA * y_C;
    stiffness_matrix[4][2] = -EA * x_C;

    // Bending stiffness with axial-bending coupling (rows 3-4, cols 3-4)
    stiffness_matrix[3][3] = bending_xx + EA * y_C * y_C;
    stiffness_matrix[3][4] = -bending_xy - EA * x_C * y_C;
    stiffness_matrix[4][3] = -bending_xy - EA * x_C * y_C;
    stiffness_matrix[4][4] = bending_yy + EA * x_C * x_C;

    // Shear-torsion coupling due to shear center offset (rows 0-1, col 5)
    stiffness_matrix[0][5] = -shear_xx * y_S - shear_xy * x_S;
    stiffness_matrix[1][5] = shear_xy * y_S + shear_yy * x_S;

    // Torsion-shear coupling due to shear center offset (row 5, cols 0-1)
    stiffness_matrix[5][0] = -shear_xx * y_S - shear_xy * x_S;
    stiffness_matrix[5][1] = shear_xy * y_S + shear_yy * x_S;

    // Torsional stiffness with shear-torsion coupling (row 5, col 5)
    stiffness_matrix[5][5] =
        GKt + shear_xx * y_S * y_S + 2. * shear_xy * x_S * y_S + shear_yy * x_S * x_S;

    return stiffness_matrix;
}

/**
 * @brief Generates a 6x6 cross-sectional mass matrix for use in beam elements.
 *
 * This function constructs the cross-sectional mass matrix that relates the generalized
 * accelerations to the generalized inertial forces in a beam cross-section. Returns the
 * mass matrix at a given point 'O' and w.r.t. given orientation axes based on the values
 * at the center of gravity 'G' and in the inertial axis frame.
 *
 * @param m Mass per unit length [Mass/Length]
 * @param I_x Mass moment of inertia about local x-axis in principal/inertial frame and at
              center of gravity G [Mass×Length]
 * @param I_y Mass moment of inertia about local y-axis in principal/inertial frame and at
              center of gravity G [Mass×Length]
 * @param I_p Polar mass moment of inertia in principal/inertial frame and at center of
              gravity G (I_x + I_y) [Mass×Length]
 * @param x_G x-coordinate of center of gravity relative to reference point i.e. distance
 *            FROM O TO G [Length]
 * @param y_G y-coordinate of center of gravity relative to reference point i.e. distance
 *            FROM O TO G [Length]
 * @param theta_i Rotation angle (around z) FROM reference axes TO principal inertial axes [radians]
 *
 * @return 6x6 cross-sectional mass matrix
 */
static std::array<std::array<double, 6>, 6> GenerateMassMatrix(
    double m, double I_x, double I_y, double I_p, double x_G = 0., double y_G = 0.,
    double theta_i = 0.
) {
    const double cos_theta_i = std::cos(theta_i);
    const double sin_theta_i = std::sin(theta_i);

    // Calculate mass moments of inertia in reference frame
    const double inertia_xx = (I_x * cos_theta_i * cos_theta_i) + (I_y * sin_theta_i * sin_theta_i);
    const double inertia_yy = (I_x * sin_theta_i * sin_theta_i) + (I_y * cos_theta_i * cos_theta_i);
    const double inertia_xy = (I_y - I_x) * sin_theta_i * cos_theta_i;

    //--------------------------------------------------------------------------
    // Assemble mass matrix by blocks
    //--------------------------------------------------------------------------
    std::array<std::array<double, 6>, 6> mass_matrix = {};  // initialized to zeros

    // Translational mass with CG coupling (rows 0-2, cols 0-2)
    mass_matrix[0][0] = m;
    mass_matrix[1][1] = m;
    mass_matrix[2][2] = m;

    // Translational-rotational coupling due to CG offset (rows 0-2, cols 3-5)
    mass_matrix[0][5] = -m * y_G;  // F_x due to α_z
    mass_matrix[1][5] = m * x_G;   // F_y due to α_z
    mass_matrix[2][3] = m * y_G;   // F_z due to α_x
    mass_matrix[2][4] = -m * x_G;  // F_z due to α_y

    // Rotational-translational coupling due to CG offset (rows 3-5, cols 0-2)
    mass_matrix[3][2] = m * y_G;   // M_x due to a_z
    mass_matrix[4][2] = -m * x_G;  // M_y due to a_z
    mass_matrix[5][0] = -m * y_G;  // M_z due to a_x
    mass_matrix[5][1] = m * x_G;   // M_z due to a_y

    // Rotational inertia with translational-rotational coupling (rows 3-4, cols 3-4)
    mass_matrix[3][3] = inertia_xx + m * y_G * y_G;
    mass_matrix[3][4] = -inertia_xy - m * x_G * y_G;
    mass_matrix[4][3] = -inertia_xy - m * x_G * y_G;
    mass_matrix[4][4] = inertia_yy + m * x_G * x_G;

    // Polar inertia with CG coupling (row 5, col 5)
    mass_matrix[5][5] = I_p + m * (x_G * x_G + y_G * y_G);

    return mass_matrix;
}

}  // namespace kynema::beams
