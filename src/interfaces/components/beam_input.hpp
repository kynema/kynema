#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace kynema::interfaces::components {

/**
 * @brief Reference axis definition for a beam
 *
 * Defines the centerline of the beam with position and twist information
 */
struct ReferenceAxis {
    /// @brief Normalized coordinate locations [0, 1] along the beam
    std::vector<double> coordinate_grid;

    /// @brief X,Y,Z coordinates of reference axis points
    std::vector<std::array<double, 3>> coordinates;

    /// @brief Normalized twist locations [0, 1] along the beam
    std::vector<double> twist_grid;

    /// @brief Structural twist values (in radians)
    std::vector<double> twist;
};

/**
 * @brief Sectional structural properties of the beam
 *
 * Defines the section properties at a given location along the beam
 */
struct Section {
    /// @brief Normalized section location [0, 1] along the beam
    double location;

    /// @brief Mass matrix (6x6) at the section
    std::array<std::array<double, 6>, 6> mass_matrix;

    /// @brief Stiffness matrix (6x6) at the section
    std::array<std::array<double, 6>, 6> stiffness_matrix;

    /**
     * @brief Construct a new Section
     * @param loc Normalized location [0,1] along the beam
     * @param m Mass matrix
     * @param k Stiffness matrix
     */
    Section(
        double loc, std::array<std::array<double, 6>, 6> m, std::array<std::array<double, 6>, 6> k
    )
        : location(loc), mass_matrix(m), stiffness_matrix(k) {
        // Check that the section location is in range [0,1]
        if (loc < 0. || loc > 1.) {
            throw std::invalid_argument("Section location must be in range [0, 1]");
        }
    }
};

/**
 * @brief Root definition for a turbine beam
 *
 * Defines the root properties of the beam
 */
struct Root {
    /// @brief Flag to use root motion as an input to the model
    bool prescribe_root_motion{false};

    /// @brief Inital position/orientation [x, y, z, qw, qx, qy, qz]
    std::array<double, 7> position{0., 0., 0., 1., 0., 0., 0.};

    /// @brief Initial translation/rotational velocity [vx, vy, vz, wx, wy, wz]
    std::array<double, 6> velocity{0., 0., 0., 0., 0., 0.};

    /// @brief Initial translation/rotational acceleration [ax, ay, az, αx, αy, αz]
    std::array<double, 6> acceleration{0., 0., 0., 0., 0., 0.};
};

/**
 * @brief Complete input specification for a beam
 *
 * Defines the input configuration for a turbine beam
 */
struct BeamInput {
    enum class QuadratureRule : std::uint8_t {
        GaussLegendre,
        GaussLobatto
    };
    /// @brief Spectral element order (num nodes - 1)
    size_t element_order{10};

    /// @brief Quadrature Rule
    QuadratureRule quadrature_rule{QuadratureRule::GaussLobatto};

    /// @brief Quadrature point refinement
    size_t section_refinement{0};

    /// @brief Structural reference axis data
    ReferenceAxis ref_axis;

    /// @brief Blade root
    Root root;

    /// @brief Section properties
    std::vector<Section> sections;
};

}  // namespace kynema::interfaces::components
