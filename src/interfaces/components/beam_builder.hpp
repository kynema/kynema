#pragma once

#include <cstdint>

#include "beam_input.hpp"

namespace kynema {
class Model;
}

namespace kynema::interfaces::components {

class Beam;

/**
 * @brief Enum to represent reference axis orientation
 */
enum class ReferenceAxisOrientation : std::uint8_t {
    X,  // X-axis
    Z   // Z-axis
};

/**
 * @brief Builder class for creating Blade objects
 *
 * @details Each of the methods for setting options returns a reference to this BeamBuilder object
 * to allow chaining calls together in a single statement.
 */
class BeamBuilder {
public:
    /**
     * @brief Sets the prder of the beam elements
     *
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetElementOrder(size_t element_order);

    /**
     * @brief Sets the Quadrature rule to use for this beam
     *
     * @details There are currently two supported Quadrature rules: Gauss Lobatto (default)
     * and Gauss Legendre.
     *
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetQuadratureRule(BeamInput::QuadratureRule rule);

    /**
     * @brief Sets the Quadrature style to use for this beam
     *
     * @details Segmented applies the quadrature rule between each of the specified
     * sectional properties.  WholeBeam has just one quadrautre applied to the enire beam.
     *
     * Combining GLL quadrature with the Segmented style wil preserve the user's specified
     * sectional properties at their exact locations.  This is expected to give better
     * performance for complex property distrubutions.
     *
     * WholeBeam will likely be the most efficient for simple property distributions.
     *
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetQuadratureStyle(BeamInput::QuadratureStyle style);

    /**
     * @brief Sets the number of section refinements to perform
     *
     * @details This adds additional quadrature points between the supplied sections and sets their
     * physical properties using linear interpolation.
     *
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetSectionRefinement(size_t section_refinement);

    /**
     * @brief Deletes all currently set reference axis points
     *
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& ClearReferenceAxisPoints();

    /**
     * @brief Adds a reference axis point with specified orientation
     * @param grid_location Location [0,1] along the blade
     * @param coordinates Coordinates (x,y,z) of the reference axis point
     * @param ref_axis Reference axis orientation ('X' or 'Z')
     * @return Reference to this builder for method chaining
     */
    BeamBuilder& AddRefAxisPoint(
        double grid_location, const std::array<double, 3>& coordinates,
        ReferenceAxisOrientation ref_axis
    );

    /**
     * @brief Adds a twist about the reference axis at a certain point
     *
     * @details This value will be linearly interpolated to the required locations along the
     * actual beam element
     *
     * @param grid_location The location where the twist is specified
     * @param twist The twist about the reference axis in degrees
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& AddRefAxisTwist(double grid_location, double twist);

    /**
     * @brief sets if this beam will have prescribed root motion or not
     *
     * @param enable If prescirbed root motion is enabled
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& PrescribedRootMotion(bool enable);

    /**
     * @brief Sets the position of this beam's root node
     *
     * @param p The root node position/orientation as a quaternion
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetRootPosition(const std::array<double, 7>& p);

    /**
     * @brief Sets the velocity of this beam's root node
     *
     * @param v The root node velocity/angular velcoity
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetRootVelocity(const std::array<double, 6>& v);

    /**
     * @brief Sets the acceleration of this beam's root node
     *
     * @param a The root node acceleration/angular acceleration
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& SetRootAcceleration(const std::array<double, 6>& a);

    /**
     * @brief Deletes all sections that have been added to the BeamBuilder
     * @return A reference to this BeamBuilder
     */
    BeamBuilder& ClearSections();

    /// @brief Adds sectional location, mass matrix, and stiffness matrix. Can handle reference
    /// axis along X or Z.
    /// @param grid_location section location [0,1] along blade
    /// @param mass_matrix sectional mass matrix
    /// @param stiffness_matrix sectional stiffness matrix
    /// @param ref_axis Reference axis ('X' or 'Z')
    /// @return Reference to this builder for method chaining
    BeamBuilder& AddSection(
        double grid_location, const std::array<std::array<double, 6>, 6>& mass_matrix,
        const std::array<std::array<double, 6>, 6>& stiffness_matrix,
        ReferenceAxisOrientation ref_axis
    );

    /**
     * @brief Get the current blade input configuration
     * @return Reference to the current blade input
     */
    [[nodiscard]] const BeamInput& Input() const;

    /**
     * @brief Build a Blade object from the current configuration
     * @param model The model to associate with this blade
     * @return A new Blade object
     */
    [[nodiscard]] Beam Build(Model& model) const;

private:
    BeamInput input;
};

}  // namespace kynema::interfaces::components
