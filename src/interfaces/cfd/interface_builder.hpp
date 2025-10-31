#pragma once

#include "interface.hpp"
#include "interface_input.hpp"

namespace kynema::interfaces::cfd {

/**
 * @brief A factory for configuring and building a CFD interface object.
 *
 * @details every method returns a reference to this InterfaceBuilder for
 * chaining multiple method calls together in a single statement.
 */
struct InterfaceBuilder {
    /**
     * @brief Sets the gravity vector for the problem
     *
     * @param gravity The gravity vector
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetGravity(const std::array<double, 3>& gravity);

    /**
     * @brief Sets the maximum number of nonlinear iterations per time step
     *
     * @param max_iter the maximum number of nonlinear iterations
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMaximumNonlinearIterations(size_t max_iter);

    /**
     * @brief Sets the time step size
     *
     * @param time_step the time step size
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetTimeStep(double time_step);

    /**
     * @brief Sets the numerical damping factor used by the generalized alpha solver
     *
     * @param rho_inf the numerical damping factor
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetDampingFactor(double rho_inf);

    /**
     * @brief Sets if the floating platform is enabled
     *
     * @param enable If the flating platform is enabled
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& EnableFloatingPlatform(bool enable);

    /**
     * @brief Sets the position of the platform
     *
     * @param p The position/orientation of the platform as a quatnernion
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetFloatingPlatformPosition(const std::array<double, 7>& p);

    /**
     * @brief Sets the velocity of the platform
     *
     * @param v The velocity/angular velocity of the platform
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetFloatingPlatformVelocity(const std::array<double, 6>& v);

    /**
     * @brief Sets the acceleration of the platform
     *
     * @param a The acceleration/angular acceleration of the platform
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetFloatingPlatformAcceleration(const std::array<double, 6>& a);

    /**
     * @brief Sets the mass matrix to represent the platform as a point mass
     *
     * @param mass_matrix The mass matrix with the mass and inertia information
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetFloatingPlatformMassMatrix(
        const std::array<std::array<double, 6>, 6>& mass_matrix
    );

    /**
     * @brief Sets the number of mooring lines and sizes the appropriate data structures
     *
     * @param number the number of mooring lines
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetNumberOfMooringLines(size_t number);

    /**
     * @brief Sets the stiffness of a given mooring line
     *
     * @param line the mooring line number to be set
     * @param stiffness the spring stiffness of the line
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineStiffness(size_t line, double stiffness);

    /**
     * @brief Sets the undeformed length of the mooring line
     *
     * @param line the mooring line number to be set
     * @param length the underformed length of the line
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineUndeformedLength(size_t line, double length);

    /**
     * @brief Sets the position of the fairlead node of the mooring line
     *
     * @param line the mooring line number to be set
     * @param p the position of the fairlead node
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineFairleadPosition(size_t line, const std::array<double, 3>& p);

    /**
     * @brief Sets the position of the anchor node of the mooring line
     *
     * @param line the mooring line number to be set
     * @param p the position of the anchor node
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineAnchorPosition(size_t line, const std::array<double, 3>& p);

    /**
     * @brief Sets the velocity of the anchor node of the mooring line
     *
     * @param line the mooring line number to be set
     * @param v the velocity of the anchor node
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineAnchorVelocity(size_t line, const std::array<double, 3>& v);

    /**
     * @brief Sets the acceleration of the anchor node of the mooring line
     *
     * @param line the mooring line number to be set
     * @param a the acceleration of the anchor node
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetMooringLineAnchorAcceleration(size_t line, const std::array<double, 3>& a);

    /**
     * @brief Sets the output file name
     *
     * @param path the output file name
     * @return A reference to this InterfaceBuilder
     */
    InterfaceBuilder& SetOutputFile(const std::string& path);

    /**
     * @brief Sets the output state prefixes to be written to the output file
     *
     * @param prefixes Vector of variable name prefixes
     * @return A reference to this InterfaceBuilder object to allow chaining
     */
    InterfaceBuilder& SetOutputStatePrefixes(const std::vector<std::string>& prefixes);

    /**
     * @brief Sets the number of output timesteps to buffer before flushing to disk
     *
     * @param buffer_size Number of timesteps to buffer
     * @return A reference to this InterfaceBuilder object to allow chaining
     */
    InterfaceBuilder& SetOutputBufferSize(std::size_t buffer_size);

    /**
     * @brief Builds the Interface based on current settings
     *
     * @return The constructed Interface
     */
    [[nodiscard]] Interface Build() const;

private:
    InterfaceInput interface_input;
};

}  // namespace kynema::interfaces::cfd
