#pragma once

#include <array>
#include <vector>

#include "interfaces/node_data.hpp"

namespace kynema {
class Model;
struct BeamSection;
}  // namespace kynema

namespace kynema::interfaces::components {

struct BeamInput;

/**
 * @brief Represents a turbine blade with nodes, elements, and constraints
 *
 * This class is responsible for creating and managing a blade based on input
 * specifications. It handles the creation of nodes, beam elements, and constraints
 * within the provided model.
 */
class Beam {
public:
    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /// @brief Placeholder node ID value for uninitialized components
    static constexpr size_t invalid_id{9999999};

    /// @brief Maximum number of points allowed in blade geometry definition
    static constexpr size_t kMaxGeometryPoints{10};

    /// @brief Beam element ID
    size_t beam_element_id{invalid_id};

    /// @brief Blade node data
    std::vector<NodeData> nodes;

    /// @brief Constraint ID of prescribed root displacement
    size_t prescribed_root_constraint_id{invalid_id};

    /// @brief Location of nodes in blade element [-1, 1]
    std::vector<double> node_xi;

    Beam() = default;

    /**
     * @brief Construct a new Blade using the provided input and model
     * @param input Configuration for the blade
     * @param model Model to which the blade elements and nodes will be added
     * @throws std::invalid_argument If the input configuration is invalid
     */
    Beam(const BeamInput& input, Model& model);

    /**
     * @brief Returns a vector of weights for distributing a point load to the nodes
     * @param s Position [0,1] along the blade
     * @return Vector of weights for each node
     */
    [[nodiscard]] std::vector<double> GetNodeWeights(double s) const;

    /**
     * @brief Adds a point load (Fx, Fy, Fz, Mx, My, Mz) to the blade at location 's' [0, 1]
     * along the material axis
     * @param s Position [0,1] along the blade material axis
     * @param loads Forces and moments (Fx, Fy, Fz, Mx, My, Mz)
     * @throws std::invalid_argument If position is outside valid range
     */
    void AddPointLoad(double s, std::array<double, 6> loads);

    /// @brief Sets blade point loads to zero
    void ClearLoads();

    /// @brief Populate node motion based on host state
    /// @param host_state Host state containing position, displacement, velocity, and acceleration
    void GetMotion(const HostState<DeviceType>& host_state);

    /// @brief Update the host state with current node forces and moments
    /// @param host_state Host state to update
    void SetLoads(HostState<DeviceType>& host_state) const;

private:
    std::vector<std::array<double, 3>> node_coordinates;  ///< Node coordinates
    std::vector<std::array<double, 3>> node_tangents;     ///< Node tangents

    /**
     * @brief Validate the input configuration
     * @param input Blade input configuration
     * @throws std::invalid_argument If configuration is invalid
     */
    static void ValidateInput(const BeamInput& input);
    /**
     * @brief Setup node locations based on input configuration
     * @param input Blade input configuration
     */
    void SetupNodeLocations(const BeamInput& input);

    /**
     * @brief Create node geometry from reference axis points
     * @param input Blade input configuration
     */
    void CreateNodeGeometry(const BeamInput& input);
    /**
     * @brief Create beam element in the model
     * @param input Blade input configuration
     * @param model Model to which the beam element will be added
     */
    void CreateBeamElement(const BeamInput& input, Model& model);
    /**
     * @brief Position the blade in space according to root properties
     * @param input Blade input configuration
     * @param model Model containing the blade
     */
    void PositionBladeInSpace(const BeamInput& input, Model& model) const;

    /**
     * @brief Setup the root node and constraints
     * @param input Blade input configuration
     * @param model Model to which constraints will be added
     */
    void SetupRootNode(const BeamInput& input, Model& model);

    /**
     * @brief Calculate tangent vectors at each node
     */
    void CalcNodeTangents();

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections(const BeamInput& input);

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections_SegmentedGLL(const BeamInput& input);

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections_SegmentedGL(const BeamInput& input);

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections_WholeBeam(
        const BeamInput& input, std::span<const double> quad_locations
    );

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections_WholeBeamGLL(const BeamInput& input);

    /**
     * @brief Create beam sections from input configuration
     * @param input Blade input configuration
     * @return Vector of beam sections
     */
    static std::vector<BeamSection> BuildBeamSections_WholeBeamGL(const BeamInput& input);
};

}  // namespace kynema::interfaces::components
