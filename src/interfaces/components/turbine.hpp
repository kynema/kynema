#pragma once

#include <span>
#include <vector>

#include "interfaces/components/beam.hpp"
#include "interfaces/components/turbine_input.hpp"
#include "interfaces/constraint_data.hpp"
#include "interfaces/host_state.hpp"
#include "interfaces/node_data.hpp"

namespace kynema {
class Model;
}

namespace kynema::interfaces::components {
/**
 * @brief Represents a turbine with nodes, elements, and constraints
 *
 * This class is responsible for creating and managing a turbine based on the
 * turbine input specifications. It handles the creation of the beam elements,
 * mass elements, nodes, and constraints within the provided model.
 *
 * --------------------------------------------------------------------------
 * Node structure
 * --------------------------------------------------------------------------
 * The turbine assembly consists of multiple interconnected nodes that
 * represent different physical components and their kinematic relationships.
 *
 * Tower Nodes
 *  - Tower nodes: Beam nodes distributed along tower height (1, 2, ..., n)
 *  - Tower base: Fixed constraint point at tower foundation (first tower node)
 *  - Tower top: Connection point to nacelle assembly (last tower node)
 *
 *     ┌─── Tower top node (connection to nacelle)
 *     │
 *     ○ <- Tower node n
 *     │
 *     ○ <- Tower node n-1
 *     |
 *     |
 *     │
 *     ○ <- Tower node 2
 *     │
 *     ○ <- Tower node 1 (Tower base node - fixed constraint)
 * -------------
 * / / / / / / /  <--  Ground / Foundation
 *
 * Nacelle/Drivetrain Nodes
 *  - Yaw bearing node: Located at tower top, allows nacelle yaw rotation
 *  - Shaft base node: Base of the main shaft within nacelle
 *  - Azimuth node: Intermediate node for rotor azimuth positioning
 *  - Hub node: Center of mass of the rotating hub assembly
 *
 *             Yaw bearing      Shaft base      Azimuth          Hub
 *                node            node           node            node
 *                 ● -------------- ● ------------ ● ------------ ● -------------
 *    yaw control  |     rigid          torque          rigid          rigid
 *     rotation    |   connection       control       connection     connection
 *   about Z-axis  |                  via revolute                    to blades
 *                                      joint
 * Blade Assembly Nodes
 *  - Blade apex nodes: Connection points between hub and blade roots (one per blade)
 *  - Blade structural nodes: Beam nodes along each blade span
 *
 *                   pitch axis,
 *                rotation control |<----- Blade nodes ----->|
 *     ● Blade apex -------------  ●  --------  ●  --------  ●
 *       node                    root                     tip
 *        │                      node                     node
 *        │ rigid connection
 *        │
 *        ● Hub node
 * --------------------------------------------------------------------------
 * Kinematic chain/constraints
 * --------------------------------------------------------------------------
 * The nodes are connected in a kinematic chain that represents the turbine's
 * degrees of freedom.
 *
 * - tower base node: Fixed boundary condition
 * - tower top node <-> yaw bearing node: Yaw rotation control
 * - yaw bearing node <-> shaft base node: Rigid joint
 * - shaft base node <-> azimuth node: Revolute joint with torque control
 * - hub <-> blade apex nodes: Rigid joint
 * - blade apex nodes <-> blade root nodes: Pitch rotation control
 */
class Turbine {
public:
    //--------------------------------------------------------------------------
    // Types and Constants
    //--------------------------------------------------------------------------

    using DeviceType =
        Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space>;

    /// @brief Placeholder node ID value for uninitialized components
    static constexpr size_t invalid_id{9999999};

    /// @brief Minimum valid hub diameter
    static constexpr double kMinHubDiameter{1e-8};

    /// @brief Tolerance for near zero comparisons
    static constexpr double kZeroTolerance{1e-12};

    //--------------------------------------------------------------------------
    // Elements
    //--------------------------------------------------------------------------

    std::vector<Beam> blades;                        //< Blades in the turbine
    Beam tower;                                      //< Tower in the turbine
    size_t hub_mass_element_id{invalid_id};          //< Hub mass element ID
    size_t yaw_bearing_mass_element_id{invalid_id};  //< Yaw bearing mass element ID

    //--------------------------------------------------------------------------
    // Nodes
    //--------------------------------------------------------------------------

    std::vector<NodeData> apex_nodes;  //< Blade root nodes
    NodeData hub_node;                 //< Hub node
    NodeData azimuth_node;             //< Azimuth node
    NodeData shaft_base_node;          //< Shaft base node
    NodeData nacelle_cm_node;          //< Nacelle center of mass node
    NodeData yaw_bearing_node;         //< Yaw bearing node

    //--------------------------------------------------------------------------
    // Constraints
    //--------------------------------------------------------------------------

    ConstraintData tower_base;                 //< Tower base constraint
    ConstraintData tower_top_to_yaw_bearing;   //< Tower top to yaw bearing constraint
    ConstraintData yaw_bearing_to_shaft_base;  //< Yaw bearing to shaft base constraint
    ConstraintData yaw_bearing_to_nacelle_cm;  //< Yaw bearing to nacelle center of mass constraint
    ConstraintData shaft_base_to_azimuth;      //< Nacelle mass to generator constraint
    ConstraintData azimuth_to_hub;             //< Azimuth to hub constraint
    std::vector<ConstraintData> blade_pitch;   //< Blade root to apex constraints
    std::vector<ConstraintData> apex_to_hub;   //< Apex to hub constraints

    //--------------------------------------------------------------------------
    // Control inputs
    //--------------------------------------------------------------------------

    std::vector<double> blade_pitch_control;  //< Blade pitch angles
    double torque_control{0.};                //< Torque control value
    double yaw_control{0.};                   //< Yaw control value

    /**
     * @brief Constructs a turbine with the specified input configuration
     *
     * @details Creates a complete turbine structure including blades, tower, hub, nacelle
     * components,and all associated nodes, elements, constraints, and control systems. The
     * turbine is positioned and oriented according to the input parameters, with proper kinematic
     * relationships established between all components.
     * --------------------------------------------------------------------------
     * Construction sequence
     * --------------------------------------------------------------------------
     *  1. Validation: Input parameters are checked for physical consistency
     *  2. Node creation + positioning: All nodes are created at the origin and
     *     assembled at their correct spatial locations/orientations
     *  3. Mass element addition: Mass and inertia properties are assigned to relevant
     *     nodes/components, e.g. lumped mass elements at hub and nacelle
     *  4. Constraint creation: Kinematic relationships are established between nodes
     *  5. Initial conditions: Initial conditions such as displacements, velocities,
     *     and accelerations are applied as required
     *
     * @param input Turbine configuration parameters
     * @param model Structural model to add turbine components to
     * @throws std::invalid_argument If input configuration is invalid
     */
    Turbine(const TurbineInput& input, Model& model);

    /**
     * @brief Populate node motion from host state
     * @param host_state Host state containing position, displacement, velocity, and acceleration
     */
    void GetMotion(const HostState<DeviceType>& host_state);

    /**
     * @brief Populate constraint loads from host constraints
     * @param host_constraints Host constraints containing constraint forces/moments
     */
    void GetLoads(const HostConstraints<DeviceType>& host_constraints);

    /**
     * @brief Update the host state with current node forces and moments
     * @param host_state Host state to update
     */
    void SetLoads(HostState<DeviceType>& host_state) const;

    /**
     * @brief Get the turbine input configuration
     * @return Turbine input configuration
     */
    [[nodiscard]] const TurbineInput& GetTurbineInput() const;

private:
    //--------------------------------------------------------------------------
    // Node ID collections of turbine components
    //--------------------------------------------------------------------------

    std::vector<size_t> tower_node_ids;       ///< All nodes on tower beam
    std::vector<size_t> drivetrain_node_ids;  ///< All nodes on drivetrain i.e. yaw bearing, shaft
                                              ///< base, azimuth, and hub nodes
    std::vector<size_t> blade_node_ids;       ///< All nodes on blade beams + blade apex nodes
    std::vector<size_t> rotor_node_ids;       ///< All blade + blade apex nodes, hub, azimuth, and
                                              ///< shaft base nodes
    std::vector<size_t> nacelle_node_ids;     ///< drivetrain_node_ids + blade_node_ids
    std::vector<size_t>
        all_turbine_node_ids;  ///< tower_node_ids + drivetrain_node_ids + blade_node_ids

    //--------------------------------------------------------------------------
    // Turbine inputs
    //--------------------------------------------------------------------------

    TurbineInput turbine_input;  ///< Turbine input configuration

    /**
     * @brief Create blades from input configuration
     * @param blade_inputs Blade input configurations
     * @param model Model to which the blades will be added
     * @return Vector of blades
     */
    [[nodiscard]] static std::vector<Beam> CreateBlades(
        std::span<const BeamInput> blade_inputs, Model& model
    );

    /**
     * @brief Validates turbine input parameters before construction
     *
     * @details Checks that all input parameters are within physically reasonable ranges
     * and that required parameters are properly specified. Throws exceptions
     * for invalid configurations that would lead to simulation errors.
     *
     * @param input Turbine configuration parameters to validate
     *
     * @throws std::invalid_argument If any parameter is outside valid range or physically
     * inconsistent
     */
    static void ValidateInput(const TurbineInput& input);

    /**
     * @brief Position all turbine nodes in the global coordinate system based on the
     * turbine input configuration
     *
     * This method performs a complex sequence of geometric transformations to position
     * all turbine components (tower, blades, hub, drivetrain) in their final reference
     * locations to achieve the correct turbine geometry.
     *
     * --------------------------------------------------------------------------
     * Transformation/assembly sequence
     * --------------------------------------------------------------------------
     *   - Tower: Rotated to align with global Z-axis (vertical)
     *   - Blades: Aligned with Z-axis -> translated to hub radius -> coned ->
     *     azimuthally positioned
     *   - Drivetrain nodes: Created at intermediate shaft positions
     *   - All rotor components: Tilted by shaft angle and translated to final hub position
     *
     * @param input Turbine configuration containing geometric parameters
     * @param model Structural model to add positioned nodes to
     */
    void PositionNodes(const TurbineInput& input, Model& model);

    /**
     * @brief Creates all intermediate nodes (drivetrain and blade apex nodes) at their initial
     * positions and populates node ID collections
     *
     * @details Creates nodes at their reference positions before any transformations are applied:
     *  - Hub node: At origin with hub CM offset
     *  - Azimuth node: At shaft length from origin
     *  - Shaft base node: At shaft length from origin
     *  - Yaw bearing node: At tower top position
     *  - Blade apex nodes: At origin (one per blade)
     *
     * @param input Turbine configuration containing geometric parameters
     * @param model Structural model to add nodes to
     */
    void CreateIntermediateNodes(const TurbineInput& input, Model& model);

    /**
     * @brief Adds mass elements to the turbine model at the yaw bearing and hub nodes
     *
     * @details Creates lumped mass elements that represent
     *  - Yaw bearing mass element: Combined nacelle system mass and yaw bearing mass
     *    with inertia tensor about tower-top
     *  - Hub mass element: Hub assembly mass and inertia properties
     *
     * @param input Turbine configuration containing inertia matrices
     * @param model Structural model to add mass elements to
     */
    void AddMassElements(const TurbineInput& input, Model& model);
    /**
     * @brief Creates all kinematic constraints and control connections for the turbine
     *
     * @details This method establishes the complete constraint system that defines the
     * turbine's degrees of freedom and control interfaces. Constraints are added in a specific
     * order to build the kinematic chain from the fixed tower base to the controllable
     * rotor and blade systems.
     *
     * --------------------------------------------------------------------------
     * Constraint schematic
     * --------------------------------------------------------------------------
     *
     *           Blade 1                 Blade 2                 Blade 3
     *     ●────────●────────●     ●────────●────────●     ●────────●────────●
     *   root      mid      tip   root      mid      tip   root      mid      tip
     *    │                       │                       │
     *    │ pitch control         │ pitch control         │ pitch control
     *    │ (rotation about       │ (rotation about       │ (rotation about
     *    │  pitch axis)          │  pitch axis)          │  pitch axis)
     *    │                       │                       │
     *    ● Apex 1                ● Apex 2                ● Apex 3
     *    │                       │                       │
     *    │ rigid connection      │ rigid connection      │ rigid connection
     *    └───────────────────────┼───────────────────────┘
     *                            │
     *                            ● Hub node
     *                            │
     *                            │ rigid connection
     *                            │
     *                            ● Azimuth node
     *                            │
     *                            │ revolute joint with
     *                            │ torque control
     *                            │ (shaft rotation)
     *                            │
     *                            ● Shaft base node
     *                            │
     *                            │ rigid connection
     *                            │
     *                            ● Yaw bearing node
     *                            │
     *                            │  yaw rotation control
     *                            │ (rotation about Z-axis)
     *                            │
     *                            ● Tower top node
     *                            │
     *                            │
     *                            ● Tower node n-1
     *                            │
     *                            ● Tower node n-2
     *                            .
     *                            .
     *                            │
     *                            ● Tower node 2
     *                            │
     *                            ● Tower node 1 (base)
     *                            │
     *                            │ fixed boundary condition
     *                            │  (all DOFs constrained)
     *                        ────┴────
     *                       / / / / / /  <- Foundation
     *
     *
     * Constraint Types:
     * - Fixed bc: Constrains all 6 DOFs (tower base to foundation)
     * - Rigid joint: Constrains relative motion between two nodes (6 DOFs)
     * - Revolute joint: Allows rotation about one axis, constrains other 5 DOFs
     * - Rotation control: Controlled rotation about specified axis
     *
     * Control Interfaces:
     * - Yaw control: Controls nacelle orientation (yaw bearing rotation)
     * - Torque control: Controls rotor speed (shaft base to azimuth rotation)
     * - Pitch control: Controls blade pitch angles (apex to root rotation)
     *
     * @param input Turbine configuration containing constraint parameters
     * @param model Structural model to add constraints to
     */
    void AddConstraints(const TurbineInput& input, Model& model);

    /**
     * @brief Set initial velocities, accelerations, and other dynamic initial
     * conditions of the nodes
     * @param input Turbine configuration parameters
     * @param model Structural model
     */
    void SetInitialConditions(const TurbineInput& input, Model& model);
    /**
     * @brief Sets initial displacements for blade pitch, nacelle yaw, and tower base positioning
     *
     * @details This method applies initial displacements to the turbine components after
     * the assembly to reference configuration is complete. These displacements represent the
     * initial state of the turbine at simulation start.
     *
     * --------------------------------------------------------------------------
     * Initial displacement sequence
     * --------------------------------------------------------------------------
     *   - Blade pitch: Rotated about pitch axes from reference (zero pitch) -> initial angle
     *   - Nacelle yaw: Rotated about tower top from reference orientation -> initial yaw
     *   - Tower base: Translated and rotated from default position -> final location
     *
     * All rotations are applied as displacements rather than reference position changes,
     * ensuring proper initial conditions for the simulation.
     *
     * @param input Turbine configuration containing initial displacement parameters
     * @param model Structural model containing the turbine nodes to be displaced
     */
    void SetInitialDisplacements(const TurbineInput& input, Model& model);

    /**
     * @brief Set initial rotational velocity about the shaft axis for rotor components
     *
     * @details This method applies initial rotational velocity to all rotor components
     * including the hub, azimuth node, blade nodes, and apex nodes. The velocity is
     * calculated based on the specified rotor speed and shaft tilt angle.
     *
     * @param input Turbine configuration containing rotor speed
     * @param model Structural model
     */
    void SetInitialRotorVelocity(const TurbineInput& input, Model& model);
};

}  // namespace kynema::interfaces::components
