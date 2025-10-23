#include "interface.hpp"

#include "floating_platform.hpp"
#include "floating_platform_input.hpp"
#include "interfaces/outputs.hpp"
#include "model/model.hpp"
#include "mooring_line_input.hpp"
#include "node_data.hpp"
#include "state/clone_state.hpp"
#include "state/copy_state_data.hpp"
#include "state/read_state_from_file.hpp"
#include "state/state.hpp"
#include "state/write_state_to_file.hpp"
#include "step/step.hpp"
#include "turbine.hpp"
#include "turbine_input.hpp"

namespace {

/**
 * @brief Populates the provided NodeData object with state information
 *
 * @param node
 * @param host_state_x
 * @param host_state_q
 * @param host_state_v
 * @param host_state_vd
 */
void GetNodeMotion(
    kynema::interfaces::cfd::NodeData& node,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_x,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_q,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_v,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_vd
) {
    for (auto component : std::views::iota(0U, 7U)) {
        node.position[component] = host_state_x(node.id, component);
        node.displacement[component] = host_state_q(node.id, component);
    }
    for (auto component : std::views::iota(0U, 6U)) {
        node.velocity[component] = host_state_v(node.id, component);
        node.acceleration[component] = host_state_vd(node.id, component);
    }
}

/**
 * @brief Adds nodes to the provided model based on the provided input configuration
 *
 * @param input The configuration for the floating platform
 * @param model The Kynema Model to be populated with mass and spring element information
 * @return A FlatingPlatform object based on the provided configuration
 */
kynema::interfaces::cfd::FloatingPlatform CreateFloatingPlatform(
    const kynema::interfaces::cfd::FloatingPlatformInput& input, kynema::Model& model
) {
    // If floating platform is not enabled, return
    if (!input.enable) {
        return {
            .active = false,                                // active
            .node = kynema::interfaces::cfd::NodeData(0U),  // platform node
            .mass_element_id = 0U,                          // mass element ID
            .mooring_lines = {},
        };
    }

    // Construct platform node and save ID
    const auto platform_node_id = model.AddNode()
                                      .SetPosition(input.position)
                                      .SetVelocity(input.velocity)
                                      .SetAcceleration(input.acceleration)
                                      .Build();

    // Add element for platform mass
    const auto mass_element_id = model.AddMassElement(platform_node_id, input.mass_matrix);

    // Instantiate platform
    kynema::interfaces::cfd::FloatingPlatform platform{
        .active = true,  // enable platform
        .node = kynema::interfaces::cfd::NodeData(platform_node_id),
        .mass_element_id = mass_element_id,
        .mooring_lines = {},
    };

    // Construct mooring lines
    std::ranges::transform(
        input.mooring_lines, std::back_inserter(platform.mooring_lines),
        [&](const kynema::interfaces::cfd::MooringLineInput& ml_input) {
            // Add fairlead node
            const auto fairlead_position = std::array{
                ml_input.fairlead_position[0],
                ml_input.fairlead_position[1],
                ml_input.fairlead_position[2],
                1.,
                0.,
                0.,
                0.
            };
            const auto fairlead_velocity = std::array{
                ml_input.fairlead_velocity[0],
                ml_input.fairlead_velocity[1],
                ml_input.fairlead_velocity[2],
                0.,
                0.,
                0.
            };
            const auto fairlead_acceleration = std::array{
                ml_input.fairlead_acceleration[0],
                ml_input.fairlead_acceleration[1],
                ml_input.fairlead_acceleration[2],
                0.,
                0.,
                0.
            };
            auto fairlead_node_id = model.AddNode()
                                        .SetPosition(fairlead_position)
                                        .SetVelocity(fairlead_velocity)
                                        .SetAcceleration(fairlead_acceleration)
                                        .Build();

            // Add anchor node
            const auto anchor_position = std::array{
                ml_input.anchor_position[0],
                ml_input.anchor_position[1],
                ml_input.anchor_position[2],
                1.,
                0.,
                0.,
                0.
            };
            const auto anchor_velocity = std::array{ml_input.anchor_velocity[0],
                                                    ml_input.anchor_velocity[1],
                                                    ml_input.anchor_velocity[2],
                                                    0.,
                                                    0.,
                                                    0.};
            const auto anchor_acceleration = std::array{
                ml_input.anchor_acceleration[0],
                ml_input.anchor_acceleration[1],
                ml_input.anchor_acceleration[2],
                0.,
                0.,
                0.
            };
            auto anchor_node_id = model.AddNode()
                                      .SetPosition(anchor_position)
                                      .SetVelocity(anchor_velocity)
                                      .SetAcceleration(anchor_acceleration)
                                      .Build();

            // Add fixed constraint to anchor node
            auto fixed_constraint_id = model.AddFixedBC3DOFs(anchor_node_id);

            // Add rigid constraint from fairlead node to platform node
            auto rigid_constraint_id =
                model.AddRigidJoint6DOFsTo3DOFs(std::array{platform.node.id, fairlead_node_id});

            // Add spring from fairlead to anchor
            auto spring_element_id = model.AddSpringElement(
                fairlead_node_id, anchor_node_id, ml_input.stiffness, ml_input.undeformed_length
            );

            // Add mooring line data to platform
            return kynema::interfaces::cfd::MooringLine{
                .fairlead_node = kynema::interfaces::cfd::NodeData(fairlead_node_id),
                .anchor_node = kynema::interfaces::cfd::NodeData(anchor_node_id),
                .fixed_constraint_id = fixed_constraint_id,
                .rigid_constraint_id = rigid_constraint_id,
                .spring_element_id = spring_element_id,
            };
        }
    );

    return platform;
}

/**
 * @brief Sets the force vector at the Platform's mass node in the host_state
 *
 * @tparam DeviceType The Kokkos Device where the State corresponding to the HostState resides
 * @param platform The platform object where the forces have been set
 * @param host_state The HostState object where the forces wil be set
 */
template <typename DeviceType>
void SetPlatformLoads(
    const kynema::interfaces::cfd::FloatingPlatform& platform,
    kynema::interfaces::HostState<DeviceType>& host_state
) {
    // Return if platform is not active
    if (!platform.active) {
        return;
    }

    // Set external loads on platform node
    for (auto component : std::views::iota(0U, 6U)) {
        host_state.f(platform.node.id, component) = platform.node.loads[component];
    }
}

/**
 * @brief Sets the state data on the provided FloatingPlatform object
 *
 * @param platform The FloatingPlatform where data will be set
 * @param host_state_x Location data
 * @param host_state_q Displacement data
 * @param host_state_v Velocity data
 * @param host_state_vd Acceleration Data
 */
void GetFloatingPlatformMotion(
    kynema::interfaces::cfd::FloatingPlatform& platform,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_x,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_q,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_v,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_vd
) {
    // If platform is not active, return
    if (!platform.active) {
        return;
    }

    // Populate platform node motion
    GetNodeMotion(platform.node, host_state_x, host_state_q, host_state_v, host_state_vd);

    // Loop through mooring lines
    for (auto& ml : platform.mooring_lines) {
        // Populate fairlead node motion
        GetNodeMotion(ml.fairlead_node, host_state_x, host_state_q, host_state_v, host_state_vd);

        // Populate anchor node motion
        GetNodeMotion(ml.anchor_node, host_state_x, host_state_q, host_state_v, host_state_vd);
    }
}

/**
 * @brief Creates the Turbine object based on the input configuration and populates
 * the provided model
 *
 * @param input The problem configuration object
 * @param model The Model to be populated with element information based on the input
 * @return The full turbine object configured based on input
 */
kynema::interfaces::cfd::Turbine CreateTurbine(
    const kynema::interfaces::cfd::TurbineInput& input, kynema::Model& model
) {
    return {
        CreateFloatingPlatform(input.floating_platform, model),
    };
}

/**
 * @brief Sets the force loads on the provided HostState and State objects
 *
 * @tparam DeviceType The Kokkos Device where the State object resides
 *
 * @param turbine The Turbine object where the forces have been set
 * @param host_state A HostState object corresponding to the provided State
 * @param state A State object where the forces will be set
 */
template <typename DeviceType>
void SetTurbineLoads(
    const kynema::interfaces::cfd::Turbine& turbine,
    kynema::interfaces::HostState<DeviceType>& host_state, kynema::State<DeviceType>& state
) {
    SetPlatformLoads(turbine.floating_platform, host_state);
    Kokkos::deep_copy(state.f, host_state.f);
}

/**
 * @brief Populates the Turbine object with the appropriate state data
 *
 * @param turbine The Turbine object to have its state data populated
 * @param host_state_x Position data
 * @param host_state_q Displacement data
 * @param host_state_v Velocity Data
 * @param host_state_vd Acceleration Data
 */
void GetTurbineMotion(
    kynema::interfaces::cfd::Turbine& turbine,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_x,
    const Kokkos::View<double* [7]>::HostMirror::const_type& host_state_q,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_v,
    const Kokkos::View<double* [6]>::HostMirror::const_type& host_state_vd
) {
    GetFloatingPlatformMotion(
        turbine.floating_platform, host_state_x, host_state_q, host_state_v, host_state_vd
    );
}

}  // namespace
namespace kynema::interfaces::cfd {

Interface::Interface(const InterfaceInput& input)
    : model(input.gravity),
      turbine(CreateTurbine(input.turbine, model)),
      state(model.CreateState<DeviceType>()),
      elements(model.CreateElements<DeviceType>()),
      constraints(model.CreateConstraints<DeviceType>()),
      parameters(true, input.max_iter, input.time_step, input.rho_inf),
      solver(CreateSolver<DeviceType>(state, elements, constraints)),
      state_save(CloneState(state)),
      host_state(state),
      outputs_(nullptr) {
    // Copy state motion members from device to host
    this->host_state.CopyFromState(this->state);
    Kokkos::deep_copy(this->host_state.f, 0.);

    // Update the turbine motion to match restored state
    GetTurbineMotion(
        this->turbine, this->host_state.x, this->host_state.q, this->host_state.v,
        this->host_state.vd
    );

    // Initialize NetCDF writer and write mesh connectivity if output path is specified
    if (!input.output_file.empty()) {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(input.output_file);

        // Initialize outputs
        this->outputs_ = std::make_unique<kynema::interfaces::Outputs>(
            input.output_file + "/cfd_interface.nc", state.num_system_nodes
        );

        // Write mesh connectivity to YAML file
        model.ExportMeshConnectivityToYAML(input.output_file + "/mesh_connectivity.yaml");
    }
}

void Interface::WriteRestart(const std::filesystem::path& filename) const {
    auto output = std::ofstream(filename);
    WriteStateToFile(output, state);
}

void Interface::ReadRestart(const std::filesystem::path& filename) {
    auto input = std::ifstream(filename);
    ReadStateFromFile(input, state);

    this->host_state.CopyFromState(state);

    GetTurbineMotion(
        this->turbine, this->host_state.x, this->host_state.q, this->host_state.v,
        this->host_state.vd
    );
}

bool Interface::Step() {
    // Transfer loads to solver
    SetTurbineLoads(this->turbine, this->host_state, this->state);

    // Solve for state at end of step
    auto converged =
        kynema::Step(this->parameters, this->solver, this->elements, this->state, this->constraints);
    if (!converged) {
        return false;
    }

    // Copy state motion members from device to host
    this->host_state.CopyFromState(this->state);

    // Update the turbine motion
    GetTurbineMotion(
        this->turbine, this->host_state.x, this->host_state.q, this->host_state.v,
        this->host_state.vd
    );

    // Write outputs and increment timestep counter
    if (this->outputs_) {
        outputs_->WriteNodeOutputsAtTimestep(this->host_state, this->current_timestep_);
    }
    this->current_timestep_++;

    return true;
}

void Interface::SaveState() {
    CopyStateData(this->state_save, this->state);
}

void Interface::RestoreState() {
    CopyStateData(this->state, this->state_save);

    // Copy state motion members from device to host
    this->host_state.CopyFromState(this->state);

    // Update the turbine motion to match restored state
    GetTurbineMotion(
        this->turbine, this->host_state.x, this->host_state.q, this->host_state.v,
        this->host_state.vd
    );
}

}  // namespace kynema::interfaces::cfd
