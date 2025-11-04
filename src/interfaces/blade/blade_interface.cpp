#include "blade_interface.hpp"

#include <filesystem>

#include "state/clone_state.hpp"
#include "state/copy_state_data.hpp"
#include "step/step.hpp"

namespace kynema::interfaces {

BladeInterface::BladeInterface(
    const components::SolutionInput& solution_input, const components::BeamInput& blade_input,
    const components::OutputsConfig& outputs_config
)
    : model(Model(solution_input.gravity)),
      blade(blade_input, model),
      state(model.CreateState<DeviceType>()),
      elements(model.CreateElements<DeviceType>()),
      constraints(model.CreateConstraints<DeviceType>()),
      parameters(
          solution_input.dynamic_solve, solution_input.max_iter, solution_input.time_step,
          solution_input.rho_inf, solution_input.absolute_error_tolerance,
          solution_input.relative_error_tolerance
      ),
      solver(CreateSolver<DeviceType>(state, elements, constraints)),
      state_save(CloneState(state)),
      host_state(state) {
    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the blade motion from state
    this->blade.GetMotion(this->host_state);

    // Initialize NetCDF writer and write mesh connectivity if output is enabled
    if (outputs_config.Enabled()) {
        // Create output directory if it doesn't exist
        std::filesystem::create_directories(outputs_config.output_file_path);

        // Write mesh connectivity to YAML file
        model.ExportMeshConnectivityToYAML(
            outputs_config.output_file_path + "/mesh_connectivity.yaml"
        );

        // Initialize outputs
        this->outputs = std::make_unique<Outputs>(
            outputs_config.output_file_path + "/blade_interface.nc", blade.nodes.size(),
            outputs_config.output_file_path + "/blade_time_series.nc",
            outputs_config.output_state_prefixes, outputs_config.buffer_size
        );
    }
}

components::Beam& BladeInterface::Blade() {
    return this->blade;
}

bool BladeInterface::Step() {
    // Transfer node loads -> state
    for (const auto& node : this->blade.nodes) {
        for (auto component : std::views::iota(0U, 6U)) {
            this->host_state.f(node.id, component) = node.loads[component];
        }
    }
    Kokkos::deep_copy(this->state.f, this->host_state.f);

    // Solve for state at end of step
    auto converged =
        kynema::Step(this->parameters, this->solver, this->elements, this->state, this->constraints);
    if (!converged) {
        return false;
    }

    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the blade motion from state
    this->blade.GetMotion(this->host_state);

    // Write outputs and increment timestep counter
    if (this->outputs) {
        outputs->WriteNodeOutputsAtTimestep(this->host_state, this->state.time_step);
    }

    return true;
}

void BladeInterface::SaveState() {
    CopyStateData(this->state_save, this->state);
}

void BladeInterface::RestoreState() {
    CopyStateData(this->state, this->state_save);

    // Update the host state with current node motion
    this->host_state.CopyFromState(this->state);

    // Update the blade motion from state
    this->blade.GetMotion(this->host_state);
}

void BladeInterface::SetRootDisplacement(const std::array<double, 7>& u) const {
    if (this->blade.prescribed_root_constraint_id == components::Beam::invalid_id) {
        throw std::runtime_error("prescribed root motion was not enabled");
    }
    this->constraints.UpdateDisplacement(this->blade.prescribed_root_constraint_id, u);
}
}  // namespace kynema::interfaces
