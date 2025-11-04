#include "blade_interface_builder.hpp"

#include "blade_interface.hpp"

namespace kynema::interfaces {

components::SolutionBuilder& BladeInterfaceBuilder::Solution() {
    return this->solution_builder;
}

components::BeamBuilder& BladeInterfaceBuilder::Blade() {
    return this->beam_builder;
}

components::OutputsBuilder& BladeInterfaceBuilder::Outputs() {
    return this->outputs_builder;
}

BladeInterface BladeInterfaceBuilder::Build() const {
    return BladeInterface(
        this->solution_builder.Input(), this->beam_builder.Input(), this->outputs_builder.Config()
    );
}

}  // namespace kynema::interfaces
