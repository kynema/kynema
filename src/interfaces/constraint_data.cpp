#include "constraint_data.hpp"

#include <ranges>

#include "host_constraints.hpp"

namespace kynema::interfaces {

void ConstraintData::GetLoads(const HostConstraints<DeviceType>& host_constraints) {
    for (auto component : std::views::iota(0U, 6U)) {
        this->loads[component] = host_constraints.loads(this->id, component);
    }
}

void ConstraintData::GetOutputs(const HostConstraints<DeviceType>& host_constraints) {
    for (auto component : std::views::iota(0U, 3U)) {
        this->outputs[component] = host_constraints.output(this->id, component);
    }
}

void ConstraintData::SetInputs(HostConstraints<DeviceType>& host_constraints) const {
    for (auto component : std::views::iota(0U, 6U)) {
        host_constraints.input(this->id, component) = this->inputs[component];
    }
}

}  // namespace kynema::interfaces
