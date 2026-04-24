#include "node_data.hpp"

#include <ranges>

#include "host_state.hpp"

namespace kynema_fmb::interfaces {
void NodeData::ClearLoads() {
    this->loads = {0., 0., 0., 0., 0., 0.};
}

void NodeData::GetMotion(const HostState<DeviceType>& host_state) {
    for (auto component : std::views::iota(0U, 7U)) {
        this->position[component] = host_state.x(this->id, component);
        this->displacement[component] = host_state.q(this->id, component);
    }
    for (auto component : std::views::iota(0U, 6U)) {
        this->velocity[component] = host_state.v(this->id, component);
        this->acceleration[component] = host_state.vd(this->id, component);
    }
}

void NodeData::SetLoads(HostState<DeviceType>& host_state) const {
    for (auto component : std::views::iota(0U, 6U)) {
        host_state.f(this->id, component) = this->loads[component];
    }
}
}  // namespace kynema_fmb::interfaces
