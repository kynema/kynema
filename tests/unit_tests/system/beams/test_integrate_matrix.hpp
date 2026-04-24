#include <Kokkos_Core.hpp>

namespace kynema_fmb::tests {

template <size_t n_elem, size_t n_nodes>
auto get_node_state_indices() {
    using IndicesView = Kokkos::View<size_t[n_elem][n_nodes]>;
    auto indices = IndicesView("node_state_indices");
    const auto host_indices = Kokkos::create_mirror_view(indices);
    for (auto element : std::views::iota(0U, n_elem)) {
        for (auto node : std::views::iota(0U, n_nodes)) {
            host_indices(element, node) = element * n_nodes + node;
        }
    }
    Kokkos::deep_copy(indices, host_indices);
    return indices;
}
}  // namespace kynema_fmb::tests
