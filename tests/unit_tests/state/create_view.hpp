#include <Kokkos_Core.hpp>

namespace kynema_fmb::tests {
template <typename ValueType, typename DataType>
typename Kokkos::View<ValueType>::const_type CreateView(
    const std::string& name, const DataType& data
) {
    const auto view = Kokkos::View<ValueType>(Kokkos::view_alloc(name, Kokkos::WithoutInitializing));
    const auto host = typename Kokkos::View<ValueType, Kokkos::HostSpace>::const_type(data.data());
    const auto mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(mirror, host);
    Kokkos::deep_copy(view, mirror);
    return view;
}

template <typename ValueType, typename DataType>
typename Kokkos::View<ValueType> CreateMutableView(const std::string& name, const DataType& data) {
    const auto view = Kokkos::View<ValueType>(Kokkos::view_alloc(name, Kokkos::WithoutInitializing));
    const auto host = typename Kokkos::View<ValueType, Kokkos::HostSpace>::const_type(data.data());
    const auto mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(mirror, host);
    Kokkos::deep_copy(view, mirror);
    return view;
}

template <typename ValueType, typename DataType>
typename Kokkos::View<ValueType, Kokkos::LayoutLeft>::const_type CreateLeftView(
    const std::string& name, const DataType& data
) {
    const auto view = Kokkos::View<ValueType, Kokkos::LayoutLeft>(
        Kokkos::view_alloc(name, Kokkos::WithoutInitializing)
    );
    const auto host = typename Kokkos::View<ValueType, Kokkos::HostSpace>::const_type(data.data());
    const auto mirror = Kokkos::create_mirror_view(Kokkos::WithoutInitializing, view);
    Kokkos::deep_copy(mirror, host);
    Kokkos::deep_copy(view, mirror);
    return view;
}

}  // namespace kynema_fmb::tests
