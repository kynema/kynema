#pragma once

#include <Kokkos_Core.hpp>
#include <umfpack.h>

#include "dss_algorithm.hpp"
#include "dss_handle_umfpack.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct NumericFunction<Handle<Algorithm::UMFPACK>, CrsMatrixType> {
    static void numeric(Handle<Algorithm::UMFPACK>& dss_handle, CrsMatrixType& A) {
        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto*& symbolic = dss_handle.get_symbolic();
        auto*& numeric = dss_handle.get_numeric();
        auto* control = dss_handle.get_control();

        umfpack_di_numeric(
            row_ptrs.data(), col_inds.data(), values.data(), symbolic, &numeric, control, nullptr
        );
    }
};

}  // namespace kynema_fmb::dss
