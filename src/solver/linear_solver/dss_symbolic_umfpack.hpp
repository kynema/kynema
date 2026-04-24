#pragma once

#include <Kokkos_Core.hpp>
#include <umfpack.h>

#include "dss_algorithm.hpp"
#include "dss_handle_umfpack.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct SymbolicFunction<Handle<Algorithm::UMFPACK>, CrsMatrixType> {
    static void symbolic(Handle<Algorithm::UMFPACK>& dss_handle, CrsMatrixType& A) {
        const auto num_rows = A.numRows();
        const auto num_cols = A.numCols();

        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto*& symbolic = dss_handle.get_symbolic();
        auto* control = dss_handle.get_control();

        umfpack_di_symbolic(
            num_rows, num_cols, row_ptrs.data(), col_inds.data(), nullptr, &symbolic, control,
            nullptr
        );
    }
};
}  // namespace kynema_fmb::dss
