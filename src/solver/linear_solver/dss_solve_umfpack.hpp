#pragma once

#include <Kokkos_Core.hpp>
#include <umfpack.h>

#include "dss_algorithm.hpp"
#include "dss_handle_umfpack.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType, typename MultiVectorType>
struct SolveFunction<Handle<Algorithm::UMFPACK>, CrsMatrixType, MultiVectorType> {
    static void solve(
        Handle<Algorithm::UMFPACK>& dss_handle, CrsMatrixType& A, MultiVectorType& b,
        MultiVectorType& x
    ) {
        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto x_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), x);
        auto b_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);

        auto*& numeric = dss_handle.get_numeric();
        auto* control = dss_handle.get_control();

        umfpack_di_solve(
            UMFPACK_At, row_ptrs.data(), col_inds.data(), values.data(), x_host.data(),
            b_host.data(), numeric, control, nullptr
        );

        Kokkos::deep_copy(x, x_host);
    }
};

}  // namespace kynema_fmb::dss
