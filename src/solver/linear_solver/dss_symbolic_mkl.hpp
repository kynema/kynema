#pragma once

#include <Kokkos_Core.hpp>
#include <mkl_dss.h>

#include "dss_algorithm.hpp"
#include "dss_handle_mkl.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct SymbolicFunction<Handle<Algorithm::MKL>, CrsMatrixType> {
    static void symbolic(Handle<Algorithm::MKL>& dss_handle, CrsMatrixType& A) {
        auto& handle = dss_handle.get_handle();
        constexpr MKL_INT opt = MKL_DSS_NON_SYMMETRIC;

        const MKL_INT num_rows = A.numRows();
        const MKL_INT num_cols = A.numCols();
        const MKL_INT num_non_zeros = A.nnz();

        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        dss_define_structure(
            handle, opt, row_ptrs.data(), num_rows, num_cols, col_inds.data(), num_non_zeros
        );

        auto& perm = dss_handle.get_perm();
        perm.resize(static_cast<size_t>(num_rows));
        constexpr MKL_INT reorder_opt = MKL_DSS_AUTO_ORDER;
        dss_reorder(handle, reorder_opt, perm.data());
    }
};
}  // namespace kynema_fmb::dss
