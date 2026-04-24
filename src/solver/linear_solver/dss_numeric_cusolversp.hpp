#pragma once

#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

#include "dss_algorithm.hpp"
#include "dss_handle_cusolversp.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct NumericFunction<Handle<Algorithm::CUSOLVER_SP>, CrsMatrixType> {
    static void numeric(Handle<Algorithm::CUSOLVER_SP>& dss_handle, CrsMatrixType& A) {
        const auto num_rows = A.numRows();
        const auto num_cols = A.numCols();
        const auto num_non_zero = A.nnz();

        auto* values = A.values.data();
        auto* row_ptrs = A.graph.row_map.data();
        auto* col_inds = A.graph.entries.data();

        auto& handle = dss_handle.get_handle();
        auto& description = dss_handle.get_description();
        auto& info = dss_handle.get_info();
        auto& buffer = dss_handle.get_buffer();

        auto size_internal = 0UL;
        auto size_chol = 0UL;

        const auto tol = 1.e-14;
        auto singularity = 0;

        cusolverSpDcsrqrBufferInfo(
            handle, num_rows, num_cols, num_non_zero, description, values, row_ptrs, col_inds, info,
            &size_internal, &size_chol
        );
        if (size_chol > buffer.extent(0)) {
            Kokkos::realloc(buffer, size_chol);
        }
        cusolverSpDcsrqrSetup(
            handle, num_rows, num_cols, num_non_zero, description, values, row_ptrs, col_inds, 0.,
            info
        );
        cusolverSpDcsrqrFactor(
            handle, num_rows, num_cols, num_non_zero, NULL, NULL, info, buffer.data()
        );
        cusolverSpDcsrqrZeroPivot(handle, info, tol, &singularity);
    }
};

}  // namespace kynema_fmb::dss
