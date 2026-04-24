#pragma once

#include <Kokkos_Core.hpp>
#include <cudss.h>

#include "dss_algorithm.hpp"
#include "dss_handle_cudss.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType, typename MultiVectorType>
struct SolveFunction<Handle<Algorithm::CUDSS>, CrsMatrixType, MultiVectorType> {
    static void solve(
        Handle<Algorithm::CUDSS>& dss_handle, CrsMatrixType& A, MultiVectorType& b,
        MultiVectorType& x
    ) {
        const auto num_rows = A.numRows();
        const auto num_cols = A.numCols();
        const auto num_non_zero = static_cast<int>(A.nnz());

        auto* values = A.values.data();
        auto* row_ptrs = A.graph.row_map.data();
        auto* col_inds = A.graph.entries.data();

        auto& handle = dss_handle.get_handle();
        auto& config = dss_handle.get_config();
        auto& data = dss_handle.get_data();

        cudssMatrix_t A_cudss;
        cudssMatrix_t x_cudss;
        cudssMatrix_t b_cudss;

        cudssMatrixCreateCsr(
            &A_cudss, num_rows, num_cols, num_non_zero, const_cast<int*>(row_ptrs), nullptr,
            col_inds, values, CUDA_R_32I, CUDA_R_64F, CUDSS_MTYPE_GENERAL, CUDSS_MVIEW_FULL,
            CUDSS_BASE_ZERO
        );
        cudssMatrixCreateDn(
            &b_cudss, num_cols, 1, num_cols, b.data(), CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR
        );
        cudssMatrixCreateDn(
            &x_cudss, num_rows, 1, num_rows, x.data(), CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR
        );

        cudssExecute(handle, CUDSS_PHASE_SOLVE, config, data, A_cudss, x_cudss, b_cudss);

        cudssMatrixDestroy(A_cudss);
        cudssMatrixDestroy(b_cudss);
        cudssMatrixDestroy(x_cudss);
    }
};

}  // namespace kynema_fmb::dss
