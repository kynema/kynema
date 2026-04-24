#pragma once

#include <Kokkos_Core.hpp>

#include "dss_algorithm.hpp"
#include "dss_handle_superlu_mt.hpp"
#include "slu_mt_ddefs.h"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct NumericFunction<DSSHandle<Algorithm::SUPERLU_MT>, CrsMatrixType> {
    static void numeric(Handle<Algorithm::SUPERLU_MT>& dss_handle, CrsMatrixType& A) {
        auto num_rows = A.numRows();
        auto num_cols = A.numCols();
        auto num_non_zero = A.nnz();

        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto& options = dss_handle.get_options();
        auto& stat = dss_handle.get_stat();
        auto& L = dss_handle.get_L();
        auto& U = dss_handle.get_U();

        SuperMatrix Amatrix;
        dCreate_CompCol_Matrix(
            &Amatrix, num_rows, num_cols, num_non_zero, values.data(), col_inds.data(),
            const_cast<int*>(row_ptrs.data()), SLU_NC, SLU_D, SLU_GE
        );
        auto info = 0;

        SuperMatrix AC;
        sp_colorder(&Amatrix, options.perm_c, &options, &AC);

        pdgstrf(&options, &AC, options.perm_r, &L, &U, &stat, &info);

        Destroy_CompCol_Permuted(&AC);

        options.refact = YES;
    }
};

}  // namespace kynema_fmb::dss
