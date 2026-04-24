#pragma once

#include <Kokkos_Core.hpp>

#include "dss_algorithm.hpp"
#include "dss_handle_superlu.hpp"
#include "slu_ddefs.h"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct SymbolicFunction<Handle<Algorithm::SUPERLU>, CrsMatrixType> {
    static void symbolic(Handle<Algorithm::SUPERLU>& dss_handle, CrsMatrixType& A) {
        auto num_rows = A.numRows();
        auto num_cols = A.numCols();
        auto num_non_zero = A.nnz();

        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto& options = dss_handle.get_options();
        options.Fact = DOFACT;
        auto& perm_r = dss_handle.get_perm_r();
        auto& perm_c = dss_handle.get_perm_c();
        auto& etree = dss_handle.get_etree();

        perm_r.resize(static_cast<size_t>(num_rows));
        perm_c.resize(static_cast<size_t>(num_cols));
        etree.resize(static_cast<size_t>(num_rows));

        SuperMatrix Amatrix;
        dCreate_CompCol_Matrix(
            &Amatrix, num_rows, num_cols, num_non_zero, values.data(), col_inds.data(),
            const_cast<int*>(row_ptrs.data()), SLU_NC, SLU_D, SLU_GE
        );
        get_perm_c(options.ColPerm, &Amatrix, perm_c.data());
    }
};

}  // namespace kynema_fmb::dss
