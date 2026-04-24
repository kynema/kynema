#pragma once

#include <Kokkos_Core.hpp>

#include "dss_algorithm.hpp"
#include "dss_handle_superlu.hpp"
#include "slu_ddefs.h"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct NumericFunction<Handle<Algorithm::SUPERLU>, CrsMatrixType> {
    static void numeric(Handle<Algorithm::SUPERLU>& dss_handle, CrsMatrixType& A) {
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
        auto& Glu = dss_handle.get_Glu();
        auto& perm_r = dss_handle.get_perm_r();
        auto& perm_c = dss_handle.get_perm_c();
        auto& etree = dss_handle.get_etree();

        SuperMatrix Amatrix;
        dCreate_CompCol_Matrix(
            &Amatrix, num_rows, num_cols, num_non_zero, values.data(), col_inds.data(),
            const_cast<int*>(row_ptrs.data()), SLU_NC, SLU_D, SLU_GE
        );
        constexpr auto relax = 1;
        constexpr auto panel_size = 1;
        auto info = 0;

        SuperMatrix AC;
        sp_preorder(&options, &Amatrix, perm_c.data(), etree.data(), &AC);
        dgstrf(
            &options, &AC, relax, panel_size, etree.data(), nullptr, 0, perm_c.data(), perm_r.data(),
            &L, &U, &Glu, &stat, &info
        );

        Destroy_CompCol_Permuted(&AC);

        options.Fact = SamePattern;
    }
};

}  // namespace kynema_fmb::dss
