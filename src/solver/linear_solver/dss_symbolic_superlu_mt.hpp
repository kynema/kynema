#pragma once

#include <Kokkos_Core.hpp>

#include "dss_algorithm.hpp"
#include "dss_handle_superlu_mt.hpp"
#include "slu_mt_ddefs.h"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct SymbolicFunction<Handle<Algorithm::SUPERLU_MT>, CrsMatrixType> {
    static void symbolic(Handle<Algorithm::SUPERLU_MT>& dss_handle, CrsMatrixType& A) {
        auto num_rows = A.numRows();
        auto num_cols = A.numCols();
        auto num_non_zero = A.nnz();

        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);
        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto& options = dss_handle.get_options();
        options.refact = NO;
        auto& stat = dss_handle.get_stat();
        auto& L = dss_handle.get_L();
        auto& U = dss_handle.get_U();
        auto& perm_r = dss_handle.get_perm_r();
        auto& perm_c = dss_handle.get_perm_c();
        auto& etree = dss_handle.get_etree();
        auto& colcnt_h = dss_handle.get_colcnt_h();
        auto& part_super_h = dss_handle.get_part_super_h();
        auto& work = dss_handle.get_work();

        StatAlloc(num_cols, options.nprocs, options.panel_size, options.relax, &stat);
        StatInit(num_cols, options.nprocs, &stat);

        perm_r.resize(static_cast<size_t>(num_rows));
        perm_c.resize(static_cast<size_t>(num_cols));
        etree.resize(static_cast<size_t>(num_rows));
        colcnt_h.resize(static_cast<size_t>(num_cols));
        part_super_h.resize(static_cast<size_t>(num_cols));
        options.perm_c = perm_c.data();
        options.perm_r = perm_r.data();
        options.etree = etree.data();
        options.colcnt_h = colcnt_h.data();
        options.part_super_h = part_super_h.data();

        SuperMatrix Amatrix;
        dCreate_CompCol_Matrix(
            &Amatrix, num_rows, num_cols, num_non_zero, values.data(), col_inds.data(),
            const_cast<int*>(row_ptrs.data()), SLU_NC, SLU_D, SLU_GE
        );
        get_perm_c(0, &Amatrix, perm_c.data());

        auto info = 0;

        SuperMatrix AC;
        sp_colorder(&Amatrix, options.perm_c, &options, &AC);

        options.lwork = -1;
        pdgstrf(&options, &AC, options.perm_r, &L, &U, &stat, &info);
        options.lwork = info;
        if (static_cast<size_t>(info) > work.size()) {
            work.resize(static_cast<size_t>(1.2 * info));
            options.work = work.data();
        }
    }
};

}  // namespace kynema_fmb::dss
