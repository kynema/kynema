#pragma once

#include <Kokkos_Core.hpp>
#include <klu.h>

#include "dss_algorithm.hpp"
#include "dss_handle_klu.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType>
struct SymbolicFunction<Handle<Algorithm::KLU>, CrsMatrixType> {
    static void symbolic(Handle<Algorithm::KLU>& dss_handle, CrsMatrixType& A) {
        const auto num_rows = A.numRows();

        auto row_ptrs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.row_map);
        auto col_inds = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.graph.entries);

        auto*& symbolic = dss_handle.get_symbolic();
        auto*& numeric = dss_handle.get_numeric();
        auto& common = dss_handle.get_common();

        if (symbolic != nullptr) {
            klu_free_symbolic(&symbolic, &common);
        }
        if (numeric != nullptr) {
            klu_free_numeric(&numeric, &common);
        }
        symbolic =
            klu_analyze(num_rows, const_cast<int*>(row_ptrs.data()), col_inds.data(), &common);
    }
};

}  // namespace kynema_fmb::dss
