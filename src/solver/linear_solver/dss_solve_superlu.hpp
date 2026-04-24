#pragma once

#include <Kokkos_Core.hpp>

#include "dss_algorithm.hpp"
#include "dss_handle_superlu.hpp"
#include "slu_ddefs.h"

namespace kynema_fmb::dss {
template <typename CrsMatrixType, typename MultiVectorType>
struct SolveFunction<Handle<Algorithm::SUPERLU>, CrsMatrixType, MultiVectorType> {
    static void solve(
        Handle<Algorithm::SUPERLU>& dss_handle, CrsMatrixType&, MultiVectorType& b,
        MultiVectorType& x
    ) {
        auto& stat = dss_handle.get_stat();
        auto& L = dss_handle.get_L();
        auto& U = dss_handle.get_U();
        auto& perm_r = dss_handle.get_perm_r();
        auto& perm_c = dss_handle.get_perm_c();

        Kokkos::deep_copy(x, b);
        auto x_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), x);

        SuperMatrix Xmatrix;
        dCreate_Dense_Matrix(
            &Xmatrix, static_cast<int>(x.extent(0)), static_cast<int>(x.extent(1)), x_host.data(),
            static_cast<int>(x.extent(0)), SLU_DN, SLU_D, SLU_GE
        );

        auto info = 0;
        dgstrs(TRANS, &L, &U, perm_c.data(), perm_r.data(), &Xmatrix, &stat, &info);

        Kokkos::deep_copy(x, x_host);
    }
};

}  // namespace kynema_fmb::dss
