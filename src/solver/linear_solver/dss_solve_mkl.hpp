#pragma once

#include <Kokkos_Core.hpp>
#include <mkl_dss.h>

#include "dss_algorithm.hpp"
#include "dss_handle_mkl.hpp"

namespace kynema_fmb::dss {
template <typename CrsMatrixType, typename MultiVectorType>
struct SolveFunction<Handle<Algorithm::MKL>, CrsMatrixType, MultiVectorType> {
    static void solve(
        Handle<Algorithm::MKL>& dss_handle, CrsMatrixType&, MultiVectorType& b, MultiVectorType& x
    ) {
        auto& handle = dss_handle.get_handle();
        constexpr MKL_INT opt = 0;

        const auto nrhs = static_cast<MKL_INT>(b.extent(1));

        auto x_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), x);
        auto b_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), b);

        dss_solve_real(handle, opt, b_host.data(), nrhs, x_host.data());

        Kokkos::deep_copy(x, x_host);
    }
};
}  // namespace kynema_fmb::dss
