#pragma once

#include <Kokkos_Core.hpp>
#include <mkl_dss.h>

#include "dss_algorithm.hpp"
#include "dss_handle_mkl.hpp"

namespace kynema_fmb::dss {

template <typename CrsMatrixType>
struct NumericFunction<Handle<Algorithm::MKL>, CrsMatrixType> {
    static void numeric(Handle<Algorithm::MKL>& dss_handle, CrsMatrixType& A) {
        auto& handle = dss_handle.get_handle();
        constexpr MKL_INT opt = MKL_DSS_INDEFINITE;

        auto values = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), A.values);

        dss_factor_real(handle, opt, values.data());
    }
};

}  // namespace kynema_fmb::dss
