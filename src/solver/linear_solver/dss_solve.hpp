#pragma once

#include "KYNEMA_FMB_config.h"

namespace kynema_fmb::dss {

template <typename DSSHandleType, typename CrsMatrixType, typename MultiVectorType>
struct SolveFunction {
    SolveFunction() = delete;
};

}  // namespace kynema_fmb::dss

#ifdef KYNEMA_FMB_ENABLE_CUSOLVERSP
#include "dss_solve_cusolversp.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_CUDSS
#include "dss_solve_cudss.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_MKL
#include "dss_solve_mkl.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_KLU
#include "dss_solve_klu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_UMFPACK
#include "dss_solve_umfpack.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU
#include "dss_solve_superlu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU_MT
#include "dss_solve_superlu_mt.hpp"
#endif

namespace kynema_fmb::dss {

template <typename DSSHandleType, typename CrsMatrixType, typename MultiVectorType>
void solve(DSSHandleType& dss_handle, CrsMatrixType& A, MultiVectorType& b, MultiVectorType& x) {
    SolveFunction<DSSHandleType, CrsMatrixType, MultiVectorType>::solve(dss_handle, A, b, x);
}

}  // namespace kynema_fmb::dss
