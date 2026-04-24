#pragma once

#include "Kynema_FMB_config.h"

namespace kynema_fmb::dss {

template <typename DSSHandleType, typename CrsMatrixType>
struct SymbolicFunction {
    SymbolicFunction() = delete;
};

}  // namespace kynema_fmb::dss

#ifdef KYNEMA_FMB_ENABLE_CUSOLVERSP
#include "dss_symbolic_cusolversp.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_CUDSS
#include "dss_symbolic_cudss.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_MKL
#include "dss_symbolic_mkl.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_KLU
#include "dss_symbolic_klu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_UMFPACK
#include "dss_symbolic_umfpack.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU
#include "dss_symbolic_superlu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU_MT
#include "dss_symbolic_superlu_mt.hpp"
#endif

namespace kynema_fmb::dss {

template <typename DSSHandleType, typename CrsMatrixType>
void symbolic_factorization(DSSHandleType& dss_handle, CrsMatrixType& A) {
    SymbolicFunction<DSSHandleType, CrsMatrixType>::symbolic(dss_handle, A);
}

}  // namespace kynema_fmb::dss
