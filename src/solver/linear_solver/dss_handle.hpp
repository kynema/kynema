#pragma once

#include "dss_algorithm.hpp"

#include "Kynema_FMB_config.h"

namespace kynema_fmb::dss {

template <Algorithm>
class Handle {
public:
    Handle() = delete;
};
}  // namespace kynema_fmb::dss

#ifdef KYNEMA_FMB_ENABLE_CUSOLVERSP
#include "dss_handle_cusolversp.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_CUDSS
#include "dss_handle_cudss.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_MKL
#include "dss_handle_mkl.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_KLU
#include "dss_handle_klu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_UMFPACK
#include "dss_handle_umfpack.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU
#include "dss_handle_superlu.hpp"
#endif

#ifdef KYNEMA_FMB_ENABLE_SUPERLU_MT
#include "dss_handle_superlu_mt.hpp"
#endif
