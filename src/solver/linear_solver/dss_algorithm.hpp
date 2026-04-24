#pragma once
#include <cstdint>

namespace kynema_fmb::dss {

enum class Algorithm : std::uint8_t {
    CUSOLVER_SP,
    CUDSS,
    KLU,
    UMFPACK,
    SUPERLU,
    SUPERLU_MT,
    MKL,
    NONE,
};

}  // namespace kynema_fmb::dss
