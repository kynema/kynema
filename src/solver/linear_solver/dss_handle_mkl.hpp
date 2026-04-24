#pragma once

#include <memory>
#include <vector>

#include <mkl_dss.h>

#include "dss_algorithm.hpp"

namespace kynema_fmb::dss {

template <>
class Handle<Algorithm::MKL> {
    struct mklDssHandleType {
        _MKL_DSS_HANDLE_t handle{};
        std::vector<MKL_INT> perm;

        mklDssHandleType() {
            constexpr MKL_INT opt = MKL_DSS_ZERO_BASED_INDEXING;
            dss_create(handle, opt);
        }

        mklDssHandleType(mklDssHandleType&) = delete;
        void operator=(mklDssHandleType&) = delete;
        mklDssHandleType(mklDssHandleType&&) = delete;
        void operator=(mklDssHandleType&&) = delete;

        ~mklDssHandleType() {
            constexpr MKL_INT opt = 0;
            dss_delete(handle, opt);
        }
    };
    std::shared_ptr<mklDssHandleType> mkl_dss_handle;

public:
    Handle() : mkl_dss_handle(std::make_shared<mklDssHandleType>()) {}

    _MKL_DSS_HANDLE_t& get_handle() { return mkl_dss_handle->handle; }

    std::vector<MKL_INT>& get_perm() { return mkl_dss_handle->perm; }
};
}  // namespace kynema_fmb::dss
