#pragma once

#include <cudss.h>

#include "dss_algorithm.hpp"

namespace kynema_fmb::dss {
template <>
class Handle<Algorithm::CUDSS> {
    struct cudssDssHandleType {
        cudssHandle_t handle;
        cudssConfig_t solverConfig;
        cudssData_t solverData;
        bool is_first_factorization;

        cudssDssHandleType() : is_first_factorization{true} {
            cudssCreate(&handle);
            cudssConfigCreate(&solverConfig);
            auto flag = CUDSS_ALG_1;
            cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG, &flag, sizeof(flag));
            cudssDataCreate(handle, &solverData);
        }

        ~cudssDssHandleType() {
            cudssDataDestroy(handle, solverData);
            cudssConfigDestroy(solverConfig);
            cudssDestroy(handle);
        }
    };
    std::shared_ptr<cudssDssHandleType> cudss_dss_handle;

public:
    Handle() : cudss_dss_handle(std::make_shared<cudssDssHandleType>()) {}

    cudssHandle_t& get_handle() { return cudss_dss_handle->handle; }

    cudssConfig_t& get_config() { return cudss_dss_handle->solverConfig; }

    cudssData_t& get_data() { return cudss_dss_handle->solverData; }

    void set_initial_factorization(bool value) { cudss_dss_handle->is_first_factorization = value; }

    bool is_initial_factorization() const { return cudss_dss_handle->is_first_factorization; }
};

}  // namespace kynema_fmb::dss
