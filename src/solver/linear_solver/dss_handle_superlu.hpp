#pragma once

#include <memory>
#include <vector>

#include "dss_algorithm.hpp"
#include "slu_ddefs.h"

namespace kynema_fmb::dss {
template <>
class Handle<Algorithm::SUPERLU> {
    struct superluDssHandleType {
        superlu_options_t options{};
        SuperLUStat_t stat{};
        SuperMatrix L{};
        SuperMatrix U{};
        GlobalLU_t Glu{};
        std::vector<int> perm_r;
        std::vector<int> perm_c;
        std::vector<int> etree;

        superluDssHandleType() {
            set_default_options(&options);
            StatInit(&stat);
        }

        superluDssHandleType(superluDssHandleType&) = delete;
        void operator=(superluDssHandleType&) = delete;
        superluDssHandleType(superluDssHandleType&&) = delete;
        void operator=(superluDssHandleType&&) = delete;

        ~superluDssHandleType() { StatFree(&stat); }
    };
    std::shared_ptr<superluDssHandleType> superlu_dss_handle;

public:
    Handle() : superlu_dss_handle(std::make_shared<superluDssHandleType>()) {}

    superlu_options_t& get_options() { return superlu_dss_handle->options; }

    SuperLUStat_t& get_stat() { return superlu_dss_handle->stat; }

    SuperMatrix& get_L() { return superlu_dss_handle->L; }

    SuperMatrix& get_U() { return superlu_dss_handle->U; }

    GlobalLU_t& get_Glu() { return superlu_dss_handle->Glu; }

    std::vector<int>& get_perm_r() { return superlu_dss_handle->perm_r; }

    std::vector<int>& get_perm_c() { return superlu_dss_handle->perm_c; }
    std::vector<int>& get_etree() { return superlu_dss_handle->etree; }
};

}  // namespace kynema_fmb::dss
