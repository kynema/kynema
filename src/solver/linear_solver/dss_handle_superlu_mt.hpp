#pragma once

#include <memory>
#include <vector>

#include "dss_algorithm.hpp"
#include "slu_mt_ddefs.h"

namespace kynema_fmb::dss {
template <>
class Handle<Algorithm::SUPERLU_MT> {
    struct superluDssHandleType {
        superlumt_options_t options{};
        Gstat_t stat{};
        SuperMatrix L{};
        SuperMatrix U{};
        std::vector<int> perm_r;
        std::vector<int> perm_c;
        std::vector<int> etree;
        std::vector<int> colcnt_h;
        std::vector<int> part_super_h;
        std::vector<char> work;

        superluDssHandleType() {
            options.nprocs = 1;
            options.fact = DOFACT;
            options.trans = TRANS;
            options.refact = NO;
            options.panel_size = 1;
            options.relax = 1;
            options.diag_pivot_thresh = 1.;
            options.usepr = NO;
            options.PrintStat = NO;
            options.drop_tol = 0.;
            options.lwork = 0;
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

    superlumt_options_t& get_options() { return superlu_dss_handle->options; }

    Gstat_t& get_stat() { return superlu_dss_handle->stat; }

    SuperMatrix& get_L() { return superlu_dss_handle->L; }

    SuperMatrix& get_U() { return superlu_dss_handle->U; }

    std::vector<int>& get_perm_r() { return superlu_dss_handle->perm_r; }
    std::vector<int>& get_perm_c() { return superlu_dss_handle->perm_c; }
    std::vector<int>& get_etree() { return superlu_dss_handle->etree; }
    std::vector<int>& get_colcnt_h() { return superlu_dss_handle->colcnt_h; }
    std::vector<int>& get_part_super_h() { return superlu_dss_handle->part_super_h; }
    std::vector<char>& get_work() { return superlu_dss_handle->work; }
};

}  // namespace kynema_fmb::dss
