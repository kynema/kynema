#pragma once

#include <memory>

#include "dss_algorithm.hpp"
#include "klu.h"

namespace kynema_fmb::dss {
template <>
class Handle<Algorithm::KLU> {
    struct kluDssHandleType {
        klu_symbolic* Symbolic = nullptr;
        klu_numeric* Numeric = nullptr;
        klu_common Common{};

        kluDssHandleType() {
            klu_defaults(&Common);
            Common.ordering = 1;
        }

        kluDssHandleType(kluDssHandleType&) = delete;
        void operator=(kluDssHandleType&) = delete;
        kluDssHandleType(kluDssHandleType&&) = delete;
        void operator=(kluDssHandleType&&) = delete;

        ~kluDssHandleType() {
            if (Symbolic != nullptr) {
                klu_free_symbolic(&Symbolic, &Common);
            }
            if (Numeric != nullptr) {
                klu_free_numeric(&Numeric, &Common);
            }
        }
    };
    std::shared_ptr<kluDssHandleType> klu_dss_handle;

public:
    Handle() : klu_dss_handle(std::make_shared<kluDssHandleType>()) {}

    klu_symbolic*& get_symbolic() { return klu_dss_handle->Symbolic; }

    klu_numeric*& get_numeric() { return klu_dss_handle->Numeric; }

    klu_common& get_common() { return klu_dss_handle->Common; }
};

}  // namespace kynema_fmb::dss
