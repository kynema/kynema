target_sources(kynema_fmb_library PRIVATE)

install(FILES
        aerodyn_inflow.hpp
        DESTINATION include/KynemaFMB/utilities/aerodynamics/
)
