target_sources(kynema_fmb_library PRIVATE)

install(FILES
        calculate_tangent_operator.hpp
        DESTINATION include/KynemaFMB/system/
)

add_subdirectory(beams)
add_subdirectory(masses)
add_subdirectory(springs)
