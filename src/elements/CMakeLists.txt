target_sources(kynema_fmb_library PRIVATE)

install(FILES
        elements.hpp
        DESTINATION include/KynemaFMB/elements
)

add_subdirectory(beams)
add_subdirectory(masses)
add_subdirectory(springs)
