target_sources(kynema_fmb_library PRIVATE)

install(FILES
        create_springs.hpp
        spring_element.hpp
        springs.hpp
        springs_input.hpp
        DESTINATION include/KynemaFMB/elements/springs/
)
