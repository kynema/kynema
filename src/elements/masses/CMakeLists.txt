target_sources(kynema_fmb_library PRIVATE)

install(FILES
        create_masses.hpp
        mass_element.hpp
        masses.hpp
        masses_input.hpp
        DESTINATION include/KynemaFMB/elements/masses/
)
