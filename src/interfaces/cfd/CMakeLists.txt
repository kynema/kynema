target_sources(kynema_fmb_library
  PRIVATE
  interface.cpp
  interface_builder.cpp
)

install(FILES
        floating_platform.hpp
        floating_platform_input.hpp
        interface.hpp
        interface_builder.hpp
        interface_input.hpp
        mooring_line.hpp
        mooring_line_input.hpp
        node_data.hpp
        turbine.hpp
        turbine_input.hpp
        DESTINATION include/KynemaFMB/interfaces/cfd/
)
