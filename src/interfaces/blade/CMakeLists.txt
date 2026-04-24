target_sources(kynema_fmb_library
  PRIVATE
  blade_interface.cpp
  blade_interface_builder.cpp
)

install(FILES
        blade_interface.hpp
        blade_interface_builder.hpp
        DESTINATION include/KynemaFMB/interfaces/blade/
)
