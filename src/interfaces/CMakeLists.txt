add_subdirectory(cfd)
add_subdirectory(blade)
add_subdirectory(components)
add_subdirectory(turbine)

target_sources(kynema_library
  PRIVATE
  node_data.cpp
  constraint_data.cpp
  outputs.cpp
)

install(FILES
	constraint_data.hpp
        constraint_data.hpp
        host_constraints.hpp
        host_state.hpp
        node_data.hpp
        outputs.hpp
        DESTINATION include/Kynema/interfaces/
)
