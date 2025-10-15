target_sources(kynema_library
  PRIVATE
  aerodynamics.cpp
  beam.cpp
  beam_builder.cpp
  inflow.cpp
  solution_builder.cpp
  turbine.cpp
  turbine_builder.cpp
)

install(FILES
        aerodynamics.hpp
        aerodynamics_builder.hpp
        aerodynamics_input.hpp
        beam.hpp
        beam_builder.hpp
        beam_input.hpp
        controller_builder.hpp
        controller_input.hpp
        inflow.hpp
        solution_builder.hpp
        solution_input.hpp
        turbine.hpp
        turbine_builder.hpp
        turbine_input.hpp
        DESTINATION include/Kynema/interfaces/components/
)
