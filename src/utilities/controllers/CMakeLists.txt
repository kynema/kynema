target_sources(kynema_fmb_library
  PRIVATE
  turbine_controller.cpp
)

target_sources(DISCON
  PRIVATE
  discon.cpp
)

target_sources(DISCON_ROTOR_TEST_CONTROLLER
  PRIVATE
  discon_rotor_test_controller.cpp
)

install(FILES
        controller_io.hpp
        discon.hpp
        discon_rotor_test_controller.hpp
        turbine_controller.hpp
        DESTINATION include/KynemaFMB/utilities/controllers/
)
