
target_sources(DISCON
  PRIVATE
  discon.cpp
)

target_sources(DISCON_ROTOR_TEST_CONTROLLER
  PRIVATE
  discon_rotor_test_controller.cpp
)

install(FILES
        discon.hpp
        discon_rotor_test_controller.hpp
        DESTINATION include/Kynema/utilities/controllers/
)
