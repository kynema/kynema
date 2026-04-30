# Specify the source files for the unit test executable
target_sources(
    kynema_fmb_unit_tests
    PRIVATE
    test_assemble_node_freedom_allocation_table.cpp
    test_compute_node_freedom_map_table.cpp
    test_create_element_freedom_table.cpp
    test_create_constraint_freedom_table.cpp
    test_freedom_signature.cpp
)
