target_sources(kynema_fmb_library PRIVATE)

install(FILES
        calculate_gravity_force.hpp
        calculate_gyroscopic_matrix.hpp
        calculate_inertial_force.hpp
        calculate_inertia_stiffness_matrix.hpp
        calculate_mass_matrix_components.hpp
        calculate_quadrature_point_values.hpp
        rotate_section_matrix.hpp
        DESTINATION include/KynemaFMB/system/masses
)
