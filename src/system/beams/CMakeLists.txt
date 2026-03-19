target_sources(kynema_library PRIVATE)

install(FILES
        calculate_D_D1.hpp
        calculate_D_D2.hpp
        calculate_force_FC.hpp
        calculate_force_FD.hpp
        calculate_force_FD1.hpp
        calculate_force_FD2.hpp
        calculate_G_D1.hpp
        calculate_G_D2.hpp
        calculate_K_D1.hpp
        calculate_K_D2.hpp
        calculate_Ouu.hpp
        calculate_P_D2.hpp
        calculate_Puu.hpp
        calculate_quadrature_point_damping_values.hpp
        calculate_quadrature_point_inertial_values.hpp
        calculate_quadrature_point_stiffness_values.hpp
        calculate_quadrature_point_values.hpp
        calculate_Quu.hpp
        calculate_strain_dot.hpp
        calculate_strain.hpp
        calculate_system_matrix.hpp
        calculate_temporary_variables.hpp
        integrate_inertia_matrix.hpp
        integrate_residual_vector.hpp
        integrate_stiffness_matrix.hpp
        interpolate_to_quadrature_point_for_damping.hpp
        interpolate_to_quadrature_point_for_inertia.hpp
        interpolate_to_quadrature_point_for_stiffness.hpp
        update_node_state.hpp
        DESTINATION include/Kynema/system/beams/
)
