target_sources(kynema_fmb_library PRIVATE)

install(FILES
        calculate_distance_components.hpp
        calculate_force_coefficients.hpp
        calculate_force_vectors.hpp
        calculate_length.hpp
        calculate_quadrature_point_values.hpp
        calculate_stiffness_matrix.hpp
        DESTINATION include/KynemaFMB/system/springs/
)
