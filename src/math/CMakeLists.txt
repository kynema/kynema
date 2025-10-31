target_sources(kynema_library PRIVATE
    least_squares_fit.cpp
)

install(FILES
        interpolation.hpp
        least_squares_fit.hpp
        matrix_operations.hpp
        project_points_to_target_polynomial.hpp
        quaternion_operations.hpp
        vector_operations.hpp
        DESTINATION include/Kynema/math
)
