target_sources(kynema_fmb_library PRIVATE)

install(FILES
        beam_element.hpp
        beam_quadrature.hpp
        beam_section.hpp
        beams.hpp
        beams_input.hpp
        calculate_jacobian.hpp
        calculate_QP_position.hpp
        create_beams.hpp
        generate_sectional_properties.hpp
        hollow_circle_properties.hpp
        interpolate_QP_position.hpp
        interpolate_QP_rotation.hpp
        interpolate_QP_state.hpp
        interpolate_QP_vector.hpp
        interpolate_to_quadrature_points.hpp
        interpolation_operations.hpp
        populate_element_views.hpp
        DESTINATION include/KynemaFMB/elements/beams/
)
