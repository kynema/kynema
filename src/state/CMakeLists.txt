target_sources(kynema_fmb_library PRIVATE)

install(FILES
        calculate_displacement.hpp
        calculate_next_state.hpp
        clone_state.hpp
        copy_state_data.hpp
        state.hpp
        update_algorithmic_acceleration.hpp
        update_dynamic_prediction.hpp
        update_global_position.hpp
        update_static_prediction.hpp
        DESTINATION include/Kynema/state/
)
