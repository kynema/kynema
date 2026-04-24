target_sources(kynema_fmb_library PRIVATE)

install(FILES
        copy_nodes_to_state.hpp
        mesh_connectivity.hpp
        model.hpp
        node.hpp
        DESTINATION include/KynemaFMB/model/
)
