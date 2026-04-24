target_sources(kynema_fmb_library PRIVATE)

add_subdirectory(linear_solver)

install(FILES
        calculate_error_sum_squares.hpp
        compute_col_inds.hpp
        compute_constraints_col_inds.hpp
        compute_constraints_row_entries.hpp
        compute_num_system_dofs.hpp
        compute_row_ptrs.hpp
        compute_system_col_inds.hpp
        compute_system_row_entries.hpp
        condition_system.hpp
        contribute_beams_to_sparse_matrix.hpp
        contribute_beams_to_vector.hpp
        contribute_constraints_system_residual_to_vector.hpp
        contribute_forces_to_vector.hpp
        contribute_lambda_to_vector.hpp
        contribute_masses_to_sparse_matrix.hpp
        contribute_masses_to_vector.hpp
        contribute_springs_to_sparse_matrix.hpp
        contribute_springs_to_vector.hpp
        copy_constraints_residual_to_vector.hpp
        copy_constraints_to_sparse_matrix.hpp
        copy_constraints_transpose_to_sparse_matrix.hpp
        create_full_matrix.hpp
        scan_row_entries.hpp
        solver.hpp
        DESTINATION include/KynemaFMB/solver/
)
