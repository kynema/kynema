target_sources(kynema_fmb_library PRIVATE)

install(FILES
        dss_algorithm.hpp
        dss_handle.hpp
        dss_handle_cusolversp.hpp
        dss_handle_cudss.hpp
        dss_handle_klu.hpp
        dss_handle_umfpack.hpp
        dss_handle_mkl.hpp
        dss_handle_superlu.hpp
        dss_handle_superlu_mt.hpp
        dss_numeric.hpp
        dss_numeric_cusolversp.hpp
        dss_numeric_cudss.hpp
        dss_numeric_klu.hpp
        dss_numeric_umfpack.hpp
        dss_numeric_mkl.hpp
        dss_numeric_superlu.hpp
        dss_numeric_superlu_mt.hpp
        dss_solve.hpp
        dss_solve_cusolversp.hpp
        dss_solve_cudss.hpp
        dss_solve_klu.hpp
        dss_solve_umfpack.hpp
        dss_solve_mkl.hpp
        dss_solve_superlu.hpp
        dss_solve_superlu_mt.hpp
        dss_symbolic.hpp
        dss_symbolic_cusolversp.hpp
        dss_symbolic_cudss.hpp
        dss_symbolic_klu.hpp
        dss_symbolic_umfpack.hpp
        dss_symbolic_mkl.hpp
        dss_symbolic_superlu.hpp
        dss_symbolic_superlu_mt.hpp
        DESTINATION include/KynemaFMB/solver/linear_solver
)
