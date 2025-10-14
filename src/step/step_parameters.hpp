#pragma once

#include <cstddef>

namespace kynema {

/**
 * @brief A Struct containing the paramters used to control the time stepping process.
 *
 * @details These paramteres
 */
struct StepParameters {
    bool is_dynamic_solve;            //< analysis type (dynamic vs static)
    size_t max_iter;                  //< maximum number of nonlinear iterations
    double h;                         //< time step size
    double alpha_m;                   //< generalized alpha solver parameter
    double alpha_f;                   //< generalized alpha solver parameter
    double gamma;                     //< generalized alpha solver parameter
    double beta;                      //< generalized alpha solver parameter
    double gamma_prime;               //< generalized alpha solver parameter
    double beta_prime;                //< generalized alpha solver parameter
    double conditioner;               //< diagonal preconditioner value
    double absolute_convergence_tol;  //< absolute convergence tolerance
    double relative_convergence_tol;  //< relative convergence tolerance
    size_t static_load_retries;  //< maximum number of load reduction attempts for static analysis

    /**
     * @brief Constructor for the StepParameters object
     *
     * @param is_dynamic_solve_ if a dynamic solve is to be performed (versus static)
     * @param max_iter_ the maximum number of nonlinear iterations to be performed
     * @param h_ the time step size
     * @param rho_inf the numerical damping pactor, used to set the generalized alpha solver
     * parameters
     * @param a_tol the absolute error convergence tolerance
     * @param r_tol the relative error convergence tolerance
     * @param num_retries the maximum number of load reduction attempts for static solves
     */
    StepParameters(
        bool is_dynamic_solve_, size_t max_iter_, double h_, double rho_inf, double a_tol = 1e-5,
        double r_tol = 1e-3, size_t num_retries = 20
    )
        : is_dynamic_solve(is_dynamic_solve_),
          max_iter(max_iter_),
          h(h_),
          alpha_m((2. * rho_inf - 1.) / (rho_inf + 1.)),
          alpha_f(rho_inf / (rho_inf + 1.)),
          gamma(0.5 + alpha_f - alpha_m),
          beta(0.25 * (gamma + 0.5) * (gamma + 0.5)),
          gamma_prime((is_dynamic_solve) ? gamma / (h * beta) : 0.),
          beta_prime((is_dynamic_solve) ? (1. - alpha_m) / (h * h * beta * (1. - alpha_f)) : 0.),
          conditioner(beta * h * h),
          absolute_convergence_tol(a_tol),
          relative_convergence_tol(r_tol),
          static_load_retries(num_retries) {}
};

}  // namespace kynema
