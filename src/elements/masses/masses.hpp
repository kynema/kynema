#pragma once

#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb {

/**
 * @brief Contains field variables for mass elements (aka, rigid bodies) to compute per-element
 * contributions to the residual vector and system/iteration matrix
 *
 * @note Mass elements consist of a single node and no quadrature points are required, hence no
 * node/qp prefix is used in the Views (as opposed to Beams)
 */
template <typename DeviceType>
struct Masses {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;

    size_t num_elems;                     //< Total number of elements
    View<size_t*> num_nodes_per_element;  //< This is always 1 for masses
    View<size_t*> state_indices;
    View<dof::FreedomSignature*> element_freedom_signature;
    View<size_t* [6]> element_freedom_table;
    View<double[3]> gravity;
    View<double* [7]> node_x0;      //< Initial position/rotation
    View<double* [6][6]> qp_Mstar;  //< Mass matrix in material frame
    View<double* [6]> residual_vector_terms;
    View<double* [6][6]> system_matrix_terms;

    explicit Masses(const size_t n_mass_elems)
        : num_elems(n_mass_elems),
          num_nodes_per_element(
              Kokkos::view_alloc("num_nodes_per_element", Kokkos::WithoutInitializing), num_elems
          ),
          state_indices(Kokkos::view_alloc("state_indices", Kokkos::WithoutInitializing), num_elems),
          element_freedom_signature(
              Kokkos::view_alloc("element_freedom_signature", Kokkos::WithoutInitializing), num_elems
          ),
          element_freedom_table(
              Kokkos::view_alloc("element_freedom_table", Kokkos::WithoutInitializing), num_elems
          ),
          gravity(Kokkos::view_alloc("gravity", Kokkos::WithoutInitializing)),
          node_x0(Kokkos::view_alloc("node_x0", Kokkos::WithoutInitializing), num_elems),
          qp_Mstar(Kokkos::view_alloc("qp_Mstar", Kokkos::WithoutInitializing), num_elems),
          residual_vector_terms(
              Kokkos::view_alloc("residual_vector_terms", Kokkos::WithoutInitializing), num_elems
          ),
          system_matrix_terms(
              Kokkos::view_alloc("system_matrix_terms", Kokkos::WithoutInitializing), num_elems
          ) {
        Kokkos::deep_copy(num_nodes_per_element, 1);  // Always 1 node per element
        Kokkos::deep_copy(element_freedom_signature, dof::FreedomSignature::AllComponents);
    }
};

}  // namespace kynema_fmb
