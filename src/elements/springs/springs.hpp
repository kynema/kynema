#pragma once

#include <Kokkos_Core.hpp>

#include "dof_management/freedom_signature.hpp"

namespace kynema_fmb {

/**
 * @brief Contains field variables for spring elements to compute per-element
 * contributions to the residual vector and system/iteration matrix
 */
template <typename DeviceType>
struct Springs {
    template <typename ValueType>
    using View = Kokkos::View<ValueType, DeviceType>;

    size_t num_elems;  //< Total number of elements

    View<size_t*> num_nodes_per_element;  //< This is always 2 for springs
    View<size_t* [2]> node_state_indices;
    View<dof::FreedomSignature* [2]> element_freedom_signature;
    View<size_t* [2][3]> element_freedom_table;  //< Only translational DOFs for springs

    View<double* [3]> x0;  //< Initial distance vector between nodes
    View<double*> l_ref;   //< Initial length of springs
    View<double*> k;       //< Spring stiffness coefficients

    View<double* [2][3]> residual_vector_terms;
    View<double* [2][2][3][3]> stiffness_matrix_terms;

    explicit Springs(const size_t n_spring_elems)
        : num_elems(n_spring_elems),
          num_nodes_per_element(
              Kokkos::view_alloc("num_nodes_per_element", Kokkos::WithoutInitializing), num_elems
          ),
          node_state_indices(
              Kokkos::view_alloc("node_state_indices", Kokkos::WithoutInitializing), num_elems
          ),
          element_freedom_signature(
              Kokkos::view_alloc("element_freedom_signature", Kokkos::WithoutInitializing), num_elems
          ),
          element_freedom_table(
              Kokkos::view_alloc("element_freedom_table", Kokkos::WithoutInitializing), num_elems
          ),
          x0(Kokkos::view_alloc("x0", Kokkos::WithoutInitializing), num_elems),
          l_ref(Kokkos::view_alloc("l_ref", Kokkos::WithoutInitializing), num_elems),
          k(Kokkos::view_alloc("k", Kokkos::WithoutInitializing), num_elems),
          residual_vector_terms(
              Kokkos::view_alloc("residual_vector_terms", Kokkos::WithoutInitializing), num_elems
          ),
          stiffness_matrix_terms(
              Kokkos::view_alloc("stiffness_matrix_terms", Kokkos::WithoutInitializing), num_elems
          ) {
        Kokkos::deep_copy(num_nodes_per_element, 2);  // Always 2 nodes per element
        Kokkos::deep_copy(
            element_freedom_signature, dof::FreedomSignature::JustPosition
        );  // Springs only have translational DOFs
    }
};

}  // namespace kynema_fmb
