#pragma once

#include <filesystem>
#include <vector>

#include <Kokkos_Core.hpp>

namespace kynema_fmb::tests {

/// Function to find the project root directory
std::filesystem::path FindProjectRoot();

/// Function to remove a directory with several retries with a delay
void RemoveDirectoryWithRetries(
    const std::filesystem::path& dir, int retries = 3, int delayMs = 1000
);

Kokkos::View<double**> create_diagonal_matrix(const std::vector<double>& values);

// Check if members of the provided 1D Kokkos view is equal to the provided expected vector
void expect_kokkos_view_1D_equal(
    const Kokkos::View<const double*>&, const std::vector<double>&, double epsilon = 1.e-6
);

// Check if members of the provided 2D Kokkos view is equal to the provided expected matrix
void expect_kokkos_view_2D_equal(
    const Kokkos::View<const double**, Kokkos::LayoutStride>&,
    const std::vector<std::vector<double>>&, double epsilon = 1.e-6
);

// Check if members of the provided 3D Kokkos view is equal to the provided expected tensor
void expect_kokkos_view_3D_equal(
    const Kokkos::View<const double***, Kokkos::LayoutStride>& view,
    const std::vector<std::vector<std::vector<double>>>& expected, double epsilon = 1.e-6
);

// Convert a 1D Kokkos view to a vector
std::vector<double> kokkos_view_1D_to_vector(const Kokkos::View<double*>& view);

// Convert a 2D Kokkos view to a vector
std::vector<std::vector<double>> kokkos_view_2D_to_vector(const Kokkos::View<double**>& view);

// Convert a 3D Kokkos view to a vector
std::vector<std::vector<std::vector<double>>> kokkos_view_3D_to_vector(
    const Kokkos::View<double***>& view
);

}  // namespace kynema_fmb::tests
