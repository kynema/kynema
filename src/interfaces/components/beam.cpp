#include "beam.hpp"

#include <numbers>
#include <stdexcept>

#include "elements/beams/beam_quadrature.hpp"
#include "interfaces/components/beam_input.hpp"
#include "math/gll_quadrature.hpp"
#include "math/interpolation.hpp"
#include "math/least_squares_fit.hpp"
#include "math/matrix_operations.hpp"
#include "math/project_points_to_target_polynomial.hpp"
#include "math/quaternion_operations.hpp"
#include "model/model.hpp"

namespace kynema::interfaces::components {

Beam::Beam(const BeamInput& input, Model& model) {
    ValidateInput(input);
    SetupNodeLocations(input);
    CreateNodeGeometry(input);
    CreateBeamElement(input, model);
    PositionBladeInSpace(input, model);
    SetupRootNode(input, model);
}

std::vector<double> Beam::GetNodeWeights(double s) const {
    std::vector<double> weights(this->node_xi.size());
    auto xi = 2. * s - 1.;
    math::LagrangePolynomialDerivWeights(xi, this->node_xi, weights);
    return weights;
}

void Beam::AddPointLoad(double s, std::array<double, 6> loads) {
    if (s < 0. || s > 1.) {
        throw std::invalid_argument("Invalid position: " + std::to_string(s));
    }

    const auto weights = this->GetNodeWeights(s);
    for (auto node : std::views::iota(0U, this->nodes.size())) {
        for (auto component : std::views::iota(0U, 6U)) {
            this->nodes[node].loads[component] += weights[node] * loads[component];
        }
    }
}

void Beam::ClearLoads() {
    for (auto& node : this->nodes) {
        node.ClearLoads();
    }
}

void Beam::GetMotion(const HostState<DeviceType>& host_state) {
    for (auto& node : this->nodes) {
        node.GetMotion(host_state);
    }
}

void Beam::SetLoads(HostState<DeviceType>& host_state) const {
    for (const auto& node : this->nodes) {
        node.SetLoads(host_state);
    }
}

void Beam::ValidateInput(const BeamInput& input) {
    if (input.ref_axis.coordinate_grid.size() < 2 || input.ref_axis.coordinates.size() < 2) {
        throw std::invalid_argument("At least two reference axis points are required");
    }
    if (input.ref_axis.coordinate_grid.size() != input.ref_axis.coordinates.size()) {
        throw std::invalid_argument("Mismatch between coordinate_grid and coordinates sizes");
    }
    if (input.sections.empty()) {
        throw std::invalid_argument("At least one section is required");
    }
    if (input.element_order < 1) {
        throw std::invalid_argument(
            "Element order must be at least 1 i.e. linear element for discretization"
        );
    }
}

void Beam::SetupNodeLocations(const BeamInput& input) {
    // Generate node locations within element [-1,1]
    this->node_xi = math::GenerateGLLPoints(input.element_order);
}

void Beam::CreateNodeGeometry(const BeamInput& input) {
    const auto n_nodes = input.element_order + 1;
    const auto n_geometry_pts =
        std::min({input.ref_axis.coordinate_grid.size(), n_nodes, kMaxGeometryPoints});

    if (n_geometry_pts < n_nodes) {
        // We need to project from n_geometry_pts -> element_order
        const auto kp_xi = math::MapGeometricLocations(input.ref_axis.coordinate_grid);
        const auto gll_points = math::GenerateGLLPoints(n_geometry_pts - 1);
        const auto phi_kn_geometry = math::ComputeShapeFunctionValues(kp_xi, gll_points);
        const auto geometry_points =
            math::PerformLeastSquaresFitting(phi_kn_geometry, input.ref_axis.coordinates);
        this->node_coordinates =
            math::ProjectPointsToTargetPolynomial(n_geometry_pts, n_nodes, geometry_points);

    } else {
        // Fit node coordinates to key points
        const auto kp_xi = math::MapGeometricLocations(input.ref_axis.coordinate_grid);
        const auto phi_kn = math::ComputeShapeFunctionValues(kp_xi, this->node_xi);
        this->node_coordinates =
            math::PerformLeastSquaresFitting(phi_kn, input.ref_axis.coordinates);
    }

    // Calculate tangent vectors at each node
    this->CalcNodeTangents();
}

void Beam::CreateBeamElement(const BeamInput& input, Model& model) {
    // Add nodes to model
    auto node_ids = std::vector<size_t>(this->node_xi.size());
    std::ranges::transform(
        std::views::iota(0U, this->node_xi.size()), std::begin(node_ids),
        [&](auto node) {
            const auto& pos = this->node_coordinates[node];
            const auto q_rot = math::TangentTwistToQuaternion(this->node_tangents[node], 0.);
            return model.AddNode()
                .SetElemLocation((this->node_xi[node] + 1.) / 2.)
                .SetPosition(pos[0], pos[1], pos[2], q_rot[0], q_rot[1], q_rot[2], q_rot[3])
                .Build();
        }
    );
    std::ranges::transform(node_ids, std::back_inserter(this->nodes), [](auto node_id) {
        return NodeData(node_id);
    });

    auto original_section_grid =
        std::vector<double>{input.sections.front().location, input.sections.back().location};
    if (input.quadrature_style == BeamInput::QuadratureStyle::Segmented) {
        original_section_grid.resize(input.sections.size());
        std::ranges::transform(
            input.sections, std::begin(original_section_grid),
            [](const auto& section) {
                return section.location;
            }
        );
    }

    // Build beam sections
    const auto sections = BuildBeamSections(input);
    auto section_grid = std::vector<double>(sections.size());
    std::ranges::transform(sections, section_grid.begin(), [](const auto& section) {
        return section.position;
    });

    // Calculate trapezoidal quadrature based on section locations
    const auto quadrature =
        (input.quadrature_rule == BeamInput::QuadratureRule::GaussLobatto)
            ? beams::CreateGaussLegendreLobattoQuadrature(
                  section_grid, original_section_grid, input.section_refinement + 1U
              )
            : beams::CreateGaussLegendreQuadrature(
                  section_grid, original_section_grid, input.section_refinement + 1U
              );

    // Add beam element and get ID
    this->beam_element_id = model.AddBeamElement(node_ids, sections, quadrature);
}

void Beam::PositionBladeInSpace(const BeamInput& input, Model& model) const {
    // Extract root location, orientation, and velocity
    const auto root_location =
        std::array{input.root.position[0], input.root.position[1], input.root.position[2]};
    const auto root_orientation = std::array{
        input.root.position[3], input.root.position[4], input.root.position[5],
        input.root.position[6]
    };
    const auto root_omega =
        std::array{input.root.velocity[3], input.root.velocity[4], input.root.velocity[5]};

    // Translate beam element to root location
    model.TranslateBeam(this->beam_element_id, root_location);

    // Rotate beam element about root location
    model.RotateBeamAboutPoint(this->beam_element_id, root_orientation, root_location);

    // Set beam velocity about root location
    model.SetBeamVelocityAboutPoint(this->beam_element_id, input.root.velocity, root_location);

    // Set beam acceleration about root location
    model.SetBeamAccelerationAboutPoint(
        this->beam_element_id, input.root.acceleration, root_omega, root_location
    );
}

void Beam::SetupRootNode(const BeamInput& input, Model& model) {
    // Add prescribed displacement constraint to root node if requested
    if (input.root.prescribe_root_motion) {
        this->prescribed_root_constraint_id = model.AddPrescribedBC(this->nodes[0].id);
    }
}

void Beam::CalcNodeTangents() {
    const auto n_nodes{this->node_coordinates.size()};

    // Calculate the derivative shape function matrix for the nodes
    const auto phi = math::ComputeShapeFunctionValues(this->node_xi, this->node_xi);
    const auto phi_prime = math::ComputeShapeFunctionDerivatives(this->node_xi, this->node_xi);

    // Calculate tangent vectors for each node
    this->node_tangents.resize(n_nodes, {0., 0., 0.});
    for (auto node_1 : std::views::iota(0U, n_nodes)) {
        for (auto node_2 : std::views::iota(0U, n_nodes)) {
            const auto shape_deriv = phi_prime[node_1][node_2];
            for (auto component : std::views::iota(0U, 3U)) {
                this->node_tangents[node_2][component] +=
                    shape_deriv * this->node_coordinates[node_1][component];
            }
        }
    }

    // Normalize tangent vectors
    std::ranges::transform(
        this->node_tangents, this->node_tangents.begin(),
        [](const std::array<double, 3>& tangent) {
            const auto normalized = Eigen::Matrix<double, 3, 1>(tangent.data()).normalized();
            return std::array{normalized(0), normalized(1), normalized(2)};
        }
    );
}

std::vector<BeamSection> Beam::BuildBeamSections(const BeamInput& input) {
    if (input.quadrature_style == BeamInput::QuadratureStyle::Segmented) {
        if (input.quadrature_rule == BeamInput::QuadratureRule::GaussLobatto) {
            return BuildBeamSections_SegmentedGLL(input);
        } else {
            return BuildBeamSections_SegmentedGL(input);
        }
    } else {
        if (input.quadrature_rule == BeamInput::QuadratureRule::GaussLobatto) {
            return BuildBeamSections_WholeBeamGLL(input);
        } else {
            return BuildBeamSections_WholeBeamGL(input);
        }
    }
}

std::vector<BeamSection> Beam::BuildBeamSections_SegmentedGLL(const BeamInput& input) {
    // Extraction section stiffness and mass matrices from blade definition
    std::vector<BeamSection> sections;

    // Add first section after rotating matrices to account for twist
    {
        const auto twist = math::LinearInterp(
            input.sections[0].location, input.ref_axis.twist_grid, input.ref_axis.twist
        );
        const auto q_twist = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(twist, Eigen::Matrix<double, 3, 1>::Unit(0))
        );
        const auto q_twist_array = std::array{q_twist.w(), q_twist.x(), q_twist.y(), q_twist.z()};
        sections.emplace_back(
            input.sections[0].location,
            math::RotateMatrix6(input.sections[0].mass_matrix, q_twist_array),
            math::RotateMatrix6(input.sections[0].stiffness_matrix, q_twist_array)
        );
    }

    // Loop through remaining section locations
    const auto quad_locations = math::GetGllLocations(input.section_refinement + 1U);
    const auto interior_nodes = input.section_refinement;
    for (auto section : std::views::iota(1U, input.sections.size())) {
        const auto section_location = input.sections[section].location;
        const auto section_mass_matrix = input.sections[section].mass_matrix;
        const auto section_stiffness_matrix = input.sections[section].stiffness_matrix;

        // Add refinement sections if requested

        for (auto refinement : std::views::iota(0U, interior_nodes)) {
            const auto left_location = input.sections[section - 1].location;
            const auto left_mass_matrix = input.sections[section - 1].mass_matrix;
            const auto left_stiffness_matrix = input.sections[section - 1].stiffness_matrix;

            // Calculate interpolation ratio between bounding sections
            const auto alpha = (quad_locations[refinement + 1] + 1.) / 2.;

            // Interpolate grid location
            const auto grid_value = (1. - alpha) * left_location + alpha * section_location;

            // Interpolate mass and stiffness matrices from bounding sections
            auto mass_matrix = std::array<std::array<double, 6>, 6>{};
            auto stiffness_matrix = std::array<std::array<double, 6>, 6>{};
            for (auto component_1 : std::views::iota(0U, 6U)) {
                for (auto component_2 : std::views::iota(0U, 6U)) {
                    mass_matrix[component_1][component_2] =
                        (1. - alpha) * left_mass_matrix[component_1][component_2] +
                        alpha * section_mass_matrix[component_1][component_2];
                    stiffness_matrix[component_1][component_2] =
                        (1. - alpha) * left_stiffness_matrix[component_1][component_2] +
                        alpha * section_stiffness_matrix[component_1][component_2];
                }
            }

            // Calculate twist at current section location via linear interpolation
            const auto twist =
                math::LinearInterp(grid_value, input.ref_axis.twist_grid, input.ref_axis.twist);

            // Add refinement section
            const auto q_twist = Eigen::Quaternion<double>(
                Eigen::AngleAxis<double>(twist, Eigen::Matrix<double, 3, 1>::Unit(0))
            );
            const auto q_twist_array =
                std::array{q_twist.w(), q_twist.x(), q_twist.y(), q_twist.z()};
            sections.emplace_back(
                grid_value, math::RotateMatrix6(mass_matrix, q_twist_array),
                math::RotateMatrix6(stiffness_matrix, q_twist_array)
            );
        }

        // Add ending section
        {
            const auto twist = math::LinearInterp(
                section_location, input.ref_axis.twist_grid, input.ref_axis.twist
            );
            const auto q_twist = Eigen::Quaternion<double>(
                Eigen::AngleAxis<double>(twist, Eigen::Matrix<double, 3, 1>::Unit(0))
            );
            const auto q_twist_array =
                std::array{q_twist.w(), q_twist.x(), q_twist.y(), q_twist.z()};
            sections.emplace_back(
                section_location, math::RotateMatrix6(section_mass_matrix, q_twist_array),
                math::RotateMatrix6(section_stiffness_matrix, q_twist_array)
            );
        }
    }

    return sections;
}

std::vector<BeamSection> Beam::BuildBeamSections_SegmentedGL(const BeamInput& input) {
    // Extraction section stiffness and mass matrices from blade definition
    std::vector<BeamSection> sections;

    // Loop through remaining section locations
    const auto quad_locations = math::GetGlLocations(input.section_refinement + 1U);
    const auto interior_nodes = input.section_refinement + 1U;
    for (auto section : std::views::iota(1U, input.sections.size())) {
        const auto section_location = input.sections[section].location;
        const auto section_mass_matrix = input.sections[section].mass_matrix;
        const auto section_stiffness_matrix = input.sections[section].stiffness_matrix;

        // Add refinement sections
        for (auto refinement : std::views::iota(0U, interior_nodes)) {
            const auto left_location = input.sections[section - 1].location;
            const auto left_mass_matrix = input.sections[section - 1].mass_matrix;
            const auto left_stiffness_matrix = input.sections[section - 1].stiffness_matrix;

            // Calculate interpolation ratio between bounding sections
            const auto alpha = (quad_locations[refinement] + 1.) / 2.;

            // Interpolate grid location
            const auto grid_value = (1. - alpha) * left_location + alpha * section_location;

            // Interpolate mass and stiffness matrices from bounding sections
            auto mass_matrix = std::array<std::array<double, 6>, 6>{};
            auto stiffness_matrix = std::array<std::array<double, 6>, 6>{};
            for (auto component_1 : std::views::iota(0U, 6U)) {
                for (auto component_2 : std::views::iota(0U, 6U)) {
                    mass_matrix[component_1][component_2] =
                        (1. - alpha) * left_mass_matrix[component_1][component_2] +
                        alpha * section_mass_matrix[component_1][component_2];
                    stiffness_matrix[component_1][component_2] =
                        (1. - alpha) * left_stiffness_matrix[component_1][component_2] +
                        alpha * section_stiffness_matrix[component_1][component_2];
                }
            }

            // Calculate twist at current section location via linear interpolation
            const auto twist =
                math::LinearInterp(grid_value, input.ref_axis.twist_grid, input.ref_axis.twist);

            // Add refinement section
            const auto q_twist = Eigen::Quaternion<double>(
                Eigen::AngleAxis<double>(twist, Eigen::Matrix<double, 3, 1>::Unit(0))
            );
            const auto q_twist_array =
                std::array{q_twist.w(), q_twist.x(), q_twist.y(), q_twist.z()};
            sections.emplace_back(
                grid_value, math::RotateMatrix6(mass_matrix, q_twist_array),
                math::RotateMatrix6(stiffness_matrix, q_twist_array)
            );
        }
    }
    return sections;
}

std::vector<BeamSection> Beam::BuildBeamSections_WholeBeam(
    const BeamInput& input, std::span<const double> quad_locations
) {
    auto sections = std::vector<BeamSection>{};
    const auto min_location = input.sections.front().location;
    const auto max_location = input.sections.back().location;

    auto grid_locations = std::vector<double>(quad_locations.size());
    std::ranges::transform(
        quad_locations, std::begin(grid_locations),
        [min_location, max_location](auto x) {
            const auto alpha = (x + 1.) / 2.;
            return (1. - alpha) * min_location + alpha * max_location;
        }
    );

    for (auto grid_value : grid_locations) {
        const auto twist =
            math::LinearInterp(grid_value, input.ref_axis.twist_grid, input.ref_axis.twist);
        const auto q_twist = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(twist, Eigen::Matrix<double, 3, 1>::Unit(0))
        );
        const auto q_twist_array = std::array{q_twist.w(), q_twist.x(), q_twist.y(), q_twist.z()};

        const auto left_bound = std::ranges::find_if(input.sections, [grid_value](const auto& s) {
            return s.location <= grid_value;
        });
        if (left_bound->location == grid_value) {
            sections.emplace_back(
                grid_value, math::RotateMatrix6(left_bound->mass_matrix, q_twist_array),
                math::RotateMatrix6(left_bound->stiffness_matrix, q_twist_array)
            );
        } else {
            const auto right_bound = std::next(left_bound);
            const auto left_location = left_bound->location;
            const auto right_location = right_bound->location;
            const auto& left_mass_matrix = left_bound->mass_matrix;
            const auto& right_mass_matrix = right_bound->mass_matrix;
            const auto& left_stiffness_matrix = left_bound->stiffness_matrix;
            const auto& right_stiffness_matrix = right_bound->stiffness_matrix;

            const auto alpha = (grid_value - left_location) / (right_location - left_location);
            auto mass_matrix = std::array<std::array<double, 6>, 6>{};
            auto stiffness_matrix = std::array<std::array<double, 6>, 6>{};
            for (auto component_1 : std::views::iota(0U, 6U)) {
                for (auto component_2 : std::views::iota(0U, 6U)) {
                    mass_matrix[component_1][component_2] =
                        (1. - alpha) * left_mass_matrix[component_1][component_2] +
                        alpha * right_mass_matrix[component_1][component_2];
                    stiffness_matrix[component_1][component_2] =
                        (1. - alpha) * left_stiffness_matrix[component_1][component_2] +
                        alpha * right_stiffness_matrix[component_1][component_2];
                }
            }

            sections.emplace_back(
                grid_value, math::RotateMatrix6(mass_matrix, q_twist_array),
                math::RotateMatrix6(stiffness_matrix, q_twist_array)
            );
        }
    }

    return sections;
}
std::vector<BeamSection> Beam::BuildBeamSections_WholeBeamGLL(const BeamInput& input) {
    return BuildBeamSections_WholeBeam(input, math::GetGllLocations(input.section_refinement + 1U));
}

std::vector<BeamSection> Beam::BuildBeamSections_WholeBeamGL(const BeamInput& input) {
    return BuildBeamSections_WholeBeam(input, math::GetGlLocations(input.section_refinement + 1U));
}

}  // namespace kynema::interfaces::components
