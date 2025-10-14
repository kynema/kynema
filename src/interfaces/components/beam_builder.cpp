#include "beam_builder.hpp"

#include <numbers>

#include "beam.hpp"
#include "math/matrix_operations.hpp"
#include "math/quaternion_operations.hpp"

namespace kynema::interfaces::components {
BeamBuilder& BeamBuilder::SetElementOrder(size_t element_order) {
    input.element_order = element_order;
    return *this;
}

BeamBuilder& BeamBuilder::SetQuadratureRule(BeamInput::QuadratureRule rule) {
    input.quadrature_rule = rule;
    return *this;
}

BeamBuilder& BeamBuilder::SetSectionRefinement(size_t section_refinement) {
    input.section_refinement = section_refinement;
    return *this;
}

BeamBuilder& BeamBuilder::ClearReferenceAxisPoints() {
    input.ref_axis.coordinate_grid.clear();
    input.ref_axis.coordinates.clear();
    return *this;
}

BeamBuilder& BeamBuilder::AddRefAxisPoint(
    double grid_location, const std::array<double, 3>& coordinates, ReferenceAxisOrientation ref_axis
) {
    input.ref_axis.coordinate_grid.emplace_back(grid_location);
    if (ref_axis == ReferenceAxisOrientation::X) {
        input.ref_axis.coordinates.emplace_back(coordinates);
        return *this;
    }
    if (ref_axis == ReferenceAxisOrientation::Z) {
        const auto q_z_to_x = math::RotationVectorToQuaternion({0., std::numbers::pi / 2., 0.});
        input.ref_axis.coordinates.emplace_back(math::RotateVectorByQuaternion(q_z_to_x, coordinates)
        );
        return *this;
    }
    throw std::invalid_argument("Invalid reference axis orientation");
}

BeamBuilder& BeamBuilder::AddRefAxisTwist(double grid_location, double twist) {
    input.ref_axis.twist_grid.emplace_back(grid_location);
    input.ref_axis.twist.emplace_back(twist);
    return *this;
}

BeamBuilder& BeamBuilder::PrescribedRootMotion(bool enable) {
    input.root.prescribe_root_motion = enable;
    return *this;
}

BeamBuilder& BeamBuilder::SetRootPosition(const std::array<double, 7>& p) {
    input.root.position = p;
    return *this;
}

BeamBuilder& BeamBuilder::SetRootVelocity(const std::array<double, 6>& v) {
    input.root.velocity = v;
    return *this;
}

BeamBuilder& BeamBuilder::SetRootAcceleration(const std::array<double, 6>& a) {
    input.root.acceleration = a;
    return *this;
}

BeamBuilder& BeamBuilder::ClearSections() {
    input.ref_axis.coordinate_grid.clear();
    input.ref_axis.coordinates.clear();
    return *this;
}

BeamBuilder& BeamBuilder::AddSection(
    double grid_location, const std::array<std::array<double, 6>, 6>& mass_matrix,
    const std::array<std::array<double, 6>, 6>& stiffness_matrix, ReferenceAxisOrientation ref_axis
) {
    if (ref_axis == ReferenceAxisOrientation::X) {
        input.sections.emplace_back(grid_location, mass_matrix, stiffness_matrix);
        return *this;
    }
    if (ref_axis == ReferenceAxisOrientation::Z) {
        const auto q_z_to_x = math::RotationVectorToQuaternion({0., std::numbers::pi / 2., 0.});
        input.sections.emplace_back(
            grid_location, math::RotateMatrix6(mass_matrix, q_z_to_x),
            math::RotateMatrix6(stiffness_matrix, q_z_to_x)
        );
        return *this;
    }
    throw std::invalid_argument("Invalid reference axis orientation");
}

const BeamInput& BeamBuilder::Input() const {
    return this->input;
}

Beam BeamBuilder::Build(Model& model) const {
    return {this->input, model};
}
}  // namespace kynema::interfaces::components
