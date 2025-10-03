#include "turbine.hpp"

#include <algorithm>
#include <array>
#include <numbers>
#include <ranges>
#include <vector>

#include "interfaces/components/turbine_input.hpp"
#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"
#include "model/model.hpp"

namespace kynema::interfaces::components {
Turbine::Turbine(const TurbineInput& input, Model& model)
    : blades(CreateBlades(input.blades, model)),
      tower(input.tower, model),
      hub_node(invalid_id),
      azimuth_node(invalid_id),
      shaft_base_node(invalid_id),
      yaw_bearing_node(invalid_id),
      tower_base{invalid_id},
      tower_top_to_yaw_bearing{invalid_id},
      yaw_bearing_to_shaft_base{invalid_id},
      shaft_base_to_azimuth{invalid_id},
      azimuth_to_hub{invalid_id},
      blade_pitch_control(input.blades.size(), input.blade_pitch_angle) {
    // Validate turbine inputs
    ValidateInput(input);

    // Add intermediate nodes and position them to reference locations
    PositionNodes(input, model);

    // Add mass elements
    AddMassElements(input, model);

    // Constraints must be added after all nodes are positioned because
    // they depend on the initial positions of the nodes
    AddConstraints(input, model);

    // Set initial conditions (velocities, accelerations etc.) after positioning
    // and adding constraints
    SetInitialConditions(input, model);
}

void Turbine::GetMotion(const HostState<DeviceType>& host_state) {
    for (auto& blade : this->blades) {
        blade.GetMotion(host_state);
    }
    this->tower.GetMotion(host_state);
}

void Turbine::SetLoads(HostState<DeviceType>& host_state) const {
    for (const auto& blade : this->blades) {
        blade.SetLoads(host_state);
    }
    this->tower.SetLoads(host_state);
}

const TurbineInput& Turbine::GetTurbineInput() const {
    return this->turbine_input;
}

std::vector<Beam> Turbine::CreateBlades(std::span<const BeamInput> blade_inputs, Model& model) {
    std::vector<Beam> blades;
    std::ranges::transform(
        blade_inputs, std::back_inserter(blades),
        [&model](const BeamInput& input) {
            return Beam(input, model);
        }
    );

    return blades;
}

void Turbine::ValidateInput(const TurbineInput& input) {
    if (input.hub_diameter <= kMinHubDiameter) {
        throw std::invalid_argument(
            "Hub diameter must be > " + std::to_string(kMinHubDiameter) +
            ", got: " + std::to_string(input.hub_diameter)
        );
    }

    if (input.rotor_speed < 0.) {
        throw std::invalid_argument(
            "rotor speed must be >= 0, got: " + std::to_string(input.rotor_speed)
        );
    }
}

void Turbine::PositionNodes(const TurbineInput& input, Model& model) {
    //--------------------------------------------------------------------------
    // Preliminary calculations
    //--------------------------------------------------------------------------

    // Define origin point for rotations
    const auto origin = std::array{0., 0., 0.};

    // Calculate angle between the blades
    const auto blade_angle_delta = 2. * std::numbers::pi / static_cast<double>(this->blades.size());

    // Calculate rotation quaternion to align a component from x-axis to z-axis (i.e.
    // rotate about y-axis by -90 degrees)
    const auto q_x_to_z = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(-std::numbers::pi / 2., Eigen::Matrix<double, 3, 1>::Unit(1))
    );
    const auto q_x_to_z_array = std::array{q_x_to_z.w(), q_x_to_z.x(), q_x_to_z.y(), q_x_to_z.z()};

    // Calculate cone angle rotation quaternion (rotate about y-axis by -cone angle)
    const auto q_cone = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(-input.cone_angle, Eigen::Matrix<double, 3, 1>::Unit(1))
    );
    const auto q_cone_array = std::array{q_cone.w(), q_cone.x(), q_cone.y(), q_cone.z()};

    //--------------------------------------------------------------------------
    // Position tower
    //--------------------------------------------------------------------------

    // Rotate tower to align with global Z-axis
    model.RotateBeamAboutPoint(this->tower.beam_element_id, q_x_to_z_array, origin);

    // Calculate hub position based on tower top and rotor apex offset
    const auto& tower_top_node = model.GetNode(this->tower.nodes.back().id);
    const auto apex_position = std::array{
        tower_top_node.x0[0] - input.tower_axis_to_rotor_apex,  // horizontal offset
        tower_top_node.x0[1],
        tower_top_node.x0[2] + input.tower_top_to_rotor_apex,  // vertical offset
    };

    //--------------------------------------------------------------------------
    // Create all intermediate nodes and populate node ID collections
    //--------------------------------------------------------------------------

    CreateIntermediateNodes(input, model);

    //--------------------------------------------------------------------------
    // Position blades
    //--------------------------------------------------------------------------

    // Loop over blades
    for (auto beam : std::views::iota(0U, this->blades.size())) {
        // Get node IDs for this blade only
        std::vector<size_t> current_blade_node_ids;
        std::ranges::transform(
            this->blades[beam].nodes, std::back_inserter(current_blade_node_ids),
            [](const NodeData& node) {
                return node.id;
            }
        );

        //----------------------------------------------------
        // Blade alignment and transformation
        //----------------------------------------------------

        // Loop through node IDs and rotate them to align with global Z-axis,
        // then translate to hub radius so blade root is at the hub radius
        for (const auto& node_id : current_blade_node_ids) {
            model.GetNode(node_id)
                .RotateAboutPoint(q_x_to_z_array, origin)  // Rotate to align with Z-axis
                .Translate(std::array{0., 0., input.hub_diameter / 2.});  // Translate to hub radius
        }

        // Add blade apex node -> current_blade_node_ids vector
        current_blade_node_ids.push_back(this->apex_nodes[beam].id);

        //----------------------------------------------------
        // Blade cone and azimuth transformations
        //----------------------------------------------------

        // Calculate azimuth angle rotation quaternion (rotate about x-axis)
        const auto q_azimuth = Eigen::Quaternion<double>(Eigen::AngleAxis<double>(
            static_cast<double>(beam) * blade_angle_delta, Eigen::Matrix<double, 3, 1>::Unit(0)
        ));
        const auto q_azimuth_array =
            std::array{q_azimuth.w(), q_azimuth.x(), q_azimuth.y(), q_azimuth.z()};

        // Rotate blade nodes (including apex node) -> cone angle -> azimuth angle
        for (const auto& node_id : current_blade_node_ids) {
            model.GetNode(node_id)
                .RotateAboutPoint(q_cone_array, origin)  // Rotate to cone angle (about y-axis)
                .RotateAboutPoint(
                    q_azimuth_array, origin
                );  // Rotate to azimuth angle (about x-axis)
        }
    }

    //--------------------------------------------------------------------------
    // Position rotor i.e. final rotor assembly
    //--------------------------------------------------------------------------

    // Create shaft rotation quaternion
    const auto q_shaft_tilt = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(input.shaft_tilt_angle, Eigen::Matrix<double, 3, 1>::Unit(1))
    );
    const auto q_shaft_tilt_array =
        std::array{q_shaft_tilt.w(), q_shaft_tilt.x(), q_shaft_tilt.y(), q_shaft_tilt.z()};

    // Rotate rotor nodes by shaft tilt angle and translate to apex position
    for (const auto& node_id : this->rotor_node_ids) {
        model.GetNode(node_id)
            .RotateAboutPoint(q_shaft_tilt_array, origin)  // Rotate about shaft tilt
            .Translate(apex_position);                     // Translate to apex position
    }
}

void Turbine::CreateIntermediateNodes(const TurbineInput& input, Model& model) {
    // Calculate shaft length from tower top to hub and overhang (adjusted for shaft tilt)
    const auto shaft_length = input.tower_axis_to_rotor_apex / cos(input.shaft_tilt_angle);

    // Get tower top node for yaw bearing positioning
    const auto& tower_top_node = model.GetNode(this->tower.nodes.back().id);

    //--------------------------------------------------------------------------
    // Create drivetrain nodes
    //--------------------------------------------------------------------------

    // Create hub node at origin with hub CM offset
    this->hub_node =
        NodeData(model.AddNode()
                     .SetPosition(std::array{-input.rotor_apex_to_hub, 0., 0., 1., 0., 0., 0.})
                     .Build());

    // Create azimuth node at shaft length from origin
    this->azimuth_node = NodeData(
        model.AddNode().SetPosition(std::array{shaft_length, 0., 0., 1., 0., 0., 0.}).Build()
    );

    // Create shaft base node at shaft length from origin
    this->shaft_base_node = NodeData(
        model.AddNode().SetPosition(std::array{shaft_length, 0., 0., 1., 0., 0., 0.}).Build()
    );

    // Create yaw bearing node at tower top position
    this->yaw_bearing_node = NodeData(model.AddNode()
                                          .SetPosition(tower_top_node.x0)
                                          .SetOrientation(std::array{1., 0., 0., 0.})
                                          .Build());

    //--------------------------------------------------------------------------
    // Create blade apex nodes
    //--------------------------------------------------------------------------

    this->apex_nodes.reserve(this->blades.size());

    // Create one apex node per blade at origin
    for ([[maybe_unused]] auto beam : std::views::iota(0U, this->blades.size())) {
        const auto apex_node_id =
            model.AddNode().SetPosition(std::array{0., 0., 0., 1., 0., 0., 0.}).Build();
        this->apex_nodes.emplace_back(apex_node_id);
    }

    //--------------------------------------------------------------------------
    // Populate node ID collections
    //--------------------------------------------------------------------------

    //----------------------------------------------------
    // Tower nodes
    //----------------------------------------------------

    // Add tower nodes to tower_node_ids vector
    this->tower_node_ids.reserve(this->tower.nodes.size());
    std::ranges::transform(
        this->tower.nodes, std::back_inserter(this->tower_node_ids),
        [](const auto& node) {
            return node.id;
        }
    );

    //----------------------------------------------------
    // Drivetrain nodes
    //----------------------------------------------------

    // Add drivetrain nodes (yaw bearing, shaft base, azimuth, hub) to drivetrain_node_ids vector
    this->drivetrain_node_ids.reserve(4);
    this->drivetrain_node_ids = {
        this->yaw_bearing_node.id, this->shaft_base_node.id, this->azimuth_node.id, this->hub_node.id
    };

    //----------------------------------------------------
    // Blade nodes
    //----------------------------------------------------

    // Populate blade_node_ids with all blade nodes and apex nodes
    this->blade_node_ids.reserve(
        this->blades.size() * (this->blades[0].nodes.size() + 1)  // +1 for apex node
    );
    for (auto beam : std::views::iota(0U, this->blades.size())) {
        // Add blade structural nodes
        std::ranges::transform(
            this->blades[beam].nodes, std::back_inserter(this->blade_node_ids),
            [](const NodeData& node) {
                return node.id;
            }
        );
        // Add apex node for this blade
        this->blade_node_ids.push_back(this->apex_nodes[beam].id);
    }

    //----------------------------------------------------
    // Rotor nodes
    //----------------------------------------------------

    // Add rotor nodes (hub, azimuth, shaft base, and all blade nodes) to rotor_node_ids vector
    this->rotor_node_ids.reserve(this->blade_node_ids.size() + 3);
    this->rotor_node_ids.insert(
        this->rotor_node_ids.end(),
        {this->hub_node.id, this->azimuth_node.id, this->shaft_base_node.id}
    );
    this->rotor_node_ids.insert(
        this->rotor_node_ids.end(), this->blade_node_ids.begin(), this->blade_node_ids.end()
    );

    //----------------------------------------------------
    // Nacelle nodes
    //----------------------------------------------------

    // Add nacelle nodes (drivetrain nodes + blade nodes) to nacelle_node_ids vector
    this->nacelle_node_ids.reserve(this->drivetrain_node_ids.size() + this->blade_node_ids.size());
    this->nacelle_node_ids.insert(
        this->nacelle_node_ids.end(), this->drivetrain_node_ids.begin(),
        this->drivetrain_node_ids.end()
    );
    this->nacelle_node_ids.insert(
        this->nacelle_node_ids.end(), this->blade_node_ids.begin(), this->blade_node_ids.end()
    );

    //----------------------------------------------------
    // All turbine nodes
    //----------------------------------------------------

    // Collect all turbine node IDs (tower + nacelle nodes)
    this->all_turbine_node_ids.reserve(this->tower_node_ids.size() + this->nacelle_node_ids.size());
    this->all_turbine_node_ids.insert(
        this->all_turbine_node_ids.end(), this->tower_node_ids.begin(), this->tower_node_ids.end()
    );
    this->all_turbine_node_ids.insert(
        this->all_turbine_node_ids.end(), this->nacelle_node_ids.begin(),
        this->nacelle_node_ids.end()
    );
}

void Turbine::AddMassElements(const TurbineInput& input, Model& model) {
    // Add mass element at yaw bearing node (nacelle mass + yaw bearing mass)
    this->yaw_bearing_mass_element_id =
        model.AddMassElement(this->yaw_bearing_node.id, input.yaw_bearing_inertia_matrix);

    // Add mass element at hub node (hub assembly mass)
    this->hub_mass_element_id = model.AddMassElement(this->hub_node.id, input.hub_inertia_matrix);
}

void Turbine::AddConstraints(const TurbineInput& input, Model& model) {
    //--------------------------------------------------------------------------
    // Blade control constraints
    //--------------------------------------------------------------------------

    // Loop through blades
    for (auto beam : std::views::iota(0U, this->blades.size())) {
        // Get the blade apex node
        const auto& apex_node = model.GetNode(this->apex_nodes[beam].id);

        // Get the blade root node
        const auto& root_node =
            model.GetNode(this->blades[beam].nodes[0].id);  // first node of blade

        // Calculate the pitch axis for the blade (from apex to root)
        const auto pitch_axis = std::array{
            root_node.x0[0] - apex_node.x0[0],
            root_node.x0[1] - apex_node.x0[1],
            root_node.x0[2] - apex_node.x0[2],
        };

        // Create pitch control constraint
        this->blade_pitch.emplace_back(ConstraintData{model.AddRotationControl(
            std::array{root_node.id, apex_node.id}, pitch_axis, &this->blade_pitch_control[beam]
        )});

        // Add rigid constraint between hub and blade apex
        this->apex_to_hub.emplace_back(ConstraintData{
            model.AddRigidJointConstraint(std::array{this->hub_node.id, apex_node.id})
        });
    }

    //--------------------------------------------------------------------------
    // Drivetrain constraints
    //--------------------------------------------------------------------------

    // Add rigid constraint between hub and azimuth node
    this->azimuth_to_hub = ConstraintData{
        model.AddRigidJointConstraint(std::array{this->azimuth_node.id, this->hub_node.id})
    };

    // Shaft axis constraint - add revolute joint between shaft base and azimuth node
    const auto shaft_axis =
        std::array{-cos(input.shaft_tilt_angle), 0., sin(input.shaft_tilt_angle)};
    this->shaft_base_to_azimuth = ConstraintData{model.AddRevoluteJointConstraint(
        std::array{this->shaft_base_node.id, this->azimuth_node.id}, shaft_axis, &torque_control
    )};

    // Add rigid constraint from yaw bearing to shaft base
    this->yaw_bearing_to_shaft_base = ConstraintData{model.AddRigidJointConstraint(
        std::array{this->yaw_bearing_node.id, this->shaft_base_node.id}
    )};

    //--------------------------------------------------------------------------
    // Nacelle control constraints
    //--------------------------------------------------------------------------

    // Add constraint from tower top to yaw bearing
    this->yaw_control = input.nacelle_yaw_angle;
    this->tower_top_to_yaw_bearing = ConstraintData{model.AddRotationControl(
        std::array{this->tower.nodes.back().id, this->yaw_bearing_node.id}, std::array{0., 0., 1.},
        &this->yaw_control
    )};

    //--------------------------------------------------------------------------
    // Tower base constraint
    //--------------------------------------------------------------------------

    // NOTE: We need to add this after we are done with setting up any initial displacements
    // to the tower base
}

void Turbine::SetInitialConditions(const TurbineInput& input, Model& model) {
    // Apply initial displacements
    SetInitialDisplacements(input, model);

    // Apply initial rotor velocity about shaft axis
    SetInitialRotorVelocity(input, model);

    // Apply initial accelerations
    // SetInitialAccelerations(input, model);
}

void Turbine::SetInitialDisplacements(const TurbineInput& input, Model& model) {
    //--------------------------------------------------------------------------
    // Apply initial blade pitch
    //--------------------------------------------------------------------------
    for (auto beam : std::views::iota(0U, this->blades.size())) {
        // Get the blade root and apex nodes
        const auto& root_node = model.GetNode(this->blades[beam].nodes.front().id);
        const auto& apex_node = model.GetNode(this->apex_nodes[beam].id);

        // Calculate the pitch axis
        const auto root_node_x0 = Eigen::Matrix<double, 3, 1>(root_node.x0.data());
        const auto apex_node_x0 = Eigen::Matrix<double, 3, 1>(apex_node.x0.data());
        const auto pitch_axis = root_node_x0 - apex_node_x0;

        // Create rotation vector about the pitch axis
        const auto pitch_axis_unit = pitch_axis.normalized();
        const auto rotation_vector = Eigen::AngleAxis(input.blade_pitch_angle, pitch_axis_unit);

        // Calculate pitch rotation quaternion
        const auto q_pitch = Eigen::Quaternion<double>(rotation_vector);
        const auto q_pitch_array = std::array{q_pitch.w(), q_pitch.x(), q_pitch.y(), q_pitch.z()};

        // Use apex node position as rotation center
        const auto pitch_center = std::array{apex_node.x0[0], apex_node.x0[1], apex_node.x0[2]};

        // Apply pitch rotation as displacement to all blade nodes and apex node
        for (const auto& blade_node : this->blades[beam].nodes) {
            model.GetNode(blade_node.id).RotateDisplacementAboutPoint(q_pitch_array, pitch_center);
        }
    }

    //--------------------------------------------------------------------------
    // Apply initial nacelle yaw rotation
    //--------------------------------------------------------------------------

    // Apply initial yaw rotation if non-zero
    if (std::abs(input.nacelle_yaw_angle) > kZeroTolerance) {
        // Get tower top node for rotation center
        const auto& tower_top_node = model.GetNode(this->tower.nodes.back().id);

        // Create yaw rotation quaternion (rotation about tower Z-axis)
        const auto q_yaw = Eigen::Quaternion<double>(
            Eigen::AngleAxis<double>(input.nacelle_yaw_angle, Eigen::Matrix<double, 3, 1>::Unit(2))
        );
        const auto q_yaw_array = std::array{q_yaw.w(), q_yaw.x(), q_yaw.y(), q_yaw.z()};

        // Rotate all nacelle components about tower top position
        for (const auto& node_id : this->nacelle_node_ids) {
            model.GetNode(node_id).RotateDisplacementAboutPoint(
                q_yaw_array,
                std::array{tower_top_node.x0[0], tower_top_node.x0[1], tower_top_node.x0[2]}
            );
        }
    }

    //--------------------------------------------------------------------------
    // Apply tower base displacement
    //--------------------------------------------------------------------------

    const auto& tower_base_node = model.GetNode(this->tower.nodes.front().id);
    const auto ref_tower_base_position =
        std::array{tower_base_node.x0[0], tower_base_node.x0[1], tower_base_node.x0[2]};

    // Calculate translation displacement from reference tower base -> input position
    const auto tower_base_displacement = std::array{
        input.tower_base_position[0] - ref_tower_base_position[0],
        input.tower_base_position[1] - ref_tower_base_position[1],
        input.tower_base_position[2] - ref_tower_base_position[2]
    };
    const auto tower_base_displacement_vec =
        Eigen::Matrix<double, 3, 1>(tower_base_displacement.data());

    // Get tower base orientation at input position
    const auto tower_base_orientation = std::array{
        input.tower_base_position[3], input.tower_base_position[4], input.tower_base_position[5],
        input.tower_base_position[6]
    };

    // Apply tower base displacement if displacement is non-zero or rotation is non-identity
    if (tower_base_displacement_vec.norm() > kZeroTolerance ||
        !math::IsIdentityQuaternion(tower_base_orientation, kZeroTolerance)) {
        // Apply displacement to all turbine nodes
        for (const auto& node_id : this->all_turbine_node_ids) {
            // first rotate about original tower base position
            model.GetNode(node_id).RotateDisplacementAboutPoint(
                tower_base_orientation, ref_tower_base_position
            );
            // then translate to new tower base position
            model.GetNode(node_id).TranslateDisplacement(tower_base_displacement);
        }
    }

    //--------------------------------------------------------------------------
    // Tower base constraint
    //--------------------------------------------------------------------------

    // Add prescribed BC constraint at the tower base with initial displacements
    this->tower_base = ConstraintData{model.AddPrescribedBC(tower_base_node.id, tower_base_node.u)};
}

void Turbine::SetInitialRotorVelocity(const TurbineInput& input, Model& model) {
    // Calculate shaft axis in current configuration
    const auto hub_position = model.GetNode(this->hub_node.id).DisplacedPosition();
    const auto shaft_base_position = model.GetNode(this->shaft_base_node.id).DisplacedPosition();
    const auto shaft_axis = (Eigen::Matrix<double, 3, 1>(hub_position.data()) -
                             Eigen::Matrix<double, 3, 1>(shaft_base_position.data()))
                                .normalized();

    // Collect all rotor node IDs (hub, azimuth, blade nodes, and apex nodes)
    std::vector<size_t> rotor_velocity_node_ids{this->hub_node.id, this->azimuth_node.id};

    // Add all blade nodes and apex nodes
    for (auto beam : std::views::iota(0U, this->blades.size())) {
        // Add blade nodes
        std::ranges::transform(
            this->blades[beam].nodes, std::back_inserter(rotor_velocity_node_ids),
            [](const auto& blade_node) {
                return blade_node.id;
            }
        );
        // Add apex node
        rotor_velocity_node_ids.push_back(this->apex_nodes[beam].id);
    }

    // Create rigid body velocity -> transl. vel = 0, angular vel. about shaft axis
    const auto rigid_body_velocity = std::array{
        0.,
        0.,
        0.,
        input.rotor_speed * shaft_axis[0],
        input.rotor_speed * shaft_axis[1],
        input.rotor_speed * shaft_axis[2]
    };

    // Apply rotational velocity to all rotor nodes about hub node
    for (const auto& node_id : rotor_velocity_node_ids) {
        model.GetNode(node_id).SetVelocityAboutPoint(
            rigid_body_velocity, std::array{hub_position[0], hub_position[1], hub_position[2]}
        );
    }
}
}  // namespace kynema::interfaces::components
