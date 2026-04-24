#include <gtest/gtest.h>

#include "regression/test_utilities.hpp"
#include "utilities/aerodynamics/aerodyn_inflow.hpp"

#include "Kynema_FMB_config.h"

namespace kynema_fmb::tests {

TEST(AerodynInflowTest, BladeInitialState_Constructor) {
    constexpr auto root = std::array{1., 2., 3., 1., 0., 0., 0.};
    constexpr auto nodes =
        std::array{std::array{4., 5., 6., 1., 0., 0., 0.}, std::array{7., 8., 9., 1., 0., 0., 0.}};
    const util::TurbineConfig::BladeInitialState blade_state{root, nodes};

    EXPECT_EQ(blade_state.root_initial_position, root);
    EXPECT_EQ(blade_state.node_initial_positions, nodes);
}

TEST(AerodynInflowTest, TurbineConfig_Constructor) {
    constexpr bool is_hawt{true};
    constexpr auto ref_pos = std::array{10.F, 20.F, 30.F};
    constexpr auto hub_pos = std::array{1., 2., 3., 1., 0., 0., 0.};
    constexpr auto nacelle_pos = std::array{4., 5., 6., 1., 0., 0., 0.};

    // Create 3 blades with 2 nodes each
    std::vector<util::TurbineConfig::BladeInitialState> blade_states;
    for (auto i : std::views::iota(0, 3)) {
        constexpr auto root = std::array{static_cast<double>(i), 0., 0., 1., 0., 0., 0.};
        constexpr auto nodes = std::array{
            std::array{static_cast<double>(i), 1., 0., 1., 0., 0., 0.},
            std::array{static_cast<double>(i), 2., 0., 1., 0., 0., 0.}
        };
        blade_states.emplace_back(root, nodes);
    }

    util::TurbineConfig turbine_config{is_hawt, ref_pos, hub_pos, nacelle_pos, blade_states};

    EXPECT_EQ(turbine_config.is_horizontal_axis, is_hawt);
    EXPECT_EQ(turbine_config.reference_position, ref_pos);
    EXPECT_EQ(turbine_config.hub_initial_position, hub_pos);
    EXPECT_EQ(turbine_config.nacelle_initial_position, nacelle_pos);
    EXPECT_EQ(turbine_config.NumberOfBlades(), 3);

    for (auto i : std::views::iota(0U, turbine_config.NumberOfBlades())) {
        EXPECT_EQ(
            turbine_config.blade_initial_states[i].root_initial_position[0], static_cast<double>(i)
        );
        EXPECT_EQ(turbine_config.blade_initial_states[i].node_initial_positions.size(), 2);
        EXPECT_EQ(
            turbine_config.blade_initial_states[i].node_initial_positions[0][0],
            static_cast<double>(i)
        );
        EXPECT_EQ(
            turbine_config.blade_initial_states[i].node_initial_positions[1][0],
            static_cast<double>(i)
        );
    }
}

TEST(AerodynInflowTest, TurbineConfig_Validate_InvalidConfiguration) {
    // Invalid configuration: No blades
    constexpr bool is_hawt{true};
    constexpr auto ref_pos = std::array{10.F, 20.F, 30.F};
    constexpr auto hub_pos = std::array{1., 2., 3., 1., 0., 0., 0.};
    constexpr auto nacelle_pos = std::array{4., 5., 6., 1., 0., 0., 0.};
    const std::vector<util::TurbineConfig::BladeInitialState> empty_blade_states;

    EXPECT_THROW(
        util::TurbineConfig(is_hawt, ref_pos, hub_pos, nacelle_pos, empty_blade_states),
        std::runtime_error
    );
}

/// Check if members of the provided array is equal to the provided expected array
template <typename T>
void ExpectArrayNear(
    std::span<const T> actual, std::array<const T> expected, T epsilon = static_cast<T>(1.e-6)
) {
    ASSERT_EQ(actual.size(), expected.size());
    for (auto i : std::views::iota(0, N)) {
        EXPECT_NEAR(actual[i], expected[i], epsilon) << "Element mismatch at index " << i;
    }
}

TEST(AerodynInflowTest, SetPositionAndOrientation) {
    constexpr auto data = std::array{1., 2., 3., 0.707107, 0.707107, 0., 0.};
    constexpr auto position = std::array{0.F, 0.F, 0.F};
    auto orientation = std::array{0., 0., 0., 0., 0., 0., 0., 0., 0.};

    util::SetPositionAndOrientation(data, position, orientation);

    ExpectArrayNear(position, std::array{1.F, 2.F, 3.F});
    ExpectArrayNear(orientation, std::array{1., 0., 0., 0., 0., -1., 0., 1., 0.});
}

TEST(AerodynInflowTest, MeshData_Constructor_NumberOfNodes) {
    const util::MeshData mesh_motion_data{1};
    EXPECT_EQ(mesh_motion_data.NumberOfMeshPoints(), 1);
    EXPECT_EQ(mesh_motion_data.position.size(), 1);
    EXPECT_EQ(mesh_motion_data.orientation.size(), 1);
    EXPECT_EQ(mesh_motion_data.velocity.size(), 1);
    EXPECT_EQ(mesh_motion_data.acceleration.size(), 1);
    EXPECT_EQ(mesh_motion_data.loads.size(), 1);
}

TEST(AerodynInflowTest, MeshData_Constructor_Data) {
    const size_t n_mesh_points{1};
    constexpr auto mesh_data = std::array{std::array{1., 2., 3., 0.707107, 0.707107, 0., 0.}};
    constexpr auto mesh_velocities = std::array{std::array{1.F, 2.F, 3.F, 4.F, 5.F, 6.F}};
    constexpr auto mesh_accelerations = std::array{std::array{7.F, 8.F, 9.F, 10.F, 11.F, 12.F}};
    constexpr auto mesh_loads = std::array{std::array{13.F, 14.F, 15.F, 16.F, 17.F, 18.F}};
    util::MeshData mesh_motion_data(
        n_mesh_points, mesh_data, mesh_velocities, mesh_accelerations, mesh_loads
    );

    EXPECT_EQ(mesh_motion_data.NumberOfMeshPoints(), n_mesh_points);
    EXPECT_EQ(mesh_motion_data.position.size(), n_mesh_points);
    EXPECT_EQ(mesh_motion_data.orientation.size(), n_mesh_points);
    EXPECT_EQ(mesh_motion_data.velocity.size(), n_mesh_points);
    EXPECT_EQ(mesh_motion_data.acceleration.size(), n_mesh_points);
    EXPECT_EQ(mesh_motion_data.loads.size(), n_mesh_points);
    ExpectArrayNear(mesh_motion_data.position[0], {1.F, 2.F, 3.F});
    ExpectArrayNear(mesh_motion_data.orientation[0], {1., 0., 0., 0., 0., -1., 0., 1., 0.});
    ExpectArrayNear(mesh_motion_data.velocity[0], {1.F, 2.F, 3.F, 4.F, 5.F, 6.F});
    ExpectArrayNear(mesh_motion_data.acceleration[0], {7.F, 8.F, 9.F, 10.F, 11.F, 12.F});
    ExpectArrayNear(mesh_motion_data.loads[0], {13.F, 14.F, 15.F, 16.F, 17.F, 18.F});
}

class MeshDataValidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a mesh data with 2 mesh points
        mesh_data = std::make_unique<util::MeshData>(2);
        mesh_data->position = {{1.F, 2.F, 3.F}, {4.F, 5.F, 6.F}};
        mesh_data->orientation = {
            {1., 0., 0., 0., 1., 0., 0., 0., 1.}, {0., 1., 0., -1., 0., 0., 0., 0., 1.}
        };
        mesh_data->velocity = {{1.F, 2.F, 3.F, 4.F, 5.F, 6.F}, {7.F, 8.F, 9.F, 10.F, 11.F, 12.F}};
        mesh_data->acceleration = {
            {1.F, 2.F, 3.F, 4.F, 5.F, 6.F}, {7.F, 8.F, 9.F, 10.F, 11.F, 12.F}
        };
        mesh_data->loads = {{1.F, 2.F, 3.F, 4.F, 5.F, 6.F}, {7.F, 8.F, 9.F, 10.F, 11.F, 12.F}};
    }

    std::unique_ptr<util::MeshData> mesh_data;
};

TEST_F(MeshDataValidationTest, InvalidNumberOfMeshPoints) {
    mesh_data->n_mesh_points = 0;  // Should be at least 1
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST_F(MeshDataValidationTest, MismatchedPositionSize) {
    mesh_data->position.pop_back();
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST_F(MeshDataValidationTest, MismatchedOrientationSize) {
    mesh_data->orientation.pop_back();
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST_F(MeshDataValidationTest, MismatchedVelocitySize) {
    mesh_data->velocity.pop_back();
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST_F(MeshDataValidationTest, MismatchedAccelerationSize) {
    mesh_data->acceleration.pop_back();
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST_F(MeshDataValidationTest, MismatchedLoadsSize) {
    mesh_data->loads.pop_back();
    EXPECT_THROW(mesh_data->Validate(), std::invalid_argument);
}

TEST(AerodynInflowTest, TurbineData_Constructor) {
    // Set up 3 blades with 2 nodes each
    std::vector<util::TurbineConfig::BladeInitialState> blade_states;
    for (auto i : std::views::iota(0, 3)) {
        const util::TurbineConfig::BladeInitialState blade_state(
            {0., 0., 90., 1., 0., 0., 0.},  // root_initial_position
            {
                {0., 5., 90., 1., 0., 0., 0.},  // node_initial_positions - 1
                {0., 10., 90., 1., 0., 0., 0.}  // node_initial_positions - 2
            }

        );
        blade_states.push_back(blade_state);
    }
    const util::TurbineConfig tc(
        true,                           // is_horizontal_axis
        {0.F, 0.F, 0.F},                // reference_position
        {0., 0., 90., 1., 0., 0., 0.},  // hub_initial_position
        {0., 0., 90., 1., 0., 0., 0.},  // nacelle_initial_position
        blade_states
    );

    util::TurbineData turbine_data(tc);

    // Check basic properties
    EXPECT_EQ(turbine_data.NumberOfBlades(), 3);
    EXPECT_EQ(turbine_data.hub.NumberOfMeshPoints(), 1);
    EXPECT_EQ(turbine_data.nacelle.NumberOfMeshPoints(), 1);
    EXPECT_EQ(turbine_data.blade_roots.NumberOfMeshPoints(), 3);
    EXPECT_EQ(turbine_data.blade_nodes.NumberOfMeshPoints(), 6);  // 3 blades * 2 nodes each

    // Check hub and nacelle positions
    ExpectArrayNear(turbine_data.hub.position[0], {0.F, 0.F, 90.F});
    ExpectArrayNear(turbine_data.nacelle.position[0], {0.F, 0.F, 90.F});

    // Check blade roots
    for (auto i : std::views::iota(0U, turbine_data.NumberOfBlades())) {
        ExpectArrayNear(turbine_data.blade_roots.position[i], {0.F, 0.F, 90.F});
    }

    // Check blade nodes
    for (auto i : std::views::iota(0U, turbine_data.blade_nodes.NumberOfMeshPoints())) {
        const float expected_y = (i % 2 == 0) ? 5.F : 10.F;
        ExpectArrayNear(turbine_data.blade_nodes.position[i], {0.F, expected_y, 90.F});
    }

    // Check blade_nodes_to_blade_num_mapping
    std::vector<int32_t> expected_blade_nodes_to_blade_num_mapping = {1, 1, 2, 2, 3, 3};
    ASSERT_EQ(
        turbine_data.blade_nodes_to_blade_num_mapping.size(),
        turbine_data.blade_nodes.NumberOfMeshPoints()
    );
    for (auto i : std::views::iota(0U, turbine_data.blade_nodes.NumberOfMeshPoints())) {
        EXPECT_EQ(
            turbine_data.blade_nodes_to_blade_num_mapping[i],
            expected_blade_nodes_to_blade_num_mapping[i]
        ) << "Mismatch at index "
          << i;
    }

    // Check node_indices_by_blade
    ASSERT_EQ(turbine_data.node_indices_by_blade.size(), turbine_data.NumberOfBlades());
    std::vector<std::vector<size_t>> expected_node_indices_by_blade = {{0, 1}, {2, 3}, {4, 5}};
    for (auto i : std::views::iota(0U, turbine_data.NumberOfBlades())) {
        ASSERT_EQ(turbine_data.node_indices_by_blade[blade].size(), 2)
            << "Incorrect number of nodes for blade " << blade;
        EXPECT_EQ(turbine_data.node_indices_by_blade[blade], expected_node_indices_by_blade[blade])
            << "Incorrect node indices for blade " << blade;
    }

    // Additional check: Verify that node_indices_by_blade correctly maps to
    // blade_nodes_to_blade_num_mapping
    for (auto blade : std::views::iota(0U, turbine_data.NumberOfBlades())) {
        for (const size_t node : turbine_data.node_indices_by_blade[blade]) {
            EXPECT_EQ(turbine_data.blade_nodes_to_blade_num_mapping[node], blade + 1)
                << "Mismatch between node_indices_by_blade and blade_nodes_to_blade_num_mapping for "
                   "blade "
                << blade;
        }
    }
}

class TurbineDataValidationTest : public ::testing::Test {
protected:
    void SetUp() override {
        blade_states.clear();
        // Set up 3 blades with 2 nodes each
        for (auto i : std::views::iota(0, 3)) {
            const util::TurbineConfig::BladeInitialState blade_state(
                {0., 0., 90., 1., 0., 0., 0.},  // root_initial_position
                {
                    {0., 5., 90., 1., 0., 0., 0.},  // node_initial_positions - 1
                    {0., 10., 90., 1., 0., 0., 0.}  // node_initial_positions - 2
                }
            );
            blade_states.push_back(blade_state);
        }
        tc = std::make_unique<util::TurbineConfig>(
            true,                                                // is_horizontal_axis
            std::array<float, 3>{0.F, 0.F, 0.F},                 // reference_position
            std::array<double, 7>{0., 0., 90., 1., 0., 0., 0.},  // hub_initial_position
            std::array<double, 7>{0., 0., 90., 1., 0., 0., 0.},  // nacelle_initial_position
            blade_states
        );
        turbine_data = std::make_unique<util::TurbineData>(*tc);
    }

    std::vector<util::TurbineConfig::BladeInitialState> blade_states;
    std::unique_ptr<util::TurbineConfig> tc;
    std::unique_ptr<util::TurbineData> turbine_data;
};

TEST_F(TurbineDataValidationTest, InvalidNumberOfBlades) {
    turbine_data->n_blades = 0;  // Should be at least 1
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, MismatchBladeRoots) {
    turbine_data->blade_roots.n_mesh_points = 2;  // Should be 3
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, MismatchBladeNodeMapping) {
    turbine_data->blade_nodes_to_blade_num_mapping.pop_back();
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, MismatchNodeIndices) {
    turbine_data->node_indices_by_blade.pop_back();
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, MismatchTotalBladeNodes) {
    turbine_data->blade_nodes.n_mesh_points = 5;  // Should be 6
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, InvalidHubMeshPoints) {
    turbine_data->hub.n_mesh_points = 2;  // Should be 1
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST_F(TurbineDataValidationTest, InvalidNacelleMeshPoints) {
    turbine_data->nacelle.n_mesh_points = 0;  // Should be 1
    EXPECT_THROW(turbine_data->Validate(), std::runtime_error);
}

TEST(AerodynInflowTest, TurbineData_SetBladeNodeValues) {
    // Set up 3 blades with 2 nodes each
    std::vector<util::TurbineConfig::BladeInitialState> blade_states;
    for (auto i : std::views::iota(0, 3)) {
        const util::TurbineConfig::BladeInitialState blade_state(
            {0., 0., 90., 1., 0., 0., 0.},  // root_initial_position
            {
                {0., 5., 90., 1., 0., 0., 0.},  // node_initial_positions - 1
                {0., 10., 90., 1., 0., 0., 0.}  // node_initial_positions - 2
            }
        );
        blade_states.push_back(blade_state);
    }
    const util::TurbineConfig tc(
        true,                           // is_horizontal_axis
        {0.F, 0.F, 0.F},                // reference_position
        {0., 0., 90., 1., 0., 0., 0.},  // hub_initial_position
        {0., 0., 90., 1., 0., 0., 0.},  // nacelle_initial_position
        blade_states
    );

    util::TurbineData turbine_data(tc);

    // Verify the current values for the first blade and node
    const size_t blade_number{1};
    const size_t node_number{0};
    const size_t node_index{turbine_data.node_indices_by_blade[blade_number][node_number]};

    ExpectArrayNear(turbine_data.blade_nodes.position[node_index], {0.F, 5.F, 90.F});
    ExpectArrayNear(
        turbine_data.blade_nodes.orientation[node_index], {1., 0., 0., 0., 1., 0., 0., 0., 1.}
    );
    ExpectArrayNear(turbine_data.blade_nodes.velocity[node_index], {0.F, 0.F, 0.F, 0.F, 0.F, 0.F});
    ExpectArrayNear(
        turbine_data.blade_nodes.acceleration[node_index], {0.F, 0.F, 0.F, 0.F, 0.F, 0.F}
    );
    ExpectArrayNear(turbine_data.blade_nodes.loads[node_index], {0.F, 0.F, 0.F, 0.F, 0.F, 0.F});

    // Define new values for the node
    const std::array<float, 3> new_position = {1.F, 2.F, 3.F};
    const std::array<double, 9> new_orientation = {1., 0., 0., 0., 1., 0., 0., 0., 1.};
    const std::array<float, 6> new_velocity = {1.F, 2.F, 3.F, 4.F, 5.F, 6.F};
    const std::array<float, 6> new_acceleration = {7.F, 8.F, 9.F, 10.F, 11.F, 12.F};
    const std::array<float, 6> new_loads = {13.F, 14.F, 15.F, 16.F, 17.F, 18.F};
    turbine_data.SetBladeNodeMotion(
        blade_number, node_number, new_position, new_orientation, new_velocity, new_acceleration,
        new_loads
    );

    // Verify that the values were set correctly
    ExpectArrayNear(turbine_data.blade_nodes.position[node_index], new_position);
    ExpectArrayNear(turbine_data.blade_nodes.orientation[node_index], new_orientation);
    ExpectArrayNear(turbine_data.blade_nodes.velocity[node_index], new_velocity);
    ExpectArrayNear(turbine_data.blade_nodes.acceleration[node_index], new_acceleration);
    ExpectArrayNear(turbine_data.blade_nodes.loads[node_index], new_loads);
}

TEST(AerodynInflowTest, SimulationControls_Default) {
    util::SimulationControls simulation_controls;
    EXPECT_EQ(simulation_controls.is_aerodyn_input_path, true);
    EXPECT_EQ(simulation_controls.is_inflowwind_input_path, true);
    EXPECT_EQ(simulation_controls.interpolation_order, 1);
    EXPECT_EQ(simulation_controls.time_step, 0.1);
    EXPECT_EQ(simulation_controls.max_time, 600.);
    EXPECT_EQ(simulation_controls.total_elapsed_time, 0.);
    EXPECT_EQ(simulation_controls.n_time_steps, 0);
    EXPECT_EQ(simulation_controls.store_HH_wind_speed, 1);
    EXPECT_EQ(simulation_controls.transpose_DCM, 1);
    EXPECT_EQ(simulation_controls.debug_level, 0);
    EXPECT_EQ(simulation_controls.output_format, 1);
    EXPECT_EQ(simulation_controls.output_time_step, 0.1);
    EXPECT_STREQ(simulation_controls.output_root_name.data(), "ADI_out");
    EXPECT_EQ(simulation_controls.n_channels, 0);
    EXPECT_STREQ(simulation_controls.channel_names_c.data(), "");
    EXPECT_STREQ(simulation_controls.channel_units_c.data(), "");
}

TEST(AerodynInflowTest, SimulationControls_Set) {
    util::SimulationControls simulation_controls;
    simulation_controls.is_aerodyn_input_path = false;
    simulation_controls.is_inflowwind_input_path = false;
    simulation_controls.interpolation_order = 2;
    simulation_controls.time_step = 0.2;
    simulation_controls.max_time = 1200.;
    simulation_controls.total_elapsed_time = 100.;
    simulation_controls.n_time_steps = 10;
    simulation_controls.store_HH_wind_speed = false;
    simulation_controls.transpose_DCM = false;
    simulation_controls.debug_level =
        static_cast<int>(util::SimulationControls::DebugLevel::kSummary);
    simulation_controls.output_format = 1;
    simulation_controls.output_time_step = 1.;
    simulation_controls.n_channels = 1;
    std::strncpy(
        simulation_controls.channel_names_c.data(), "test_channel",
        simulation_controls.channel_names_c.size() - 1
    );
    simulation_controls.channel_names_c[simulation_controls.channel_names_c.size() - 1] = '\0';
    std::strncpy(
        simulation_controls.channel_units_c.data(), "test_unit",
        simulation_controls.channel_units_c.size() - 1
    );
    simulation_controls.channel_units_c[simulation_controls.channel_units_c.size() - 1] = '\0';

    EXPECT_EQ(simulation_controls.is_aerodyn_input_path, 0);
    EXPECT_EQ(simulation_controls.is_inflowwind_input_path, 0);
    EXPECT_EQ(simulation_controls.interpolation_order, 2);
    EXPECT_EQ(simulation_controls.time_step, 0.2);
    EXPECT_EQ(simulation_controls.max_time, 1200.);
    EXPECT_EQ(simulation_controls.total_elapsed_time, 100.);
    EXPECT_EQ(simulation_controls.n_time_steps, 10);
    EXPECT_EQ(simulation_controls.store_HH_wind_speed, 0);
    EXPECT_EQ(simulation_controls.transpose_DCM, 0);
    EXPECT_EQ(simulation_controls.debug_level, 1);
    EXPECT_EQ(simulation_controls.output_format, 1);
    EXPECT_EQ(simulation_controls.output_time_step, 1.);
    EXPECT_STREQ(simulation_controls.output_root_name.data(), "ADI_out");
    EXPECT_EQ(simulation_controls.n_channels, 1);
    EXPECT_STREQ(simulation_controls.channel_names_c.data(), "test_channel");
    EXPECT_STREQ(simulation_controls.channel_units_c.data(), "test_unit");
}

TEST(AerodynInflowTest, ErrorHandling_NoThrow) {
    const util::ErrorHandling error_handling;
    EXPECT_NO_THROW(error_handling.CheckError());
}

TEST(AerodynInflowTest, ErrorHandling_Throw) {
    util::ErrorHandling error_handling;
    error_handling.error_status = 1;  // Set error status to 1 to trigger an error
    EXPECT_THROW(error_handling.CheckError(), std::runtime_error);
}

TEST(AerodynInflowTest, FluidProperties_Default) {
    const util::FluidProperties fluid_properties;
    EXPECT_NEAR(fluid_properties.density, 1.225F, 1e-6F);
    EXPECT_NEAR(fluid_properties.kinematic_viscosity, 1.464E-5F, 1e-6F);
    EXPECT_NEAR(fluid_properties.sound_speed, 335.F, 1e-6F);
    EXPECT_NEAR(fluid_properties.vapor_pressure, 1700.F, 1e-6F);
}

TEST(AerodynInflowTest, FluidProperties_Set) {
    util::FluidProperties fluid_properties;
    fluid_properties.density = 1.1F;
    fluid_properties.kinematic_viscosity = 1.5E-5F;
    fluid_properties.sound_speed = 340.F;
    fluid_properties.vapor_pressure = 1800.F;

    EXPECT_NEAR(fluid_properties.density, 1.1F, 1e-6F);
    EXPECT_NEAR(fluid_properties.kinematic_viscosity, 1.5E-5F, 1e-6F);
    EXPECT_NEAR(fluid_properties.sound_speed, 340.F, 1e-6F);
    EXPECT_NEAR(fluid_properties.vapor_pressure, 1800.F, 1e-6F);
}

TEST(AerodynInflowTest, EnvironmentalConditions_Default) {
    const util::EnvironmentalConditions environmental_conditions;
    EXPECT_NEAR(environmental_conditions.gravity, 9.80665F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.atm_pressure, 103500.F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.water_depth, 0.F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.msl_offset, 0.F, 1e-6F);
}

TEST(AerodynInflowTest, EnvironmentalConditions_Set) {
    util::EnvironmentalConditions environmental_conditions;
    environmental_conditions.gravity = 9.79665F;
    environmental_conditions.atm_pressure = 103000.F;
    environmental_conditions.water_depth = 100.F;
    environmental_conditions.msl_offset = 10.F;

    EXPECT_NEAR(environmental_conditions.gravity, 9.79665F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.atm_pressure, 103000.F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.water_depth, 100.F, 1e-6F);
    EXPECT_NEAR(environmental_conditions.msl_offset, 10.F, 1e-6F);
}

TEST(AerodynInflowTest, VTKSettings_Default) {
    const util::VTKSettings vtk_settings;
    EXPECT_EQ(vtk_settings.write_vtk, false);
    EXPECT_EQ(vtk_settings.vtk_type, 1);
    ExpectArrayNear(vtk_settings.vtk_nacelle_dimensions, {-2.5F, -2.5F, 0.F, 10.F, 5.F, 5.F});
    EXPECT_EQ(vtk_settings.vtk_hub_radius, 1.5F);
}

TEST(AerodynInflowTest, VTKSettings_Set) {
    util::VTKSettings vtk_settings;
    vtk_settings.write_vtk = true;
    vtk_settings.vtk_type = 2;
    vtk_settings.vtk_nacelle_dimensions = {-1.5F, -1.5F, 0.F, 5.F, 2.5F, 2.5F};
    vtk_settings.vtk_hub_radius = 1.F;

    EXPECT_EQ(vtk_settings.write_vtk, true);
    EXPECT_EQ(vtk_settings.vtk_type, 2);
    ExpectArrayNear(vtk_settings.vtk_nacelle_dimensions, {-1.5F, -1.5F, 0.F, 5.F, 2.5F, 2.5F});
    EXPECT_EQ(vtk_settings.vtk_hub_radius, 1.F);
}

/// Helper function to get the shared library path
std::string GetSharedLibraryPath() {
    return std::string(KYNEMA_FMB_ADI_LIBRARY);
}

TEST(AerodynInflowTest, AeroDynInflowLibrary_DefaultConstructor) {
    // Load the shared library
    const std::string path = GetSharedLibraryPath();
    const util::AeroDynInflowLibrary aerodyn_inflow_library(path);

    // Check initial error handling state
    EXPECT_EQ(aerodyn_inflow_library.GetErrorHandling().error_status, 0);
    EXPECT_STREQ(aerodyn_inflow_library.GetErrorHandling().error_message.data(), "");

    // Check default values for other important members
    EXPECT_EQ(aerodyn_inflow_library.GetFluidProperties().density, 1.225F);
    EXPECT_EQ(aerodyn_inflow_library.GetEnvironmentalConditions().gravity, 9.80665F);
    EXPECT_EQ(aerodyn_inflow_library.GetSimulationControls().debug_level, 0);
    EXPECT_EQ(aerodyn_inflow_library.GetSimulationControls().transpose_DCM, 1);
    EXPECT_EQ(aerodyn_inflow_library.GetVTKSettings().write_vtk, 0);
}

TEST(AerodynInflowTest, AeroDynInflowLibrary_FullLoopSimulation) {
    const std::filesystem::path project_root = FindProjectRoot();
    const std::filesystem::path input_path = project_root / "tests/regression_tests/external/";

    // Set up simulation parameters
    util::SimulationControls sim_controls;
    sim_controls.is_aerodyn_input_path = true;
    sim_controls.is_inflowwind_input_path = true;
    sim_controls.aerodyn_input = (input_path / "ad_primary.dat").string();
    sim_controls.inflowwind_input = (input_path / "ifw_primary.dat").string();
    sim_controls.time_step = 0.0125;
    sim_controls.max_time = 10.;
    sim_controls.interpolation_order = 2;
    sim_controls.store_HH_wind_speed = false;
    sim_controls.transpose_DCM = true;
    sim_controls.debug_level = static_cast<int>(util::SimulationControls::DebugLevel::kNone);

    // Set up fluid properties
    util::FluidProperties fluid_props;
    fluid_props.density = 1.225F;
    fluid_props.kinematic_viscosity = 1.464E-05F;
    fluid_props.sound_speed = 335.F;
    fluid_props.vapor_pressure = 1700.F;

    // Set up environmental conditions
    util::EnvironmentalConditions env_conditions;
    env_conditions.gravity = 9.80665F;
    env_conditions.atm_pressure = 103500.F;
    env_conditions.water_depth = 0.F;
    env_conditions.msl_offset = 0.F;

    // Set up VTK settings
    const util::VTKSettings vtk_settings{};

    // Set up turbine configuration
    const std::array<double, 7> hub_pos = {0., 0., 90., 1., 0., 0., 0.};
    const std::array<double, 7> nacelle_pos = {0., 0., 90., 1., 0., 0., 0.};
    std::vector<util::TurbineConfig::BladeInitialState> blade_states;
    for (auto i : std::views::iota(0, 3)) {
        const util::TurbineConfig::BladeInitialState blade_state(
            {0., 0., 90., 1., 0., 0., 0.},  // root_initial_position
            {
                {0., 5., 90., 1., 0., 0., 0.},   // node_initial_positions - 1
                {0., 10., 90., 1., 0., 0., 0.},  // node_initial_positions - 2
            }
        );
        blade_states.push_back(blade_state);
    }
    const util::TurbineConfig turbine_config(
        true,             // is_horizontal_axis
        {0.F, 0.F, 0.F},  // reference_position
        hub_pos,          // hub_initial_position
        nacelle_pos,      // nacelle_initial_position
        blade_states
    );

    // Load the shared library and initialize AeroDynInflowLibrary
    const std::string path = GetSharedLibraryPath();
    util::AeroDynInflowLibrary aerodyn_inflow_library(
        path, util::ErrorHandling{}, fluid_props, env_conditions, sim_controls, vtk_settings
    );

    // Initialize with turbine configuration
    const std::vector<util::TurbineConfig> turbine_configs = {turbine_config};
    EXPECT_NO_THROW(aerodyn_inflow_library.Initialize(turbine_configs));

    // Simulate for 10 time steps
    double current_time = 0.;
    double next_time = sim_controls.time_step;
    std::vector<float> output_channel_values;

    for (auto i : std::views::iota(0, 10)) {
        // Update motion data for each time step (if needed)
        EXPECT_NO_THROW(aerodyn_inflow_library.SetRotorMotion());

        EXPECT_NO_THROW(aerodyn_inflow_library.UpdateStates(current_time, next_time));
        EXPECT_NO_THROW(aerodyn_inflow_library.CalculateOutput(next_time, output_channel_values));
        EXPECT_NO_THROW(aerodyn_inflow_library.GetRotorAerodynamicLoads());

        // Assert loads on blade nodes - they don't change since we're not updating the motion
        auto expected_loads = std::vector<std::array<float, 6>>{
            {11132.2F, -1938.43F, 0.F, 44472.6F, 323391.F, 50444.8F},  // Blade 1, Node 1
            {0.F, 0.F, 0.F, 0.F, 0.F, 0.F},                            // Blade 1, Node 2
            {11132.2F, -1938.43F, 0.F, 44472.6F, 323391.F, 50444.8F},  // Blade 2, Node 1
            {0.F, 0.F, 0.F, 0.F, 0.F, 0.F},                            // Blade 2, Node 2
            {11132.2F, -1938.43F, 0.F, 44472.6F, 323391.F, 50444.8F},  // Blade 3, Node 1
            {0.F, 0.F, 0.F, 0.F, 0.F, 0.F}                             // Blade 3, Node 2
        };

        const auto& turbine = aerodyn_inflow_library.GetTurbines()[0];
        for (auto ii : std::views::iota(0U, turbine.NumberOfBlades())) {
            for (auto jj : std::views::iota(0U, turbine.node_indices_by_blade[ii].size())) {
                auto node_index = turbine.node_indices_by_blade[ii][jj];
                const auto& node_loads = turbine.blade_nodes.loads[node_index];
                const auto& e_loads = expected_loads[ii * 2 + jj];
                for (auto k : std::views::iota(0, 3)) {
                    EXPECT_NEAR(node_loads[k], e_loads[k], 1e-1F);
                }
            }
        }
        current_time = next_time;
        next_time += sim_controls.time_step;
    }

    // End simulation
    EXPECT_NO_THROW(aerodyn_inflow_library.Finalize());
}

}  // namespace kynema_fmb::tests
