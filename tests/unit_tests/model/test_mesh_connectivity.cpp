#include <array>
#include <cstddef>
#include <filesystem>
#include <vector>

#include <gtest/gtest.h>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/node/parse.h>

#include "model/mesh_connectivity.hpp"

namespace kynema_fmb::model::tests {

class MeshConnectivityTest : public ::testing::Test {
protected:
    MeshConnectivity mesh_connectivity;

    void SetUp() override {
        // Beam connectivity
        mesh_connectivity.AddBeamElementConnectivity(1, std::array{4UL, 5UL, 6UL, 7UL, 8UL});
        mesh_connectivity.AddBeamElementConnectivity(2, std::array{9UL, 10UL, 11UL, 12UL, 13UL});

        // Mass connectivity
        mesh_connectivity.AddMassElementConnectivity(1, 0);

        // Spring connectivity
        mesh_connectivity.AddSpringElementConnectivity(1, std::array{1UL, 2UL});

        // Constraint connectivity
        mesh_connectivity.AddConstraintConnectivity(1, std::array{0UL, 1UL});
        mesh_connectivity.AddConstraintConnectivity(2, std::array{2UL, 4UL});
        mesh_connectivity.AddConstraintConnectivity(3, std::array{8UL, 9UL});
    }

    void TearDown() override {
        if (std::filesystem::exists("test_connectivity.yaml")) {
            std::filesystem::remove("test_connectivity.yaml");
        }
    }
};

TEST_F(MeshConnectivityTest, BeamElementConnectivity) {
    // Retrieve existing beam element connectivity
    auto nodes = mesh_connectivity.GetBeamElementConnectivity(1);
    ASSERT_EQ(nodes.size(), 5);
    EXPECT_EQ(nodes[0], 4);
    EXPECT_EQ(nodes[1], 5);
    EXPECT_EQ(nodes[2], 6);
    EXPECT_EQ(nodes[3], 7);
    EXPECT_EQ(nodes[4], 8);

    // Add new beam element connectivity
    mesh_connectivity.AddBeamElementConnectivity(3, std::array{14UL, 15UL, 16UL});
    nodes = mesh_connectivity.GetBeamElementConnectivity(3);
    ASSERT_EQ(nodes.size(), 3);
    EXPECT_EQ(nodes[0], 14);
    EXPECT_EQ(nodes[1], 15);
    EXPECT_EQ(nodes[2], 16);

    // Overwrite existing beam element connectivity
    mesh_connectivity.AddBeamElementConnectivity(1, std::array{20UL, 21UL, 22UL});
    nodes = mesh_connectivity.GetBeamElementConnectivity(1);
    ASSERT_EQ(nodes.size(), 3);
    EXPECT_EQ(nodes[0], 20);
    EXPECT_EQ(nodes[1], 21);
    EXPECT_EQ(nodes[2], 22);
}

TEST_F(MeshConnectivityTest, MassElementConnectivity) {
    // Retrieve existing mass element connectivity
    auto nodes = mesh_connectivity.GetMassElementConnectivity(1);
    ASSERT_EQ(nodes.size(), 1);
    EXPECT_EQ(nodes[0], 0);

    // Add new mass element connectivity
    mesh_connectivity.AddMassElementConnectivity(2, 5);
    nodes = mesh_connectivity.GetMassElementConnectivity(2);
    ASSERT_EQ(nodes.size(), 1);
    EXPECT_EQ(nodes[0], 5);
}

TEST_F(MeshConnectivityTest, SpringElementConnectivity) {
    // Retrieve existing spring element connectivity
    auto nodes = mesh_connectivity.GetSpringElementConnectivity(1);
    ASSERT_EQ(nodes.size(), 2);
    EXPECT_EQ(nodes[0], 1);
    EXPECT_EQ(nodes[1], 2);

    // Add new spring element connectivity
    mesh_connectivity.AddSpringElementConnectivity(2, std::array{3UL, 4UL});
    nodes = mesh_connectivity.GetSpringElementConnectivity(2);
    ASSERT_EQ(nodes.size(), 2);
    EXPECT_EQ(nodes[0], 3);
    EXPECT_EQ(nodes[1], 4);
}

TEST_F(MeshConnectivityTest, ConstraintConnectivity) {
    // Retrieve existing constraint element connectivity
    auto nodes = mesh_connectivity.GetConstraintConnectivity(1);
    ASSERT_EQ(nodes.size(), 2);
    EXPECT_EQ(nodes[0], 0);
    EXPECT_EQ(nodes[1], 1);

    nodes = mesh_connectivity.GetConstraintConnectivity(3);
    ASSERT_EQ(nodes.size(), 2);
    EXPECT_EQ(nodes[0], 8);
    EXPECT_EQ(nodes[1], 9);
}

TEST_F(MeshConnectivityTest, ExportToYAML) {
    auto yaml_stream = std::ostringstream();
    mesh_connectivity.ExportToYAML(yaml_stream);
    YAML::Node root = YAML::Load(yaml_stream.str());

    // Mass elements
    ASSERT_TRUE(root["masses"]);
    ASSERT_TRUE(root["masses"]["1"]);
    auto mass_nodes = root["masses"]["1"].as<std::vector<size_t>>();
    ASSERT_EQ(mass_nodes.size(), 1);
    EXPECT_EQ(mass_nodes[0], 0);

    // Spring elements
    ASSERT_TRUE(root["springs"]);
    ASSERT_TRUE(root["springs"]["1"]);
    auto spring_nodes = root["springs"]["1"].as<std::vector<size_t>>();
    ASSERT_EQ(spring_nodes.size(), 2);
    EXPECT_EQ(spring_nodes[0], 1);
    EXPECT_EQ(spring_nodes[1], 2);

    // Beam elements
    ASSERT_TRUE(root["beams"]);
    ASSERT_TRUE(root["beams"]["1"]);
    auto beam_nodes = root["beams"]["1"].as<std::vector<size_t>>();
    ASSERT_EQ(beam_nodes.size(), 5);
    EXPECT_EQ(beam_nodes[0], 4);
    EXPECT_EQ(beam_nodes[4], 8);

    // Constraints
    ASSERT_TRUE(root["constraints"]);
    ASSERT_TRUE(root["constraints"]["3"]);
    auto constraint_nodes = root["constraints"]["3"].as<std::vector<size_t>>();
    ASSERT_EQ(constraint_nodes.size(), 2);
    EXPECT_EQ(constraint_nodes[0], 8);
    EXPECT_EQ(constraint_nodes[1], 9);
}

TEST_F(MeshConnectivityTest, ImportFromYAML) {
    auto yaml_stream = std::ostringstream();
    mesh_connectivity.ExportToYAML(yaml_stream);

    MeshConnectivity imported_mesh;
    imported_mesh.ImportFromYAML(YAML::Load(yaml_stream.str()));

    // Beam elements
    ASSERT_EQ(imported_mesh.GetBeamElementConnectivity(1).size(), 5);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(1)[0], 4);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(1)[1], 5);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(1)[2], 6);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(1)[3], 7);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(1)[4], 8);

    ASSERT_EQ(imported_mesh.GetBeamElementConnectivity(2).size(), 5);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(2)[0], 9);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(2)[1], 10);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(2)[2], 11);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(2)[3], 12);
    EXPECT_EQ(imported_mesh.GetBeamElementConnectivity(2)[4], 13);

    // Mass elements
    ASSERT_EQ(imported_mesh.GetMassElementConnectivity(1).size(), 1);
    EXPECT_EQ(imported_mesh.GetMassElementConnectivity(1)[0], 0);

    // Spring elements
    ASSERT_EQ(imported_mesh.GetSpringElementConnectivity(1).size(), 2);
    EXPECT_EQ(imported_mesh.GetSpringElementConnectivity(1)[0], 1);
    EXPECT_EQ(imported_mesh.GetSpringElementConnectivity(1)[1], 2);

    // Constraints
    ASSERT_EQ(imported_mesh.GetConstraintConnectivity(1).size(), 2);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(1)[0], 0);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(1)[1], 1);

    ASSERT_EQ(imported_mesh.GetConstraintConnectivity(2).size(), 2);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(2)[0], 2);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(2)[1], 4);

    ASSERT_EQ(imported_mesh.GetConstraintConnectivity(3).size(), 2);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(3)[0], 8);
    EXPECT_EQ(imported_mesh.GetConstraintConnectivity(3)[1], 9);
}

}  // namespace kynema_fmb::model::tests
