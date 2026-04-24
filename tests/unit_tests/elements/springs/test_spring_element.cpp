#include <gtest/gtest.h>

#include "elements/springs/spring_element.hpp"
#include "model/model.hpp"

namespace kynema_fmb::tests {

TEST(SpringElement, CreateSpringElement_ZeroUndeformedLength) {
    Model model;

    // Create two nodes
    auto node1_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    auto node2_id = model.AddNode().SetPosition(1., 0., 0., 1., 0., 0., 0.).Build();

    constexpr double k = 100.;  // stiffness
    constexpr double l0 = 0.;   // undeformed length
    SpringElement spring(0U, {node1_id, node2_id}, k, l0);

    EXPECT_EQ(spring.stiffness, k);
    EXPECT_EQ(spring.undeformed_length, l0);
    EXPECT_EQ(spring.node_ids[0], 0U);
    EXPECT_EQ(spring.node_ids[1], 1U);
}

TEST(SpringElement, CreateSpringElement_NonZeroUndeformedLength) {
    Model model;

    // Create two nodes
    auto node1_id = model.AddNode().SetPosition(0., 0., 0., 1., 0., 0., 0.).Build();
    auto node2_id = model.AddNode().SetPosition(2., 2., 1., 1., 0., 0., 0.).Build();

    constexpr double k = 100.;  // stiffness
    constexpr double l0 = 3.;   // undeformed length
    SpringElement spring(0U, {node1_id, node2_id}, k, l0);

    EXPECT_EQ(spring.stiffness, k);
    EXPECT_EQ(spring.undeformed_length, l0);
    EXPECT_EQ(model.GetNode(spring.node_ids[0]).x0[0], 0.);
    EXPECT_EQ(model.GetNode(spring.node_ids[1]).x0[0], 2.);
}

}  // namespace kynema_fmb::tests
