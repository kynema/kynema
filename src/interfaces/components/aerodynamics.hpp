#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <ranges>
#include <span>
#include <vector>

#include "aerodynamics_input.hpp"
#include "interfaces/host_state.hpp"
#include "math/interpolation.hpp"
#include "math/quaternion_operations.hpp"
#include "math/vector_operations.hpp"
#include "model/node.hpp"

namespace kynema::interfaces::components {

double CalculateAngleOfAttack(std::span<const double, 3> v_rel);

std::array<double, 6> CalculateAerodynamicLoad(
    std::span<double, 3> ref_axis_moment, std::span<const double, 3> v_inflow,
    std::span<const double, 3> v_motion, std::span<const double> aoa_polar,
    std::span<const double> cl_polar, std::span<const double> cd_polar,
    std::span<const double> cm_polar, double chord, double delta_s, double fluid_density,
    std::span<const double, 3> con_force, std::span<const double, 4> qqr
);

std::array<double, 3> CalculateConMotionVector(
    double ac_to_ref_axis_horizontal, double chord_to_ref_axis_vertical
);

std::vector<double> CalculateJacobianXi(std::span<const double> aero_node_xi);

std::vector<double> CalculateAeroNodeWidths(
    std::span<const double> jacobian_xi, std::span<const double> jacobian_integration_matrix,
    std::span<const double> node_x
);

class AerodynamicBody {
public:
    size_t id;

    std::vector<size_t> node_ids;
    std::vector<std::array<double, 7>> node_u;
    std::vector<std::array<double, 6>> node_v;
    std::vector<std::array<double, 6>> node_f;

    std::vector<std::array<double, 7>> xr_motion_map;
    std::vector<std::array<double, 7>> u_motion_map;
    std::vector<std::array<double, 6>> v_motion_map;
    std::vector<std::array<double, 4>> qqr_motion_map;
    std::vector<std::array<double, 3>> con_motion;
    std::vector<std::array<double, 3>> x_motion;
    std::vector<std::array<double, 3>> v_motion;

    std::vector<std::array<double, 3>> con_force;
    std::vector<std::array<double, 6>> loads;
    std::vector<std::array<double, 3>> ref_axis_moments;

    std::vector<double> jacobian_xi;
    std::vector<std::array<double, 3>> v_inflow;
    std::vector<std::array<double, 3>> v_rel;
    std::vector<double> twist;
    std::vector<double> chord;
    std::vector<double> delta_s;
    std::vector<size_t> polar_size;
    std::vector<std::vector<double>> aoa;
    std::vector<std::vector<double>> cl;
    std::vector<std::vector<double>> cd;
    std::vector<std::vector<double>> cm;

    std::vector<double> motion_interp;
    std::vector<double> shape_deriv_jac;

private:
    static std::vector<double> ExtractSectionXi(const AerodynamicBodyInput& input);

    static std::vector<double> ExtractBeamNodeXi(
        const AerodynamicBodyInput& input, std::span<const Node> nodes
    );

    static std::vector<double> ComputeMotionInterp(
        std::span<const double> section_xi, std::span<const double> beam_node_xi
    );

    static std::vector<std::array<double, 7>> ExtractNodeX(
        const AerodynamicBodyInput& input, std::span<const Node> nodes
    );

    static void InterpolateQuaternionFromNodesToSections(
        std::span<std::array<double, 7>> xr, std::span<const std::array<double, 7>> node_x,
        std::span<const double> interp
    );

    static std::vector<std::array<double, 7>> InterpolateNodePositionsToSections(
        const AerodynamicBodyInput& input, std::span<const std::array<double, 7>> node_x,
        std::span<const double> interp, std::span<const double> section_xi,
        std::span<const double> beam_node_xi
    );

    static std::vector<double> ComputeShapeDerivNode(
        std::span<const double> section_xi, std::span<const double> beam_node_xi
    );

    static void AddTwistToReferenceLocation(
        std::vector<std::array<double, 7>>& xr, std::span<const std::array<double, 7>> node_x,
        const AerodynamicBodyInput& input, std::span<const double> shape_deriv_node
    );

    static std::vector<std::array<double, 3>> ComputeConMotion(const AerodynamicBodyInput& input);

    static std::vector<double> ComputeShapeDerivJacobian(
        std::span<const double> jacobian_xi, std::span<const double> beam_node_xi
    );

    static std::vector<double> ComputeDeltaS(
        std::span<const std::array<double, 7>> node_x, std::span<const double> jacobian_xi,
        std::span<const double> shape_deriv_jac
    );

    template <typename T>
    static std::vector<std::vector<double>> ExtractPolar(size_t n_sections, T polar_extractor) {
        auto output = std::vector<std::vector<double>>(n_sections);
        for (auto section : std::views::iota(0U, n_sections)) {
            const auto& polar_data = polar_extractor(section);
            const auto n_polar_points = polar_data.size();
            output[section].resize(n_polar_points);
            for (auto polar : std::views::iota(0U, n_polar_points)) {
                output[section][polar] = polar_data[polar];
            }
        }
        return output;
    }

    static std::vector<std::array<double, 3>> InitializeConForce(
        std::span<const std::array<double, 3>> con_motion
    );

public:
    AerodynamicBody(const AerodynamicBodyInput& input, std::span<const Node> nodes);

    template <typename DeviceType>
    void CalculateMotion(const HostState<DeviceType>& state) {
        // Copy beam node displacements from state
        for (auto node : std::views::iota(0U, node_u.size())) {
            for (auto component : std::views::iota(0U, 7U)) {
                node_u[node][component] = state.q(node_ids[node], component);
            }
        }

        // Copy beam node velocities from state
        for (auto node : std::views::iota(0U, node_v.size())) {
            for (auto component : std::views::iota(0U, 6U)) {
                node_v[node][component] = state.v(node_ids[node], component);
            }
        }

        InterpolateQuaternionFromNodesToSections(u_motion_map, node_u, motion_interp);

        // Interpolate beam node velocities to aerodynamic sections on the reference axis
        for (auto i : std::views::iota(0U, v_motion_map.size())) {
            for (auto component : std::views::iota(0U, 6U)) {
                v_motion_map[i][component] = 0.;
            }
            for (auto j : std::views::iota(0U, node_v.size())) {
                const auto coeff = motion_interp[i * node_v.size() + j];
                for (auto component : std::views::iota(0U, 6U)) {
                    v_motion_map[i][component] += coeff * node_v[j][component];
                }
            }
        }

        // Calculate global rotation of each section
        for (auto i : std::views::iota(0U, qqr_motion_map.size())) {
            const auto xr = std::array{
                xr_motion_map[i][3], xr_motion_map[i][4], xr_motion_map[i][5], xr_motion_map[i][6]
            };
            const auto u = std::array{
                u_motion_map[i][3], u_motion_map[i][4], u_motion_map[i][5], u_motion_map[i][6]
            };
            const auto qqr = math::QuaternionCompose(u, xr);
            for (auto component = 0U; component < 4U; ++component) {
                qqr_motion_map[i][component] = qqr[component];
            }
        }

        // Calculate motion of aerodynamic centers in global coordinates
        for (auto i : std::views::iota(0U, x_motion.size())) {
            const auto qqr_con = math::RotateVectorByQuaternion(qqr_motion_map[i], con_motion[i]);

            for (auto component : std::views::iota(0U, 3U)) {
                x_motion[i][component] =
                    xr_motion_map[i][component] + u_motion_map[i][component] + qqr_con[component];
            }

            const auto omega =
                std::array{v_motion_map[i][3], v_motion_map[i][4], v_motion_map[i][5]};
            const auto omega_qqr_con = math::CrossProduct(omega, qqr_con);
            for (auto component : std::views::iota(0U, 3U)) {
                v_motion[i][component] = v_motion_map[i][component] + omega_qqr_con[component];
            }
        }

        auto node_x = std::vector<double>(node_u.size() * 3U);
        for (auto node : std::views::iota(0U, node_f.size())) {
            for (auto component : std::views::iota(0U, 3U)) {
                node_x[component * node_u.size() + node] = state.x(node_ids[node], component);
            }
        }

        delta_s = CalculateAeroNodeWidths(jacobian_xi, shape_deriv_jac, node_x);
    }

    void SetInflowFromVector(std::span<const std::array<double, 3>> inflow_velocity);

    template <typename T>
    void SetInflowFromFunction(const T& inflow_velocity_function) {
        for (auto node = 0U; node < v_inflow.size(); ++node) {
            const auto inflow_velocity = inflow_velocity_function(x_motion[node]);
            for (auto component = 0U; component < 3U; ++component) {
                v_inflow[node][component] = inflow_velocity[component];
            }
        }
    }

    void SetAerodynamicLoads(std::span<const std::array<double, 6>> aerodynamic_loads);

    void CalculateAerodynamicLoads(double fluid_density);

    void CalculateNodalLoads();

    template <typename DeviceType>
    void AddNodalLoadsToState(HostState<DeviceType>& state) {
        for (auto node : std::views::iota(0U, node_f.size())) {
            for (auto component : std::views::iota(0U, 6U)) {
                state.f(node_ids[node], component) += node_f[node][component];
            }
        }
    }
};

class Aerodynamics {
public:
    std::vector<AerodynamicBody> bodies;

    Aerodynamics(std::span<const AerodynamicBodyInput> inputs, std::span<const Node> nodes);

    template <typename DeviceType>
    void CalculateMotion(HostState<DeviceType>& state) {
        for (auto& body : bodies) {
            body.CalculateMotion(state);
        }
    }

    void SetInflowFromVector(
        std::span<const std::vector<std::array<double, 3>>> body_inflow_velocities
    );

    template <typename T>
    void SetInflowFromFunction(const T& body_inflow_velocity_function) {
        for (auto& body : bodies) {
            body.SetInflowFromFunction(body_inflow_velocity_function);
        }
    }

    void SetAerodynamicLoads(std::span<const std::vector<std::array<double, 6>>> body_aero_loads);

    void CalculateAerodynamicLoads(double fluid_density);

    void CalculateNodalLoads();

    template <typename DeviceType>
    void AddNodalLoadsToState(HostState<DeviceType>& state) {
        for (auto& body : bodies) {
            body.AddNodalLoadsToState(state);
        }
    }
};
}  // namespace kynema::interfaces::components
