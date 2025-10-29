#include <array>
#include <initializer_list>
#include <vector>

#include <gtest/gtest.h>

#include "model/model.hpp"
#include "solver/solver.hpp"
#include "step/predict_next_state.hpp"
#include "step/reset_constraints.hpp"
#include "step/step.hpp"
#include "step/step_parameters.hpp"
#include "step/update_constraint_variables.hpp"
#include "test_utilities.hpp"

namespace kynema::tests {

inline void SetUpSolverAndAssemble() {
    // Mass matrix for uniform composite beam section
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},   std::array{0., 8.538e-2, 0., 0., 0., 0.},
        std::array{0., 0., 8.538e-2, 0., 0., 0.},   std::array{0., 0., 0., 1.4433e-2, 0., 0.},
        std::array{0., 0., 0., 0., 0.40972e-2, 0.}, std::array{0., 0., 0., 0., 0., 1.0336e-2},
    };

    // Stiffness matrix for uniform composite beam section
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
        std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
        std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
    };

    // Create model for adding nodes and constraints
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    const auto x0_root = Eigen::Matrix<double, 3, 1>(2., 0., 0.);

    // Node locations (GLL quadrature)
    constexpr auto node_s = std::array{
        0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242, 1.
    };

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    constexpr auto omega = 0.1;
    std::vector<size_t> beam_node_ids;
    std::ranges::transform(node_s, std::back_inserter(beam_node_ids), [&](auto s) {
        const auto x = 10 * s + x0_root[0];
        return model.AddNode()
            .SetElemLocation(s)
            .SetPosition(x, 0., 0., 1., 0., 0., 0.)
            .SetVelocity(0., x * omega, 0., 0., 0., omega)
            .Build();
    });

    // Add beam element
    model.AddBeamElement(
        beam_node_ids,
        std::array{
            BeamSection(0., mass_matrix, stiffness_matrix),
            BeamSection(1., mass_matrix, stiffness_matrix),
        },
        std::array{
            std::array{-0.9491079123427585, 0.1294849661688697},
            std::array{-0.7415311855993943, 0.27970539148927664},
            std::array{-0.40584515137739696, 0.3818300505051189},
            std::array{6.123233995736766e-17, 0.4179591836734694},
            std::array{0.4058451513773971, 0.3818300505051189},
            std::array{0.7415311855993945, 0.27970539148927664},
            std::array{0.9491079123427585, 0.1294849661688697},
        }
    );

    // Add prescribed BC on beam root node
    model.AddPrescribedBC(beam_node_ids[0]);

    // Solution parameters
    constexpr auto max_iter = 10U;
    constexpr auto is_dynamic_solve = true;
    constexpr auto step_size = 0.01;  // seconds
    constexpr auto rho_inf = 0.9;
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    const auto u_rot = Eigen::Quaternion<double>(
        Eigen::AngleAxis(omega * step_size, Eigen::Matrix<double, 3, 1>::Unit(2))
    );
    const auto x_root = u_rot._transformVector(x0_root);
    const auto u_trans = x_root - x0_root;
    const auto displacement =
        std::array{u_trans(0), u_trans(1), u_trans(2), u_rot.w(), u_rot.x(), u_rot.y(), u_rot.z()};
    constraints.UpdateDisplacement(0, displacement);

    // Predict the next state for the solver
    step::PredictNextState(parameters, state);
    step::ResetConstraints(constraints);
    step::ResetSolver(solver);

    // Update beam elements state from solvers
    step::UpdateSystemVariables(parameters, elements, state);
    step::AssembleSystemMatrix(parameters, solver, elements);
    step::AssembleSystemResidual(solver, elements, state);

    step::UpdateConstraintVariables(state, constraints);
    step::AssembleConstraintsMatrix(solver, constraints);
    step::AssembleConstraintsResidual(solver, constraints);

    expect_kokkos_view_2D_equal(constraints.lambda, {{0., 0., 0., 0., 0., 0.}});
    expect_kokkos_view_2D_equal(
        state.q_prev,
        {
            {0., 0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0., 0.},
            {0., 0., 0., 1., 0., 0., 0.},
        }
    );
    expect_kokkos_view_2D_equal(
        state.a,
        {
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
        }
    );
    expect_kokkos_view_2D_equal(
        state.v,
        {
            {0, 0.20000000000000000, 0, 0, 0, 0.1},
            {0, 0.31747233803526764, 0, 0, 0, 0.1},  // NOLINT
            {0, 0.55738424175967749, 0, 0, 0, 0.1},
            {0, 0.84261575824032242, 0, 0, 0, 0.1},
            {0, 1.08252766196473240, 0, 0, 0, 0.1},
            {0, 1.20000000000000000, 0, 0, 0, 0.1},
        }
    );
    expect_kokkos_view_2D_equal(
        state.vd,
        {
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
            {0., 0., 0., 0., 0., 0.},
        }
    );
    expect_kokkos_view_2D_equal(
        state.q_delta,
        {
            {0., 0.20000000000000001, 0., 0., 0., 0.10000000000000001},
            {0., 0.31747233803526764, 0., 0., 0., 0.10000000000000001},  // NOLINT
            {0., 0.55738424175967749, 0., 0., 0., 0.10000000000000001},
            {0., 0.84261575824032242, 0., 0., 0., 0.10000000000000001},
            {0., 1.08252766196473240, 0., 0., 0., 0.10000000000000001},
            {0., 1.20000000000000020, 0., 0., 0., 0.10000000000000001},
        }
    );
    expect_kokkos_view_2D_equal(
        state.q,
        {
            {0, 0.0020000000000000, 0, 0.999999875, 0, 0, 0.0004999999},
            {0, 0.0031747233803526, 0, 0.999999875, 0, 0, 0.0004999999},
            {0, 0.0055738424175967, 0, 0.999999875, 0, 0, 0.0004999999},
            {0, 0.0084261575824032, 0, 0.999999875, 0, 0, 0.0004999999},
            {0, 0.0108252766196473, 0, 0.999999875, 0, 0, 0.0004999999},
            {0, 0.0120000000000000, 0, 0.999999875, 0, 0, 0.0004999999},
        }
    );
    expect_kokkos_view_2D_equal(
        constraints.residual_terms, {{
                                        9.9999991642896191E-7,
                                        3.3333331667800836E-10,
                                        0.,
                                        0.,
                                        0.,
                                        0.,
                                    }}
    );
    expect_kokkos_view_1D_equal(
        Kokkos::subview(solver.b, Kokkos::ALL, 0),
        {
            -0.68408451644565105,
            -0.00065456473269251652,
            0,
            1.4663215017154278E-17,
            1.5487939846687054E-17,
            0.0000098399278916272914,
            -4.9960036108132044E-16,
            7.2063733561229804E-15,
            0,
            -9.2033949200707357E-18,
            -9.7210350349908583E-18,
            0.000055862494393220151,
            -8.1878948066105295E-16,
            -3.4811400723838704E-14,
            0,
            -4.27839293851558E-18,
            -4.5190289029179014E-18,
            0.000081896496420578814,
            1.7208456881689926E-15,
            1.0463827612543219E-13,
            0,
            1.991039683802277E-17,
            2.1030246653970917E-17,
            0.000081896496635873786,
            -5.1070259132757201E-15,
            -3.455288355783126E-13,
            0,
            -7.9977507519330913E-17,
            -8.4475800436550718E-17,
            0.000055862494238845577,
            0.68408451644565516,
            0.00065456473296101171,
            0,
            5.8885683522740175E-17,
            6.2197677873801522E-17,
            0.0000098399277591661938,
            9.9999991642896191E-7,
            3.3333331667800836E-10,
            0,
            0,
            0,
            0,
        }
    );
}

TEST(NewSolverTest, SolverPredictNextState) {
    SetUpSolverAndAssemble();
}

inline void SetupAndTakeNoSteps() {
    // Mass matrix for uniform composite beam section
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},   std::array{0., 8.538e-2, 0., 0., 0., 0.},
        std::array{0., 0., 8.538e-2, 0., 0., 0.},   std::array{0., 0., 0., 1.4433e-2, 0., 0.},
        std::array{0., 0., 0., 0., 0.40972e-2, 0.}, std::array{0., 0., 0., 0., 0., 1.0336e-2},
    };

    // Stiffness matrix for uniform composite beam section
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
        std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
        std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
    };

    // Create model for adding nodes and constraints
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Node locations (GLL quadrature)
    constexpr auto node_s = std::array{
        0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242, 1.
    };

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    constexpr auto omega = 0.1;
    const auto x0_root = Eigen::Matrix<double, 3, 1>(2., 0., 0.);
    std::vector<size_t> beam_node_ids;
    std::ranges::transform(node_s, std::back_inserter(beam_node_ids), [&](auto s) {
        const auto x = 10 * s + x0_root(0);
        return model.AddNode()
            .SetElemLocation(s)
            .SetPosition(x, 0., 0., 1., 0., 0., 0.)
            .SetVelocity(0., x * omega, 0., 0., 0., omega)
            .Build();
    });

    // Add beam element
    model.AddBeamElement(
        beam_node_ids,
        std::array{
            BeamSection(0., mass_matrix, stiffness_matrix),
            BeamSection(1., mass_matrix, stiffness_matrix),
        },
        std::array{
            std::array{-0.9491079123427585, 0.1294849661688697},
            std::array{-0.7415311855993943, 0.27970539148927664},
            std::array{-0.40584515137739696, 0.3818300505051189},
            std::array{6.123233995736766e-17, 0.4179591836734694},
            std::array{0.4058451513773971, 0.3818300505051189},
            std::array{0.7415311855993945, 0.27970539148927664},
            std::array{0.9491079123427585, 0.1294849661688697},
        }
    );

    // Add prescribed boundary condition to first node
    model.AddPrescribedBC(beam_node_ids[0]);

    // Solution parameters
    constexpr auto max_iter = 1U;
    constexpr auto is_dynamic_solve = true;
    constexpr auto step_size = 0.01;  // seconds
    constexpr auto rho_inf = 0.9;
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    const auto u_rot = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(omega * step_size, Eigen::Matrix<double, 3, 1>::Unit(2))
    );
    const auto x_root = u_rot._transformVector(x0_root);
    const auto u_trans = x_root - x0_root;
    auto displacement =
        std::array{u_trans(0), u_trans(1), u_trans(2), u_rot.w(), u_rot.x(), u_rot.y(), u_rot.z()};
    constraints.UpdateDisplacement(0, displacement);

    Step(parameters, solver, elements, state, constraints);

    const auto x = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), solver.x);
    expect_kokkos_view_1D_equal(
        Kokkos::subview(solver.x, Kokkos::ALL, 0),
        {
            -9.9999991642894645E-7,    -3.3333331667800779E-10, -2.7693137528041648E-28,
            -6.5954208494030506E-30,   1.3190840599775882E-26,  7.7898219222917266E-24,
            -0.0000014964806473087194, -7.295038225205836E-10,  -4.5434716051775647E-14,
            -3.1920684921905924E-12,   -4.0723290197473225E-14, -2.2269841387499837E-10,
            -0.0000025281897479172782, -1.9662810249790896E-9,  -3.0976187800764434E-14,
            -4.8228594588737896E-12,   -9.9091066510552565E-14, -4.1568756453239842E-10,
            -0.0000038030054024772931, -3.5935658864740807E-9,  6.4411049537235008E-14,
            -3.6024743778286134E-12,   -1.1579932606573532E-13, -4.3409555957226201E-10,
            -0.000004933638516780269,  -4.9475925730408304E-9,  2.5185545987603492E-13,
            -2.7436857930422583E-12,   -1.6107798285714521E-13, -4.2116615327927008E-10,
            -0.0000055119072875186336, -5.625802843908016E-9,   4.2749628020387186E-13,
            -2.6313864405287534E-12,   -1.7915982090983506E-13, -4.1964251580055361E-10,
            0.10816660597819647,       0.000095455157310304377, -2.1220160418770112E-9,
            2.0316113427113777E-8,     2.1788291524533811E-8,   -0.000033726153743119725,
        }
    );
    EXPECT_NEAR(solver.convergence_err[0], 1669.0477021697936, 1.e-6);
    expect_kokkos_view_2D_equal(
        state.q_delta,
        {
            {-0.000099999991642894646, 0.19999996666666833, -2.7693137528041646E-26,
             -6.5954208494030504E-28, 1.3190840599775882E-24, 0.10000000000000001},
            {-0.00014964806473087192, 0.31747226508488541, -4.5434716051775645E-12,  // NOLINT
             -3.1920684921905925E-10, -4.072329019747322E-12, 0.099999977730158617},
            {-0.00025281897479172784, 0.55738404513157502, -3.0976187800764433E-12,
             -4.8228594588737892E-10, -9.9091066510552555E-12, 0.099999958431243555},
            {-0.00038030054024772928, 0.84261539888373382, 6.4411049537235007E-12,
             -3.6024743778286132E-10, -1.1579932606573532E-11, 0.099999956590444042},
            {-0.00049336385167802685, 1.082527167205475, 2.5185545987603493E-11,
             -2.743685793042258E-10, -1.6107798285714522E-11, 0.099999957883384671},
            {-0.00055119072875186334, 1.1999994374197158, 4.2749628020387184E-11,
             -2.6313864405287535E-10, -1.7915982090983507E-11, 0.099999958035748431},
        }
    );
    expect_kokkos_view_2D_equal(
        state.v,
        {
            {-0.00019949998332757481, 0.19999993350000334, -5.5247809368443091E-26,
             -1.3157864594559085E-27, 2.6315726996552886E-24, 0.10000000000000001},
            {-0.00029854788913808952, 0.31747219249925507, -9.0642258523292409E-12,  // NOLINT
             -6.3681766419202313E-10, -8.1242963943959083E-12, 0.099999955571666437},
            {-0.00050437385470949699, 0.55738384948661301, -6.1797494662525043E-12,
             -9.62160462045321E-10, -1.9768667768855238E-11, 0.099999917070330888},
            {-0.00075869957779421996, 0.84261504132392806, 1.2850004382678384E-11,
             -7.186936383768084E-10, -2.3101965550114195E-11, 0.099999913397935874},
            {-0.00098426088409766361, 1.082526674920014, 5.0245164245268966E-11,
             -5.4736531571193053E-10, -3.2135057580000468E-11, 0.09999991597735243},
            {-0.0010996255038599673, 1.1999988776523327, 8.5285507900672438E-11,
             -5.2496159488548634E-10, -3.5742384271512094E-11, 0.099999916281318101},
        }
    );
    expect_kokkos_view_2D_equal(
        state.vd,
        {
            {-0.039709996681393474, -0.000013236666005283692, -1.0996944912385341E-23,
             -2.619041619297952E-25, 5.2380828021710043E-22, 3.0933382853420454E-19},
            {-0.059425246504629256, -0.000028968596792292381, -1.8042125744160112E-9,
             -1.2675703982488845E-7, -1.617121853741662E-9, -0.0000088433540149761862},
            {-0.10039441488979514, -0.000078081019501919665, -1.2300644175683558E-9,
             -1.9151574911187822E-7, -3.9349062511340432E-9, -0.000016506953187581543},
            {-0.15101734453237334, -0.00014270050135188576, 2.5577627771236027E-9,
             -1.4305425754357427E-7, -4.5983912380703506E-9, -0.000017237934670614527},
            {-0.19591478550134453, -0.00019646890107545141, 1.0001180311677349E-8,
             -1.0895176284170809E-7, -6.3964066992572377E-9, -0.000016724507946719819},
            {-0.21887783838736499, -0.00022340063093158737, 1.6975877286895753E-8,
             -1.0449235555339682E-7, -7.1144364883295513E-9, -0.000016664004302439988},
        }
    );
    expect_kokkos_view_2D_equal(
        state.q,
        {
            {-9.9999991642894645E-7, 0.0019999996666666834, -2.7693137528041648E-28,
             0.99999987500000264, -3.2977102872969255E-30, 6.5954200250787654E-27,
             0.0004999999791666669},
            {-0.0000014964806473087194, 0.003174722650848854, -4.5434716051775647E-14,
             0.99999987500005827, -1.5960341795939002E-12, -2.0361644250335121E-14,
             0.00049999986781747402},
            {-0.0000025281897479172786, 0.0055738404513157504, -3.0976187800764434E-14,
             0.99999987500010656, -2.4114296289607408E-12, -4.9545531190880803E-14,
             0.0004999997713229107},
            {-0.0000038030054024772931, 0.0084261539888373388, 6.4411049537235008E-14,
             0.99999987500011112, -1.8012371138628233E-12, -5.7899660620383827E-14,
             0.0004999997621189143},
            {-0.000004933638516780269, 0.010825271672054751, 2.5185545987603492E-13,
             0.99999987500010789, -1.3718428393610572E-12, -8.0538988072784165E-14,
             0.00049999976858361656},
            {-0.0000055119072875186336, 0.011999994374197159, 4.2749628020387186E-13,
             0.99999987500010756, -1.3156931654438724E-12, -8.9579906722424441E-14,
             0.00049999976934543534},
        }
    );
    expect_kokkos_view_2D_equal(
        constraints.lambda, {{
                                0.10816660597819647,
                                0.000095455157310304377,
                                2.1220160418770112E-9,
                                -2.0316113427113777E-8,
                                -2.1788291524533811E-8,
                                -0.000033726153743119725,
                            }}
    );
}

TEST(SolverStep1Test, SolutionVector) {
    SetupAndTakeNoSteps();
}

inline auto SetupAndTakeTwoSteps() {
    // Mass matrix for uniform composite beam section
    constexpr auto mass_matrix = std::array{
        std::array{8.538e-2, 0., 0., 0., 0., 0.},   std::array{0., 8.538e-2, 0., 0., 0., 0.},
        std::array{0., 0., 8.538e-2, 0., 0., 0.},   std::array{0., 0., 0., 1.4433e-2, 0., 0.},
        std::array{0., 0., 0., 0., 0.40972e-2, 0.}, std::array{0., 0., 0., 0., 0., 1.0336e-2},
    };

    // Stiffness matrix for uniform composite beam section
    constexpr auto stiffness_matrix = std::array{
        std::array{1368.17e3, 0., 0., 0., 0., 0.},
        std::array{0., 88.56e3, 0., 0., 0., 0.},
        std::array{0., 0., 38.78e3, 0., 0., 0.},
        std::array{0., 0., 0., 16.9600e3, 17.6100e3, -0.3510e3},
        std::array{0., 0., 0., 17.6100e3, 59.1200e3, -0.3700e3},
        std::array{0., 0., 0., -0.3510e3, -0.3700e3, 141.470e3},
    };

    // Create model for adding nodes and constraints
    auto model = Model();

    // Gravity vector
    model.SetGravity(0., 0., 0.);

    // Node locations (GLL quadrature)
    constexpr auto node_s = std::array{
        0., 0.11747233803526763, 0.35738424175967748, 0.64261575824032247, 0.88252766196473242, 1.
    };

    // Build vector of nodes (straight along x axis, no rotation)
    // Calculate displacement, velocity, acceleration assuming a
    // 0.1 rad/s angular velocity around the z axis
    constexpr auto omega = 0.1;
    std::vector<size_t> beam_node_ids;
    std::ranges::transform(node_s, std::back_inserter(beam_node_ids), [&](auto s) {
        const auto x = 10 * s;
        return model.AddNode()
            .SetElemLocation(s)
            .SetPosition(x, 0., 0., 1., 0., 0., 0.)
            .SetVelocity(0., x * omega, 0., 0., 0., omega)
            .Build();
    });

    // Add beam element
    model.AddBeamElement(
        beam_node_ids,
        std::array{
            BeamSection(0., mass_matrix, stiffness_matrix),
            BeamSection(1., mass_matrix, stiffness_matrix),
        },
        std::array{
            std::array{-0.9491079123427585, 0.1294849661688697},
            std::array{-0.7415311855993943, 0.27970539148927664},
            std::array{-0.40584515137739696, 0.3818300505051189},
            std::array{6.123233995736766e-17, 0.4179591836734694},
            std::array{0.4058451513773971, 0.3818300505051189},
            std::array{0.7415311855993945, 0.27970539148927664},
            std::array{0.9491079123427585, 0.1294849661688697},
        }
    );

    // Constraint inputs
    model.AddPrescribedBC(beam_node_ids[0]);

    // Solution parameters
    constexpr auto max_iter = 2U;
    constexpr auto is_dynamic_solve = true;
    constexpr auto step_size = 0.01;  // seconds
    constexpr auto rho_inf = 0.9;
    auto parameters = StepParameters(is_dynamic_solve, max_iter, step_size, rho_inf);

    // Create solver, elements, constraints, and state
    auto [state, elements, constraints] = model.CreateSystem();
    auto solver = CreateSolver<>(state, elements, constraints);

    const auto q = Eigen::Quaternion<double>(
        Eigen::AngleAxis<double>(omega * step_size, Eigen::Matrix<double, 3, 1>::Unit(2))
    );
    auto displacement = std::array{0., 0., 0., q.w(), q.x(), q.y(), q.z()};
    constraints.UpdateDisplacement(0, displacement);

    Step(parameters, solver, elements, state, constraints);

    expect_kokkos_view_2D_equal(
        constraints.residual_terms, {{
                                        0.0,
                                        0.0,
                                        -2.769313752804165e-28,
                                        -1.31908395004359e-29,
                                        1.3190835103592412e-26,
                                        0.0,
                                    }}
    );
    expect_kokkos_view_1D_equal(
        Kokkos::subview(solver.b, Kokkos::ALL, 0),
        {
            -1.6976609407260757e-15,
            -5.288375554255456E-9,
            1.1889571903333178E-13,
            -1.1240653642481715E-12,
            -1.2020880297284837E-12,
            1.8684840013084998E-9,
            -3.6875789625611664E-15,
            -4.1819420587222298E-15,
            3.7352345516767378E-15,
            -1.1877649082466982E-15,
            -1.0506791059026386E-14,
            -5.8939891699170042E-15,
            -2.8792767262448265E-16,
            -1.6256730932454705E-15,
            -3.7325943099788989E-16,
            -5.3808016138201058E-16,
            -1.6282824212620447E-14,
            -1.1016714995625436E-14,
            5.2084529364118357E-16,
            -1.1638779207501438E-16,
            -1.4881053281336617E-15,
            9.5313284912395902E-17,
            -1.0413534631705933E-14,
            -1.3446531440415119E-14,
            2.7455938415381741E-15,
            -2.5729773964724179E-16,
            -3.7154086685782557E-16,
            1.5306103158449729E-16,
            -5.1065389766195233E-15,
            -9.9765062644864507E-15,
            -9.9677609771866959E-16,
            7.3693341584230349E-15,
            -2.8348323425684356E-15,
            3.4513177152054705E-17,
            -8.3796234226909556E-16,
            -1.8000463717889823E-15,
            0,
            0,
            -2.7693137528041648E-28,
            -1.3190839500435901E-29,
            1.3190835103592412E-26,
            0,
        }
    );
    const auto x = Kokkos::create_mirror_view(solver.x);
    Kokkos::deep_copy(x, solver.x);
    expect_kokkos_view_1D_equal(
        Kokkos::subview(solver.x, Kokkos::ALL, 0),
        {
            2.0701786547467534E-22,  -9.7804349489591981E-26, 2.7272832797488856E-28,
            2.0547315819961114E-29,  -1.4712947016254835E-26, 5.4331648183819847E-25,
            3.8348324508389868E-17,  9.4969920440321861E-16,  -1.3580769206195316E-17,
            -5.2335075093141956E-18, 3.1855954047525528E-15,  9.2846344878465848E-16,
            -6.0860076883771766E-17, 8.028169690742347E-16,   -1.479979763477676E-17,
            -2.7808660671715477E-17, 4.8007090441525481E-15,  8.6051705604834966E-16,
            -2.179748742218857E-16,  -2.7745453151127729E-15, 1.1652103769994212E-17,
            -3.0548898515545171E-17, 3.569536108741796E-15,   -4.1743489803882303E-17,
            -3.0721911895997357E-16, -9.1411516388420519E-15, 7.8855513035796911E-17,
            -6.4825821213015845E-17, 2.6903185401857071E-15,  -1.1004628463528681E-15,
            -2.8438123968141702E-16, -1.3829871537755664E-14, 1.4293886089416713E-16,
            -7.9353195214248362E-17, 2.5707926448736107E-15,  -1.3780235858853489E-15,
            -1.179058346509202e-11,  8.6709875777896587e-11,  -4.2445828198259645E-9,
            4.0632386057053241E-8,   4.3576751058673394E-8,   1.9904737297204642e-10,
        }
    );
}

TEST(SolverStep2Test, ConstraintResidualVector) {
    SetupAndTakeTwoSteps();
}

}  // namespace kynema::tests
