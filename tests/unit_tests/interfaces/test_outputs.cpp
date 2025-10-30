#include <array>
#include <filesystem>
#include <vector>

#include <gtest/gtest.h>

#include "interfaces/host_state.hpp"
#include "interfaces/outputs.hpp"
#include "state/state.hpp"
#include "utilities/netcdf/netcdf_file.hpp"

namespace kynema::tests {

using kynema::interfaces::Outputs;

TEST(OutputsTest, ConstructWithAndWithoutTimeSeriesFile) {
    const ::testing::TestInfo* const test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string test_name = std::string(test_info->test_case_name()) + "_" + test_info->name();
    const std::string node_state_file = "node_states_outputs_" + test_name + ".nc";
    const std::string time_series_file = "time_series_outputs_" + test_name + ".nc";
    std::filesystem::remove(node_state_file);
    std::filesystem::remove(time_series_file);

    // With time series
    {
        Outputs outputs(
            node_state_file, /*num_nodes=*/1, time_series_file, {"u"}, /*buffer_size=*/0
        );
        ASSERT_NE(outputs.GetOutputWriter(), nullptr);
        ASSERT_NE(outputs.GetTimeSeriesWriter(), nullptr);
    }

    // Without time series
    {
        Outputs outputs(
            node_state_file, /*num_nodes=*/1, /*time_series_file=*/"", {"u"}, /*buffer_size=*/0
        );
        ASSERT_NE(outputs.GetOutputWriter(), nullptr);
        ASSERT_EQ(outputs.GetTimeSeriesWriter(), nullptr);
    }

    std::filesystem::remove(node_state_file);
    std::filesystem::remove(time_series_file);
}

TEST(OutputsTest, WriteTimeSeriesAndReopen) {
    const ::testing::TestInfo* const test_info =
        ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string test_name = std::string(test_info->test_case_name()) + "_" + test_info->name();
    const std::string node_state_file = "node_states_outputs_" + test_name + ".nc";
    const std::string time_series_file = "time_series_outputs_" + test_name + ".nc";
    std::filesystem::remove(node_state_file);
    std::filesystem::remove(time_series_file);

    Outputs outputs(
        node_state_file, /*num_nodes=*/1, time_series_file, /*enabled_state_prefixes=*/{"u"},
        /*buffer_size=*/0
    );

    // Write values at two timesteps, close/reopen in between
    outputs.WriteValueAtTimestep(0, "azimuth_angle", 12.34);
    outputs.Close();

    // Check files actually closed
    EXPECT_EQ(outputs.GetOutputWriter()->GetFile().GetNetCDFId(), -1);
    EXPECT_EQ(outputs.GetTimeSeriesWriter()->GetFile().GetNetCDFId(), -1);

    // Reopen the files
    outputs.Open();
    outputs.WriteValueAtTimestep(1, "azimuth_angle", 56.78);

    // Validate time-series file contents
    const util::NetCdfFile time_series_output(time_series_file, /*create=*/false);
    const std::vector<size_t> shape = time_series_output.GetShape("azimuth_angle");
    ASSERT_EQ(shape.size(), 2U);  // time, value_dim
    ASSERT_EQ(shape[1], 1U);      // scalar time series

    std::vector<double> read_values(1);
    {
        std::array<size_t, 2> start = {0, 0};
        const std::array<size_t, 2> count = {1, 1};
        time_series_output.ReadVariableAt(
            "azimuth_angle", std::span<const size_t>(start), std::span<const size_t>(count),
            read_values.data()
        );
        EXPECT_NEAR(read_values[0], 12.34, 1e-14);

        start[0] = 1;
        time_series_output.ReadVariableAt(
            "azimuth_angle", std::span<const size_t>(start), std::span<const size_t>(count),
            read_values.data()
        );
        EXPECT_NEAR(read_values[0], 56.78, 1e-14);
    }

    std::filesystem::remove(node_state_file);
    std::filesystem::remove(time_series_file);
}

}  // namespace kynema::tests
