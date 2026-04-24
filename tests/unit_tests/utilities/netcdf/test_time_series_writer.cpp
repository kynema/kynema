#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "utilities/netcdf/netcdf_file.hpp"
#include "utilities/netcdf/time_series_writer.hpp"

namespace kynema_fmb::tests {

class TimeSeriesWriterTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a unique filename using the test info to avoid race condition
        const testing::TestInfo* const test_info =
            testing::UnitTest::GetInstance()->current_test_info();
        test_file =
            std::string(test_info->test_case_name()) + "_" + std::string(test_info->name()) + ".nc";
        std::filesystem::remove(test_file);
    }

    void TearDown() override { std::filesystem::remove(test_file); }

    std::string test_file;
};

TEST_F(TimeSeriesWriterTest, ConstructorCreatesTimeDimension) {
    const util::TimeSeriesWriter writer(test_file);
    const auto& file = writer.GetFile();

    EXPECT_NO_THROW({
        const int time_dim_id = file.GetDimensionId("time");
        EXPECT_GE(time_dim_id, 0);
    });
}

TEST_F(TimeSeriesWriterTest, WriteValuesCreatesVariableAndDimension) {
    util::TimeSeriesWriter writer(test_file);
    const std::vector<double> values = {100., 200., 300.};

    EXPECT_NO_THROW({ writer.WriteValuesAtTimestep("hub_height_wind_speed", 0, values); });

    const auto& file = writer.GetFile();
    EXPECT_GE(file.GetVariableId("hub_height_wind_speed"), 0);
    EXPECT_GE(file.GetDimensionId("hub_height_wind_speed_dimension"), 0);

    std::vector<double> read_data(values.size());
    file.ReadVariable("hub_height_wind_speed", read_data.data());
    EXPECT_EQ(read_data, values);
}

TEST_F(TimeSeriesWriterTest, WriteValuesSavesCorrectData) {
    util::TimeSeriesWriter writer(test_file);
    const std::vector<double> values1 = {100., 200., 300.};  // dimensions = 3
    const std::vector<double> values2 = {400., 500., 600.};  // dimensions = 3

    EXPECT_NO_THROW({
        writer.WriteValuesAtTimestep("rotor_torque", 0, values1);
        writer.WriteValuesAtTimestep("rotor_torque", 1, values2);
    });

    const auto& file = writer.GetFile();
    std::vector<double> read_data(values1.size());

    const std::vector<size_t> start_1 = {0, 0};
    const std::vector<size_t> count = {1, values1.size()};
    file.ReadVariableAt("rotor_torque", start_1, count, read_data.data());
    EXPECT_EQ(read_data, values1);

    const std::vector<size_t> start_2 = {1, 0};
    file.ReadVariableAt("rotor_torque", start_2, count, read_data.data());
    EXPECT_EQ(read_data, values2);
}

TEST_F(TimeSeriesWriterTest, WriteValueCreatesVariableAndWritesSingleValue) {
    util::TimeSeriesWriter writer(test_file);

    EXPECT_NO_THROW({ writer.WriteValueAtTimestep("rotor_power", 0, 100.); });

    const auto& file = writer.GetFile();
    EXPECT_GE(file.GetVariableId("rotor_power"), 0);
    EXPECT_GE(file.GetDimensionId("rotor_power_dimension"), 0);

    std::vector<double> read_data(1);
    const std::vector<size_t> start = {0, 0};
    const std::vector<size_t> count = {1, 1};
    file.ReadVariableAt("rotor_power", start, count, read_data.data());
    EXPECT_EQ(read_data[0], 100.);
}

}  // namespace kynema_fmb::tests
