#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <netcdf.h>

#include "utilities/netcdf/netcdf_file.hpp"

namespace kynema::tests {

class NetCDFFileTest : public ::testing::Test {
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

TEST_F(NetCDFFileTest, CreateAndCloseFile) {
    EXPECT_NO_THROW({
        const util::NetCDFFile file(test_file);
        EXPECT_NE(file.GetNetCDFId(), -1);
    });
    EXPECT_TRUE(std::filesystem::exists(test_file));
}

TEST_F(NetCDFFileTest, AddDimension) {
    const util::NetCDFFile file(test_file);
    EXPECT_NO_THROW({
        const int dim_id = file.AddDimension("test_dim", 10);
        EXPECT_GE(dim_id, 0);
        EXPECT_EQ(file.GetDimensionId("test_dim"), dim_id);
    });
}

TEST_F(NetCDFFileTest, AddAndWriteVariableWithTypeDouble) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<double>("position", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<double> data = {1.1, 2.2, 3.3, 4.4, 5.5};
    EXPECT_NO_THROW(file.WriteVariable("position", data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableWithTypeFloat) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<float>("velocity", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<float> data = {1.F, 2.F, 3.F, 4.F, 5.F};
    EXPECT_NO_THROW(file.WriteVariable("velocity", data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableWithTypeInt) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<int>("count", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<int> data = {1, 2, 3, 4, 5};
    EXPECT_NO_THROW(file.WriteVariable("count", data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableWithTypeString) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<std::string>("labels", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<std::string> data = {"one", "two", "three", "four", "five"};
    EXPECT_NO_THROW(file.WriteVariable("labels", data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableAtWithTypeDouble) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<double>("position", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<double> data = {1.1, 2.2, 3.3};
    const std::vector<size_t> start = {2};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.WriteVariableAt("position", start, count, data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableAtWithTypeFloat) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<float>("velocity", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<float> data = {1.F, 2.F, 3.F};
    const std::vector<size_t> start = {4};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.WriteVariableAt("velocity", start, count, data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableAtWithTypeInt) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<int>("count", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<int> data = {1, 2, 3};
    const std::vector<size_t> start = {6};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.WriteVariableAt("count", start, count, data));
}

TEST_F(NetCDFFileTest, AddAndWriteVariableAtWithTypeString) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};

    const int var_id = file.AddVariable<std::string>("labels", dim_ids);
    EXPECT_GE(var_id, 0);

    const std::vector<std::string> data = {"one", "two", "three"};
    const std::vector<size_t> start = {0};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.WriteVariableAt("labels", start, count, data));
}

TEST_F(NetCDFFileTest, AddAttribute) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};
    EXPECT_NO_THROW({
        const int var_id = file.AddVariable<double>("position", dim_ids);
        EXPECT_GE(var_id, 0);
    });

    EXPECT_NO_THROW(file.AddAttribute("position", "units", std::string("meters")));
}

TEST_F(NetCDFFileTest, OpenExistingFile) {
    {
        const util::NetCDFFile file(test_file);
        const int time_dim = file.AddDimension("time", 5);
        const int position_var = file.AddVariable<double>("position", std::array{time_dim});
        EXPECT_GE(position_var, 0);
    }

    EXPECT_NO_THROW({
        const util::NetCDFFile file(test_file, false);
        EXPECT_NE(file.GetNetCDFId(), -1);
        EXPECT_NO_THROW({
            const int time_dim = file.GetDimensionId("time");
            EXPECT_GE(time_dim, 0);
        });
        EXPECT_NO_THROW({
            const int position_var = file.GetVariableId("position");
            EXPECT_GE(position_var, 0);
        });
    });
}

TEST_F(NetCDFFileTest, GetNumberOfDimensions) {
    const util::NetCDFFile file(test_file);

    const int time_dim = file.AddDimension("time", 5);
    const int space_dim = file.AddDimension("nodes", 3);
    const std::vector<int> dim_ids = {time_dim, space_dim};
    const int position_1D_var = file.AddVariable<double>("position_1D", std::array{time_dim});
    const int position_2D_var = file.AddVariable<double>("position_2D", dim_ids);
    EXPECT_GE(position_1D_var, 0);
    EXPECT_GE(position_2D_var, 0);

    EXPECT_EQ(file.GetNumberOfDimensions("position_1D"), 1);
    EXPECT_EQ(file.GetNumberOfDimensions("position_2D"), 2);
}

TEST_F(NetCDFFileTest, GetDimensionLength) {
    const util::NetCDFFile file(test_file);

    const int time_dim = file.AddDimension("time", 5);
    const int space_dim = file.AddDimension("nodes", 3);

    EXPECT_EQ(file.GetDimensionLength(time_dim), 5);
    EXPECT_EQ(file.GetDimensionLength(space_dim), 3);

    EXPECT_EQ(file.GetDimensionLength("time"), 5);
    EXPECT_EQ(file.GetDimensionLength("nodes"), 3);
}

TEST_F(NetCDFFileTest, GetShape) {
    const util::NetCDFFile file(test_file);

    const int time_dim = file.AddDimension("time", 5);
    const int space_dim = file.AddDimension("nodes", 3);

    const int position_1D_var = file.AddVariable<double>("position_1D", std::array{time_dim});
    EXPECT_GE(position_1D_var, 0);
    const std::vector<size_t> expected_shape_1D = {5};
    EXPECT_EQ(file.GetShape("position_1D"), expected_shape_1D);

    const std::vector<int> dim_ids_2D = {time_dim, space_dim};
    const int position_2D_var = file.AddVariable<double>("position_2D", dim_ids_2D);
    EXPECT_GE(position_2D_var, 0);
    const std::vector<size_t> expected_shape_2D = {5, 3};
    EXPECT_EQ(file.GetShape("position_2D"), expected_shape_2D);
}

TEST_F(NetCDFFileTest, ReadVariableWithTypeDouble) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int position_var = file.AddVariable<double>("position", dim_ids);
    EXPECT_GE(position_var, 0);
    const std::vector<double> write_data = {1.1, 2.2, 3.3, 4.4, 5.5};
    file.WriteVariable("position", write_data);

    std::vector<double> read_data(5);
    EXPECT_NO_THROW(file.ReadVariable("position", read_data.data()));
    EXPECT_EQ(read_data, write_data);
}

TEST_F(NetCDFFileTest, ReadVariableWithTypeFloat) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int velocity_var = file.AddVariable<float>("velocity", dim_ids);
    EXPECT_GE(velocity_var, 0);
    const std::vector<float> write_data = {1.F, 2.F, 3.F, 4.F, 5.F};
    file.WriteVariable("velocity", write_data);

    std::vector<float> read_data(5);
    EXPECT_NO_THROW(file.ReadVariable("velocity", read_data.data()));
    EXPECT_EQ(read_data, write_data);
}

TEST_F(NetCDFFileTest, ReadVariableWithTypeInt) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 5);
    const std::vector<int> dim_ids = {dim_id};

    const int count_var = file.AddVariable<int>("count", dim_ids);
    EXPECT_GE(count_var, 0);
    const std::vector<int> write_data = {1, 2, 3, 4, 5};
    file.WriteVariable("count", write_data);

    std::vector<int> read_data(5);
    EXPECT_NO_THROW(file.ReadVariable("count", read_data.data()));
    EXPECT_EQ(read_data, write_data);
}

TEST_F(NetCDFFileTest, ReadVariableAtWithTypeDouble) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int position_var = file.AddVariable<double>("position", dim_ids);
    EXPECT_GE(position_var, 0);

    // Write full data
    const std::vector<double> write_data = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.};
    file.WriteVariable("position", write_data);

    // Read partial data
    std::vector<double> read_data(3);
    const std::vector<size_t> start = {2};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.ReadVariableAt("position", start, count, read_data.data()));
    EXPECT_EQ(read_data, std::vector<double>({3.3, 4.4, 5.5}));
}

TEST_F(NetCDFFileTest, ReadVariableAtWithTypeFloat) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int velocity_var = file.AddVariable<float>("velocity", dim_ids);
    EXPECT_GE(velocity_var, 0);

    // Write full data
    const std::vector<float> write_data = {1.F, 2.F, 3.F, 4.F, 5.F, 6.F, 7.F, 8.F, 9.F, 10.F};
    file.WriteVariable("velocity", write_data);

    // Read partial data
    std::vector<float> read_data(3);
    const std::vector<size_t> start = {4};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.ReadVariableAt("velocity", start, count, read_data.data()));
    EXPECT_EQ(read_data, std::vector<float>({5.F, 6.F, 7.F}));
}

TEST_F(NetCDFFileTest, ReadVariableAtWithTypeInt) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int count_var = file.AddVariable<int>("count", dim_ids);
    EXPECT_GE(count_var, 0);

    // Write full data
    const std::vector<int> write_data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    file.WriteVariable("count", write_data);

    // Read partial data
    std::vector<int> read_data(3);
    const std::vector<size_t> start = {6};
    const std::vector<size_t> count = {3};
    EXPECT_NO_THROW(file.ReadVariableAt("count", start, count, read_data.data()));
    EXPECT_EQ(read_data, std::vector<int>({7, 8, 9}));
}

TEST_F(NetCDFFileTest, ReadVariableWithStrideTypeDouble) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int position_var = file.AddVariable<double>("position", dim_ids);
    EXPECT_GE(position_var, 0);

    // Write full data
    const std::vector<double> write_data = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.};
    file.WriteVariable("position", write_data);

    // Read every other element starting from index 1
    std::vector<double> read_data(3);
    const std::vector<size_t> start = {1};
    const std::vector<size_t> count = {3};
    const std::vector<ptrdiff_t> stride = {2};
    EXPECT_NO_THROW(file.ReadVariableWithStride("position", start, count, stride, read_data.data()));
    EXPECT_EQ(read_data, std::vector<double>({2.2, 4.4, 6.6}));
}

TEST_F(NetCDFFileTest, ReadVariableWithStrideTypeFloat) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int velocity_var = file.AddVariable<float>("velocity", dim_ids);
    EXPECT_GE(velocity_var, 0);

    // Write full data
    const std::vector<float> write_data = {1.F, 2.F, 3.F, 4.F, 5.F, 6.F, 7.F, 8.F, 9.F, 10.F};
    file.WriteVariable("velocity", write_data);

    // Read every third element starting from index 0
    std::vector<float> read_data(3);
    const std::vector<size_t> start = {0};
    const std::vector<size_t> count = {3};
    const std::vector<ptrdiff_t> stride = {3};
    EXPECT_NO_THROW(file.ReadVariableWithStride("velocity", start, count, stride, read_data.data()));
    EXPECT_EQ(read_data, std::vector<float>({1.F, 4.F, 7.F}));
}

TEST_F(NetCDFFileTest, ReadVariableWithStrideTypeInt) {
    const util::NetCDFFile file(test_file);
    const int dim_id = file.AddDimension("time", 10);
    const std::vector<int> dim_ids = {dim_id};
    const int count_var = file.AddVariable<int>("count", dim_ids);
    EXPECT_GE(count_var, 0);

    // Write full data
    const std::vector<int> write_data = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    file.WriteVariable("count", write_data);

    // Read every fourth element starting from index 2
    std::vector<int> read_data(2);
    const std::vector<size_t> start = {2};
    const std::vector<size_t> count = {2};
    const std::vector<ptrdiff_t> stride = {4};
    EXPECT_NO_THROW(file.ReadVariableWithStride("count", start, count, stride, read_data.data()));
    EXPECT_EQ(read_data, std::vector<int>({3, 7}));
}

TEST_F(NetCDFFileTest, SetChunking_1D) {
    const util::NetCDFFile file(test_file);

    const int time_dim = file.AddDimension("time", NC_UNLIMITED);
    const int var_id = file.AddVariable<double>("state_1d", std::array{time_dim});
    EXPECT_GE(var_id, 0);

    // Set chunking for the variable to 4 timesteps per chunk
    const std::array<size_t, 1> expected_chunk = {4};  // 4 timesteps per chunk
    EXPECT_NO_THROW(file.SetChunking("state_1d", expected_chunk));

    // Verify that NetCDF reports the chunking storage for the variable as expected
    int storage = NC_CONTIGUOUS;          // default storage is contiguous
    std::array<size_t, 1> chunk_sizes{};  // chunk sizes for each dimension
    ASSERT_EQ(
        nc_inq_var_chunking(
            file.GetNetCDFId(), file.GetVariableId("state_1d"), &storage, chunk_sizes.data()
        ),
        NC_NOERR
    );
    EXPECT_EQ(storage, NC_CHUNKED);                // chunked storage is expected
    EXPECT_EQ(chunk_sizes[0], expected_chunk[0]);  // chunk size for the time dimension is 4
}

TEST_F(NetCDFFileTest, SetChunking_2D) {
    const util::NetCDFFile file(test_file);

    const int time_dim = file.AddDimension("time", NC_UNLIMITED);
    const int node_dim = file.AddDimension("nodes", 8);
    const std::vector<int> dim_ids = {time_dim, node_dim};

    const int var_id = file.AddVariable<double>("state_2d", dim_ids);
    EXPECT_GE(var_id, 0);

    // Set chunking for the variable to 3 timesteps per chunk, all 8 nodes
    const std::array<size_t, 2> expected_chunk = {3, 8};  // 3 timesteps, all nodes
    EXPECT_NO_THROW(file.SetChunking("state_2d", expected_chunk));

    // Verify that NetCDF reports the chunking storage for the variable as expected
    int storage = NC_CONTIGUOUS;          // default storage is contiguous
    std::array<size_t, 2> chunk_sizes{};  // chunk sizes for each dimension
    ASSERT_EQ(
        nc_inq_var_chunking(
            file.GetNetCDFId(), file.GetVariableId("state_2d"), &storage, chunk_sizes.data()
        ),
        NC_NOERR
    );
    EXPECT_EQ(storage, NC_CHUNKED);                // chunked storage is expected
    EXPECT_EQ(chunk_sizes[0], expected_chunk[0]);  // chunk size for the time dimension is 3
    EXPECT_EQ(chunk_sizes[1], expected_chunk[1]);  // chunk size for the nodes dimension is 8
}

}  // namespace kynema::tests
