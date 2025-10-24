#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "utilities/netcdf/netcdf_file.hpp"
#include "utilities/netcdf/node_state_writer.hpp"

namespace kynema::tests {

class NodeStateWriterTest : public ::testing::Test {
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
    size_t num_nodes{3};
    size_t no_buffering{0};
};

TEST_F(NodeStateWriterTest, ConstructorCreatesExpectedDimensionsAndVariables) {
    util::NodeStateWriter writer(test_file, true, num_nodes, no_buffering);
    const auto& file = writer.GetFile();

    EXPECT_EQ(writer.GetNumNodes(), num_nodes);

    EXPECT_NO_THROW({
        // dimensions
        EXPECT_GE(file.GetDimensionId("time"), 0);
        EXPECT_GE(file.GetDimensionId("nodes"), 0);

        // position, x
        EXPECT_GE(file.GetVariableId("x_x"), 0);
        EXPECT_GE(file.GetVariableId("x_y"), 0);
        EXPECT_GE(file.GetVariableId("x_z"), 0);
        EXPECT_GE(file.GetVariableId("x_i"), 0);
        EXPECT_GE(file.GetVariableId("x_j"), 0);
        EXPECT_GE(file.GetVariableId("x_k"), 0);
        EXPECT_GE(file.GetVariableId("x_w"), 0);

        // displacement, u
        EXPECT_GE(file.GetVariableId("u_x"), 0);
        EXPECT_GE(file.GetVariableId("u_y"), 0);
        EXPECT_GE(file.GetVariableId("u_z"), 0);
        EXPECT_GE(file.GetVariableId("u_i"), 0);
        EXPECT_GE(file.GetVariableId("u_j"), 0);
        EXPECT_GE(file.GetVariableId("u_k"), 0);
        EXPECT_GE(file.GetVariableId("u_w"), 0);

        // velocity, v
        EXPECT_GE(file.GetVariableId("v_x"), 0);
        EXPECT_GE(file.GetVariableId("v_y"), 0);
        EXPECT_GE(file.GetVariableId("v_z"), 0);
        EXPECT_GE(file.GetVariableId("v_i"), 0);
        EXPECT_GE(file.GetVariableId("v_j"), 0);
        EXPECT_GE(file.GetVariableId("v_k"), 0);

        // acceleration, a
        EXPECT_GE(file.GetVariableId("a_x"), 0);
        EXPECT_GE(file.GetVariableId("a_y"), 0);
        EXPECT_GE(file.GetVariableId("a_z"), 0);
        EXPECT_GE(file.GetVariableId("a_i"), 0);
        EXPECT_GE(file.GetVariableId("a_j"), 0);
        EXPECT_GE(file.GetVariableId("a_k"), 0);

        // force, f
        EXPECT_GE(file.GetVariableId("f_x"), 0);
        EXPECT_GE(file.GetVariableId("f_y"), 0);
        EXPECT_GE(file.GetVariableId("f_z"), 0);
        EXPECT_GE(file.GetVariableId("f_i"), 0);
        EXPECT_GE(file.GetVariableId("f_j"), 0);
        EXPECT_GE(file.GetVariableId("f_k"), 0);
    });
}

TEST_F(NodeStateWriterTest, WriteStateDataAtTimestepForPosition) {
    util::NodeStateWriter writer(test_file, true, num_nodes, no_buffering);

    const std::vector<double> x = {1., 2., 3.};
    const std::vector<double> y = {4., 5., 6.};
    const std::vector<double> z = {7., 8., 9.};
    const std::vector<double> i = {0.1, 0.2, 0.3};
    const std::vector<double> j = {0.4, 0.5, 0.6};
    const std::vector<double> k = {0.7, 0.8, 0.9};
    const std::vector<double> w = {1., 1., 1.};

    EXPECT_NO_THROW(writer.WriteStateDataAtTimestep(0, "x", x, y, z, i, j, k, w));

    const auto& file = writer.GetFile();
    std::vector<double> read_data(num_nodes);

    file.ReadVariable("x_x", read_data.data());
    EXPECT_EQ(read_data, x);

    file.ReadVariable("x_y", read_data.data());
    EXPECT_EQ(read_data, y);

    file.ReadVariable("x_z", read_data.data());
    EXPECT_EQ(read_data, z);

    file.ReadVariable("x_i", read_data.data());
    EXPECT_EQ(read_data, i);

    file.ReadVariable("x_j", read_data.data());
    EXPECT_EQ(read_data, j);

    file.ReadVariable("x_k", read_data.data());
    EXPECT_EQ(read_data, k);

    file.ReadVariable("x_w", read_data.data());
    EXPECT_EQ(read_data, w);
}

TEST_F(NodeStateWriterTest, WriteStateDataAtTimestepForVelocity) {
    util::NodeStateWriter writer(test_file, true, num_nodes, no_buffering);

    const std::vector<double> x = {1., 2., 3.};
    const std::vector<double> y = {4., 5., 6.};
    const std::vector<double> z = {7., 8., 9.};
    const std::vector<double> i = {0.1, 0.2, 0.3};
    const std::vector<double> j = {0.4, 0.5, 0.6};
    const std::vector<double> k = {0.7, 0.8, 0.9};

    EXPECT_NO_THROW(writer.WriteStateDataAtTimestep(0, "v", x, y, z, i, j, k));

    const auto& file = writer.GetFile();
    std::vector<double> read_data(num_nodes);
    const std::vector<size_t> start = {0, 0};
    const std::vector<size_t> count = {1, num_nodes};

    file.ReadVariableAt("v_x", start, count, read_data.data());
    EXPECT_EQ(read_data, x);

    file.ReadVariableAt("v_y", start, count, read_data.data());
    EXPECT_EQ(read_data, y);

    file.ReadVariableAt("v_z", start, count, read_data.data());
    EXPECT_EQ(read_data, z);

    file.ReadVariableAt("v_i", start, count, read_data.data());
    EXPECT_EQ(read_data, i);

    file.ReadVariableAt("v_j", start, count, read_data.data());
    EXPECT_EQ(read_data, j);

    file.ReadVariableAt("v_k", start, count, read_data.data());
    EXPECT_EQ(read_data, k);
}

TEST_F(NodeStateWriterTest, ThrowsOnInvalidComponentPrefix) {
    util::NodeStateWriter writer(test_file, true, num_nodes, no_buffering);
    const std::vector<double> data(num_nodes, 1.);

    EXPECT_THROW(
        writer.WriteStateDataAtTimestep(0, "invalid_prefix", data, data, data, data, data, data),
        std::invalid_argument
    );
}

TEST_F(NodeStateWriterTest, ThrowsOnMismatchedVectorSizes) {
    util::NodeStateWriter writer(test_file, true, num_nodes, no_buffering);

    const std::vector<double> correct_size(num_nodes, 1.);    // write data to 3 nodes
    const std::vector<double> wrong_size(num_nodes + 1, 1.);  // write data to 4 nodes

    EXPECT_THROW(
        writer.WriteStateDataAtTimestep(
            0, "position", correct_size, wrong_size, correct_size, correct_size, correct_size,
            correct_size, correct_size
        ),
        std::invalid_argument
    );
}

TEST_F(NodeStateWriterTest, SetsChunkingForVariables) {
    const size_t buffer_size{3};
    util::NodeStateWriter writer(test_file, true, num_nodes, buffer_size);
    const auto& file = writer.GetFile();

    const std::array<size_t, 2> expected_chunk = {buffer_size, num_nodes};

    auto expect_chunk = [&](const std::string& var) {
        int storage = NC_CONTIGUOUS;
        std::array<size_t, 2> created_chunk{};
        ASSERT_EQ(
            nc_inq_var_chunking(
                file.GetNetCDFId(), file.GetVariableId(var), &storage, created_chunk.data()
            ),
            NC_NOERR
        );
        EXPECT_EQ(storage, NC_CHUNKED);                  // chunked storage is expected
        EXPECT_EQ(created_chunk[0], expected_chunk[0]);  // chunk size for the time dimension is 3
        EXPECT_EQ(created_chunk[1], expected_chunk[1]);  // chunk size for the nodes dimension is 3
    };

    // position, x
    expect_chunk("x_x");
    expect_chunk("x_y");
    expect_chunk("x_z");
    expect_chunk("x_i");
    expect_chunk("x_j");
    expect_chunk("x_k");
    expect_chunk("x_w");

    // displacement, u
    expect_chunk("u_x");
    expect_chunk("u_y");
    expect_chunk("u_z");
    expect_chunk("u_w");
    expect_chunk("u_i");
    expect_chunk("u_j");
    expect_chunk("u_k");

    // velocity, v
    expect_chunk("v_x");
    expect_chunk("v_y");
    expect_chunk("v_z");
    expect_chunk("v_i");
    expect_chunk("v_j");
    expect_chunk("v_k");

    // acceleration, a
    expect_chunk("a_x");
    expect_chunk("a_y");
    expect_chunk("a_z");
    expect_chunk("a_i");
    expect_chunk("a_j");
    expect_chunk("a_k");

    // force, f
    expect_chunk("f_x");
    expect_chunk("f_y");
    expect_chunk("f_z");
    expect_chunk("f_i");
    expect_chunk("f_j");
    expect_chunk("f_k");

    // deformation, x
    expect_chunk("deformation_x");
    expect_chunk("deformation_y");
    expect_chunk("deformation_z");
}

TEST_F(NodeStateWriterTest, BuffersAndFlushesOnCapacity_State) {
    const size_t buffer_size{2};  // buffer size is 2 timesteps
    util::NodeStateWriter writer(test_file, true, num_nodes, buffer_size);

    // Data for timestep 0
    const std::vector<double> x0 = {1., 2., 3.};
    const std::vector<double> y0 = {4., 5., 6.};
    const std::vector<double> z0 = {7., 8., 9.};
    const std::vector<double> i0 = {0.1, 0.2, 0.3};
    const std::vector<double> j0 = {0.4, 0.5, 0.6};
    const std::vector<double> k0 = {0.7, 0.8, 0.9};

    // After first write (buffer not full) -> no data written yet
    ASSERT_NO_THROW(writer.WriteStateDataAtTimestep(0, "v", x0, y0, z0, i0, j0, k0));
    EXPECT_EQ(writer.GetFile().GetDimensionLength("time"), 0U);

    // Data for timestep 1
    const std::vector<double> x1 = {10., 20., 30.};
    const std::vector<double> y1 = {40., 50., 60.};
    const std::vector<double> z1 = {70., 80., 90.};
    const std::vector<double> i1 = {1.1, 1.2, 1.3};
    const std::vector<double> j1 = {1.4, 1.5, 1.6};
    const std::vector<double> k1 = {1.7, 1.8, 1.9};

    // After second write (flush triggered) -> two timesteps written
    ASSERT_NO_THROW(writer.WriteStateDataAtTimestep(1, "v", x1, y1, z1, i1, j1, k1));
    EXPECT_EQ(writer.GetFile().GetDimensionLength("time"), 2U);

    const auto& file = writer.GetFile();
    std::vector<double> read(num_nodes);
    const std::vector<size_t> start0 = {0, 0};
    const std::vector<size_t> start1 = {1, 0};
    const std::vector<size_t> count = {1, num_nodes};

    file.ReadVariableAt("v_x", start0, count, read.data());
    EXPECT_EQ(read, x0);
    file.ReadVariableAt("v_y", start0, count, read.data());
    EXPECT_EQ(read, y0);
    file.ReadVariableAt("v_z", start0, count, read.data());
    EXPECT_EQ(read, z0);
    file.ReadVariableAt("v_i", start0, count, read.data());
    EXPECT_EQ(read, i0);
    file.ReadVariableAt("v_j", start0, count, read.data());
    EXPECT_EQ(read, j0);
    file.ReadVariableAt("v_k", start0, count, read.data());
    EXPECT_EQ(read, k0);

    file.ReadVariableAt("v_x", start1, count, read.data());
    EXPECT_EQ(read, x1);
    file.ReadVariableAt("v_y", start1, count, read.data());
    EXPECT_EQ(read, y1);
    file.ReadVariableAt("v_z", start1, count, read.data());
    EXPECT_EQ(read, z1);
    file.ReadVariableAt("v_i", start1, count, read.data());
    EXPECT_EQ(read, i1);
    file.ReadVariableAt("v_j", start1, count, read.data());
    EXPECT_EQ(read, j1);
    file.ReadVariableAt("v_k", start1, count, read.data());
    EXPECT_EQ(read, k1);
}

TEST_F(NodeStateWriterTest, BuffersAndFlushesOnCapacity_Deformation) {
    const size_t buffer_size{2};  // buffer size is 2 timesteps
    util::NodeStateWriter writer(test_file, true, num_nodes, buffer_size);

    // Data for timestep 0
    const std::vector<double> x0 = {1., 2., 3.};
    const std::vector<double> y0 = {4., 5., 6.};
    const std::vector<double> z0 = {7., 8., 9.};

    // After first write (buffer not full) -> no data written yet
    ASSERT_NO_THROW(writer.WriteDeformationDataAtTimestep(0, x0, y0, z0));
    EXPECT_EQ(writer.GetFile().GetDimensionLength("time"), 0U);

    // Data for timestep 1
    const std::vector<double> x1 = {10., 20., 30.};
    const std::vector<double> y1 = {40., 50., 60.};
    const std::vector<double> z1 = {70., 80., 90.};

    // After second write (flush triggered) -> two timesteps written
    ASSERT_NO_THROW(writer.WriteDeformationDataAtTimestep(1, x1, y1, z1));
    EXPECT_EQ(writer.GetFile().GetDimensionLength("time"), 2U);

    const auto& file = writer.GetFile();
    std::vector<double> read(num_nodes);
    const std::vector<size_t> start0 = {0, 0};
    const std::vector<size_t> start1 = {1, 0};
    const std::vector<size_t> count = {1, num_nodes};

    file.ReadVariableAt("deformation_x", start0, count, read.data());
    EXPECT_EQ(read, x0);
    file.ReadVariableAt("deformation_y", start0, count, read.data());
    EXPECT_EQ(read, y0);
    file.ReadVariableAt("deformation_z", start0, count, read.data());
    EXPECT_EQ(read, z0);

    file.ReadVariableAt("deformation_x", start1, count, read.data());
    EXPECT_EQ(read, x1);
    file.ReadVariableAt("deformation_y", start1, count, read.data());
    EXPECT_EQ(read, y1);
    file.ReadVariableAt("deformation_z", start1, count, read.data());
    EXPECT_EQ(read, z1);
}

TEST_F(NodeStateWriterTest, DestructorFlushesRemainingData_State) {
    // Scope the writer object to ensure it goes out of scope and flushes remaining data in
    // destructor
    {
        const size_t buffer_size{3};  // buffer size is 3 timesteps
        util::NodeStateWriter writer(test_file, true, num_nodes, buffer_size);

        // Data for timestep 0
        const std::vector<double> x0 = {1., 2., 3.};
        const std::vector<double> y0 = {4., 5., 6.};
        const std::vector<double> z0 = {7., 8., 9.};
        const std::vector<double> i0 = {0.1, 0.2, 0.3};
        const std::vector<double> j0 = {0.4, 0.5, 0.6};
        const std::vector<double> k0 = {0.7, 0.8, 0.9};
        const std::vector<double> w0 = {1., 1., 1.};

        // Data for timestep 1
        const std::vector<double> x1 = {10., 20., 30.};
        const std::vector<double> y1 = {40., 50., 60.};
        const std::vector<double> z1 = {70., 80., 90.};
        const std::vector<double> i1 = {1.1, 1.2, 1.3};
        const std::vector<double> j1 = {1.4, 1.5, 1.6};
        const std::vector<double> k1 = {1.7, 1.8, 1.9};
        const std::vector<double> w1 = {2., 2., 2.};

        ASSERT_NO_THROW(writer.WriteStateDataAtTimestep(0, "x", x0, y0, z0, i0, j0, k0, w0));
        ASSERT_NO_THROW(writer.WriteStateDataAtTimestep(1, "x", x1, y1, z1, i1, j1, k1, w1));

        // Flush has not been triggered yet, so there should be no data in file
        EXPECT_EQ(writer.GetFile().GetDimensionLength("time"), 0U);
    }

    // Verify that the data has been written to the file in the correct order
    const util::NetCDFFile file(test_file, false);
    std::vector<double> read(num_nodes);
    const std::vector<size_t> start0 = {0, 0};
    const std::vector<size_t> start1 = {1, 0};
    const std::vector<size_t> count = {1, num_nodes};

    file.ReadVariableAt("x_x", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{1., 2., 3.}));
    file.ReadVariableAt("x_y", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{4., 5., 6.}));
    file.ReadVariableAt("x_z", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{7., 8., 9.}));
    file.ReadVariableAt("x_i", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{0.1, 0.2, 0.3}));
    file.ReadVariableAt("x_j", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{0.4, 0.5, 0.6}));
    file.ReadVariableAt("x_k", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{0.7, 0.8, 0.9}));
    file.ReadVariableAt("x_w", start0, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{1., 1., 1.}));

    file.ReadVariableAt("x_x", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{10., 20., 30.}));
    file.ReadVariableAt("x_y", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{40., 50., 60.}));
    file.ReadVariableAt("x_z", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{70., 80., 90.}));
    file.ReadVariableAt("x_i", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{1.1, 1.2, 1.3}));
    file.ReadVariableAt("x_j", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{1.4, 1.5, 1.6}));
    file.ReadVariableAt("x_k", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{1.7, 1.8, 1.9}));
    file.ReadVariableAt("x_w", start1, count, read.data());
    EXPECT_EQ(read, (std::vector<double>{2., 2., 2.}));
}

// Disabled by default due to environment variability. Run locally to inspect timings.
TEST_F(NodeStateWriterTest, DISABLED_BufferedChunkingPerformance) {
    using clock = std::chrono::steady_clock;
    using ms = std::chrono::milliseconds;

    const std::string unbuffered_file{test_file + "_unbuffered.nc"};
    const std::string buffered_file{test_file + "_buffered.nc"};
    std::filesystem::remove(unbuffered_file);
    std::filesystem::remove(buffered_file);

    // Node and time step sizes
    const size_t nodes{100};
    const size_t steps{1000};
    auto make_vector = [&](double base) {
        std::vector<double> vector(nodes);
        for (size_t i = 0; i < nodes; ++i)
            vector[i] = base + static_cast<double>(i) * 0.001;
        return vector;
    };

    // Pre-create per-timestep data (keeps timing focused on NetCDF IO)
    std::vector<std::array<std::vector<double>, 6>> data;
    data.reserve(steps);
    for (size_t t = 0; t < steps; ++t) {
        data.push_back(
            {make_vector(1. + static_cast<double>(t)), make_vector(2. + static_cast<double>(t)),
             make_vector(3. + static_cast<double>(t)),  // x, y, z
             make_vector(4. + static_cast<double>(t)), make_vector(5. + static_cast<double>(t)),
             make_vector(6. + static_cast<double>(t))}  // i, j, k
        );
    }

    // Unbuffered (buffer_size = 0) — chunk size along time will be 1
    auto start_time_unbuffered = clock::now();
    {
        util::NodeStateWriter writer(unbuffered_file, true, nodes, /*buffer_size=*/0);
        for (size_t t = 0; t < steps; ++t) {
            const auto& rec = data[t];
            writer.WriteStateDataAtTimestep(t, "v", rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]);
        }
    }
    auto end_time_unbuffered = clock::now();

    // Buffered (buffer_size = 32) — chunk size along time equals buffer size
    auto start_time_buffered = clock::now();
    {
        util::NodeStateWriter writer(buffered_file, true, nodes, /*buffer_size=*/32);
        for (size_t t = 0; t < steps; ++t) {
            const auto& rec = data[t];
            writer.WriteStateDataAtTimestep(t, "v", rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]);
        }
    }
    auto end_time_buffered = clock::now();

    const auto unbuffered_time =
        std::chrono::duration_cast<ms>(end_time_unbuffered - start_time_unbuffered).count();
    const auto buffered_time =
        std::chrono::duration_cast<ms>(end_time_buffered - start_time_buffered).count();

    std::cout << "Performance of NodeStateWriter with " << nodes << " nodes and " << steps
              << " steps: " << "\n"
              << " unbuffered : " << unbuffered_time << " ms" << "\n"
              << " buffered : " << buffered_time << " ms" << "\n";
    EXPECT_LT(buffered_time, unbuffered_time);  // buffered should be faster

    std::filesystem::remove(unbuffered_file);
    std::filesystem::remove(buffered_file);
}

}  // namespace kynema::tests
