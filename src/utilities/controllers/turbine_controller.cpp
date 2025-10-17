#include "turbine_controller.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace kynema::util {

TurbineController::TurbineController(
    std::string shared_lib_path, std::string controller_function_name, std::string input_file_path,
    std::string output_file_path, bool yaw_control_enabled
)
    : io{},
      input_file_path_(std::move(input_file_path)),
      output_file_path_(std::move(output_file_path)),
      shared_lib_path_(std::move(shared_lib_path)),
      controller_function_name_(std::move(controller_function_name)),
      lib_(shared_lib_path_, util::dylib::no_filename_decorations) {
    // Make sure we have a valid shared library path + controller function name
    try {
        lib_.get_function<void(
            float* avrSWAP, int* aviFAIL, const char* accINFILE, const char* avcOUTNAME,
            const char* avcMSG
        )>(this->controller_function_name_);
    } catch (const util::dylib::load_error& e) {
        throw std::runtime_error("Failed to load shared library: " + shared_lib_path_);
    } catch (const util::dylib::symbol_error& e) {
        throw std::runtime_error("Failed to get function: " + controller_function_name_);
    }

    // Store the controller function from the shared lib in a function pointer for later use
    this->controller_function_ =
        lib_.get_function<void(float*, int*, const char* const, char* const, char* const)>(
            this->controller_function_name_
        );

    // Initialize some values required for calling the controller function
    // std::fill(&this->swap_array_[0], &this->swap_array_[kSwapArraySize], 0.F);
    // this->status_ = 0;                        // Status of the controller function call
    // this->message_ = std::string(1024, ' ');  // 1024 characters for message

    // Map swap array to ControllerIO structure for easier access
    // this->io = reinterpret_cast<ControllerIO*>(this->swap_array_);
    this->io.infile_array_size = input_file_path_.size();
    this->io.outname_array_size = output_file_path_.size();
    this->io.message_array_size = 1024U;
}

void TurbineController::CallController() {
    auto swap_array = std::array<float, kSwapArraySize>{};
    int status{};
    auto message = std::string(1024, ' ');
    io.CopyToSwapArray(swap_array);
    this->controller_function_(
        swap_array.data(), &status, this->input_file_path_.c_str(), this->output_file_path_.data(),
        message.data()
    );
    io.CopyFromSwapArray(swap_array);
    if (status < 0) {
        throw std::runtime_error("Error raised in controller: " + message);
    } else if (status > 0) {
        std::cout << "Warning from controller: " << message << "\n";
    }

    // Integrate yaw angle command from yaw rate command
    if (this->yaw_control_enabled_) {
        this->yaw_angle_command_ += this->io.nacelle_yaw_rate_command * this->io.dt;
    }
}

}  // namespace kynema::util
