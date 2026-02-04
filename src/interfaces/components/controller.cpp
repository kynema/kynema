#include "controller.hpp"

#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace kynema::interfaces::components {

Controller::Controller(const ControllerInput& input)
    : io{},
      pitch_control_enabled_(input.pitch_control_enabled),
      torque_control_enabled_(input.torque_control_enabled),
      yaw_control_enabled_(input.yaw_control_enabled),
      gearbox_ratio_(input.gearbox_ratio),
      torque_command_(0.0),
      pitch_angle_command_(input.pitch_angle),
      yaw_angle_command_(input.yaw_angle),
      input_file_path_(input.input_file_path),
      output_file_path_(input.output_file_path),
      shared_lib_path_(input.shared_lib_path),
      controller_function_name_(input.function_name),
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

    // Map swap array to ControllerIO structure for easier access
    // this->io = reinterpret_cast<ControllerIO*>(this->swap_array_);
    this->io.infile_array_size = input_file_path_.size();
    this->io.outname_array_size = output_file_path_.size();
    this->io.message_array_size = 1024U;
}

void Controller::CallController() {
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

}  // namespace kynema::interfaces::components
