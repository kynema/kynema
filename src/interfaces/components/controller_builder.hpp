#pragma once

#include <string>
#include <string_view>

#include "controller_input.hpp"

namespace kynema::interfaces::components {

class ControllerBuilder {
public:
    ControllerInput& Input() { return input; }

    ControllerBuilder& SetLibraryPath(std::string_view lib_path) {
        input.shared_lib_path = std::string(lib_path);
        return *this;
    }

    ControllerBuilder& SetFunctionName(std::string_view func_name) {
        input.function_name = std::string(func_name);
        return *this;
    }

    ControllerBuilder& SetInputFilePath(std::string_view inp_file_path) {
        input.input_file_path = std::string(inp_file_path);
        return *this;
    }

    ControllerBuilder& SetControllerInput(std::string_view sim_name) {
        input.simulation_name = std::string(sim_name);
        return *this;
    }

private:
    ControllerInput input;
};

}  // namespace kynema::interfaces::components
