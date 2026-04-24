#pragma once

#include <vector>

namespace kynema_fmb::interfaces::components {

struct AerodynamicSection {
    size_t id{};
    double s{};
    double chord{};
    double section_offset_x{};
    double section_offset_y{};
    double aerodynamic_center{};
    double twist{};
    std::vector<double> aoa;
    std::vector<double> cl;
    std::vector<double> cd;
    std::vector<double> cm;
};

struct AerodynamicBodyInput {
    size_t id{};
    std::vector<size_t> beam_node_ids;
    std::vector<AerodynamicSection> aero_sections;
};

class AerodynamicsInput {
public:
    bool is_enabled = false;
    std::vector<std::vector<AerodynamicSection>> aero_inputs;
    std::vector<size_t> airfoil_map;
};
}  // namespace kynema_fmb::interfaces::components
