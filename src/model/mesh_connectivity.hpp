#pragma once

#include <algorithm>
#include <fstream>
#include <span>
#include <unordered_map>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace kynema::model {

/**
 * @brief Class to manage element-to-node connectivity information for a mesh
 *
 * This class stores and manages the relationships between elements/constraints
 * and the nodes they connect to, providing YAML import/export functionality.
 */
class MeshConnectivity {
public:
    MeshConnectivity() = default;

    /**
     * @brief Adds a beam element's node connectivity
     * @param elem_id The beam element ID
     * @param node_ids The vector of node IDs connected to this beam
     */
    void AddBeamElementConnectivity(size_t elem_id, std::span<const size_t> node_ids) {
        beams_[elem_id].assign(std::begin(node_ids), std::end(node_ids));
    }

    /**
     * @brief Adds a mass element's node connectivity
     * @param elem_id The mass element ID
     * @param node_id The node ID connected to this mass
     */
    void AddMassElementConnectivity(size_t elem_id, size_t node_id) {
        masses_[elem_id] = std::vector<size_t>{node_id};
    }

    /**
     * @brief Adds a spring element's node connectivity
     * @param elem_id The spring element ID
     * @param node_ids The array of node IDs connected to this spring
     */
    void AddSpringElementConnectivity(size_t elem_id, const std::array<size_t, 2>& node_ids) {
        springs_[elem_id] = std::vector<size_t>{node_ids[0], node_ids[1]};
    }

    /**
     * @brief Adds a constraint's node connectivity
     * @param constraint_id The constraint ID
     * @param node_ids The vector of node IDs connected to this constraint
     */
    void AddConstraintConnectivity(size_t constraint_id, std::span<const size_t> node_ids) {
        constraints_[constraint_id].assign(std::begin(node_ids), std::end(node_ids));
    }

    /**
     * @brief Get nodes connected to a specific beam element
     * @param elem_id The beam element ID
     * @return Vector of node IDs
     */
    [[nodiscard]] const std::vector<size_t>& GetBeamElementConnectivity(size_t elem_id) const {
        return beams_.at(elem_id);
    }

    /**
     * @brief Get nodes connected to a specific mass element
     * @param elem_id The mass element ID
     * @return Vector of node IDs
     */
    [[nodiscard]] const std::vector<size_t>& GetMassElementConnectivity(size_t elem_id) const {
        return masses_.at(elem_id);
    }

    /**
     * @brief Get nodes connected to a specific spring element
     * @param elem_id The spring element ID
     * @return Vector of node IDs
     */
    [[nodiscard]] const std::vector<size_t>& GetSpringElementConnectivity(size_t elem_id) const {
        return springs_.at(elem_id);
    }

    /**
     * @brief Get nodes connected to a specific constraint
     * @param constraint_id The constraint ID
     * @return Vector of node IDs
     */
    [[nodiscard]] const std::vector<size_t>& GetConstraintConnectivity(size_t constraint_id) const {
        return constraints_.at(constraint_id);
    }

    /**
     * @brief Export mesh connectivity inforation to a YAML file
     * @param file Stream to output YAML file
     */
    void ExportToYAML(std::ostream& file) const {
        YAML::Node root;

        ExportElementTypeToYAML(root, "beams", beams_);
        ExportElementTypeToYAML(root, "masses", masses_);
        ExportElementTypeToYAML(root, "springs", springs_);
        ExportElementTypeToYAML(root, "constraints", constraints_);

        file << root;
    }

    /**
     * @brief Export mesh connectivity information to a YAML file
     * @param filename Path to the output YAML file
     */
    void ExportToYAML(const std::string& filename) const {
        std::ofstream file(filename);
        ExportToYAML(file);
    }

    /**
     * @brief Import mesh connectivity information from a YAML file
     * @param root YAML node with the input YAML file loaded
     */
    void ImportFromYAML(const YAML::Node& root) {
        masses_.clear();
        springs_.clear();
        beams_.clear();
        constraints_.clear();

        // Import masses
        if (root["masses"]) {
            for (const auto& entry : root["masses"]) {
                const size_t id = std::stoul(entry.first.as<std::string>());
                masses_[id] = entry.second.as<std::vector<size_t>>();
            }
        }

        // Import springs
        if (root["springs"]) {
            for (const auto& entry : root["springs"]) {
                const size_t id = std::stoul(entry.first.as<std::string>());
                springs_[id] = entry.second.as<std::vector<size_t>>();
            }
        }

        // Import beams
        if (root["beams"]) {
            for (const auto& entry : root["beams"]) {
                const size_t id = std::stoul(entry.first.as<std::string>());
                beams_[id] = entry.second.as<std::vector<size_t>>();
            }
        }

        // Import constraints
        if (root["constraints"]) {
            for (const auto& entry : root["constraints"]) {
                const size_t id = std::stoul(entry.first.as<std::string>());
                constraints_[id] = entry.second.as<std::vector<size_t>>();
            }
        }
    }

    /**
     * @brief Import mesh connectivity information from a YAML file
     * @param filename Path to the input YAML file
     */
    void ImportFromYAML(const std::string& filename) { ImportFromYAML(YAML::LoadFile(filename)); }

private:
    std::unordered_map<size_t, std::vector<size_t>> beams_;
    std::unordered_map<size_t, std::vector<size_t>> masses_;
    std::unordered_map<size_t, std::vector<size_t>> springs_;
    std::unordered_map<size_t, std::vector<size_t>> constraints_;

    /**
     * @brief Helper function to export a specific element type with ordered IDs
     * @param root The root YAML node
     * @param element_type The type of element (beams, masses, etc.)
     * @param map The map containing the element data
     */
    template <typename MapType>
    static void ExportElementTypeToYAML(
        YAML::Node& root, const std::string& element_type, const MapType& map
    ) {
        if (map.empty()) {
            return;
        }

        YAML::Node element_node;

        // Get all keys from the map and sort them to maintain consistent order
        std::vector<size_t> keys;
        keys.reserve(map.size());
        for (const auto& [id, _] : map) {
            keys.push_back(id);
        }
        std::ranges::sort(keys);

        // Add each element in sorted order of id
        for (const auto& id : keys) {
            YAML::Node array_node;
            std::for_each(map.at(id).begin(), map.at(id).end(), [&array_node](const auto& node_id) {
                array_node.push_back(node_id);
            });

            array_node.SetStyle(YAML::EmitterStyle::Flow);
            element_node[std::to_string(id)] = array_node;
        }

        root[element_type] = element_node;
    }
};

}  // namespace kynema::model
