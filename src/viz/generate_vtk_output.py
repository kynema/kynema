import argparse
import numpy as np
import os
import vtk
import yaml

from netCDF4 import Dataset
from typing import Dict, List, Optional

# -------------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------------


def quaternion_to_rotation_matrix(quaternion: List[float]) -> np.ndarray:
    """Converts a 4x1 quaternion [w, i, j, k] to a 3x3 rotation matrix.

    The rotation matrix is computed using the following formula (used in Kynema):
    R = | 1 - 2(j^2 + k^2)  2(i*j - w*k)  2(i*k + w*j) |
        | 2(i*j + w*k)  1 - 2(i^2 + k^2)  2(j*k - w*i) |
        | 2(i*k - w*j)  2(j*k + w*i)  1 - 2(i^2 + j^2) |

    Args:
        quaternion (list[float]): A list of 4 floats representing the quaternion

    Returns:
        np.ndarray: A 3x3 rotation matrix
    """
    # Unpack components from quaternion
    w, i, j, k = quaternion

    # Compute intermediate products
    ii, jj, kk = i * i, j * j, k * k
    ij, ik, jk = i * j, i * k, j * k
    wi, wj, wk = w * i, w * j, w * k

    # Create rotation matrix
    R = np.zeros((3, 3))

    R[0, 0] = 1.0 - 2.0 * (jj + kk)
    R[0, 1] = 2.0 * (ij - wk)
    R[0, 2] = 2.0 * (ik + wj)

    R[1, 0] = 2.0 * (ij + wk)
    R[1, 1] = 1.0 - 2.0 * (ii + kk)
    R[1, 2] = 2.0 * (jk - wi)

    R[2, 0] = 2.0 * (ik - wj)
    R[2, 1] = 2.0 * (jk + wi)
    R[2, 2] = 1.0 - 2.0 * (ii + jj)

    return R


def create_vector_array(name: str, num_components: int = 3):
    """Creates a VTK double array with the given name and number of components.

    Args:
        name (str): Name of the array
        num_components (int): Number of components in the array

    Returns:
        vtk.vtkDoubleArray: A VTK double array with the given name and number of components
    """
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(num_components)
    array.SetName(name)

    return array


# -------------------------------------------------------------------------------
# Core VTK output logic
# -------------------------------------------------------------------------------


class VTKOutput:
    """Class to generate VTK files from Kynema (NetCDF-based) outputs and mesh connectivity (YAML-based)."""

    def __init__(self, netcdf_path: str, connectivity_path: str):
        """Initializes the visualizer with the path to the NetCDF file and mesh connectivity file.

        Args:
            netcdf_path (str): Path to the NetCDF output file
            connectivity_path (str): Path to the mesh connectivity YAML file
        """
        self.netcdf_path = netcdf_path
        self.data = Dataset(netcdf_path, "r")

        # Initialize mesh connectivity data
        self.mesh_connectivity = {
            "beams": {},
            "masses": {},
            "springs": {},
            "constraints": {},
        }
        self._load_mesh_connectivity(connectivity_path)

        # Get dimensions from the NetCDF file
        self.num_nodes = len(self.data.dimensions["nodes"])
        self.num_timesteps = len(self.data.dimensions["time"])

        print(
            f"Loaded output data from {self.netcdf_path} with {self.num_nodes} nodes and {self.num_timesteps} timesteps"
        )

    def _load_mesh_connectivity(self, connectivity_path: str):
        """Loads mesh connectivity information from a YAML file.

        Args:
            connectivity_path (str): Path to the connectivity YAML file
        """
        try:
            with open(connectivity_path, "r") as file:
                data = yaml.safe_load(file)

                # Process all element types
                for elem_type in ["beams", "masses", "springs", "constraints"]:
                    if elem_type in data:
                        for elem_id, node_ids in data[elem_type].items():
                            self.mesh_connectivity[elem_type][int(elem_id)] = node_ids

            print(f"Loaded connectivity data from {connectivity_path}")
            print(f"  Beams: {len(self.mesh_connectivity['beams'])}")
            print(f"  Masses: {len(self.mesh_connectivity['masses'])}")
            print(f"  Springs: {len(self.mesh_connectivity['springs'])}")
            print(f"  Constraints: {len(self.mesh_connectivity['constraints'])}")

        except Exception as e:
            print(f"Error loading connectivity data: {e}")

    def _extract_node_data_at_timestep(
        self, timestep: int, node_indices: Optional[List[int]] = None
    ) -> List[Dict[str, List[float]]]:
        """Extracts node data for a specific timestep and returns a list of Kynema's NodeData-like structures.

        Node data contains the following components:
        - Position: x_x, x_y, x_z, x_w, x_i, x_j, x_k
        - Velocity: v_x, v_y, v_z, v_i, v_j, v_k
        - Acceleration: a_x, a_y, a_z, a_i, a_j, a_k
        - Force (if available): f_x, f_y, f_z
        - Moment (if available): f_i, f_j, f_k
        - Deformation (if available): deformation_x, deformation_y, deformation_z

        Args:
            timestep (int): The timestep to extract node data from
            node_indices (Optional[List[int]]):
                If provided, only extract data for these specific node indices.
                If None, extract data for all nodes.

        Returns:
            list[dict]: A list of NodeData-like structures (dicts)
        """
        if timestep >= self.num_timesteps:
            raise ValueError(
                f"Timestep {timestep} out of range (max: {self.num_timesteps-1})"
            )

        # Validate that all requested indices are within range
        if node_indices is not None:
            for idx in node_indices:
                if idx < 0 or idx >= self.num_nodes:
                    raise ValueError(
                        f"Node index {idx} out of range (max: {self.num_nodes-1})"
                    )
            indices_to_extract = node_indices

        # If no specific nodes requested, use all nodes
        if node_indices is None:
            indices_to_extract = range(self.num_nodes)

        # ------------------------------------------------------------
        # Extract data that are always available
        # ------------------------------------------------------------

        # Read all data from NetCDF file at once for the requested timestep
        position = [
            self.data.variables["x_x"][timestep, :],
            self.data.variables["x_y"][timestep, :],
            self.data.variables["x_z"][timestep, :],
            self.data.variables["x_w"][timestep, :],
            self.data.variables["x_i"][timestep, :],
            self.data.variables["x_j"][timestep, :],
            self.data.variables["x_k"][timestep, :],
        ]
        velocity = [
            self.data.variables["v_x"][timestep, :],
            self.data.variables["v_y"][timestep, :],
            self.data.variables["v_z"][timestep, :],
            self.data.variables["v_i"][timestep, :],
            self.data.variables["v_j"][timestep, :],
            self.data.variables["v_k"][timestep, :],
        ]
        acceleration = [
            self.data.variables["a_x"][timestep, :],
            self.data.variables["a_y"][timestep, :],
            self.data.variables["a_z"][timestep, :],
            self.data.variables["a_i"][timestep, :],
            self.data.variables["a_j"][timestep, :],
            self.data.variables["a_k"][timestep, :],
        ]

        # ------------------------------------------------------------
        # Read optional force/moment and deformation data if available
        # ------------------------------------------------------------

        force = None
        moment = None
        deformation = None

        force_available = all(
            f"f_{comp}" in self.data.variables for comp in ["x", "y", "z"]
        )
        if force_available:
            force = [
                self.data.variables["f_x"][timestep, :],
                self.data.variables["f_y"][timestep, :],
                self.data.variables["f_z"][timestep, :],
            ]
        moment_available = all(
            f"f_{comp}" in self.data.variables for comp in ["i", "j", "k"]
        )
        if moment_available:
            moment = [
                self.data.variables["f_i"][timestep, :],
                self.data.variables["f_j"][timestep, :],
                self.data.variables["f_k"][timestep, :],
            ]
        deformation_available = all(
            f"deformation_{comp}" in self.data.variables for comp in ["x", "y", "z"]
        )
        if deformation_available:
            deformation = [
                self.data.variables["deformation_x"][timestep, :],
                self.data.variables["deformation_y"][timestep, :],
                self.data.variables["deformation_z"][timestep, :],
            ]

        # Create node data dictionaries only for the requested indices
        nodes = []
        for i_node in indices_to_extract:
            node = {
                "position": [position[j_comp][i_node] for j_comp in range(7)],
                "velocity": [velocity[j_comp][i_node] for j_comp in range(6)],
                "acceleration": [acceleration[j_comp][i_node] for j_comp in range(6)],
                "force": (
                    [force[j_comp][i_node] for j_comp in range(3)]
                    if force is not None
                    else None
                ),
                "moment": (
                    [moment[j_comp][i_node] for j_comp in range(3)]
                    if moment is not None
                    else None
                ),
                "deformation": (
                    [deformation[j_comp][i_node] for j_comp in range(3)]
                    if deformation is not None
                    else None
                ),
            }
            nodes.append(node)

        return nodes

    def _add_node_data_to_vtk_object(
        self, vtk_object: vtk.vtkObject, nodes: List[Dict[str, List[float]]]
    ):
        """Adds common node data to a VTK object (polydata or unstructured grid).

        Following data are added to the VTK object:
        - Orientation (as 3x3 rotation matrix)
        - Velocity (as 6x1 vector)
        - Acceleration (as 6x1 vector)

        Args:
            vtk_object: VTK object to add data to
            nodes: List of node data dictionaries
        """
        # Add orientation data
        orientation_arrays = {
            axis: create_vector_array(f"Orientation{axis}") for axis in ["X", "Y", "Z"]
        }

        for node in nodes:
            # Convert quaternion to rotation matrix
            quaternion = [
                node["position"][3],  # w
                node["position"][4],  # i
                node["position"][5],  # j
                node["position"][6],  # k
            ]
            R = quaternion_to_rotation_matrix(quaternion)

            orientation_arrays["X"].InsertNextTuple3(R[0, 0], R[1, 0], R[2, 0])
            orientation_arrays["Y"].InsertNextTuple3(R[0, 1], R[1, 1], R[2, 1])
            orientation_arrays["Z"].InsertNextTuple3(R[0, 2], R[1, 2], R[2, 2])

        for axis in ["X", "Y", "Z"]:
            vtk_object.GetPointData().AddArray(orientation_arrays[axis])

        # Add velocity data
        translation_velocity = create_vector_array("TranslationalVelocity")
        rotation_velocity = create_vector_array("RotationalVelocity")

        for node in nodes:
            translation_velocity.InsertNextTuple3(*node["velocity"][0:3])
            rotation_velocity.InsertNextTuple3(*node["velocity"][3:6])

        vtk_object.GetPointData().AddArray(translation_velocity)
        vtk_object.GetPointData().AddArray(rotation_velocity)

        # Add acceleration data
        trans_accel = create_vector_array("TranslationalAcceleration")
        rot_accel = create_vector_array("RotationalAcceleration")

        for node in nodes:
            trans_accel.InsertNextTuple3(*node["acceleration"][0:3])
            rot_accel.InsertNextTuple3(*node["acceleration"][3:6])

        vtk_object.GetPointData().AddArray(trans_accel)
        vtk_object.GetPointData().AddArray(rot_accel)

        # Add load data
        force = create_vector_array("Force")
        moment = create_vector_array("Moment")

        for node in nodes:
            force.InsertNextTuple3(*node["force"])
            moment.InsertNextTuple3(*node["moment"])

        vtk_object.GetPointData().AddArray(force)
        vtk_object.GetPointData().AddArray(moment)

    def generate_visualization(self, timestep: int, output_dir: str):
        """Generates visualization for the specified timestep based on mesh connectivity.

        This method automatically determines what elements to create based on the
        available connectivity data.

        Args:
            timestep (int): Timestep to visualize
            output_dir (str): Directory to save the output files
        """
        os.makedirs(output_dir, exist_ok=True)

        nodes = self._extract_node_data_at_timestep(timestep)
        grid = vtk.vtkUnstructuredGrid()

        # Create points and add node IDs
        points = vtk.vtkPoints()
        node_id_array = create_vector_array("NodeID", 1)
        for i in range(len(nodes)):
            position = nodes[i]["position"][0:3]
            points.InsertNextPoint(position)
            node_id_array.InsertNextValue(i)  # Add node ID

        grid.SetPoints(points)
        grid.GetPointData().AddArray(node_id_array)  # Add node IDs to point data

        # Add node data to the grid
        self._add_node_data_to_vtk_object(grid, nodes)

        # Process each element type based on connectivity
        cell_types = {}
        element_ids = {}
        element_type_names = {}

        # Add beam elements to the grid
        if "beams" in self.mesh_connectivity:
            self._add_beams_to_grid(grid, cell_types, element_ids, element_type_names)

        # Add mass nodes to the grid
        if "masses" in self.mesh_connectivity:
            self._add_masses_to_grid(grid, cell_types, element_ids, element_type_names)

        # Add spring elements to the grid
        if "springs" in self.mesh_connectivity:
            self._add_springs_to_grid(grid, cell_types, element_ids, element_type_names)

        # Add constraints to the grid
        if "constraints" in self.mesh_connectivity:
            self._add_constraints_to_grid(
                grid, cell_types, element_ids, element_type_names
            )

        # Add element type array
        type_array = create_vector_array("ElementType", 1)
        type_name_array = vtk.vtkStringArray()
        type_name_array.SetName("ElementTypeName")

        # Add element ID array
        element_id_array = create_vector_array("ElementID", 1)
        for cell_id in range(grid.GetNumberOfCells()):
            type_array.InsertNextValue(cell_types.get(cell_id, 0))
            type_name_array.InsertNextValue(element_type_names.get(cell_id, "Unknown"))
            element_id_array.InsertNextValue(element_ids.get(cell_id, -1))

        grid.GetCellData().AddArray(type_array)
        grid.GetCellData().AddArray(type_name_array)
        grid.GetCellData().AddArray(element_id_array)

        # Write the file
        filename = os.path.join(output_dir, f"timestep_{timestep:04d}.vtu")
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(grid)
        writer.Write()

        print(f"Wrote visualization to {filename}")

        return filename

    def _add_beams_to_grid(self, grid, cell_types, element_ids, element_type_names):
        """Adds beam elements to the unstructured grid. Represented as polyline cells.

        Args:
            grid (vtk.vtkUnstructuredGrid): Grid to add beams to
            cell_types (Dict[int, int]): Dictionary to track cell types
            element_ids (Dict[int, int]): Dictionary to track original element IDs
            element_type_names (Dict[int, str]): Dictionary to track element type names
        """
        beam_type_id = 1  # ID for beam elements -> 1

        for beam_id, node_ids in self.mesh_connectivity["beams"].items():
            # Beams should have at least two nodes i.e. be linear element
            if len(node_ids) < 2:
                continue

            # Create polyline cell
            line = vtk.vtkPolyLine()
            line.GetPointIds().SetNumberOfIds(len(node_ids))

            for i, node_id in enumerate(node_ids):
                line.GetPointIds().SetId(i, node_id)

            cell_id = grid.InsertNextCell(line.GetCellType(), line.GetPointIds())
            cell_types[cell_id] = beam_type_id
            element_ids[cell_id] = beam_id
            element_type_names[cell_id] = "Beam"

    def _add_masses_to_grid(self, grid, cell_types, element_ids, element_type_names):
        """Adds mass elements to the unstructured grid. Represented as vertices.

        Args:
            grid (vtk.vtkUnstructuredGrid): Grid to add masses to
            cell_types (Dict[int, int]): Dictionary to track cell types
            element_ids (Dict[int, int]): Dictionary to track original element IDs
            element_type_names (Dict[int, str]): Dictionary to track element type names
        """
        mass_type_id = 2  # ID for mass elements -> 2

        for mass_id, node_ids in self.mesh_connectivity["masses"].items():
            # Masses should have exactly one node i.e. be point element
            if len(node_ids) != 1:
                continue

            # Create vertex cell
            vertex = vtk.vtkVertex()
            vertex.GetPointIds().SetId(0, node_ids[0])

            cell_id = grid.InsertNextCell(vertex.GetCellType(), vertex.GetPointIds())
            cell_types[cell_id] = mass_type_id
            element_ids[cell_id] = mass_id
            element_type_names[cell_id] = "Mass"

    def _add_springs_to_grid(self, grid, cell_types, element_ids, element_type_names):
        """Adds spring elements to the unstructured grid. Represented as line cells.

        Args:
            grid (vtk.vtkUnstructuredGrid): Grid to add springs to
            cell_types (Dict[int, int]): Dictionary to track cell types
            element_ids (Dict[int, int]): Dictionary to track original element IDs
            element_type_names (Dict[int, str]): Dictionary to track element type names
        """
        spring_type_id = 3  # ID for spring elements -> 3
        for spring_id, node_ids in self.mesh_connectivity["springs"].items():
            # Springs should have exactly two nodes i.e. be linear element
            if len(node_ids) != 2:
                continue

            # Create line cell
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, node_ids[0])
            line.GetPointIds().SetId(1, node_ids[1])

            cell_id = grid.InsertNextCell(line.GetCellType(), line.GetPointIds())
            cell_types[cell_id] = spring_type_id
            element_ids[cell_id] = spring_id
            element_type_names[cell_id] = "Spring"

    def _add_constraints_to_grid(
        self, grid, cell_types, element_ids, element_type_names
    ):
        """Adds constraint elements to the unstructured grid. Represented as line cells.

        Args:
            grid (vtk.vtkUnstructuredGrid): Grid to add constraints to
            cell_types (Dict[int, int]): Dictionary to track cell types
            element_ids (Dict[int, int]): Dictionary to track original element IDs
            element_type_names (Dict[int, str]): Dictionary to track element type names
        """
        constraint_type_id = 4  # ID for constraint elements -> 4

        for constraint_id, node_ids in self.mesh_connectivity["constraints"].items():
            # NOTE: We are not considering single node constraints for viz
            if len(node_ids) != 2:
                continue

            # Create line cell
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, node_ids[0])
            line.GetPointIds().SetId(1, node_ids[1])

            cell_id = grid.InsertNextCell(line.GetCellType(), line.GetPointIds())
            cell_types[cell_id] = constraint_type_id
            element_ids[cell_id] = constraint_id
            element_type_names[cell_id] = "Constraint"

    def visualize_all_timesteps(self, output_dir: str):
        """Generates visualization for all timesteps.

        Args:
            output_dir (str): Directory to save the output files
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Create PVD file to group all timesteps
        pvd_filename = os.path.join(output_dir, "simulation.pvd")
        with open(pvd_filename, "w") as pvd_file:
            pvd_file.write('<?xml version="1.0"?>\n')
            pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
            pvd_file.write("  <Collection>\n")

            # Get time values
            times = np.arange(self.num_timesteps, dtype=float)

            # Generate visualization for each timestep
            for timestep in range(self.num_timesteps):
                vtu_file = self.generate_visualization(timestep, output_dir)
                vtu_basename = os.path.basename(vtu_file)

                # Add to collection
                time_value = times[timestep]
                pvd_file.write(
                    f'    <DataSet timestep="{time_value}" file="{vtu_basename}"/>\n'
                )

            pvd_file.write("  </Collection>\n")
            pvd_file.write("</VTKFile>\n")

        print(f"Wrote PVD file to {pvd_filename}")


# -------------------------------------------------------------------------------
# Main function
# -------------------------------------------------------------------------------


def main():
    """Main function to parse arguments and generate VTK files.

    Example usage (run from the build directory):
        python ../src/viz/generate_vtk_output.py \
            tests/regression_tests/TurbineInterfaceTest.IEA15/turbine_interface.nc \
            tests/regression_tests/TurbineInterfaceTest.IEA15/mesh_connectivity.yaml \
            --output_dir tests/regression_tests/TurbineInterfaceTest.IEA15/vtk_output \
            --start-timestep 0 \
            --end-timestep 5

    NOTE: Files are overwritten in the output directory if they already exist.
    """
    parser = argparse.ArgumentParser(
        description="Generate VTK files from Kynema NetCDF output"
    )

    # ------------------------------------------------------------
    # Input arguments
    # ------------------------------------------------------------

    # NetCDF input file -- required argument
    parser.add_argument(
        "netcdf_file",
        type=str,
        help="Path to Kynema NetCDF output file e.g. blade_interface.nc",
    )

    # Mesh connectivity file -- required argument
    parser.add_argument(
        "connectivity_file",
        type=str,
        help="Path to mesh connectivity YAML file e.g. mesh_connectivity.yaml",
    )

    # Output directory -- optional argument
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="vtk_output",
        help="Directory for writing the vtk output files",
    )

    # Timestep range to visualize -- optional arguments
    parser.add_argument(
        "--start-timestep",
        "-s",
        type=int,
        default=None,
        help="Starting timestep to visualize (default: 0)",
    )

    parser.add_argument(
        "--end-timestep",
        "-e",
        type=int,
        default=None,
        help="Ending timestep to visualize (default: last timestep)",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Generate VTK output
    # ------------------------------------------------------------

    # Issue warning before overwriting files
    if os.path.exists(args.output_dir):
        print(
            f"* Warning: {args.output_dir} already exists -- files may be overwritten"
        )
        os.makedirs(args.output_dir, exist_ok=True)

    # Create the VTK output object
    vtk_output = VTKOutput(args.netcdf_file, args.connectivity_file)

    # Generate visualization for timestep range or all timesteps
    if args.start_timestep is not None or args.end_timestep is not None:
        start = args.start_timestep if args.start_timestep is not None else 0
        end = (
            args.end_timestep
            if args.end_timestep is not None
            else vtk_output.num_timesteps - 1
        )

        # Validate range
        if start < 0 or start >= vtk_output.num_timesteps:
            raise ValueError(
                f"Start timestep {start} out of range (0-{vtk_output.num_timesteps-1})"
            )
        if end < 0 or end >= vtk_output.num_timesteps:
            raise ValueError(
                f"End timestep {end} out of range (0-{vtk_output.num_timesteps-1})"
            )
        if start > end:
            raise ValueError(
                f"Start timestep {start} cannot be greater than end timestep {end}"
            )

        # Generate visualization for the specified range
        print(f"Generating visualization for timesteps {start} to {end}")
        for timestep in range(start, end + 1):
            vtk_output.generate_visualization(timestep, args.output_dir)
    else:
        # Generate visualization for all timesteps
        print("Generating visualization for all timesteps")
        vtk_output.visualize_all_timesteps(args.output_dir)


if __name__ == "__main__":
    main()
