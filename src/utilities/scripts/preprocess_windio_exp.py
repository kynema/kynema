import numpy as np
from sympy import false
import yaml
from scipy.interpolate import PchipInterpolator, make_splprep

# Change these settings to modify the interpolation
N_AIRFOIL_POINTS = 60
N_POLAR_POINTS = 180

INPUT_FILE = "NREL-5MW.yaml"
OUTPUT_FILE = "NREL-5MW-aero.yaml"
N_AERO_NODES = 50

# Spanwise locations of aerodynamic nodes (from 0 to 1)
aero_locs = np.linspace(0.0, 1.0, N_AERO_NODES)


def normalize_airfoil_coordinates(coords, n_points=N_AIRFOIL_POINTS):
    """Interpolate airfoil coordinates so they have the specified number of points along upper and lower surfaces."""

    # Get index of the minimum x-coordinate to split the airfoil into upper and lower parts
    i_min_x = np.argmin(coords["x"])

    # Split the coordinates into upper and lower parts
    upper_coords = np.array([coords["x"][i_min_x:], coords["y"][i_min_x:]])
    lower_coords = np.array([coords["x"][: i_min_x + 1], coords["y"][: i_min_x + 1]])

    # Sort the upper and lower coordinates by x-coordinate
    upper_coords = upper_coords[:, np.argsort(upper_coords[0])]
    lower_coords = lower_coords[:, np.argsort(lower_coords[0])]

    # Ensure that the upper and lower coordinates are ordered correctly
    if upper_coords[1][1] < lower_coords[1][1]:
        upper_coords, lower_coords = lower_coords, upper_coords

    # Create a cubic parametric spline representation of the upper and lower airfoil shapes
    spl_upper = make_splprep(upper_coords, s=0, k=3)[0]
    spl_lower = make_splprep(lower_coords, s=0, k=3)[0]

    # Evaluate the spline at equally spaced points along the airfoil
    upper_coords = spl_upper(np.linspace(0, 1, n_points // 2 + 1))
    lower_coords = spl_lower(np.linspace(0, 1, n_points // 2))

    # Merge upper and lower coordinates into a single array and update the original coords dictionary
    coords_x = np.concatenate((upper_coords[0][::-1], lower_coords[0][1:]))
    coords_y = np.concatenate((upper_coords[1][::-1], lower_coords[1][1:]))

    return (coords_x, coords_y)


data = yaml.load(open(INPUT_FILE), Loader=yaml.FullLoader)
blade = data["components"]["blade"]["outer_shape"]

# Get airfoils in span-wise order
af_map = {airfoil["name"]: airfoil for airfoil in data["airfoils"]}
base_airfoils = [af_map[af["name"]] for af in blade["airfoils"]]

# Extract locations of base airfoils along the span
af_grid = [af["spanwise_position"] for af in blade["airfoils"]]

# Interpolate blade data to aerodynamic node locations
blade_data = {
    p: PchipInterpolator(af_grid, [af[p] for af in base_airfoils])(aero_locs).tolist()
    for p in ["aerodynamic_center", "rthick"]
}

# Set the spanwise position
blade_data["spanwise_position"] = aero_locs

# Use PCHIP interpolation to get values at aerodynamic node positions
for param in ["twist", "chord", "section_offset_x", "section_offset_y"]:
    if param not in blade:
        blade_data[param] = np.zeros(N_AERO_NODES)
    else:
        p = blade[param]
        blade_data[param] = PchipInterpolator(p["grid"], p["values"])(aero_locs)

# Loop through airfoils
polar_grid_points = []
for af in base_airfoils:
    for cfg in af["polars"]:
        for re_set in cfg["re_sets"]:
            for p in ["cl", "cd", "cm"]:
                polar_grid_points.extend(re_set[p]["grid"])

# Divide polar grid points into quantiles to create a uniform grid
quantiles = np.quantile(polar_grid_points, [0, 0.20, 0.40, 0.60, 0.80, 1])
polar_grid = np.concatenate(
    [
        np.linspace(quantiles[0], quantiles[1], N_POLAR_POINTS // 5, endpoint=False),
        np.linspace(quantiles[1], quantiles[2], N_POLAR_POINTS // 5, endpoint=False),
        np.linspace(quantiles[2], quantiles[3], N_POLAR_POINTS // 5, endpoint=False),
        np.linspace(quantiles[3], quantiles[4], N_POLAR_POINTS // 5, endpoint=False),
        np.linspace(quantiles[4], quantiles[5], N_POLAR_POINTS // 5, endpoint=True),
    ]
)

# Initialize polar data structure to store data from base airfoils
# This will be used for the interpolation of polar data
polar_data = {
    p: np.zeros([len(base_airfoils), N_POLAR_POINTS]) for p in ["cl", "cd", "cm"]
}

# Initialize coordinate data structure to store airfoil coordinates
coordinate_data = {
    p: np.zeros([len(base_airfoils), N_AIRFOIL_POINTS]) for p in ["x", "y"]
}

# Loop through airfoils
for i, af in enumerate(base_airfoils):

    # Normalize the airfoil coordinates so they have a consistent number of points
    coordinate_data["x"][i, :], coordinate_data["y"][i, :] = (
        normalize_airfoil_coordinates(af["coordinates"])
    )

    # Store the polar data in the polar_data structure
    # This assumes that there is one polar configuration and one re_set per airfoil
    # TODO: Handle multiple configurations and re_sets if necessary
    for param in polar_data.keys():
        p = af["polars"][0]["re_sets"][0][param]
        polar_data[param][i, :] = np.interp(polar_grid, p["grid"], p["values"])

    # Apply 3D correction to the polars
    # p = Polar(alpha=polar_grid, cl=polar_data["cl"][i, :], cd=polar_data["cd"][i, :], cm=polar_data["cm"][i, :])

    # p.correction3D(r_over_R=blade_data["spanwise_position"][i], Re=1e6, Mach=0.0)


# Interpolate the polar data to the aerodynamic node positions using PCHIP
polar_data_interp = {
    p: PchipInterpolator(af_grid, v)(aero_locs) for p, v in polar_data.items()
}

# Interpolate the coordinates to the aerodynamic node positions
coordinate_data_interp = {
    p: PchipInterpolator(af_grid, coordinate_data[p])(aero_locs) for p in ["x", "y"]
}

# Update file with new airfoil data
data["airfoils"] = [
    {
        "name": f"AF{i+1:02d}",
        "aerodynamic_center": blade_data["aerodynamic_center"][i],
        "coordinates": {
            "x": coordinate_data_interp["x"][i, :].tolist(),
            "y": coordinate_data_interp["y"][i, :].tolist(),
        },
        "rthick": float(blade_data["rthick"][i]),
        "twist": float(blade_data["twist"][i]),
        "chord": float(blade_data["chord"][i]),
        "section_offset_x": float(blade_data["section_offset_x"][i]),
        "section_offset_y": float(blade_data["section_offset_y"][i]),
        "spanwise_position": float(blade_data["spanwise_position"][i]),
        "polars": [
            {
                "configuration": "default",
                "re_sets": [
                    {
                        p: {
                            "grid": polar_grid.tolist(),
                            "values": polar_data_interp[p][i, :].tolist(),
                        }
                        for p in polar_data_interp.keys()
                    }
                ],
            }
        ],
    }
    for i in range(N_AERO_NODES)
]

# Update the blade airfoil list
blade["airfoils"] = [
    {
        "name": af["name"],
        "spanwise_position": float(loc),
        "configuration": ["default"],
        "weight": [1.0],
    }
    for af, loc in zip(data["airfoils"], aero_locs)
]

# Update the blade data with the interpolated values
for param in ["chord", "section_offset_x", "section_offset_y"]:
    blade[param] = {
        "grid": aero_locs.tolist(),
        "values": blade_data[param].tolist(),
    }

# ------------------------------------------------------------------------------
# Missing blade elastic_properties fields
# ------------------------------------------------------------------------------

ep = data["components"]["blade"]["structure"]["elastic_properties"]
for v in ["cm_x", "cm_y", "i_cp"]:
    if v not in ep["inertia_matrix"]:
        ep["inertia_matrix"][v] = [0.0] * len(ep["inertia_matrix"]["grid"])

for v in [
    "K12",
    "K13",
    "K14",
    "K15",
    "K16",
    "K23",
    "K24",
    "K25",
    "K26",
    "K34",
    "K35",
    "K36",
    "K45",
    "K46",
    "K56",
]:
    if v not in ep["stiffness_matrix"]:
        ep["stiffness_matrix"][v] = [0.0] * len(ep["stiffness_matrix"]["grid"])


# ------------------------------------------------------------------------------
# Common blade reference axis grid
# ------------------------------------------------------------------------------

grids = [data["components"]["blade"]["reference_axis"][c]["grid"] for c in "xyz"]
new_grid = sorted(list(set([item for sublist in grids for item in sublist])))

for c in "xyz":
    values = data["components"]["blade"]["reference_axis"][c]["values"]
    grid = data["components"]["blade"]["reference_axis"][c]["grid"]
    data["components"]["blade"]["reference_axis"][c]["grid"] = new_grid
    data["components"]["blade"]["reference_axis"][c]["values"] = PchipInterpolator(
        grid, values
    )(new_grid).tolist()

# ------------------------------------------------------------------------------
# Common tower reference axis grid
# ------------------------------------------------------------------------------

grids = [data["components"]["tower"]["reference_axis"][c]["grid"] for c in "xyz"]
new_grid = sorted(list(set([item for sublist in grids for item in sublist])))

for c in "xyz":
    values = data["components"]["tower"]["reference_axis"][c]["values"]
    grid = data["components"]["tower"]["reference_axis"][c]["grid"]
    data["components"]["tower"]["reference_axis"][c]["grid"] = new_grid
    data["components"]["tower"]["reference_axis"][c]["values"] = PchipInterpolator(
        grid, values
    )(new_grid).tolist()

# -------------------------------------------------------------------------------
# Output new file
# -------------------------------------------------------------------------------

# Save the processed data to a new YAML file
yaml.dump(data, open(OUTPUT_FILE, "w"))
