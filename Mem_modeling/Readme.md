cryo-EM membrane modeling
==============================


This is a tool to directly sample and calculate membrane curvature from cryo-EM density maps, as well as build lipids bilayers into your membrane density maps.

Features
--------------

With cryo-MEM you can:

- Sample your membrane density with PO4 molecules and convert it into pdb files, where you can check or visualize it in your own way.
- Derive 2D curvature contour maps.
- Build decent lipids pdb files.


Installation
--------------

cryo-MEM dependencies are mrcfile, numpy, matplotlib, pandas, coot (version: 0.9.x or 1.1), matplotlib and scipy.

These modules are available via `pip`:

```
pip install xxxx
```

Usage
--------------

This is a quick example on how to run cryo-MEM:

```Python

map_path = "~/example_data/" 
membrane_width = 55
output_path = script_dir + '/output/'

# Load and sample your density map into two PO4 molecules pdb files.
df_lower, df_upper = load_and_sample_maps(map_path, membrane_width,save_pdb=True) 
po4_lower, po4_upper = replace_into_po4(df_lower, leaflet="lower"), replace_into_po4(df_upper, leaflet="upper")

# To better represent the curavture of your membrane
rotated_lower, rotated_upper = rotated_and_visualize_curvature(df_lower, df_upper)
plot_contour_old(rotated_lower)

# If you want to improve your accuracy. Then, perform coot refine

mmbr_pdb = output_path+"po4.pdb"
protein_pdb = example_path+"mito_protein.pdb"
mmbr_id, map_id, protein_id = coot_input(mmbr_pdb, map_path, protein_pdb)

num_patches=10
overlap_factor=0.15
patch_ranges = calculate_patch_dimensions(num_patches, po4_upper, overlap_factor)
patch_mol_ids_coot = extract_atoms_in_patch(po4_upper, patch_ranges)

coot_init(map_id)
coot_refine(mmbr_id, patch_mol_ids_coot)
print(protein_id)
coot_refine_end(protein_id)

```
Function
--------------
# WarnOnce Class

A class to issue a warning about atoms out of grid boundaries. It is designed to issue a full warning message only once, and a generic message for subsequent occurrences.

## Attributes

- `msg` (str): The full warning message to be displayed on the first occurrence.
- `msg_multiple` (str): The generic warning message for subsequent occurrences.
- `warned` (bool): Indicates whether the full warning has been issued at least once.
- `warned_multiple` (bool): Indicates whether the generic warning has been issued at least once.

## Methods

### `__init__(self, msg, msg_multiple)`

Initializes the WarnOnce instance with specific messages.

#### Parameters

- `msg` (str): The full warning message to be displayed on the first occurrence.
- `msg_multiple` (str): The generic warning message for subsequent occurrences.

### `warn(self, *args)`

Issues the appropriate warning message based on the instance's state. It issues a full warning message on the first call and a generic warning on subsequent calls.

#### Parameters

- `*args`: Arguments to format the warning message.

## Examples

```python
warn_once_instance = WarnOnce("Full warning message.", "Generic warning.")
warn_once_instance.warn()
```

# load_density_map Function

Read density maps that should be gaussian-ed by Chimerax.

## Parameters

- `file_path` (str): The path of the density map.

## Returns

- `density_map` (numpy.ndarray): The density map data.
- `voxel_size` (numpy.ndarray): The size of the voxel unit. An array with (3,1) shape.
- `origin` (numpy.ndarray): The origin of the density map. It should be [0,0,0] in most cases. An array with (3,1) shape.
- `dimensions` (tuple): The dimensions of the map (x_dim, y_dim, z_dim).

## Examples

```python
file_path = "path/to/density_map.mrc"
density_map, voxel_size, origin, dimensions = load_density_map(file_path)
```

# identify_surfaces Function

Binarizes the density map. For voxels where density values are larger than the threshold, they are recorded as 1; otherwise, they are recorded as 0.

## Parameters

- `density_map` (numpy.ndarray): The density map, which is the first output of `load_density_map()`.
- `threshold` (float, optional): The threshold for binarization. Considered as the level parameter shown in the volume viewer in ChimeraX. You can pick a threshold from the sidebar when visualizing your density map in ChimeraX. Default is 0.02.

## Returns

- `lower_surface` (numpy.ndarray): The lower surface of the binarized density map. The shape is the same as the density map array.
- `upper_surface` (numpy.ndarray): The upper surface of the binarized density map. The shape is the same as the density map array.

## Examples

```python
density_map, _, _, _ = load_density_map("path/to/density_map.mrc")
lower_surface, upper_surface = identify_surfaces(density_map, 0.02)
```

# surface_to_coords_without_meditated Function

Convert the binarized map's array into a 3D array with [x, y, z] coordinates. This function is used to convert surface data into spatial coordinates.

## Parameters

- `surface` (numpy.ndarray): An array representing a surface, typically obtained from the output of `identify_surfaces()`.

## Returns

- `coords` (numpy.ndarray): A 3D array where each element is a coordinate in the form [x, y, z]. The z-coordinate represents the density value of that voxel. NaN values are used to indicate missing or undefined data points.

## Examples

```python
density_map, _, _, _ = load_density_map("path/to/density_map.mrc")
lower_surface, upper_surface = identify_surfaces(density_map, 0.02)
lower_coords = surface_to_coords_without_meditated(lower_surface)
```

# physical_coords Function

Converts the [x, y, z] coordinates into coordinates with physical meaning, representing the coordinates in real space.

## Parameters

- `coords` (numpy.ndarray): An array of [x, y, z] coordinates, typically obtained from `surface_to_coords_without_meditated()`. This array represents coordinates in the map's grid space.
- `voxel_size` (dict): A dictionary containing the voxel size for each dimension ('x', 'y', 'z'). Each entry should be a float representing the size of a voxel in that dimension.
- `origin` (dict): A dictionary specifying the origin of the density map in real space. It typically contains 'x', 'y', and 'z' keys with float values.

## Returns

- `physical_coords` (numpy.ndarray): A numpy array of [x, y, z] coordinates that represent the positions in real space.

## Examples

```python
coords = surface_to_coords_without_meditated(surface)
voxel_size = {'x': 1.0, 'y': 1.0, 'z': 1.0}
origin = {'x': 0.0, 'y': 0.0, 'z': 0.0}
real_space_coords = physical_coords(coords, voxel_size, origin)
```

# subsample_coordinates Function

Subsamples coordinates obtained from a density map at a specified rate. This is useful for reducing the number of points for further processing or visualization.

## Parameters

- `coords` (numpy.ndarray): A 3D array of coordinates, typically the output from `physical_coords()`. The shape of the array should be (x_dim, y_dim, 3), representing [x, y, z] coordinates.
- `rate` (float): The subsample rate, specifying the step size for slicing the coordinates. A rate of 1.0 means no subsampling, while a rate of 0.5 will select every second coordinate.

## Returns

- `subsampled_coords` (numpy.ndarray): A subsampled array of coordinates, reduced according to the specified rate.

## Examples

```python
coords = physical_coords(some_surface, voxel_size, origin)
subsampled_coords = subsample_coordinates(coords, 0.5)
```

# surface_to_coords Function

Converts the [x, y, z] array into PDB coordinates. This function combines several operations: `surface_to_coords_without_meditated`, `physical_coords`, and optionally `subsample_coordinates`.

## Parameters

- `lower_surface` (numpy.ndarray): The lower surface array obtained from `identify_surfaces()`.
- `upper_surface` (numpy.ndarray): The upper surface array obtained from `identify_surfaces()`.
- `membrane_width` (float): The distance between the upper and lower leaflets of the membrane. Used to correct errors from direct density map sampling.
- `voxel_size` (dict): A dictionary containing the voxel sizes for each dimension ('x', 'y', 'z').
- `origin` (dict): A dictionary specifying the origin of the density map in real space.
- `subsample` (bool, optional): If True, the coordinates are subsampled. Default is True.

## Returns

- `coordinates` (tuple of numpy.ndarray): A tuple containing the coordinates arrays. If `subsample` is True, it returns four arrays: meditated lower and upper coordinates, and subsampled lower and upper coordinates. If `subsample` is False, it returns only the meditated lower and upper coordinates.

## Examples

```python
lower_surface, upper_surface = identify_surfaces(density_map, threshold)
voxel_size = {'x': 1.0, 'y': 1.0, 'z': 1.0}
origin = {'x': 0.0, 'y': 0.0, 'z': 0.0}
coords = surface_to_coords(lower_surface, upper_surface, 55, voxel_size, origin, subsample=True)
```

# coords_to_dataframe Function

Converts coordinates of upper and lower leaflets into a Pandas DataFrame.

## Parameters

- `coords` (numpy.ndarray): A 3D array of [x, y, z] coordinates, typically the output from functions like `surface_to_coords_without_meditated` or `subsample_coordinates`.

## Returns

- `dataframe` (pandas.DataFrame): A DataFrame containing columns ['atom_id', 'res_id', 'chain_id', 'atom_type', 'coordinates']. Each row represents a point in the membrane with its corresponding attributes.

## Examples

```python
lower_surface, upper_surface = identify_surfaces(density_map, threshold)
lower_coords = surface_to_coords_without_meditated(lower_surface)
df_lower = coords_to_dataframe(lower_coords)
```

# write_pdb Function

Outputs membrane coordinates to a PDB file. The format starts with 'HETATM' and includes only one chain. The residue IDs are encoded in the hybrid36 format.

## Parameters

- `coords` (numpy.ndarray): A 3D array of coordinates. Each coordinate corresponds to a point in the membrane. These coordinates should be in the format suitable for a PDB file.
- `output_file` (str): Path to the output PDB file where the coordinates will be saved.

## Notes

The PDB file format is used for the storage and exchange of 3D structures of large biological molecules, such as proteins and nucleic acids. This function is specifically tailored for membrane structure data.

## Examples

```python
coords = surface_to_coords(lower_surface, upper_surface, membrane_width, voxel_size, origin)
write_pdb(coords, 'path/to/output.pdb')
```

# rotated_points Function

Rotates the membrane coordinates to depict the true curvature of the membrane. This involves finding the fitting plane and then rotating all points to make them parallel with the XoY plane. The rotation method uses quaternions to avoid the gimbal lock problem.

## Parameters

- `lower_leaflet` (pandas.DataFrame): A DataFrame containing the membrane coordinates for the lower leaflet. It should include 'coordinates' column with [x, y, z] data.
- `upper_leaflet` (pandas.DataFrame): A DataFrame containing the membrane coordinates for the upper leaflet. Similar structure to the `lower_leaflet`.

## Returns

- tuple: A tuple containing two elements, each is a list: [Rotation matrix, rotation vector, rotated points]. The first element corresponds to the upper leaflet and the second to the lower leaflet.

## Examples

```python
df_lower, df_upper = load_and_sample_maps(map_path, membrane_width, save_pdb=True)
rotated_lower, rotated_upper = rotated_points(df_lower, df_upper)
```
# fitting_plane Function

Calculates the coefficients of the best-fitting plane for a given set of 3D coordinates. This function is typically used in structural analysis to find the plane that best represents a set of points in 3D space.

## Parameters

- `coords` (numpy.ndarray): A 2D array where each row represents a 3D coordinate (x, y, z).

## Returns

- `A_coeff` (float): The coefficient 'A' in the plane equation Ax + By + Cz + D = 0.
- `B_coeff` (float): The coefficient 'B' in the plane equation.
- `D_coeff` (float): The coefficient 'D' in the plane equation.

## Notes

- The function uses least squares fitting to calculate the plane coefficients. It does not require GPU acceleration and uses numpy for computations.

## Examples

```python
points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
A_coeff, B_coeff, D_coeff = fitting_plane(points)
```


# quaternion_to_matrix Function

Converts a quaternion into a rotation matrix. This is used in rotational transformations where quaternions provide a more robust solution against gimbal lock.

## Parameters

- `quaternion` (list or numpy.ndarray): A quaternion represented as [w, x, y, z], where w is the scalar component, and x, y, z are the vector components.

## Returns

- numpy.ndarray: A 3x3 rotation matrix derived from the input quaternion.

## Notes

Quaternions are a system of complex numbers used to represent rotations in three-dimensional space. They are particularly useful in applications like 3D graphics, aerospace, and robotics.

## Examples

```python
quaternion = [1, 0, 0, 0]  # Represents no rotation
rotation_matrix = quaternion_to_matrix(quaternion)
```

# rotate_operation Function

Rotates a point cloud to align its best-fitting plane with a target vector. This is used for aligning the membrane plane with a reference plane, typically the XoY plane.

## Parameters

- `plane_coeff` (list or numpy.ndarray): Coefficients [A_coeff, B_coeff, D_coeff] of the best-fitting plane to the point cloud.
- `target_vector` (list or numpy.ndarray): The target direction vector for the plane's normal after rotation. Typically [0, 0, -1] for alignment with the XoY plane.
- `points` (numpy.ndarray): The coordinates of the point cloud to be rotated, in the format of a 2D array.

## Returns

- list: A list containing the rotation matrix, the rotated normal vector, and the rotated points as numpy.ndarray.

## Notes

The rotation is performed using quaternion-based methods to ensure stable and smooth rotation without the gimbal lock issue.

## Examples

```python
plane_coeff = fitting_plane(points)
target_vector = [0, 0, -1]
rotated_result = rotate_operation(plane_coeff, target_vector, points)
```

# write_pdb_for_rotated Function

Writes a PDB file for the rotated points. This function is specifically designed to output rotated membrane coordinates into a standard PDB format.

## Parameters

- `coords` (numpy.ndarray): A 2D array of [x, y, z] coordinates, representing the rotated points of a membrane or molecular structure.
- `output_file` (str): The path to the output PDB file where the coordinates will be saved.

## Notes

The function formats each coordinate point according to the PDB file standards, using 'HETATM' records. It is tailored to handle the output from membrane rotation operations.

## Examples

```python
rotated_lower, _ = rotated_points(df_lower, df_upper)
write_pdb_for_rotated(rotated_lower[-1], 'path/to/rotated_lower_leaflet.pdb')
```

# load_and_sample_maps Function

Loads and samples density maps, processing them for both lower and upper leaflets. This function is responsible for handling the entire workflow of loading, processing, and optionally saving the sampled membrane data as PDB files.

## Parameters

- `map_path` (str): The absolute path to the density map file.
- `membrane_width` (float): The thickness of the membrane, used in processing the density map. Unit is typically in Angstroms.
- `save_pdb` (bool, optional): Whether to save the sampled points as PDB files. Default is True.
- `output_path` (str, optional): The output directory path where PDB files will be saved if `save_pdb` is True. Default is '/output/'.

## Returns

- `tuple` of `pandas.DataFrame`: A tuple containing two DataFrames. The first DataFrame corresponds to the lower leaflet, and the second to the upper leaflet. Each DataFrame contains columns ['atom_id', 'res_id', 'chain_id', 'atom_type', 'coordinates'].

## Examples

```python
map_path = "/path/to/density_map.mrc"
membrane_width = 55
df_lower, df_upper = load_and_sample_maps(map_path, membrane_width, save_pdb=True)
```

# rotated_and_visualize_curvature Function

Rotates and visualizes the curvature of the membrane based on the provided lower and upper leaflet data. It handles the rotation of points to align with a reference plane and optionally saves the rotated structure as PDB files.

## Parameters

- `df_lower` (pandas.DataFrame): A DataFrame containing the lower leaflet membrane coordinates. Expected to have a column named 'coordinates' with [x, y, z] data.
- `df_upper` (pandas.DataFrame): A DataFrame containing the upper leaflet membrane coordinates. Structured similarly to `df_lower`.
- `save_pdb` (bool, optional): Whether to save the rotated coordinates as PDB files. Default is False.
- `output_path` (str, optional): The directory path for saving output PDB files if `save_pdb` is set to True. Default is '/output/'.

## Returns

- `tuple`: A tuple containing two elements, each a list: [Rotation matrix, rotation vector, rotated points]. The first element corresponds to the rotated lower leaflet, and the second to the rotated upper leaflet.

## Examples

```python
df_lower, df_upper = load_and_sample_maps(map_path, membrane_width, save_pdb=True)
rotated_lower, rotated_upper = rotated_and_visualize_curvature(df_lower, df_upper, save_pdb=True)
```

# plot_contour_old Function

Uses the griddata method to smooth/regress the original image and visualize it with contour plotting. This function is specifically designed for visualizing the curvature of rotated membrane points.

## Parameters

- `rotated_points_list` (list): A list containing the rotated points as [Rotation matrix, rotation vector, [x, y, z]].
- `interpolate_method` (str, optional): The interpolation method to use. Options include 'linear', 'nearest', 'cubic', and 'bilinear'. Default is 'bilinear'.
- `level` (str, optional): The level of oversampling for contouring. Default is '2'.
- `my_color_map` (str, optional): The color palette to use for contouring. Default is 'RdBu'.
- `distance_threshold` (float, optional): The threshold distance to consider for contouring. Default is 6.5.

## Returns

- None: The function does not return any value but produces a contour plot.

## Examples

```python
rotated_lower, _ = rotated_points(df_lower, df_upper)
plot_contour_old(rotated_lower, interpolate_method="linear", level="3", my_color_map='viridis')
```

# gaussian_smoothing Function

Applies Gaussian smoothing to an image using a specified kernel size and standard deviation (sigma).
This function is primarily used for smoothing the surface data of the membrane.

## Parameters

- `image` (numpy.ndarray): The input image or 2D data array to be smoothed.
- `kernel_size` (int, optional): The size of the Gaussian kernel. Default is 5.
- `sigma` (float, optional): The standard deviation of the Gaussian kernel. Default is 1.0.

## Returns

- `numpy.ndarray`: The smoothed image as a 2D array.

## Notes

The function computes a Gaussian kernel, applies padding to the input image, and then performs convolution to achieve Gaussian smoothing. It handles NaN values appropriately, ensuring they do not affect the smoothing of neighboring data points.

## Examples

```python
image = np.random.rand(100, 100)  # Randomly generated image
smoothed_image = gaussian_smoothing(image, kernel_size=7, sigma=2.0)
```

# batched_cdist_min Function

Computes the minimum Euclidean distances between each point in array X and any point in array Y, 
in a batched manner to optimize memory usage.

## Parameters

- `X` (numpy.ndarray): A 2D array where each row represents a point in space (e.g., [x, y] or [x, y, z]).
- `Y` (numpy.ndarray): Another 2D array of points in the same space as X.
- `batch_size` (int, optional): The number of points from X to process in each batch. This helps to manage memory usage when dealing with large arrays. Default is 1000.

## Returns

- `numpy.ndarray`: A 1D array where each element is the minimum distance from a point in X to any point in Y.

## Examples

```python
X = np.array([[1, 2], [3, 4], [5, 6]])
Y = np.array([[7, 8], [9, 10]])
min_distances = batched_cdist_min(X, Y)
```

# replace_into_po4 Function

Replaces membrane sampling points with PO4 (phosphate) molecules. This is useful for creating a more realistic representation of the membrane in molecular simulations or visualizations.

## Parameters

- `target_positions_df` (pandas.DataFrame): A DataFrame containing membrane sampling points information. It should include a 'coordinates' column with [x, y, z] data.
- `save_pdb` (bool, optional): If True, the resulting PO4 coordinates will be saved to a PDB file. Default is False.
- `leaflet` (str, optional): Specifies which leaflet (upper or lower) the PDB file represents, used in the output file name. Relevant only if `save_pdb` is True.

## Returns

- `pandas.DataFrame`: A DataFrame similar in structure to the input but with PO4 molecules replacing the original points. The columns include ['atom_id', 'res_id', 'chain_id', 'atom_type', 'coordinates'].

## Examples

```python
df_lower = coords_to_dataframe(lower_coords)
po4_lower = replace_into_po4(df_lower, save_pdb=True, leaflet="lower")
```

# calculate_patch_dimensions Function

Calculates the dimensions of patches for a given number of patches and points array, considering an overlap factor. This is typically used in segmenting a membrane or surface for detailed analysis or processing in patches.

## Parameters

- `num_patches` (int): The total number of patches into which the points array is to be divided.
- `points` (pandas.DataFrame): A DataFrame containing points data, specifically membrane coordinates. The DataFrame should include a 'coordinates' column with [x, y, z] data.
- `overlap_factor` (float, optional): The factor by which adjacent patches should overlap each other, expressed as a fraction of the patch size. Default is 0.15 (15%).

## Returns

- `dict`: A dictionary where keys are patch IDs and values are tuples specifying the x and y range of each patch as ((x_start, x_end), (y_start, y_end)).

## Examples

```python
points_df = coords_to_dataframe(some_coords)
patch_ranges = calculate_patch_dimensions(10, points_df)
```

# extract_atoms_in_patch Function

Extracts atoms within specified patch ranges from a given points dataset. This function is useful in segmenting a membrane or surface data into smaller patches for detailed analysis.

## Parameters

- `points` (pandas.DataFrame): A DataFrame containing points data. It should include 'coordinates' as well as other molecular attributes.
- `patch_ranges` (dict): A dictionary with patch IDs as keys and tuples specifying the x and y range of each patch as values, formatted as ((x_start, x_end), (y_start, y_end)).

## Returns

- `dict`: A dictionary with patch IDs as keys and lists of molecular IDs (residue IDs) within each patch as values.

## Examples

```python
points_df = coords_to_dataframe(some_coords)
patch_ranges = calculate_patch_dimensions(10, points_df)
patch_atoms = extract_atoms_in_patch(points_df, patch_ranges)
```

# coot_input Function

Prepares the input for the COOT (Crystallographic Object-Oriented Toolkit) program by reading PDB and map files. This function is designed to be used in a structural biology context for model building and validation.

## Parameters

- `mmbr_pdb` (str): The file path to the membrane PDB file.
- `density_map` (str): The file path to the density map, typically in CCP4 or MRC format.
- `protein_pdb` (str): The file path to the protein PDB file.

## Returns

- `tuple`: A tuple containing the molecule IDs for the membrane, density map, and protein, respectively, as used in COOT.

## Examples

```python
mmbr_pdb = 'path/to/membrane.pdb'
density_map = 'path/to/density.map'
protein_pdb = 'path/to/protein.pdb'
mmbr_id, map_id, protein_id = coot_input(mmbr_pdb, density_map, protein_pdb)
```

# coot_init Function

Initializes the refinement process in COOT by setting the reference density map. This function is part of the workflow for molecular modeling and refinement in structural biology.

## Parameters

- `map_id` (int): The molecule ID of the density map as read by COOT. This map will be used as the reference during refinement.

## Notes

- This function should be called after loading all the necessary molecular and map files into COOT using `coot_input`. It sets the refinement map and configures the refinement parameters.

## Examples

```python
mmbr_id, map_id, protein_id = coot_input(mmbr_pdb, density_map, protein_pdb)
coot_init(map_id)
```

# coot_refine Function

Performs rigid body refinement on specified molecular patches in COOT. This function is used for refining membrane or protein structures against a density map.

## Parameters

- `mmbr_id` (int): The molecule ID of the membrane or protein structure as loaded in COOT.
- `patch_mol_ids_coot` (dict): A dictionary with patch IDs as keys and lists of molecular IDs (residue IDs) for each patch as values. These IDs are formatted for COOT's refinement process.

## Notes

- Before calling this function, ensure that the necessary structures and density maps are loaded into COOT, and the refinement map is set using `coot_init`.

## Examples

```python
mmbr_id, map_id, protein_id = coot_input(mmbr_pdb, density_map, protein_pdb)
coot_init(map_id)
patch_mol_ids_coot = extract_atoms_in_patch(points_df, patch_ranges)
coot_refine(mmbr_id, patch_mol_ids_coot)
```

# coot_refine_end Function

Finalizes the refinement process in COOT, specifically converting the graphics representation of the protein to C-alpha (CA) representation. This function is part of the workflow for molecular modeling and refinement in structural biology.

## Parameters

- `protein_id` (int): The molecule ID of the protein structure as loaded in COOT.

## Notes

- This function should be called after performing all necessary refinements using `coot_refine`. It changes the visual representation of the protein to focus on the alpha carbon atoms, which is a common practice in protein structure analysis.

## Examples

```python
protein_id = coot.read_pdb(protein_pdb_file)
coot_refine_end(protein_id)
```

