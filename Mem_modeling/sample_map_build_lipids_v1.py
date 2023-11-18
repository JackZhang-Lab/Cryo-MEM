import os
import sys
import mrcfile
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import warnings
import coot
# import cupy as cp # uncomment it when you have GPU
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist

# Obtain current paths
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
module_path = os.path.join(script_dir, 'modules/')
example_path = os.path.join(script_dir, 'example_data/')

sys.path.append(module_path)
# Encode residue ID in pdb files into the hybrid_36 format
import hybrid_36

class WarnOnce:
    """
    Class to warn atoms out of grid boundaries only once with full message.
    After the first ocurrance, message will be generic.
    """
    def __init__(self, msg, msg_multiple) -> None:
        self.msg = msg
        self.warned = False
        self.warned_multiple = False
        self.msg_multiple = msg_multiple

    def warn(self, *args):
        if not self.warned:
            self.warned = True
            warnings.warn(self.msg.format(*args))
            # logger.warning(self.msg.format(*args))
        elif not self.warned_multiple:
            warnings.warn(self.msg_multiple)

def load_density_map(file_path):
    '''
    Read density maps. The density map should be gaussian-ed by Chimerax.

    Parameters
    ----------
        file_path: the path of your density map. 
    Returns
    ----------
        density_map
        voxel_size: the size of voxel unit. An array with (3,1) shape
        origin: The origin of the density map. It should be [0,0,0] in most of the time. An array with (3,1) shape
        (x_dim, y_dim, z_dim): Dimension of the map.
    '''
    with mrcfile.open(file_path, mode='r') as mrc:
        density_map = mrc.data
        voxel_size = mrc.voxel_size 
        origin = mrc.header['origin']
        x_dim, y_dim, z_dim = density_map.shape
    return density_map, voxel_size, origin, (x_dim, y_dim, z_dim)

def identify_surfaces(density_map, threshold: float =0.02):
    '''
    Binarized the density map. For the voxels where density values larger than threshold record as 1, the other record as 0.

    Parameters
    ----------
        density map: This is the 1st output of load_density_map()
        threshold: The threshold of binarization. Considered as the level parameters shows in volume viewer in ChimeraX.
                    You can pick a threshold from the side bar when visualizing your density map in ChimeraX.
                    float
    Returns
    ----------
        lower_surface, upper_surface: . numpy array, shape is same as the (array of density map)

    '''
    binary_map = np.where(density_map >= threshold, 1, 0)
    # Initialization
    lower_surface = np.full((density_map.shape[1], density_map.shape[2]), np.nan, dtype=float)
    upper_surface = np.full((density_map.shape[1], density_map.shape[2]), np.nan, dtype=float)
    
    # for y in range(density_map.shape[1]):
    #     for z in range(density_map.shape[2]):
    #         mask = binary_map[:, y, z] == 1
    #         if mask.sum() > 0:
    #             lower_surface[y, z] = np.argmax(mask)
    #             upper_surface[y, z] = len(mask) - 1 - np.argmax(mask[::-1])
    #         else:
    #             lower_surface[y, z] = np.nan  
    #             upper_surface[y, z] = np.nan
    # return lower_surface, upper_surface

    for x in range(density_map.shape[1]):
        for y in range(density_map.shape[2]):
            mask = binary_map[x, y, :] == 1
            if mask.sum() > 0:
                lower_surface[x, y] = np.argmax(mask)
                upper_surface[x, y] = len(mask) - 1 - np.argmax(mask[::-1])
            else:
                lower_surface[x, y] = np.nan  
                upper_surface[x, y] = np.nan
    return lower_surface, upper_surface


def surface_to_coords_without_meditated(surface):
    '''
    Convert the binarized map's array into [x,y,z] array
    Parameters
    ----------
        Output from the identify_surfaces()
    Returns
    ----------
        coords: [x,y,z] array. Z represents the density value of that voxel.
    
    '''
    coords = np.zeros((surface.shape[0], surface.shape[1], 3), dtype=float)
    for i in range(surface.shape[0]):
        for j in range(surface.shape[1]):
            coords[i, j, 0] = i
            coords[i, j, 1] = j
            coords[i, j, 2] = surface[i, j]
            if np.isnan(surface[i, j]):
                coords[i, j, :] = [np.nan, np.nan, np.nan]
    return coords

def physical_coords(coords, voxel_size, origin):
    '''
    Convert the [x,y,z] into coordinates with physcial meaning. In other words, making the [x,y,z] coordinates represent the coords in real space.
    Parameters
    ----------
        coords: Output from the surface_to_coords_without_meditated(). array, shape is (x_dim, y_dim, 3)
        voxel_size, origin: Output 
    Returns
    ----------
        physical_coords: [x,y,z] array. array, shape is (x_dim, y_dim, 3)
    
    '''

    # print("physical coords input shape: ", coords.shape)
    voxel_size_array = np.array([voxel_size['x'], voxel_size['y'], voxel_size['z']])
    origin_array = np.array([origin['x'], origin['y'], origin['z']])
    physical_coords = coords.copy()
    for i in range(3):  # Loop over the x, y, z dimensions
        physical_coords[..., i] = physical_coords[..., i] * voxel_size_array[i] + origin_array[i]
    # print("physical coords shape at the output point", physical_coords.shape)
    return physical_coords

def subsample_coordinates(coords, rate: float):
    '''
    Subsample your coordinates which obtained from density map.
    That will be useful when combined with coot refine. Highly recommend using this coordinates as the input of coot_input()
    Parameters
    ----------
        coords: The 1st output from the physical_coords()
        rate: subsample rate. the step size of slice the coords.
    Returns
    ----------
        coords: [x,y,z] array.
    
    '''

    return coords[::rate, ::rate, :]

def surface_to_coords(lower_surface, upper_surface, membrane_width, voxel_size, origin, subsample=True):
    '''
    Convert the [x,y,z] array into pdb coordinates. 
    This function actually combines the surface_to_coords_without_meditated(), physical_coords(), subsample_coordinates() (denpend on your need)
    Parameters
    ----------
        lower_surface: upper_surface: Output from the identify_surfaces()
        membrane_width: the distance between upper leaflet and lower leaflet. Rectify the errors from directly sampling density map.
        voxel_size, origin: input for physical_coords(). refer to physical_coords()
        subsample: True or False
    Returns
    ----------
        coords: [x,y,z] array. This one can be used for output pdb
    
    '''

    membrane_width = float(membrane_width)
    subsample_rate = 5

    lower_coords = surface_to_coords_without_meditated(lower_surface)
    upper_coords = surface_to_coords_without_meditated(upper_surface)
    middle_coords = np.divide((lower_coords + upper_coords), 2)

    middle_coords_physical = physical_coords(middle_coords, voxel_size, origin)
    lower_coords_meditated = middle_coords_physical.copy()
    upper_coords_meditated = middle_coords_physical.copy()

    lower_coords_meditated[..., 2] -= (membrane_width / 2)
    upper_coords_meditated[..., 2] += (membrane_width / 2)
    if subsample == False:
        return lower_coords_meditated, upper_coords_meditated
    else:
        # Subsample the coordinates
        subsampled_lower_coords = subsample_coordinates(lower_coords_meditated, subsample_rate)
        subsampled_upper_coords = subsample_coordinates(upper_coords_meditated, subsample_rate)
        return lower_coords_meditated, upper_coords_meditated, subsampled_lower_coords, subsampled_upper_coords

def coords_to_dataframe(coords):
    '''
    Convert coords of upper_leaflet and lower_leaflet into relatively dataframe
    Parameters
    -----------
        coords, [x, y, z]
    
    Returns
    -----------
        dataframe ["atom_id", "res_id", "chain_id", "atom_type", "coordinates"]
        df = pd.DataFrame({
        "atom_id": atom_id_list,
        "res_id": res_id_list,
        "chain_id": ["A"] * len(atom_id_list), 
        "atom_type": ["P"] * len(atom_id_list),
        "coordinates": coordinates_list
    })

    '''
    atom_id_list = []
    res_id_list = []
    coordinates_list = []

    atom_counter = 1
    for i in range(coords.shape[0]):
        for j in range(coords.shape[1]):
            if not np.isnan(coords[i, j]).any() and not np.array_equal(coords[i, j], [0,0,0]):
                x, y, z = coords[i,j][0], coords[i,j][1], coords[i,j][2]
                x = round(x, 3)
                y = round(y, 3)
                z = round(z, 3)
                # Encode atom_counter in hybrid-36 format for use in columns 23-26
                encoded_res_id = hybrid_36.hy36encode(4, atom_counter)
                atom_id_list.append(atom_counter)
                res_id_list.append(encoded_res_id)
                coordinates_list.append([y, x, z])
                atom_counter += 1

    df = pd.DataFrame({
        "atom_id": atom_id_list,
        "res_id": res_id_list,
        "chain_id": ["A"] * len(atom_id_list), 
        "atom_type": ["P"] * len(atom_id_list),
        "coordinates": coordinates_list
    })
    return df

def write_pdb(coords, output_file):
    '''
    Output your membrane coordinates into pdb file.
    Start with HETATM.
    Only one chain in this pdb.
    Encoded in the hybrid36 way.
    '''
    lines = []
    atom_counter = 1
    for i in range(coords.shape[0]):
        for j in range(coords.shape[1]):
            if not np.isnan(coords[i, j]).any() and not np.array_equal(coords[i, j], [0,0,0]):
                x, y, z = coords[i,j][0], coords[i,j][1], coords[i,j][2]
                encoded_res_id = hybrid_36.hy36encode(4, atom_counter)
                # Write PDB format line
                lines.append(
                    f"ATOM  {atom_counter:5d}  P   PO4 A{encoded_res_id:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           P\n"
                )
                atom_counter += 1

    with open(output_file, 'w') as file:
        for line in lines:
            file.write(line)
        file.write("END")


def rotated_points(lower_leaflet, upper_leaflet):
    '''
    Rotate the membrane in order to depict the true curvature of the membrane.
        Start from finding the fitting plane. Then rotate all the points to make them parallel with the XoY plane. 
        Rotation method is querternion. Avoiding auto-lock problem
    Parameters
    ------------
        the dataframe which cotains membrane coordinates.
    
    Returns
        A list: [Rotation matrix, rotation vector, rotated points]
    '''
    p_in_lower_leaflet = lower_leaflet["coordinates"].to_numpy()
    p_in_upper_leaflet = upper_leaflet["coordinates"].to_numpy()
    p_in_lower_leaflet = np.stack(p_in_lower_leaflet)
    p_in_upper_leaflet = np.stack(p_in_upper_leaflet)

    def fitting_plane(coords):
        x = coords[:,0]
        y = coords[:,1]
        z = coords[:,2]
        # print(x)
        # Load numpy array into cupy array
        # x_gpu = cp.array(x)
        # y_gpu = cp.array(y)
        # z_gpu = cp.array(z)

        # Setup the coefficient matrix A
        # A_gpu = cp.vstack((x_gpu, y_gpu, cp.ones_like(x_gpu))).T
        A = np.vstack((x, y, np.ones_like(x))).T

        # coeff_gpu, _, _, _ = cp.linalg.lstsq(A_gpu, z_gpu, rcond=-0.01)
        coeff, _, _, _ = np.linalg.lstsq(A, z, rcond=-0.01)

        # Extract the coefficients A, B, and D
        # coeff = coeff_gpu.get()
        [A_coeff, B_coeff, D_coeff] = coeff
        return A_coeff, B_coeff, D_coeff
    
    # return upper_leaflet_plane, lower_leaflet_plane
    # The normal to the plane is given by [A_coeff, B_coeff, -1]

    def quaternion_to_matrix(quaternion):
        '''
        Convert a quaternion into a rotation matrix.
        quaternion q = (q0, q1, q2, q3)
        ''' 
        w, x, y, z = quaternion
        return np.array([
            [1 - 2*y*y - 2*z*z, 2*x*y - 2*z*w, 2*x*z + 2*y*w],
            [2*x*y + 2*z*w, 1 - 2*x*x - 2*z*z, 2*y*z - 2*x*w],
            [2*x*z - 2*y*w, 2*y*z + 2*x*w, 1 - 2*x*x - 2*y*y]
        ])

    target_vector = [0, 0, -1]
    def rotate_operation(plane_coeff, target_vector, points): 
        '''
        Function:
            Rotate the points of cloud following the rotation of its best fitting plane
        Input:
            plane_coeff: a list which contains three elements, [A_coeff, B_coeff, D_coeff]
            target_vector: [0, 0, -1]
            
        '''
        # v is represent the normal of the plane of coords
        v = np.array([plane_coeff[0], plane_coeff[1], -1])
        v = v / np.linalg.norm(v)
        t = target_vector / np.linalg.norm(target_vector)
        
        rotation_axis = np.cross(v, t)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        
        cos_theta = np.dot(v, t)
        theta = np.arccos(cos_theta) / 2.0
        
        w = np.cos(theta)
        x = rotation_axis[0] * np.sin(theta)
        y = rotation_axis[1] * np.sin(theta)
        z = rotation_axis[2] * np.sin(theta)
        
        R = quaternion_to_matrix([w, x, y, z])
        # points_rotated = np.dot(points, R)
        points_rotated = np.dot(R, points.T).T
        v_rotated = np.dot(R, v)
        # print(v.shape)
        rotated_result = [R, v_rotated, points_rotated]
        return rotated_result

    lower_leaflet_plane = fitting_plane(p_in_lower_leaflet)
    upper_leaflet_plane = fitting_plane(p_in_upper_leaflet)
    # print("lower_leaflet_plane", lower_leaflet_plane)

    upper_rotated = rotate_operation(upper_leaflet_plane,
                                     target_vector,
                                     p_in_upper_leaflet
                                      )
    lower_rotated = rotate_operation(lower_leaflet_plane,
                                     target_vector,
                                     p_in_lower_leaflet)
    upper_rotated[-1] = np.round(upper_rotated[-1], 3)
    lower_rotated[-1] = np.round(lower_rotated[-1], 3)
    return upper_rotated, lower_rotated

def write_pdb_for_rotated(coords, output_file):
    '''
    Function: 
        Write a pdb file for the rotated_points
    Input:
        the [x,y,z] coordinates of rotated points
    '''
    lines = []
    atom_counter = 1
    for line in coords:
        x = line[0]
        y = line[1]
        z = line[2]
        encoded_res_id = hybrid_36.hy36encode(4, atom_counter)
        # Write PDB format line
        lines.append(
            f"ATOM  {atom_counter:5d}  P   PO4 A{encoded_res_id:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           P\n"
        )
        atom_counter += 1
    with open(output_file, 'w') as file:
        for line in lines:
            file.write(line)
        file.write("END")


def load_and_sample_maps(map_path,
                         membrane_width,
                         save_pdb=True,
                         output_path=(script_dir + '/output/')):
    '''
    # Function name: 
        load_and_sample_maps()

    # Define:
        This funtion is in charge of loading and sampling density maps.

    # Input:
        map_path: absolute path, ; type: string
        membrane_width: take A as the unit, float
        the rate of subsample: rate range 0-1.0, float 
        original_membrane_sample_points_file: True /False , type: bool

    # Outout:
        lower_surface 
            The result of sampling the lower leaflet of a density map. Returned when 'original_membrane_sample_points_file==True'
            Type: PDB
        upper_surface 
            The result of sampling the upper leaflet of a density map. Returned when 'original_membrane_sample_points_file==True'
            Type: PDB
        lower_leaflet 
            The result of sampling the lower leaflet of a density map. Returns a dataframe['atom_id', 'res_id', 'chain_id', 'atom_type', 'coordinates[x,y,z]']
            Type: dataframe
        upper_leaflet 
            The result of sampling the upper leaflet of a density map. Returns a dataframe['atom_id', 'res_id', 'chain_id', 'atom_type', 'coordinates[x,y,z]']
            Type: dataframe
    '''
    density_map, voxel_size, origin, _ = load_density_map(map_path) # correct
    lower_surface, upper_surface = identify_surfaces(density_map) # correct
    _, _, lower_leaflet, upper_leaflet = surface_to_coords(lower_surface, upper_surface, membrane_width, voxel_size, origin) 
    df_lower = coords_to_dataframe(lower_leaflet)
    df_upper = coords_to_dataframe(upper_leaflet)

    if save_pdb == True:
        write_pdb(lower_leaflet, output_path+"vsv_lower_leaflets_subsample_coords2"+".pdb")
        write_pdb(upper_leaflet, output_path+"vsv_upper_leaflets_subsample_coords2"+".pdb")
    else:
        pass
    return df_lower, df_upper

def rotated_and_visualize_curvature(df_lower, df_upper,
                                    save_pdb=False,
                                    output_path=(script_dir + '/output/')):
    rotated_lower, rotated_upper = rotated_points(df_lower, df_upper)
    if save_pdb == True:
        write_pdb_for_rotated(rotated_lower[-1], output_path+"rotated_lower_leaflets_coords2"+".pdb")
        write_pdb_for_rotated(rotated_upper[-1], output_path+"rotated_upper_leaflets_coords2"+".pdb")
    else:
        pass
    return rotated_lower, rotated_upper

def plot_contour_old(rotated_points_list,
                     interpolate_method="bilinear",
                     level: str = 2,
                     my_color_map: str ='RdBu',
                     distance_threshold = 6.5
):
    '''
    Function:
        Using griddata method to smooth/regress the original image.
        Cut the hole!
    Input:
        rotated_points_list: A list contains [R, vector, [x,y,z]]
        interpolate_method = linear, nearest, cubic
        level: the level of oversampling
        my_color_map: specify the color plattee
        output_name: 
    Output:
        An image
    '''
    output_path = script_dir + '/output/'
    coords = rotated_points_list[-1]
    x = np.copy(coords[:, 0])
    y = np.copy(coords[:, 1])
    z = np.copy(coords[:, 2])
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    # Creating a grid
    
    grid_x, grid_y = np.mgrid[x_min: x_max: 400j, y_min: y_max:400j]
    grid_z = griddata((x, y), z, (grid_x, grid_y), method="cubic")
    grid_x = grid_x.astype(np.float32)
    grid_y = grid_y.astype(np.float32)
    grid_z = grid_z.astype(np.float32)

    plt.figure(figsize=(10, 8))
    plt.imshow(grid_z, extent=(x_min+50, x_max+50, y_min, y_max), 
               origin='lower', cmap=my_color_map, aspect='auto', interpolation=interpolate_method)
    plt.colorbar(label='Height')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.grid(True)
    plt.show()
    plt.savefig()

    def gaussian_smoothing(image, kernel_size: int = 5, sigma: float = 1.0):
    # Directly computing the Gaussian kernel outside the loop
        base_kernel = np.fromfunction(
            lambda x, y: (1 / (2 * np.pi * sigma ** 2)) * 
                        np.exp(-((x - (kernel_size - 1) / 2) ** 2 + (y - (kernel_size - 1) / 2) ** 2) / (2 * sigma ** 2)),
            (kernel_size, kernel_size)
        )
        
        center = (kernel_size - 1) / 2
        distances = np.fromfunction(lambda x, y: np.sqrt((x - center) ** 2 + (y - center) ** 2), (kernel_size, kernel_size))
        
        # Padding and Convolution
        pad_size = kernel_size // 2
        padded_image = np.pad(image, ((pad_size, pad_size), (pad_size, pad_size)), mode='constant')
        
        output = np.zeros_like(image)
        for i in range(image.shape[0]):
            for j in range(image.shape[1]):
                patch = padded_image[i:i+kernel_size, j:j+kernel_size]
                
                # Compute the weights for the patch
                weights = np.where(np.isnan(patch), 0, patch) * distances
                # Compute the weighted kernel
                weighted_kernel = base_kernel * weights
                weighted_kernel = weighted_kernel / (np.sum(weighted_kernel)+ 1e-10)
                
                # Convolve the patch with the weighted kernel
                output[i, j] = np.sum(patch * weighted_kernel)
        
        return output

    def batched_cdist_min(X, Y, batch_size=1000):
        mins = []
        for i in range(0, len(X), batch_size):
            batch = X[i:i+batch_size]
            dists = cdist(batch, Y)
            mins.extend(dists.min(axis=1))
        return np.array(mins)

    distances_to_nearest_point = batched_cdist_min(np.column_stack([grid_x.ravel(), grid_y.ravel()]), coords[:, :2])
    distances_grid = distances_to_nearest_point.reshape(grid_x.shape)
    distance_mask = distances_grid > distance_threshold
    grid_z_mask = np.where(distance_mask, np.NAN, grid_z)
    grid_z_mask_v = np.copy(grid_z_mask)
    kernel_size = 8
    grid_z_filter_mask = gaussian_smoothing(grid_z_mask_v, 22, kernel_size)
    plt.figure(figsize=(12, 10))
    contour = plt.contourf(grid_x, grid_y, grid_z_filter_mask,10, cmap='RdBu_r',norm="linear")
    plt.contour(grid_x, grid_y, grid_z_filter_mask,10, linewidths=1, colors='#27374D', linestyles="dashed")
    plt.colorbar(contour)
    plt.grid(False)
    # plt.show()
    plt.savefig(output_path+"contour_final.png", dpi=500, transparent=True)

def replace_into_po4(target_positions_df, save_pdb=False, leaflet=''):
    '''
    Change membrane sampling points into PO4 molecules.
    Input:
        target_positions_df: dataframe which contains the membrane sampling points information 
        save_pdb: True of False. You call.
        leaflet: if you choose saving pdb files, then you need to specify which leaflet your pdb should be. 
    Output:
        dataframe, similar structure with the sampling points.
    '''
    output_path = script_dir + '/output/'
    target_positions= target_positions_df["coordinates"].to_numpy()
    target_positions = np.vstack(target_positions)
    po4_coords_info = np.array([['P', '0.0', '-0.0001', '0.0'], 
       ['O', '1.5234', '-0.0688', '-0.1977'],
       ['O', '-0.6406', '0.718', '-1.1991'],
       ['O', '-0.5679', '-1.4248', '0.1073'],
       ['O', '-0.3149', '0.7757', '1.2895']], dtype='<U32') # po4 3-, PubChem
    po4_info_df = pd.DataFrame(po4_coords_info, columns=['Element', 'X', 'Y', 'Z'])
    current_p_position = po4_info_df[po4_info_df['Element'] == 'P'][['X', 'Y', 'Z']].values[0] # extract the position of P atoms and convert it into array
    current_p_position = current_p_position.astype(np.float16)
    displacements = target_positions - current_p_position # the shape of target_positions is (n, 3). the shape of current_p_position is (5,3). Broadcasting
    new_po4_df = []
    atom_ids = 0
    for i, displacement in enumerate(displacements):
        p_position = po4_info_df[['X', 'Y', 'Z']].astype(np.float16)
        moved_po4 = p_position + displacement
        new_po4 = pd.DataFrame({
            "atom_id": [(atom_ids + i) for i in range(1, len(moved_po4)+1)],
            "res_id": [target_positions_df.iloc[i]["res_id"]] * len(moved_po4),
            "chain_id": [target_positions_df.iloc[i]["chain_id"]] * len(moved_po4),
            "atom_type": po4_info_df['Element'].to_list(),
            "coordinates": moved_po4.values.tolist()
        })
        new_po4_df.append(new_po4)
        atom_ids += 5
    combined_df = pd.concat(new_po4_df, ignore_index=True)
    if save_pdb == True:
        # output pdb
        lines = []
        for index, row in combined_df.iterrows():
            x, y, z = row["coordinates"]
            atom_id = row["atom_id"]
            atom_type = row["atom_type"]
            res_id = row["res_id"]
            encoded_res_id = hybrid_36.hy36encode(4, int(res_id))
            chain_id = row['chain_id']
            line = f"ATOM  {atom_id:5d}  {atom_type:1}   PO4 {chain_id}{encoded_res_id:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom_type:1}\n"
            lines.append(line)
        with open(output_path+leaflet+"_po4.pdb", "w") as file:
            for i in lines:
                file.write(i)
            file.write("END")
    else:
        pass
    return combined_df

'''
Automatically refine your membrane. Restriant is the density map.
When refineing, cut the membrane into patches with some overlaps.
'''

def calculate_patch_dimensions(num_patches, points, overlap_factor=0.15):
    '''
    mmbr should be a dataframe, which contains... 
    we only take the x,y part of membrane
    '''
    mmbr_xoy_coordinates = points["coordinates"].to_numpy()
    mmbr_xoy_coordinates = np.stack(mmbr_xoy_coordinates)

    x_min, x_max = min(mmbr_xoy_coordinates[:,0]), max(mmbr_xoy_coordinates[:,0])
    y_min, y_max = min(mmbr_xoy_coordinates[:,1]), max(mmbr_xoy_coordinates[:,1])
    aspect_ratio = (y_max - y_min) / (x_max - x_min)
    num_patches_y = int(round(np.sqrt(num_patches * aspect_ratio)))
    num_patches_x = int(num_patches / num_patches_y)

    # Calculate patch width and height without and with overlap
    base_patch_width = (x_max - x_min) / num_patches_x
    base_patch_height = (y_max - y_min) / num_patches_y

    overlap_width = base_patch_width * overlap_factor
    overlap_height = base_patch_height * overlap_factor

    adjusted_patch_width = base_patch_width + overlap_width
    adjusted_patch_height = base_patch_height + overlap_height

    patch_ranges = {}
    patch_id = 0
    for i in range(num_patches_y):  # y direction
        for j in range(num_patches_x):  # x direction
            x_start = x_min + j * (adjusted_patch_width - overlap_width)
            x_end = x_start + adjusted_patch_width
            y_start = y_min + i * (adjusted_patch_height - overlap_height)
            y_end = y_start + adjusted_patch_height
            x_start, x_end, y_start, y_end = [math.ceil(x) for x in (x_start, x_end, y_start, y_end)]
            patch_id += 1
            patch_ranges[patch_id] = [(x_start, x_end), (y_start, y_end)]
    # print(patch_ranges)
    return patch_ranges

def extract_atoms_in_patch(points, patch_ranges):
    patch_mol_ids_coot = {}
    patch_mol_ids = {}
    for i, key in enumerate(patch_ranges):
        ranges = patch_ranges[key]
        x_range, y_range = ranges

    # Find atoms within the patch range
        coords = points["coordinates"].to_numpy() 
        coords = np.stack(coords)
        x_coord, y_coord = coords[:, 0], coords[:, 1]
        in_x_range_indices = np.where((x_range[0] <= x_coord) & (x_coord <= x_range[1]))
        in_y_range_indices = np.where((y_range[0] <= y_coord) & (y_coord <= y_range[1]))
        in_both_ranges_indices = np.intersect1d(in_x_range_indices, in_y_range_indices).tolist()
        res_ids_in_patch = []
        chains_id = []
        for row_index in in_both_ranges_indices:
            res_id_in_patch = points.iloc[row_index, points.columns.get_loc("res_id")] 
            res_ids_in_patch.append(res_id_in_patch)
            res_ids_in_patch = list(set(res_ids_in_patch))

            chain_id = points.iloc[row_index, points.columns.get_loc("chain_id")]
            chains_id.append(chain_id)
            chains_id = list(set(chains_id))
        patch_mol_ids[key] = res_ids_in_patch # build the dict key:path_id, value:res_ids_in_patch (list)
    
        mol_id_coot = [] # store the residue selection in a coot format
        for j in res_ids_in_patch:
            mol_id_list = []
            row_index = points.index[points["res_id"] == j].tolist()
            row_index = row_index[0]
            chain_id = points.loc[row_index, "chain_id"]
            j_numerical_format = hybrid_36.hy36decode(4, j)
            # mol_id_list = ["{}".format(chain_id), "{}".format(j_numerical_format), ""]
            mol_id_list = ["{}".format(chain_id), "{}".format(j_numerical_format), ""]

            mol_id_coot.append(mol_id_list)
        patch_mol_ids_coot[key] = mol_id_coot
    return patch_mol_ids_coot

def coot_input(mmbr_pdb, density_map, protein_pdb):
    mmbr_id = coot.read_pdb(mmbr_pdb)
    protein_id = coot.read_pdb(protein_pdb)
    map_id = coot.read_ccp4_map(density_map, 0) # do not show the negative value of a difference map
    return mmbr_id, map_id, protein_id

def coot_init(map_id):
    '''
    Function:
        Designate the reference density map
    '''
    coot.set_imol_refinement_map(map_id)
    # Do not ask about accept/reject refinement. Set it before refinement
    coot.set_refinement_immediate_replacement(1) 

def coot_refine(mmbr_id ,patch_mol_ids_coot):
    # ipatches = []
    for patch_id, mol_id in patch_mol_ids_coot.items():
        # print(mol_id)
        # ipatch = coot.new_molecule_by_residue_specs_py(mmbr_id, mol_id)
        # ipatches.append(ipatch)
        print("mmbr_id", mmbr_id)
        print(mol_id)
        coot.rigid_body_refine_by_residue_ranges_py(mmbr_id, mol_id)
        # coot.accept_moving_atoms_py()

def coot_refine_end(protein_id):
    coot.graphics_to_ca_representation(protein_id)

def main():
    '''
    Input paras: 
        - map_path, the path of density map; type, string
        - membrane_width
        - save_pdb, if you want to see the sample result? ; type, bool; default_value=False
        - output_path,
        - threshold
        - subsampe_rate (optional)
    '''
    map_path = "/Users/boyuan/Desktop/cryosparc_P29_VSV_sym_ave_R19-and-L20_onGrid_P29_J428.mrc" # input paras
    membrane_width = 55 # input paras
    output_path = script_dir + '/output/' # input paras

    # rotated_lower, rotated_upper = rotated_and_visualize_curvature(lower_coords)
    df_lower, df_upper = load_and_sample_maps(map_path, membrane_width,save_pdb=True) 
    po4_lower, po4_upper = replace_into_po4(df_lower, leaflet="lower"), replace_into_po4(df_upper, leaflet="upper")
    rotated_lower, rotated_upper = rotated_and_visualize_curvature(df_lower, df_upper)
    plot_contour_old(rotated_lower)

    # For coot refine
    num_patches=10
    overlap_factor=0.15
    patch_ranges = calculate_patch_dimensions(num_patches, po4_upper, overlap_factor)
    patch_mol_ids_coot = extract_atoms_in_patch(po4_upper, patch_ranges)
    
    # mmbr_pdb = output_path+"upper_leaflets_subsample_coords.pdb"
    mmbr_pdb = output_path+"po4.pdb"
    protein_pdb = example_path+"mito_protein.pdb"
    mmbr_id, map_id, protein_id = coot_input(mmbr_pdb, map_path, protein_pdb)
    coot_init(map_id)
    coot_refine(mmbr_id, patch_mol_ids_coot)
    print(protein_id)
    coot_refine_end(protein_id)


if __name__ == "__main__":
    main()
