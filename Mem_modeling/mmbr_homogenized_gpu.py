import cupy as cp
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.neighbors
from sklearn.neighbors import KDTree
from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interactive
from scipy.spatial import distance_matrix

def parse_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    coords = {}
    for line in lines:
        if line[0: 6] == 'HETATM':
            atom_type = line[13:15].strip()
            molecule_number = line[23: 28].strip()
            x = float(line[31: 38])
            y = float(line[39: 46])
            z = float(line[47: 54])
            if atom_type in ['P', 'O1', 'O2', 'O3', 'O4']:
                if molecule_number not in coords:
                    coords[molecule_number] = {}
                coords[molecule_number][atom_type] = cp.array([x, y, z])
    return coords

def remove_duplicates(coords, threshold=0.1, cluster_radius=6.0):
    # tree = KDTree(cp.asnumpy([molecule['P'] for molecule in coords.values()]))
    tree = KDTree([molecule['P'].get() for molecule in coords.values()])

    duplicates = set()
    for molecule_number, coord in coords.items():
        if molecule_number not in duplicates:
            sphere_heart_coord = cp.array([coord['P']])
            indices = tree.query_radius(cp.asnumpy(sphere_heart_coord), r=cluster_radius)
            neighbor_molecule_numbers = [list(coords.keys())[index] for index in indices[0]]
            for neighbor_molecule_number in neighbor_molecule_numbers:
                if neighbor_molecule_number != molecule_number:
                    neighbor_P_coord = coords[neighbor_molecule_number]['P']
                    distance = cp.linalg.norm(sphere_heart_coord - neighbor_P_coord)
                    if distance < threshold:
                        duplicates.add(neighbor_molecule_number)
    for molecule_number in duplicates:
        del coords[molecule_number]
    return coords

def calculate_force(coord1, coord2, min_distance, max_distance):
    distance = cp.linalg.norm(coord1 - coord2)
    if distance == 0:
        force = cp.array([0, 0, 0])
    elif distance < min_distance:
        force = (min_distance - distance) * (coord1 - coord2) / distance
    elif distance > max_distance:
        force = (max_distance - distance) * (coord1 - coord2) / distance
    else:
        force = cp.array([0, 0, 0])
    return force

def update_position(molecule_coords, force, time_step, max_z_movement):
    delta = force * time_step
    delta[2] = cp.clip(delta[2], -max_z_movement, max_z_movement)
    for atom_type in molecule_coords.keys():
        molecule_coords[atom_type] += delta
    return molecule_coords

def write_pdb(coords, template_pdb_file, output_pdb_file):
    with open(template_pdb_file, 'r') as f:
        lines = f.readlines()
    new_lines = []
    a = 0
    for line in lines:
        a += 1
        if line[0: 6] == 'HETATM':
            atom_type = line[13:15].strip()
            molecule_number = line[23: 28].strip()
            if molecule_number in coords.keys():
                x, y, z = coords[molecule_number][atom_type].tolist()
                new_line = f"{line[:6]}{a:>5}  {atom_type:<3} PO4 B {int(molecule_number):>4}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{line[54:60]}{float(line[60:66]):>6.2f} {atom_type[0]:>12} \n"
                new_lines.append(new_line)
    with open(output_pdb_file, 'w') as f:
        f.writelines(new_lines)
        f.writelines('END')

def run_optimization_and_plot(min_distance, max_distance, time_step, max_iterations, pdb_file, max_z_movement=0.1):
    for iteration in range(max_iterations+4):
        coords = parse_pdb(pdb_file)
        coords = remove_duplicates(coords)
        molecule_numbers = list(coords.keys())
        tree = KDTree([molecule['P'].get() for molecule in coords.values()])

        for iteration in range(max_iterations):
            forces = {}
            coords = coords
            molecule_searched_number = []
            for i in range(len(list(coords.keys()))):
                coords = coords
                molecule_number, molecule = list(coords.keys())[i], list(coords.values())[i]
                if molecule_number not in molecule_searched_number:
                    indices = tree.query_radius([molecule['P'].get()], r=5)[0]
                    for index in indices:
                        neighbor_molecule_number = molecule_numbers[index]
                        molecule_searched_number.append(neighbor_molecule_number)
                        if neighbor_molecule_number != molecule_number:
                            molecule_in_cluster_coord = coords[neighbor_molecule_number]
                            force = calculate_force(molecule['P'], molecule_in_cluster_coord['P'], min_distance, max_distance)
                            update_position(coords[neighbor_molecule_number], force, time_step, max_z_movement)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        x = [molecule['P'].get()[0] for molecule in coords.values()]
        y = [molecule['P'].get()[1] for molecule in coords.values()]
        z = [molecule['P'].get()[2] for molecule in coords.values()]
        ax.scatter(x, y, z)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title(f'Final positions of P atoms (min_distance={min_distance}, max_distance={max_distance}, time_step={time_step}, max_iterations={max_iterations})')

        plot_filename = f"/home/boyuan/Desktop/mmbr_fit/chatGPT/plot_min_distance={min_distance}_max_distance={max_distance}_time_step={time_step}_max_iterations={max_iterations}.png"
        plt.savefig(plot_filename)
        plt.close(fig)

        pdb_output_filename = f"/home/boyuan/Desktop/mmbr_fit/chatGPT/nineiters_output_min_dist={min_distance}_max_dist={max_distance}_time_step={time_step}_max_iters={max_iterations}.pdb"
        write_pdb(coords, pdb_file, pdb_output_filename)

# Further reduce the parameter ranges for faster execution
min_distance_range = cp.arange(2, 3, 1) # Step size of 2 to further reduce the number of iterations
max_distance_range = cp.arange(9, 10, 1) # Step size of 2 to further reduce the number of iterations
time_step_range = cp.arange(0.4, 0.5, 0.1) # Step size of 0.3 to further reduce the number of iterations
max_iterations_range = cp.arange(4, 5, 1) # Step size of 1 to cover all possible values

parameter_grid = [(int(min_distance), int(max_distance), float(time_step), int(max_iterations))
                  for min_distance in min_distance_range
                  for max_distance in max_distance_range
                  for time_step in time_step_range
                  for max_iterations in max_iterations_range]

print(len(parameter_grid)) # Number of combinations


for i, (min_distance, max_distance, time_step, max_iterations) in enumerate(parameter_grid):
    run_optimization_and_plot(min_distance, max_distance, time_step, max_iterations, pdb_file="/home/boyuan/Desktop/mmbr_fit/input/mmbr_fit_original.pdb")
    print("Current paras are: ", (min_distance, max_distance, time_step, max_iterations))
print("All optimizations, plots, and pdb files completed.")
