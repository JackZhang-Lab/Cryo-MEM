import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.neighbors
from sklearn.neighbors import KDTree
from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interactive
from scipy.spatial import distance_matrix

def parse_pdb(pdb_file):
    """
    Parse a PDB file and return a dictionary with molecule numbers as keys and coordinates as values.
    """
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
                coords[molecule_number][atom_type] = np.array([x, y, z])
    print("Number of atoms before removing duplicates", len(list(coords.keys())))
    return coords

def remove_duplicates(coords, threshold=0.1, cluster_radius=6.0):
    """
    Remove duplicate atoms from the coords dictionary.
    Two atoms are considered duplicates if their Euclidean distance is smaller than a threshold.
    """
    # Build the KDTree
    tree = KDTree([molecule['P'] for molecule in coords.values()])
    duplicates = set()
    for molecule_number, coord in coords.items(): # return molucule_number = dict.key(), coord = dict.values()
        if molecule_number not in duplicates:
            # Query the tree for neighbors within the cluster radius
            sphere_heart_coord = np.array([coord['P']])
            indices = tree.query_radius(sphere_heart_coord, r=cluster_radius)
            # Convert indices to molecule numbers
            neighbor_molecule_numbers = [list(coords.keys())[index] for index in indices[0]]
            # Find duplicates within the neighbors
            for neighbor_molecule_number in neighbor_molecule_numbers:
                if neighbor_molecule_number != molecule_number:
                    neighbor_P_coord = coords[neighbor_molecule_number]['P']
                    distance = np.linalg.norm(sphere_heart_coord - neighbor_P_coord)
                    if distance < threshold:
                        duplicates.add(neighbor_molecule_number)
    # Remove duplicates from coords
    for molecule_number in duplicates:
        del coords[molecule_number]
    print(f"Number of remaining atoms after removing duplicates: {len(coords)}")
    return coords

def calculate_force(coord1, coord2, min_distance, max_distance):
    """
    Calculate the force between two P atoms using Hooke's law.
    """
    distance = np.linalg.norm(coord1 - coord2)
    if distance == 0:
        force = np.array([0, 0, 0])
    elif distance < min_distance:
        # If the distance is too small, the force is repulsive
        force = (min_distance - distance) * (coord1 - coord2) / distance
    elif distance > max_distance:
        # If the distance is too large, the force is attractive
        force = (max_distance - distance) * (coord1 - coord2) / distance
    else:
        # If the distance is in the right range, there is no force
        force = np.array([0, 0, 0])
    return force

def update_position(molecule_coords, force, time_step, max_z_movement):
    """
    Update the position of a P atom and its associated O atoms using Newton's second law.
    """
    # The mass and the acceleration are assumed to be 1 for simplicity
    # The type of delta is numpy array. It is defined in the for loop.
    delta = force * time_step
    # Limit the z axis movement and update the z value of delta np.array
    delta[2] = np.clip(delta[2], -max_z_movement, max_z_movement)
    for atom_type in molecule_coords.keys():
        molecule_coords[atom_type] += delta
    # print(molecule_coords)
    return molecule_coords
#
def write_pdb(coords, template_pdb_file, output_pdb_file):
    """
    Write a PDB file using new coordinates and a template PDB file.
    """
    with open(template_pdb_file, 'r') as f:
        lines = f.readlines()
    new_lines = []
    a = 0
    for line in lines:
        if line[0: 6] == 'HETATM':
            atom_type = line[13:15].strip()
            molecule_number = line[23: 28].strip()
            if molecule_number in coords.keys():
                a += 1
                x, y, z = coords[molecule_number][atom_type]
                # >: line to the right;
                # 8: means taking up 8 positions included space
                # .3f: keep three numbers after the dot
                new_line = f"{line[:6]}{a:>5}  {atom_type:<3} PO4 B {int(molecule_number):>4}   {x:>8.3f}{y:>8.3f}{z:>8.3f}{line[54:60]}{float(line[60:66]):>6.2f} {atom_type[0]:>12} \n"
                new_lines.append(new_line)

    with open(output_pdb_file, 'w') as f:
        f.writelines(new_lines)
        f.writelines('END')


path = '/Users/boyuan/Desktop/Yale/mmbr_fit/'
coords = parse_pdb(path + "input_data/mmbr_fit_jul_11th_1228_adjusted.pdb")
# Parameters
min_distance = 2
max_distance = 7
time_step = 0.5
max_z_movement = 0.01
max_iterations = 4

# Remove duplicate atoms
coords = remove_duplicates(coords)
molecule_numbers = list(coords.keys())
# Build the KDTree with P atom coordinates. More efficiently store and abstract the position information P atoms
tree = KDTree([molecule['P'] for molecule in coords.values()])

# Implementation aboved functions
for iteration in range(max_iterations):
    forces = {}
    # distances = distance_matrix()
    # for i in range(len(coords.keys())):
    #     for j in range((i+1, len(coords.keys()))):
    #         if distances[i, j] < min_distance or distances[i, j] > max_distance:

    # New coords dictionary without duplicates will be read in the following steps
    molecule_searched_number = []
    for molecule_number, molecule in coords.items():
        if molecule_number not in molecule_searched_number:
        # Query the tree for the neighbors within max_distance
            indices = tree.query_radius([molecule['P']], r=9)[0]
        # 在该搜索半径内所有的分子
            for index in indices:
                # Convert index to molecule_number
                neighbor_molecule_number = molecule_numbers[index]
                # 确保对已经移动距离的分子不进行重复操作
                molecule_searched_number.append(neighbor_molecule_number)
                if neighbor_molecule_number != molecule_number:
                    molecule_in_cluster_coord = coords[neighbor_molecule_number] # 该PO4分子 所有原子的coordinates
                    # 根据P原子的距离，确定吸引力/排斥力的数值
                    force = calculate_force(molecule['P'], molecule_in_cluster_coord['P'], min_distance, max_distance)
                    update_position(coords[neighbor_molecule_number], force, time_step, max_z_movement)
            # test = int(indices[0])
            # # print(test, type(test))
            # print("end loop updated: ", coords[test]['O1'], "\n")
        # Plot the positions of the P atoms after each iteration
    fig = plt.figure(figsize=(10, 10))  # increase the size of the figure
    ax = fig.add_subplot(111, projection='3d')
    x = [molecule['O1'][0] for molecule in coords.values()]
    y = [molecule['O1'][1] for molecule in coords.values()]
    z = [molecule['O1'][2] for molecule in coords.values()]
    ax.scatter(x, y, z)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title(f'Iteration {iteration}')
    plt.show()

# Use the function to write a new PDB file

write_pdb(coords, path+"input_data/mmbr_fit_jul_11th_1228_adjusted.pdb", path+"output_data/new_mmbr_fit_test_jul_18th.pdb")

# # if __name__ == '__main__':




