#!/usr/bin/env python
# coding: utf-8

# ## Import Libs

# In[ ]:


import pythreejs, ipywidgets
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal
import sklearn.neighbors
from sklearn.neighbors import KDTree

from scipy.spatial.distance import cdist
from scipy.spatial import distance
from scipy.interpolate import griddata
from scipy import ndimage

# In[ ]:


# Loading my data
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
                coords[molecule_number][atom_type] = np.array([x, y, z])
    return coords

pdb_file="/home/boyuan/Desktop/mmbr_fit/chatGPT/downleaflet/final.pdb"

# Dictionary of the molecule coordinates
test_data = parse_pdb(pdb_file)

# Abstract P atoms' coordinates
p_coordinates = []
for i in test_data.keys():
    a = test_data[i]["P"]
    p_coordinates.append(a)

# Convert P coordinates to cupy.ndarray
p_coordinates = np.array(p_coordinates)
p_coordinates = np.float32(p_coordinates)
p_coordinates.shape

# In[ ]:


x = np.copy(p_coordinates[:, 0])
y = np.copy(p_coordinates[:, 1])
z = np.copy(p_coordinates[:, 2])

# ## Find the density center (Optional)

# In[ ]:


from scipy.spatial.distance import pdist, squareform

# Compute the pairwise distances between all points
distance_matrix = squareform(pdist(p_coordinates))

distance_matrix.shape

weights = 1 / np.sum(distance_matrix, axis=1)
weights.shape

density_center_3d = np.average(p_coordinates, weights=weights, axis=0)

centroid = np.mean(p_coordinates, axis=0)

# ## 1. Find the best fitting plane
# Using cupy to perform least square fitting

# In[ ]:


import time
import cupy as cp

x_gpu = cp.array(x)
y_gpu = cp.array(y)
z_gpu = cp.array(z)

# Setup the coefficient matrix A
A_gpu = cp.vstack((x_gpu, y_gpu, cp.ones_like(x_gpu))).T

%time coeff_gpu, _, _, _ = cp.linalg.lstsq(A_gpu, z_gpu, rcond=None)

# Extract the coefficients A, B, and D
coeff = coeff_gpu.get()
A_coeff, B_coeff, D_coeff = coeff

# The normal to the plane is given by [A_coeff, B_coeff, -1]
normal = np.array([A_coeff, B_coeff, -1])

# In[ ]:


A_coeff, B_coeff, D_coeff, normal

# In[ ]:


centroid = [np.mean(x), np.mean(y), np.mean(z)]
centroid

# ### 1.2 Plot the fitting plane

# In[ ]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Define the plane using the coefficients
def plane(x, y, A_coeff, B_coeff, D_coeff):
    return A_coeff * x + B_coeff * y + D_coeff

# Create a meshgrid for the fitting plane
x_plane = np.linspace(min(x) - 1, max(x) + 1, 100)
y_plane = np.linspace(min(y) - 1, max(y) + 1, 100)
X_plane, Y_plane = np.meshgrid(x_plane, y_plane)
Z_plane = plane(X_plane, Y_plane, A_coeff, B_coeff, D_coeff)

from ipywidgets import interactive

def plot_3d_plane_and_points(angle1=30, angle2=30):
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Rotate the axes according to the slider values
    ax.view_init(angle1, angle2)

    # Plot the plane
    ax.plot_surface(X_plane, Y_plane, Z_plane, alpha=0.5, rstride=100, cstride=100, color='grey', label='Fitting Plane')

    # Plot the scattered points
    ax.scatter(x, y, z, color='lightpink', label='Scattered Points')

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Least Squares Fitting Plane and Normal Vector')

    plt.show()


interactive_plot = interactive(plot_3d_plane_and_points, angle1=(0, 360), angle2=(0, 360))
output = interactive_plot.children[-1]
output.layout.height = '550px'
interactive_plot


# ## 2. Find the rotation matrix

# In[ ]:


def rotate_operation(v, target_vector, D_coeff, points):

    print("Your normal is: ", v)
    print("Your target_v is:", target_v)

    v = v / np.linalg.norm(v)
    t = target_vector

    rotation_axis = np.cross(v, t)

    cos_theta = np.dot(v, t) / (np.linalg.norm(v) * np.linalg.norm(t))
    theta = np.arccos(cos_theta)

    # Rodrigues' rotation formula
    K = np.array([
        [0, -rotation_axis[2], rotation_axis[1]],
        [rotation_axis[2], 0, -rotation_axis[0]],
        [-rotation_axis[1], rotation_axis[0], 0]
    ])

    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * np.dot(K, K)

    points_rotated = np.dot(points, rotation_matrix.T)

    v_rotated = np.dot(v, rotation_matrix.T)
    print("Done")
    return R, v_rotated, points_rotated


# In[ ]:


target_v = np.array([0, 0, 1])
points = np.copy(p_coordinates)
points[1:5]
R1, normal_rotated, points_rotated= rotate_operation(normal, target_v, D_coeff, points)

# In[ ]:


x_rotated = np.copy(points_rotated[:, 0])
y_rotated = np.copy(points_rotated[:, 1])
z_rotated = np.copy(points_rotated[:, 2])

def plot_rotated_points(angle1=30, angle2=30):


    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Rotate the axes according to the slider values
    ax.view_init(angle1, angle2)

    # Plot the scattered points
    ax.scatter(x, y, z, color='lightpink', label='Scattered Points')
    ax.scatter(x_rotated, y_rotated, z_rotated, color="mediumpurple")

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    plt.show()

# Call the function with ipywidgets' interactive for sliders
interactive_plot = interactive(plot_rotated_points, angle1=(0, 360), angle2=(0, 360))
output = interactive_plot.children[-1]
output.layout.height = '550px'
interactive_plot


# ## 3. Output the new coordinates after rotation
# Write the rotated points into txt

# In[ ]:


# Convert coordinates to text format
formatted_coords = "\n".join([f"{x:.5f},{y:.5f},{z:.5f}" for x, y, z in points_rotated])

# Write to a txt file
file_path = "/home/boyuan/Desktop/mmbr_fit/input/class01/rotated_points_class01.txt"
with open(file_path, 'w') as file:
    file.write(formatted_coords)

file_path
