#!/usr/bin/env python
# coding: utf-8

# ## Load un-curved mmbr

# In[17]:


import os
import numpy as np

p_coordinates = []

with open("/home/boyuan/Desktop/mmbr_fit/input/double_mmbr_square.pdb", "r") as file:
    contents = file.readlines()
    
p_lines = [line for line in contents if line[0:6]=="HETATM" and line[-2]=='P']


# Displaying the first few lines for inspection
p_lines[:10]

for line in p_lines:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:55])
        p_coordinates.append((x, y, z))
    except ValueError:
        # In case of any issues with a line, print it for inspection
        print(line)

p_coord_lipid = np.array(p_coordinates)

x_lipid = np.copy(p_coord_lipid[:,0])
y_lipid = np.copy(p_coord_lipid[:,1])
z_lipid = np.copy(p_coord_lipid[:,2])

# In[18]:


# extract one leaflet
mask = z_lipid < 0
x_lipid = x_lipid[mask]
y_lipid = y_lipid[mask]
z_lipid = z_lipid[mask]

x_lipid = x_lipid +201
y_lipid = y_lipid +202

# In[19]:


x_lipid.min(), x_lipid.max()

# In[ ]:


z_lipid_fit = griddata((x_lipid, y_lipid), z_lipid, (target_x, target_y), method='linear')

target_coords = np.copy(target_coords_raw)
target_coords = target_coords[~np.isnan(target_coords[:, 2])]

target_coords[:,0].shape

target_x = np.copy(target_coords[:,0])
target_y = np.copy(target_coords[:,1])
target_z = np.copy(target_coords[:,2])

# In[ ]:


# Rbf Regression
from scipy.interpolate import Rbf

rbf_model = Rbf(target_x, target_y, target_z, function='gaussian',epsilon=8)
z_original_rbf = rbf_model(x_lipid, y_lipid)

# Plot the Rbf regression result.
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')
# ax.view_init(30, 30)
ax.scatter(x_lipid, y_lipid, z_original_rbf, color='red', s=1, label='P atoms Points')
# ax.scatter(target_x, target_y, target_z, color='red', s=1, label='P atoms Points')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()

# Convert coordinates to text format
x_lipid = list(x_lipid)
y_lipid = list(y_lipid)
z_original_rbf = list(z_original_rbf)


formatted_coords = "\n".join([f"{x:.3f},{y:.3f},{z:.3f}" for x,y,z in zip(x_lipid, y_lipid, z_original_rbf)])

# Write to a txt file
file_path = "/home/boyuan/Desktop/mmbr_fit/process/lipids_down_p_target_coordinates.txt"
with open(file_path, 'w') as file:
    file.write(formatted_coords)

file_path
