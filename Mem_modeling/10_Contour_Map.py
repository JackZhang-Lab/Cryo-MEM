#!/usr/bin/env python
# coding: utf-8

# ## 0 Loading
# ### 0.0 Import Libs

# In[2]:


import pythreejs, ipywidgets
import cupy as cp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal
import sklearn.neighbors
from sklearn.neighbors import KDTree

from scipy.spatial import distance
from scipy import ndimage

from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interactive

# ### 0.1 Load Data

# In[5]:


import os

def load_coords(path):

    with open(path, 'r') as f:
        lines = f.readlines()
        matrix = []
        for line in lines:
            x, y, z = map(float, line.strip().split(","))
            matrix.append([x, y, z])
        coords = np.copy(matrix).astype(np.float32)
        x = np.copy(coords[:, 0])
        y = np.copy(coords[:, 1])
        z = np.copy(coords[:, 2])
    return x, y, z, coords

filepath = "/home/boyuan/Desktop/mmbr_fit/input/class01/rotated_points_class01.txt"
x, y, z, coords = load_coords(filepath)
coords.shape

# In[13]:


z.min()

# ## 1 Grid data 

# In[35]:


import time
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()
print(y_max)
# Creating a grid
%time grid_x, grid_y = np.mgrid[x_min: x_max: 400j, y_min: y_max:400j]
%time grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

# In[36]:


grid_x = grid_x.astype(np.float32)
grid_y = grid_y.astype(np.float32)
grid_z = grid_z.astype(np.float32)
grid_x.shape

# ### 1.1 Plot

# In[42]:


plt.figure(figsize=(10, 8))
plt.imshow(grid_z, extent=(x_min+50, x_max+50, y_min, y_max), origin='lower', cmap='RdBu', aspect='auto', interpolation='none')
plt.colorbar(label='Height')
# plt.scatter(x, y, color='green', label='P atoms Points')
plt.title('P atoms Distribution')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
# plt.legend()
plt.grid(True)
plt.show()

# In[48]:


from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

ax = plt.figure().add_subplot(projection='3d')
# X, Y, Z = axes3d.get_test_data(0.05)

# Plot the 3D surface
ax.plot_surface(grid_x, grid_y, grid_z, edgecolor='royalblue', lw=0.1, rstride=8, cstride=8,
                alpha=0.5)

# Plot projections of the contours for each dimension.  By choosing offsets
# that match the appropriate axes limits, the projected contours will sit on
# the 'walls' of the graph
# ax.contour(grid_x, grid_y, grid_z, zdir='z', offset=100, cmap='coolwarm')
# ax.contour(grid_x, grid_y, grid_z, zdir='x', offset=0, cmap='coolwarm')
ax.contour(grid_x, grid_y, grid_z, zdir='y', offset=0, cmap='coolwarm', levels=3)

ax.set(xlim=(0, 400), ylim=(0, 400), zlim=(100, 200),
       xlabel='X', ylabel='Y', zlabel='Z')

plt.show()

# ## 2 Gaussian Filter
# ### 2.0 Def gaussian function
# 

# In[5]:


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

# ### 2.1 Apply Mask

# In[6]:


from scipy.spatial.distance import cdist
import numpy as np

def batched_cdist_min(X, Y, batch_size=1000):
    mins = []
    for i in range(0, len(X), batch_size):
        batch = X[i:i+batch_size]
        dists = cdist(batch, Y)
        mins.extend(dists.min(axis=1))
    return np.array(mins)

distances_to_nearest_point = batched_cdist_min(np.column_stack([grid_x.ravel(), grid_y.ravel()]), coords[:, :2])

distance_threshold = 6.5
distances_grid = distances_to_nearest_point.reshape(grid_x.shape)
distance_mask = distances_grid > distance_threshold

# In[7]:


grid_z_mask = np.where(distance_mask, np.NAN, grid_z)
grid_z_mask_v = np.copy(grid_z_mask)
grid_z_mask.shape

# In[8]:


kernel_size = 8
grid_z_filter_mask = gaussian_smoothing(grid_z_mask, 22, kernel_size)

# In[39]:


grid_z_filter_mask.shape

# ## 3 Contour map

# In[22]:


font_axis_properties = {'family':'monospace', 'size':20}
font_ticks_properties = {'family':'monospace', 'weight':'bold', 'size':13}
plt.figure(figsize=(12, 10))
contour = plt.contourf(grid_x, grid_y, grid_z_filter_mask,10, cmap='RdBu_r',norm="linear")
plt.contour(grid_x, grid_y, grid_z_filter_mask,10, linewidths=1, colors='#27374D', linestyles="dashed")
cbar = plt.colorbar(contour)
cbar.set_label('Curvature',fontsize=15, fontname='monospace',fontweight='bold')
cbar.ax.tick_params(labelsize=12, labelrotation=45)
for label in cbar.ax.yaxis.get_ticklabels():
    label.set_fontname('monospace')
    label.set_weight('bold')
# plt.colorbar.set_label(**font_ticks_properties)
# plt.scatter(x, y, color='green', label='P Points', zorder=5)
# plt.title('Contour Map of P atoms Distribution')
plt.xlabel('X (Å)', **font_axis_properties)
plt.ylabel('Y (Å)', **font_axis_properties)
plt.xticks(**font_ticks_properties)
plt.yticks(**font_ticks_properties)
# plt.legend()
plt.grid(False)
# plt.show()
plt.savefig("/home/boyuan/Desktop/mmbr_fit/output/Pic/contour_final.png", dpi=500, transparent=True)
