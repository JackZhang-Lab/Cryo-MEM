#!/usr/bin/env python
# coding: utf-8

# # try to find analytical equations

# In[ ]:


import numpy as np
from sklearn.linear_model import ElasticNet, BayesianRidge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# In[3]:


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

# In[15]:


x_grid, y_grid = np.meshgrid(x, y)


# In[16]:


X_samples = np.column_stack((x_grid.ravel(), y_grid.ravel()))

# In[17]:


X_samples.shape

# In[19]:


poly_features = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly_features.fit_transform(X_samples)
assert X_poly.shape[0] == z.ravel().shape[0]

scaler = StandardScaler()
elastic_net = make_pipeline(scaler, ElasticNet(alpha=0.1, l1_ratio=0.5))
elastic_net.fit(X_poly, z.ravel())
z_elastic_net = elastic_net.predict(poly_features.transform(X_samples)).reshape(x_grid.shape)

# In[ ]:


z_elastic_net.shape

# In[13]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ipywidgets import interactive
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plotting the surface
ax.plot_surface(x_grid, y_grid, z_elastic_net, alpha=0.5, cmap='viridis')

# Plotting the scatter points
ax.scatter(x, y, z, c='r', marker='o', s=5)

ax.set_title(title)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# In[ ]:


X = 
