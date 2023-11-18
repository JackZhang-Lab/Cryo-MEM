#!/usr/bin/env python
# coding: utf-8

# In[1]:


!pip install matplotlib plotly mayavi


# In[11]:


# Import necessary libraries
import numpy as np
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.neighbors import KNeighborsRegressor
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import plotly.graph_objects as go

# In[6]:


# Load the data
# data = pd.read_csv('/Users/boyuan/Desktop/Yale/mmbr_fit/nineiters_output_min_dist=2_max_dist=9_time_step=0.4_max_iters=4 (1).pdb', delim_whitespace=True, header=None)
# Define the column positions (indices) and names
colspecs = [(0, 6), (6, 11), (12, 16), (17, 20), (21, 22), (22, 26), 
            (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78)]

column_names = ["record_name", "atom_number", "atom_name", "residue_name", 
                "chain_id", "residue_number", "x", "y", "z", 
                "occupancy", "temperature_factor", "element_symbol"]

# Load the PDB file with the correct column names and positions
data = pd.read_fwf('/Users/boyuan/Desktop/Yale/mmbr_fit/nineiters_output_min_dist=2_max_dist=9_time_step=0.4_max_iters=4 (1).pdb', 
                   colspecs=colspecs, 
                   header=None, 
                   names=column_names)

# Convert the x, y, z coordinates to numeric data type
data[['x', 'y', 'z']] = data[['x', 'y', 'z']].apply(pd.to_numeric, errors='coerce')

# Remove rows with NaN values in the x, y, or z columns
data = data.dropna(subset=['x', 'y', 'z'])

# Extract the coordinates
x = data['x']
y = data['y']
z = data['z']


# Fit a k-NN model to the existing data
X = np.array([x, y]).T
knn = KNeighborsRegressor(n_neighbors=3)
knn.fit(X, z)

# Create a grid of points where we want to evaluate the model
x_grid, y_grid = np.meshgrid(np.linspace(min(x), max(x), 10), np.linspace(min(y), max(y), 10))

# Use the k-NN model to predict the z values on the entire grid (including the hole)
X_grid = np.array([x_grid.ravel(), y_grid.ravel()]).T
z_grid_filled = knn.predict(X_grid).reshape(x_grid.shape)

# Fit the polynomial model to the entire grid data (including the filled hole)
model = make_pipeline(PolynomialFeatures(2), LinearRegression())
model.fit(X_grid, z_grid_filled.ravel())

# Use the model to predict the z values on the entire grid
z_grid_predicted = model.predict(X_grid).reshape(x_grid.shape)

# Plot the original points and the surface fitted by the polynomial model
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, s=5, color='b')
ax.plot_surface(x_grid, y_grid, z_grid_predicted, color='r', alpha=0.3)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# In[15]:


# Create a 3D plot
fig = go.Figure(data=[go.Surface(z=z_grid_predicted, x=x_grid, y=y_grid)])

# Update layout options
fig.update_layout(title='Your Surface', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))

# Display the plot
fig.show()

# In[16]:


# import plotly.graph_objects as go

# Create a scatter plot of the original points
scatter = go.Scatter3d(
    x=x, y=y, z=z,
    mode='markers',
    marker=dict(
        size=2,
        color='blue',
    )
)

# Create a surface plot of the generated surface
surface = go.Surface(z=z_grid_predicted, x=x_grid, y=y_grid, colorscale='Reds')

# Add both plots to the figure
fig = go.Figure(data=[scatter, surface])

# Update layout options
fig.update_layout(title='Your Surface', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))

# Display the plot
fig.show()


# In[14]:


# Extract the polynomial features and linear regression model from the pipeline
poly_features = model.named_steps['polynomialfeatures']
linear_regression = model.named_steps['linearregression']

# Get the coefficients and intercept
coefficients = linear_regression.coef_
intercept = linear_regression.intercept_

# The coefficients correspond to the terms in the polynomial in the following order:
# 1, x, y, x^2, x*y, y^2
# We rearrange them to match the order in the equation: x^2, y^2, x*y, x, y, 1
coefficients = [coefficients[i] for i in [3, 5, 4, 1, 2, 0]]

# Print the equation
equation = "z = "
for i, coeff in enumerate(coefficients):
    if i == 0:
        equation += f"{coeff:.5f}x^2 "
    elif i == 1:
        equation += f"+ {coeff:.5f}y^2 "
    elif i == 2:
        equation += f"+ {coeff:.5f}xy "
    elif i == 3:
        equation += f"+ {coeff:.5f}x "
    elif i == 4:
        equation += f"+ {coeff:.5f}y "
    elif i == 5:
        equation += f"+ {coeff:.5f}"

print(equation)


# In[ ]:




# In[ ]:



