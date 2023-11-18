import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Define the path to your data files
file1 = 'smoothed_x_points.cor'
file2 = 'scaled_prot.cor'

# Define the X range (minimum and maximum X values)
x_min_map = 20  # Change this to your desired minimum X
x_min_prot = 50  # Change this to your desired minimum X
x_max_map = 110  # Change this to your desired maximum X
x_max_prot = 80  # Change this to your desired maximum X

# Load points from file 1 (assuming the file contains one [x, y, z] coordinate per line)
points_file1 = np.loadtxt(file1)
x_values_file1 = points_file1[:, 0]
y_values_file1 = points_file1[:, 1]
z_values_file1 = points_file1[:, 2]

# Filter points in file 1 based on the X range
mask_file1 = (x_values_file1 >= x_min_map) & (x_values_file1 <= x_max_map)
x_values_file1 = x_values_file1[mask_file1]
y_values_file1 = y_values_file1[mask_file1]
z_values_file1 = z_values_file1[mask_file1]

# Create a contour plot for file 1
plt.figure(figsize=(12, 8))

# Create a grid of Y and Z values using meshgrid
y_range = np.linspace(np.min(y_values_file1), np.max(y_values_file1), 100)
z_range = np.linspace(np.min(z_values_file1), np.max(z_values_file1), 100)
Y, Z = np.meshgrid(y_range, z_range)

# Interpolate the X values onto the YZ grid
from scipy.interpolate import griddata
points_file1 = np.column_stack((y_values_file1, z_values_file1))
X_interpolated = griddata(points_file1, x_values_file1, (Y, Z), method='linear')

# Create the contour plot
contour = plt.contourf(Y, Z, X_interpolated, cmap='coolwarm', levels=15)
plt.colorbar(contour, label='intermembrane space side membrane surface altitude')
plt.xlabel('Y')
plt.ylabel('Z')
plt.title('Contour plot of membrane around the type O supercomplex\n', fontsize=24)
# Circle each 2D area with dashed lines
## Circle each 2D area with dashed lines
#for c in contour.collections:
#    c.set_edgecolor("black")  # Set edge color to black for dashed lines
#    c.set_linestyle('dashed')
#    plt.plot(*zip(*c.get_paths()[0].vertices), color='black', linewidth=1, linestyle='dashed')
#

# Load points from file 2 (assuming the file contains one [x, y, z] coordinate per line)
points_file2 = np.loadtxt(file2)
x_values_file2 = points_file2[:, 0]
y_values_file2 = points_file2[:, 1]
z_values_file2 = points_file2[:, 2]

# Filter points in file 2 based on the X range
mask_file2 = (x_values_file2 >= x_min_prot) & (x_values_file2 <= x_max_prot)
x_values_file2 = x_values_file2[mask_file2]
y_values_file2 = y_values_file2[mask_file2]
z_values_file2 = z_values_file2[mask_file2]

# Plot points from file 2 on top of file 1
plt.scatter(y_values_file2, z_values_file2, c='white', marker='.', label='')

# Create white squares around each point in file 2
square_size = 6  # Size of the white square
for x, y, z in zip(x_values_file2, y_values_file2, z_values_file2):
    rect = Rectangle((y - square_size / 2, z - square_size / 2), square_size, square_size, color='white', alpha=0.5)
    plt.gca().add_patch(rect)

plt.legend()
plt.show()


