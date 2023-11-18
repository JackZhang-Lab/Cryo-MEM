import mrcfile
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the MRC file path and the subsampling factor
mrc_file_path = 'typeO_membrane.mrc'
subsample_factor = 2  # Adjust as needed

# Read the MRC file
with mrcfile.open(mrc_file_path, mode='r') as mrc:
    density_map = mrc.data
    x_dim, y_dim, z_dim = density_map.shape

# Set the threshold for binarization
threshold = 0.03  # Adjust as needed

# Create a binary map based on the threshold
binary_map = np.where(density_map >= threshold, 1, 0)

# Generate points in 3D space and subsample based on the factor
points = []
for x in range(x_dim):
    for y in range(y_dim):
        for z in range(z_dim):
            if binary_map[x, y, z] == 1:
                # Only add the point if it passes the subsampling check
                if (x % subsample_factor == 0) and (y % subsample_factor == 0) and (z % subsample_factor == 0):
                    points.append((x, y, z))

# Define the output file path
output_file_path = 'mem.cor'

# Write the coordinates of the selected points to the output file
with open(output_file_path, 'w') as output_file:
    for point in points:
        x, y, z = point
        output_file.write(f"{x} {y} {z}\n")

print(f"Coordinates of {len(points)} points have been saved to {output_file_path}.")

## Create a 3D plot of the selected points
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
## Extract x, y, and z values from the selected points
#x_values = [point[0] for point in points]
#y_values = [point[1] for point in points]
#z_values = [point[2] for point in points]
#
## Plot the points
#ax.scatter(x_values, y_values, z_values, c='b', marker='o', label='Selected Points')
#
## Set labels and title
#ax.set_xlabel('X')
#ax.set_ylabel('Y')
#ax.set_zlabel('Z')
#plt.title('3D Plot of Selected Points')
#
## Show the plot
#plt.show()
#
