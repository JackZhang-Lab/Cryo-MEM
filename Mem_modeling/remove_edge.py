import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# Step 1: Read 3D points from a file (replace 'your_points_file.txt' with your file path)
points_3d = np.loadtxt('surface.cor')

# Define the number of iterations (N)
N = 95  # Change this to the desired number of iterations

for iteration in range(N):
    # Step 2: Project 3D points onto the YZ plane to get 2D points (y, z)
    points_2d = points_3d[:, 1:3]  # Assuming your file has columns x, y, z

    # Step 3: Use the Convex Hull algorithm to find the border line
    hull = ConvexHull(points_2d)

    # Step 4: Get the vertices of the convex hull (border points)
    border_points_2d = points_2d[hull.vertices]

    # Step 5: Print the number of points before and after removing border points
    num_points_before = len(points_3d)
    points_3d = np.delete(points_3d, hull.vertices, axis=0)
    num_points_after = len(points_3d)
    points_removed = num_points_before - num_points_after

    print(f"Iteration {iteration+1}: Removed {points_removed} points")

# Step 6: After N iterations, you have the final set of points with the border points removed

# Step 7: Find the corresponding 3D points for the final 2D points
corresponding_3d_points = []

for point_2d in points_2d:
    # Find the nearest 2D point in the original set
    nearest_idx = np.argmin(np.linalg.norm(points_2d - point_2d, axis=1))
    
    # Check if nearest_idx is within the valid range
    if nearest_idx < len(points_3d):
        corresponding_3d_points.append(points_3d[nearest_idx])

# The list 'corresponding_3d_points' contains the corresponding 3D points for the final 2D points

# Step 8: Save the corresponding 3D points to a file
np.savetxt('edge_removed.cor', corresponding_3d_points)

# Step 9: Plot the results after N rounds are finished
# Plot the original points and the border line
plt.scatter(points_2d[:, 0], points_2d[:, 1], c='b', label='Original Points')
plt.plot(border_points_2d[:, 0], border_points_2d[:, 1], c='r', lw=2, label='Border Line')
plt.xlabel('Y')
plt.ylabel('Z')
plt.legend()
plt.title(f'Border Line Enclosing Points (After {N} Rounds)')
plt.show()
