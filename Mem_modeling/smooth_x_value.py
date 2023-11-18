import numpy as np

# Define the valid range for y and z coordinates
y_min, y_max = 0, 600
z_min, z_max = 0, 600

# Read points from the file
points = []
try:
    with open('edge_removed.cor', 'r') as file:
        for line in file:
            x, y, z = map(float, line.strip().split())
            # Check if y and z coordinates are within the valid range
            if y_min <= y <= y_max and z_min <= z <= z_max:
                points.append((x, y, z))
            else:
                print(f"Warning: Point ({x}, {y}, {z}) is outside the valid range. Skipping this point.")
except FileNotFoundError:
    print("Error: File not found. Please make sure 'points.txt' exists in the current directory.")
    exit()
except ValueError:
    print("Error: Invalid data format in 'points.txt'. Please ensure each line contains three numbers.")
    exit()

# Define the number of grids (N)
N = 70  # You can adjust this value based on your requirements

# Calculate grid boundaries
y_grid_size = (y_max - y_min) / N
z_grid_size = (z_max - z_min) / N

# Calculate average x-values within each grid
average_x_values = np.zeros((N, N))
counts = np.zeros((N, N))

for x, y, z in points:
    try:
        grid_y = int(np.clip((y - y_min) / y_grid_size, 0, N - 1))
        grid_z = int(np.clip((z - z_min) / z_grid_size, 0, N - 1))
        average_x_values[grid_y][grid_z] += x
        counts[grid_y][grid_z] += 1
    except IndexError:
        print(f"Error: Point ({x}, {y}, {z}) falls outside the grid. Skipping this point.")

for i in range(N):
    for j in range(N):
        if counts[i][j] != 0:
            average_x_values[i][j] /= counts[i][j]
        else:
            average_x_values[i][j] = None  # Mark grid cells with no points as None

# Generate new points for non-empty grid cells
new_points = []
for i in range(N):
    for j in range(N):
        if average_x_values[i][j] is not None:
            grid_center_y = y_min + (i + 0.5) * y_grid_size
            grid_center_z = z_min + (j + 0.5) * z_grid_size
            new_x = average_x_values[i][j]
            new_points.append((new_x, grid_center_y, grid_center_z))

# Save new points to a file
with open('smoothed_x_points.cor', 'w') as output_file:
    for x, y, z in new_points:
        output_file.write(f'{x} {y} {z}\n')

print("New points generated and saved successfully.")


