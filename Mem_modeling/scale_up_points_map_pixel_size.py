scale = 1.17
import numpy as np
points_file_path = 'prot.cor'


with open(points_file_path, 'r') as points_file:
    points_lines = points_file.readlines()

scaled_up_coor = []
for line in points_lines:
    xyz = list(line.strip().split())
    scaled_xyz = [(float(i) * scale) for i in xyz ]
    scaled_up_coor.append(scaled_xyz)
with open("scaled_points.coor", 'w') as output_file:
    for point in scaled_up_coor:
        output_file.write(f'{point[0]} {point[1]} {point[2]}\n')

#np_array= np.array(points_data)
#new_array = np_array* scale
#print(np_array)
#print(new_array)
