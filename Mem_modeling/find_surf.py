def read_points_from_file(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                x, y, z = map(float, parts)
                points.append((x, y, z))
    return points

def filter_largest_x(points, threshold=0.1): 
    filtered_points = {}
    
    for x, y, z in points:
        key = (round(y, 2), round(z, 2))  # Round y and z to 2 decimal places for comparison
        
        if key not in filtered_points or x > filtered_points[key][0]:
            filtered_points[key] = (x, y, z)
    
    return list(filtered_points.values())

def main(input_filename, output_filename):
    points = read_points_from_file(input_filename)
    
    if not points:
        print("No valid points found in the input file.")
        return
    
    filtered_points = filter_largest_x(points)
    
    with open(output_filename, 'w') as output_file:
        for x, y, z in filtered_points:
            output_file.write(f"{x} {y} {z}\n")
    
    print(f"Filtered points have been saved to '{output_filename}'.")

if __name__ == "__main__":
    input_file = "scaled_mem.cor"  # Replace with your input file name
    output_file = "surface.cor"  # Replace with your desired output file name
    main(input_file, output_file)

