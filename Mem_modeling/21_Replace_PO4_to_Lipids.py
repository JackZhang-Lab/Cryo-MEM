#!/usr/bin/env python
# coding: utf-8

# In[723]:


import os
import numpy as np

def load_coords(path):
    data_list = []
    with open(path, 'r') as f:
        lines = f.readlines()
        matrix = []
        for line in lines:
            x, y, z = map(float, line.strip().split(","))
            matrix.append([x, y, z])
        coords = np.copy(matrix).astype(np.float32)
        x = coords[:, 0]-200
        y = coords[:, 1]-200
        z = coords[:, 2]
        
        for xi, yi, zi in zip(x, y, z):
            data_list.append({
                        'X': float(xi),
                        'Y': float(yi),
                        'Z': float(zi)
                    })
        
    df = pd.DataFrame(data_list)
    df.insert(0, "Lipid_ID", range(1, len(df)+1))
    df["Lipid_ID"] = df["Lipid_ID"].astype(float)
    df['X'] = df['X'].round(3)
    df['Y'] = df['Y'].round(3)
    df['Z'] = df['Z'].round(3)
    return x, y, z, df

filepath = "/home/boyuan/Desktop/mmbr_fit/process/lipids_up_p_target_coordinates.txt"
x, y, z, target_coords = load_coords(filepath)

x_min, x_max = x.min(), x.max()
y_min, y_max = y.min(), y.max()
print(y_max)

# In[720]:


target_p_coords = target_coords["Lipid_ID"]

# In[695]:


# Target Coords
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')
# ax.view_init(30, 30)
ax.scatter(x, y, z, color='green', s=1, label='P atoms Points')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
plt.show()

# In[725]:


# Load original coords
import pandas as pd

def parse_pdb(file_name):
    data_list = []
    with open(file_name, 'r') as file:
        for line in file:
            if line[0:6] == "HETATM":
                x1 = line[30:38].strip()
                y1 = line[38:46].strip()
                z1 = line[46:54].strip()
                atom_type = line[12:16].strip()
                lipids_id = int(line[22:26].strip())-2343

                data_list.append({
                    'Lipid_ID': lipids_id,
                    'Atom_Type': atom_type,
                    'X': float(x1),
                    'Y': float(y1),
                    'Z': float(z1)
                })
    df = pd.DataFrame(data_list)
    lipids_to_remove = df[(df['Atom_Type'] == 'P') & (df['Z'].astype(float) > 0)]['Lipid_ID'].unique()
    df = df[~df['Lipid_ID'].isin(lipids_to_remove)]
    df = df.reset_index(drop=True)
    return df

# In[726]:


filename = '/home/boyuan/Desktop/mmbr_fit/input/double_mmbr_square.pdb'
original_coords = parse_pdb(filename)

# In[729]:


def update_lipid_position(original_coords, target_coords):
        
    translated_df_new = pd.DataFrame(columns=["Lipid_ID", "Atom_Type", "X", "Y", "Z"])
    for lipid_id in original_coords["Lipid_ID"].unique():
        original_p_coords = original_coords[(original_coords["Lipid_ID"] == lipid_id) & (original_coords["Atom_Type"] == "P")][["X", "Y", "Z"]].values[0]
        target_p_coords = target_coords[(target_coords["Lipid_ID"] == lipid_id)][["X","Y","Z"]].values[0]
        translation_vector =  original_p_coords - target_p_coords
        sub_df = original_coords[original_coords["Lipid_ID"] == lipid_id].copy()
        sub_df["X"] += translation_vector[0]
        sub_df["Y"] += translation_vector[1]
        sub_df["Z"] += translation_vector[2]

        translated_df_new = pd.concat([translated_df_new, sub_df])
        translated_df_new['X'] = translated_df_new['X'].round(3)
        translated_df_new['Y'] = translated_df_new['Y'].round(3)
        translated_df_new['Z'] = translated_df_new['Z'].round(3)

    return translated_df_new

updated_coords = update_lipid_position(original_coords, target_coords)

# In[730]:


template = "HETATM{:<5} {:<4} POPC {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00           {}"

# In[647]:


def int_to_corrected_hybrid36(num):
    if num < 0:
        raise ValueError("Negative numbers are not supported for hybrid36 encoding.")
    
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # For numbers less than 100000, return the number as a string
    if num < 100000:
        return str(num)
    
    # For numbers greater than or equal to 100000
    num -= 100000
    result = []
    
    # Extract the prefix letter
    num, remainder = divmod(num, 36**4)
    result.append(digits[num])
    
    # Extract remaining digits
    for _ in range(4):
        num, remainder = divmod(remainder, 36)
        result.append(digits[remainder])
    
    return ''.join(result)


# In[731]:


# Modify PDB format into hybrid-36 encoding
pdb_output_corrected_hybrid36 = "\n".join(template.format(int_to_corrected_hybrid36(index+1),row['Atom_Type'], row['Lipid_ID'], row['X'], row['Y'], row['Z'], row['Atom_Type'][0]) for index, row in updated_coords.iterrows())

# In[732]:


# Output to a file
filename = "/home/boyuan/Desktop/mmbr_fit/output/pdb/lipid_down.pdb"
with open(filename, 'w') as file:
    file.write(pdb_output_corrected_hybrid36)
