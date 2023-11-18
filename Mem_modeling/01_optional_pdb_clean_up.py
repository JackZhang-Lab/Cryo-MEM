#!/usr/bin/env python
# coding: utf-8

# # Lipids Clean-Up
# This script can be used after discarding lipids occupied protein's position.  
# To discard lipids, use select zone command in Chimearx and set the discarding range as u like.

# In[117]:


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
                lipids_id = int(line[22:26].strip())

                data_list.append({
                    'Lipid_ID': lipids_id,
                    'Atom_Type': atom_type,
                    'X': float(x1),
                    'Y': float(y1),
                    'Z': float(z1)
                })
    df = pd.DataFrame(data_list)
    df = df.reset_index(drop=True)
    return df

# In[118]:


filename = '/home/boyuan/Desktop/mmbr_fit/input/class01/class01_down_lipids.pdb'
original_coords = parse_pdb(filename)

# In[ ]:


# fetch the lipid
lip_id = original_coords["Lipid_ID"].unique()

# remove non-intact lipids
cleaned_coords = original_coords.copy()
cleaned_coords = cleaned_coords.groupby("Lipid_ID").filter(lambda x: len(x) == 52)

# update lipids id
unique_lip_ids = cleaned_coords["Lipid_ID"].unique()
new_lip_ids = range(1, len(unique_lip_ids) + 1)
lipid_id_mapping = dict(zip(unique_lip_ids, new_lip_ids))
cleaned_coords["Lipid_ID"] = cleaned_coords["Lipid_ID"].replace(lipid_id_mapping)
cleaned_coords = cleaned_coords.reset_index(drop=True)

cleaned_coords.index

# In[119]:


template = "HETATM{:>5} {:<4} POPC {:>4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00           {}"
# Modify PDB format into hybrid-36 encoding
pdb_output = "\n".join(
    template.format(index+1, row['Atom_Type'], row['Lipid_ID'], row['X'], row['Y'], row['Z'], row['Atom_Type'][0]) 
    for index, row in cleaned_coords.iterrows()
)
filename = "/home/boyuan/Desktop/mmbr_fit/input/class01/cleaned_lipids_down.pdb"

# In[122]:


with open(filename, 'w') as file:
    file.write(pdb_output)

# # PO4 Clean Up
# I use this script for removing non-intact PO4 mols. It may help to perform further refinement in COOT. 

# In[137]:


# Load PO4 original coords
import pandas as pd

def parse_pdb(file_name):
    data_list = []
    with open(file_name, 'r') as file:
        for line in file:
            if line[0:6] == "HETATM":
                atom_type = line[13:16].strip()
                po4_id = int(line[22:26].strip())
                x1 = line[30:38].strip()
                y1 = line[38:46].strip()
                z1 = line[46:54].strip()
                data_list.append({
                    'PO4_ID': po4_id,
                    'Atom_Type': atom_type,
                    'X': float(x1),
                    'Y': float(y1),
                    'Z': float(z1)
                })
    df = pd.DataFrame(data_list)
    df = df.reset_index(drop=True)
    return df

# In[144]:


import os

directory_path = '/home/boyuan/Desktop/mmbr_fit/input/class00/'
all_files_wait4clean = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f)) and f.startswith("trim_fitted_up_")]

for file_name in all_files_wait4clean:
    filepath = os.path.join(directory_path,file_name)
    data_list = []
    with open(filepath, "r") as f:
        for line in f:
            if line[0:6] == "HETATM":
                atom_type = line[13:16].strip()
                po4_id = int(line[22:26].strip())
                x1 = line[31:38].strip()
                y1 = line[39:46].strip()
                z1 = line[47:54].strip()
                data_list.append({
                    'PO4_ID': po4_id,
                    'Atom_Type': atom_type,
                    'X': float(x1),
                    'Y': float(y1),
                    'Z': float(z1)
                })
    ori_po4_coords = pd.DataFrame(data_list)
    p_id = ori_po4_coords["PO4_ID"].unique()

    # remove non-intact PO4 mols
    cleaned_po4_coords = ori_po4_coords.copy()
    cleaned_po4_coords = cleaned_po4_coords.groupby("PO4_ID").filter(lambda x: len(x)==5)

    cleaned_po4_coords = cleaned_po4_coords.reset_index(drop=True)

    unique_p_ids = cleaned_po4_coords["PO4_ID"].unique()
    new_p_ids = range(1, len(unique_p_ids) + 1)
    p_id_mapping = dict(zip(unique_p_ids, new_p_ids))
    cleaned_po4_coords["PO4_ID"] = cleaned_po4_coords["PO4_ID"].replace(p_id_mapping)
    
    # write_updated files
    template = "HETATM {:>4}  {:<2}  PO4 B{:<4}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  30.00           {}"
    po4_output = "\n".join(
        template.format(index+1, row['Atom_Type'], row['PO4_ID'], row['X'], row['Y'], row['Z'], row['Atom_Type'][0]) 
        for index, row in cleaned_po4_coords.iterrows()
    )
    newfilename = os.path.join(directory_path,'cleaned_t_f_up'+file_name[15:])
    with open(newfilename, 'w') as file:
        file.write(po4_output)
        file.write("\n"+"END")
