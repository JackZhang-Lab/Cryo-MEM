from cryosparc.dataset import Dataset
import numpy as np
import glob
import json
import os
import mrcfile
from subtraction import MembraneSubtract
import time
import multiprocessing

folder_path = 'J713/extract'
# make a new dictionary
subtracted_folder_path = folder_path.replace('extract', 'subtracted')
os.makedirs(subtracted_folder_path,exist_ok=True)

dset = Dataset.load('/data2/J729/particles_selected.cs')
particle_filenames_list = np.array(dset['blob/path'])
particle_idx_list = np.array(dset['blob/idx'])
psi_list = np.array(dset['alignments2D/pose'])
pixel_size_list = np.array(dset['blob/psize_A'])
shift_list = np.array(dset['alignments2D/shift'])
class_list = np.array(dset['alignments2D/class'])

# read 'control_points.json' file into a dictionary
with open('control_points.json', 'r') as f:
    control_points = json.load(f)
name_pattern = '*.mrc'
file_name_list = glob.glob(f'{folder_path}/{name_pattern}')


def membrane_subtract(particle_filename):
    start_time = time.time()
    with mrcfile.open(particle_filename, permissive=True) as mrc:
        particle_stack = mrc.data
    subtracted_particle_stack = particle_stack.copy()
    mask = (particle_filenames_list == particle_filename)
    particle_idxes = particle_idx_list[mask]
    psis = psi_list[mask]
    pixel_sizes = pixel_size_list[mask]
    shifts = shift_list[mask]
    classes = class_list[mask]
    for particle_idx, psi, pixel_size, shift, class_ in zip(particle_idxes, psis, pixel_sizes, shifts, classes):
        # class_得根据control_points.json找到control_points
        subtractor = MembraneSubtract(class_, particle_stack[particle_idx], psi, shift[0], shift[1], pixel_size).subtract()
        subtracted_particle_stack[particle_idx] = subtractor.subtract()
    with mrcfile.new(particle_filename.replace('/extract/', '/subtracted/'), overwrite=True) as mrc:
        mrc.set_data(subtracted_particle_stack)
    end_time = time.time()
    print(f'{particle_filename} finished in {end_time - start_time} seconds')

def multiprocess_membrane_subtract(file_name_list, num_cpu):
    pool = multiprocessing.Pool(processes=num_cpu)
    pool.map(membrane_subtract, file_name_list)
    pool.close()
    pool.join()

