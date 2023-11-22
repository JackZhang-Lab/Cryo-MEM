from cryosparc.dataset import Dataset
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import zoom
import glob
import json
import os
import mrcfile
import time
import multiprocessing
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--particle', type=str, default='./particles_selected.cs')
parser.add_argument('--template', type=str, default='./templates_selected.cs')
parser.add_argument('--control_points', type=str, default='./control_points.json')
parser.add_argument('--points_step', type=float, default=0.005)
parser.add_argument('--physical_membrane_dist', type=int, default=35)

args = parser.parse_args()

particle_dset = Dataset.load(args.particle)
particle_filenames_list = np.array(particle_dset['blob/path'])
particle_idx_list = np.array(particle_dset['blob/idx'])
psi_list = np.array(particle_dset['alignments2D/pose'])
pixel_size_list = np.array(particle_dset['blob/psize_A'])
shift_list = np.array(particle_dset['alignments2D/shift'])
class_list = np.array(particle_dset['alignments2D/class'])

template_dset = Dataset.load(args.template)
template_shape = np.array(template_dset['blob/shape'])[0]
particle_shape = np.array(particle_dset['blob/shape'])[0]


def get_folder_path_list(particle_filenames_list):
    folder_path_list = []
    for particle_filename in particle_filenames_list:
        folder_path_list.append(os.path.dirname(particle_filename))
    return folder_path_list

def mkdir_subtracted_folder_path(folder_path_list):
    folder_path_list_norepeat = list(set(folder_path_list))
    subtracted_folder_path_list = []
    for folder_path in folder_path_list_norepeat:
        subtracted_folder_path_list.append(folder_path.replace('extract', 'subtracted'))
    
    for subtracted_folder_path in subtracted_folder_path_list:
        os.makedirs(subtracted_folder_path,exist_ok=True)

mkdir_subtracted_folder_path(get_folder_path_list(particle_filenames_list))
# read 'control_points.json' file into a dictionary
with open(args.control_points, 'r') as f:
    control_points_dict = json.load(f)
# name_pattern = '*.mrc'
# file_name_list = glob.glob(f'{folder_path}/{name_pattern}')

def membrane_subtract(particle_filename):
    from ..lib.subtraction import MembraneSubtract
    start_time = time.time()
    with mrcfile.open(particle_filename, permissive=True) as mrc:
        particle_stack = mrc.data
    particle_stack = cp.array(particle_stack)
    # particle_stack = particle_stack.get()
    subtracted_particle_stack = particle_stack.copy()
    mask = (particle_filenames_list == particle_filename)
    particle_idxes = particle_idx_list[mask]
    psis = psi_list[mask]
    pixel_sizes = pixel_size_list[mask]
    shifts = shift_list[mask]
    classes = class_list[mask]
    for particle_idx, psi, pixel_size, shift, class_ in zip(particle_idxes, psis, pixel_sizes, shifts, classes):
        # class_得根据control_points.json找到control_points
        # if str(class_) in control_points_dict:
        control_points = np.array(control_points_dict[str(class_)])
        # print(control_points)
        # print(f'particle_filename:{particle_filename}, particle_index:{particle_idx}, psi: {psi}, pixel_size: {pixel_size}, shift: {shift}, class: {class_}')
        if particle_stack.ndim == 2:
            subtractor = MembraneSubtract(control_points, particle_stack, psi, shift[0], shift[1], pixel_size, args.points_step, args.physical_membrane_dist)
            subtracted_particle_stack = subtractor.mem_subtract()
        elif particle_stack.ndim == 3:
            subtractor = MembraneSubtract(control_points, particle_stack[particle_idx], psi, shift[0], shift[1], pixel_size, args.points_step, args.physical_membrane_dist)
            subtracted_particle_stack[particle_idx] = subtractor.mem_subtract()
        # print(f'{particle_idx}@{particle_filename} finished')
    # subtracted_particle_stack = cp.array(subtracted_particle_stack)
    # subtracted_particle_stack = subtracted_particle_stack.get()
    with mrcfile.new(particle_filename.replace('/extract/', '/subtracted/'), overwrite=True) as mrc:
        mrc.set_data(subtracted_particle_stack.get())
    end_time = time.time()
    print(f'{particle_filename}, {len(particle_idxes)} particles finished in {end_time - start_time} seconds')
    del subtracted_particle_stack
    del particle_stack
    del MembraneSubtract
    del mrc
    # release memory
    cp.cuda.MemoryPool().free_all_blocks()

def remove_duplicates_preserve_order(seq):
    return list(dict.fromkeys(seq))

def process_membrane_subtract(file_name_list):
    file_name_list = remove_duplicates_preserve_order(file_name_list)
    for file_name in file_name_list:
        membrane_subtract(file_name)


process_membrane_subtract(particle_filenames_list)

# def multiprocess_membrane_subtract(file_name_list, num_cpu):
#     pool = multiprocessing.Pool(processes=num_cpu)
#     pool.map(membrane_subtract, file_name_list)
#     pool.close()
#     pool.join()

