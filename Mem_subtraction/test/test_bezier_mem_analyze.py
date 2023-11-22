from memxterminator.bezierfit.bin import mem_analyze_main

from cryosparc.dataset import Dataset

import numpy as np

particle_csfile = '/data2/J729/particles_selected.cs'
passthrough_csfile = '/data3/Zhen/VSV/membrane_subtraction/VSV/J419/J419_passthrough_particles_selected.cs'
particle_dset = Dataset.load(particle_csfile)
pass_dset = Dataset.load(passthrough_csfile)
particle_filenames_list = np.array(particle_dset['blob/path'])
particle_idx_list = np.array(particle_dset['blob/idx'])
psi_list = np.array(particle_dset['alignments2D/pose'])
pixel_size_list = np.array(particle_dset['blob/psize_A'])
shift_list = np.array(particle_dset['alignments2D/shift'])
class_list = np.array(particle_dset['alignments2D/class'])

mask = (particle_filenames_list == 'J713/extract/2023-04-29_18.02.45_mito4_28_00025_X-1Y-1-1_patch_aligned_doseweighted_particles.mrc')
particle_idxes = particle_idx_list[mask]
psis = psi_list[mask]
pixel_sizes = pixel_size_list[mask]
shifts = shift_list[mask]
classes = class_list[mask]

print(pass_dset)