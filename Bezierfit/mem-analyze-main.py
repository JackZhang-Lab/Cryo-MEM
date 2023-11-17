from cryosparc.dataset import Dataset
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import zoom
import mrcfile
from bezierfit import Coarsefit, GA_Refine
import json

template_dset = Dataset.load('/data2/J729/templates_selected.cs')
particle_dset = Dataset.load('/data2/J729/particles_selected.cs')
template_filenames_list = np.array(template_dset['blob/path'])
template_idx_list = np.array(template_dset['blob/idx'])
template_pixel_szs = np.array(template_dset['blob/psize_A'])
template_shape = np.array(template_dset['blob/shape'])[0]
particle_shape = np.array(particle_dset['blob/shape'])[0]
zoom_factor = particle_shape[0] / template_shape[0]

control_points_with_idx = {}
for template_filename, template_idx, template_pixel_sz in zip(template_filenames_list, template_idx_list, template_pixel_szs):
    with mrcfile.open(template_filename, permissive=True) as mrc:
        template_original = mrc.data[template_idx]
    template_original = cp.array(template_original)
    template_image = zoom(template_original, zoom_factor)
    template_image_cpu = template_image.get()
    coarsefit = Coarsefit(template_image_cpu, 600, 3, 300, 20)
    initial_control_points = coarsefit()
    ga_refine = GA_Refine(template_image_cpu, template_pixel_sz, 0.05, 50, 700, 18)
    refined_control_points = ga_refine(initial_control_points, template_image_cpu)
    # save refined control points and template_idx into a dictionary and save it as a JSON file
    control_points_with_idx[template_idx] = refined_control_points.tolist()
    print(f'Template {template_filename} {template_idx} done, control points: {refined_control_points}')

with open('control_points.json', 'a') as f:
    json.dump(control_points_with_idx, f)
    f.write('\n')
