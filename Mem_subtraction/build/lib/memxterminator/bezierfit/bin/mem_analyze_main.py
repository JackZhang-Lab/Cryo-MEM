from cryosparc.dataset import Dataset
import numpy as np
import cupy as cp
from cupyx.scipy.ndimage import zoom
import mrcfile
from ..lib.bezierfit import Coarsefit, GA_Refine
import json
import argparse
import multiprocessing

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--template', type=str, default='./templates_selected.cs')
    parser.add_argument('--particle', type=str, default='./particles_selected.cs')
    parser.add_argument('--output', type=str, default='./control_points.json')
    parser.add_argument('--num_points', type=int, default=600)
    parser.add_argument('--degree', type=int, default=3)
    parser.add_argument('--coarsefit_iter', type=int, default=300)
    parser.add_argument('--coarsefit_cpus', type=int, default=20)
    parser.add_argument('--cur_penalty_thr', type=float, default=0.05)
    parser.add_argument('--dithering_range', type=int, default=50)
    parser.add_argument('--refine_iter', type=int, default=700)
    parser.add_argument('--refine_cpus', type=int, default=12)
    parser.add_argument('--physical_membrane_dist', type=int, default=35)
    args = parser.parse_args()

    template_dset = Dataset.load(args.template)
    particle_dset = Dataset.load(args.particle)

    template_filenames_list = np.array(template_dset['blob/path'])
    template_idx_list = np.array(template_dset['blob/idx'])
    particle_pixel_szs = np.array(particle_dset['blob/psize_A'])
    template_pixel_szs = np.array(template_dset['blob/psize_A'])
    template_shape = np.array(template_dset['blob/shape'])[0]
    particle_shape = np.array(particle_dset['blob/shape'])[0]
    zoom_factor = particle_shape[0] / template_shape[0]

    template_pixel_szs = particle_pixel_szs

    control_points_with_idx = {}

    multiprocessing.set_start_method('spawn', force=True)
    for template_filename, template_idx, template_pixel_sz in zip(template_filenames_list, template_idx_list, template_pixel_szs):
        with mrcfile.open(template_filename, permissive=True) as mrc:
            template_original = mrc.data[template_idx]
            template_original = cp.array(template_original)
            template_image = zoom(template_original, zoom_factor)
            template_image_cpu = template_image.get()
            # template_image_cpu = mrc.data[template_idx]
        coarsefit = Coarsefit(template_image_cpu, args.num_points, args.degree, args.coarsefit_iter, args.coarsefit_cpus)
        initial_control_points = coarsefit()
        ga_refine = GA_Refine(template_image_cpu, template_pixel_sz, args.cur_penalty_thr, args.dithering_range, args.refine_iter, args.refine_cpus, args.physical_membrane_dist)
        refined_control_points = ga_refine(initial_control_points, template_image_cpu)
        # save refined control points and template_idx into a dictionary and save it as a JSON file
        control_points_with_idx[str(template_idx)] = refined_control_points.tolist()
        print(f'Template {template_filename} {template_idx} done, control points: {refined_control_points}')
        with open(args.output, 'w') as f:
            json.dump(control_points_with_idx, f)
            f.write('\n')
        cp.cuda.MemoryPool().free_all_blocks()
