import numpy as np
import matplotlib.pyplot as plt
import json
import mrcfile
from cryosparc.dataset import Dataset
from .bezierfit import generate_curve_within_boundaries
from scipy.ndimage import zoom

def getzoomfactor(template_csfile, particle_csfile):
    template_dset = Dataset.load(template_csfile)
    particle_dset = Dataset.load(particle_csfile)
    template_shape = np.array(template_dset['blob/shape'])[0]
    particle_shape = np.array(particle_dset['blob/shape'])[0]
    zoom_factor = particle_shape[0] / template_shape[0]
    return zoom_factor

def readimage(template_csfile, json_file):
    with open(json_file, 'r') as f:
        control_points_with_idx = json.load(f)
    template_dset = Dataset.load(template_csfile)
    template_filenames_list = np.array(template_dset['blob/path'])
    template_filenames_list_norepeat = np.unique(template_filenames_list)
    if len(template_filenames_list_norepeat) > 1:
        raise ValueError('Multiple templates in the template cs file')
    idx_with_control_points_and_template = {}
    for idx, control_points in control_points_with_idx.items():
        idx_with_control_points_and_template[idx] = {'control_points': np.array(control_points), 'template': template_filenames_list_norepeat[0]}
    return idx_with_control_points_and_template

def visualize(control_points, image, step=0.01):
    fitted_curve_points, t_values = generate_curve_within_boundaries(control_points, image.shape, step)
    plt.imshow(image, cmap='gray')
    plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'r-')
    plt.plot(control_points[:, 0], control_points[:, 1], 'g.')
    plt.show()

def main(template_csfile, particle_csfile, json_file, idx):
    template_dset = Dataset.load(template_csfile)
    template_idx_list = np.array(template_dset['blob/idx'])
    zoom_factor = getzoomfactor(template_csfile, particle_csfile)
    if idx in template_idx_list:
        idx_with_control_points_and_template = readimage(template_csfile, json_file)
        control_points = idx_with_control_points_and_template[str(idx)]['control_points']
        template = idx_with_control_points_and_template[str(idx)]['template']
        with mrcfile.open(template, permissive=True) as mrc:
            template_original = mrc.data[idx]
        template_original = np.array(template_original)
        template_image = zoom(template_original, zoom_factor)
        print(control_points)
        visualize(control_points, template_image)
    else:
        raise ValueError('Invalid idx')

main('/data2/J729/templates_selected.cs', '/data2/J729/particles_selected.cs', '/data2/control_points.json', 29)