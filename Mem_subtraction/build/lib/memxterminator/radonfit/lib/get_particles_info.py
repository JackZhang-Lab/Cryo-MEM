import pandas as pd
import numpy as np
import cupy as cp
import starfile
import mrcfile
from ._utils import *
from cupyx.scipy.ndimage import zoom
# from mem_subtract_old import *

df_membrane_analysis_star = starfile.read('mem_analysis.star')
average_raw_series = df_membrane_analysis_star['rln2DAverageimageName'].apply(lambda x: x.split('@')[1]), df_membrane_analysis_star['rln2DAverageimageName'].apply(lambda x: x.split('@')[0])
average_raw_lst = list(map(lambda x: (x[0], int(x[1])), list(zip(*average_raw_series))))
print(average_raw_lst)
average_raw = readmrc('neuron_templates.mrc', section=4, mode='gpu')
average_raw = zoom(average_raw, 256/64)
average_2d = readmrc('neuron_templates_averaged_4.mrc', section=0, mode='gpu')
membrane_mask = readmrc('neuron_templates_mask_4.mrc', section=0, mode='gpu')
x0 = df_membrane_analysis_star.loc[0, 'rlnCenterX']
y0 = df_membrane_analysis_star.loc[0, 'rlnCenterY']


starfile_name = 'particles_selected.star'
df_star = starfile.read(starfile_name)
df_star_selected = df_star[df_star['rlnClassNumber'] == 96]
rawimage_series = df_star_selected['rlnImageName'].apply(lambda x: x.split('@')[1]), df_star_selected['rlnImageName'].apply(lambda x: x.split('@')[0]), df_star_selected['rlnAnglePsi'], df_star_selected['rlnOriginX'], df_star_selected['rlnOriginY']
rawimage_lst = list(map(lambda x: (x[0], int(x[1]), x[2], x[3], x[4]), list(zip(*rawimage_series))))
# rawimage_lst = rawimage_lst[20:40]
# print([item[0] for item in rawimage_lst])
filenames = np.array([item[0] for item in rawimage_lst])
sections = np.array([item[1] for item in rawimage_lst])
values = cp.asarray([list(item[2:]) for item in rawimage_lst])

rawimage_dict = {}
for item in rawimage_lst:
    if item[0] not in rawimage_dict:
        rawimage_dict[item[0]] = [item[1]]
    else:
        rawimage_dict[item[0]].append(item[1])

# def get_raw_params(filename, section):
#     mask = (filenames == filename) & (sections == section)
#     matched_values = values[mask]
#     return matched_values[0].tolist() if matched_values.size > 0 else None
# for key, section_lst in rawimage_dict.items():
#     image_stack = cp.asarray(mrcfile.open(key).data[list(map(lambda x: x-1, section_lst))].copy())
#     subtracted_image_stack = cp.zeros_like(image_stack)
#     print(key)
#     print(section_lst)
#     i = 0
#     for section in section_lst:
#         print('section', section)
#         raw_params = get_raw_params(key, section)
#         print('raw_params', raw_params)
#         rawimage = image_stack[i]
#         psi = raw_params[0]
#         dx = raw_params[1]
#         dy = raw_params[2]
#         get_to_raw = get2raw(average_raw, average_2d, rawimage, membrane_mask, x0, y0, psi, dx, dy)
#         average_2d_s, average_raw_s, mask = get_to_raw.rotate_average_to_raw()
#         subtracted_image = get_to_raw.raw_membrane_average_subtract()
#         subtracted_image_stack[i] = subtracted_image
#         i += 1
#     # savemrc(cp.asnumpy(subtracted_image_stack), key.replace('.mrc', '_subtracted.mrc'))