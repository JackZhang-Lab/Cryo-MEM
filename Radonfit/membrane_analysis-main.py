import numpy as np
import cupy as cp
from _utils import *
from radonanalyser import *
from template_centerfitting import *
from calculate_curve import *
from generate_membrane_mask import *
from mem_average import *
from mem_subtract_old import *
from cupyx.scipy.ndimage import zoom
import pandas as pd
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--templates_starfile_name', '-ts', type=str, default='templates_selected.star', help='star file name of templates')
parser.add_argument('--output_filename', '-o', type=str, default='mem_analysis.star', help='output star file name')
parser.add_argument('--particle_starfile_name', '-ps', type=str, default='particles_selected.star', help='star file name of particles')
parser.add_argument('--kappa_template', '-k', type=int, default=False, help='whether to generate kappa templates and which template to use')
parser.add_argument('--kappanum', '-kn', type=int, default=40, help='number of kappa templates')
parser.add_argument('--kappastart', '-ks', type=float, default=-0.008, help='start value of kappa templates')
parser.add_argument('--kappaend', '-ke', type=float, default=0.008, help='end value of kappa templates')
parser.add_argument('--info_json', '-j', type=str, default='radonanalysis_info.json', help='json file name of radon analysis info')
parser.add_argument('--sigma1', '-s1', type=float, default=6, help='sigma1 of template center fitting')
parser.add_argument('--sigma2', '-s2', type=float, default=6, help='sigma2 of template center fitting')
parser.add_argument('--template_size', type=int, default=64, help='template size of template center fitting')
parser.add_argument('--sigma_range', type=int, default=3, help='sigma range of template center fitting')
parser.add_argument('--sigma_step', type=float, default=0.5, help='sigma step of template center fitting')
parser.add_argument('--curve_kappa_start', type=float, default=-0.01, help='start value of kappa in curve fitting')
parser.add_argument('--curve_kappa_end', type=float, default=0.01, help='end value of kappa in curve fitting')
parser.add_argument('--curve_kappa_step', type=float, default=0.0002, help='step of kappa in curve fitting')
parser.add_argument('--edge_sigma', type=float, default=3, help='sigma of edge detection in membrane mask generation')
parser.add_argument('--extra_mem_dist', type=int, default=15, help='extra membrane distance in membrane average')
parser.add_argument('--mem_edge_sigma', type=float, default=5, help='sigma of membrane average')
args = parser.parse_args()

# templates_starfile_name = 'templates_selected.star'
# output_filename = 'mem_analysis.star'
# particle_starfile_name = 'particles_selected.star'

templates_starfile_name = args.templates_starfile_name
output_filename = args.output_filename
particle_starfile_name = args.particle_starfile_name
select_kappa_template = args.kappa_template
kappa_num = args.kappanum
kappa_start = args.kappastart
kappa_end = args.kappaend
info_json = args.info_json
initial_sigma1 = args.sigma1
initial_sigma2 = args.sigma2
template_size = args.template_size
sigma_range = args.sigma_range
sigma_step = args.sigma_step
curve_kappa_start = args.curve_kappa_start
curve_kappa_end = args.curve_kappa_end
curve_kappa_step = args.curve_kappa_step
edge_sigma = args.edge_sigma
extra_mem_dist = args.extra_mem_dist
mem_edge_sigma = args.mem_edge_sigma


def get_zoom_factor(particle_starfile_name, templates_starfile_name):
    def get_image_size(starfile_name):
        df_star = starfile.read(starfile_name)
        image_sample = df_star['rlnImageName'][0].split('@')[1], int(df_star['rlnImageName'][0].split('@')[0])
        image_sample_array = readmrc(image_sample[0], section=image_sample[1]-1, mode='gpu')
        image_size = image_sample_array.shape[0]
        return image_size
    particle_image_size = get_image_size(particle_starfile_name)
    templates_image_size = get_image_size(templates_starfile_name)
    zoom_factor = particle_image_size / templates_image_size
    return zoom_factor

def get_parameters_for_section(i):
    with open(info_json, 'r') as file:
        data = json.load(file)
        section_key = str(i)
        if section_key in data:
            return data[section_key]['crop_rate'], data[section_key]['thr'], data[section_key]['theta_start'], data[section_key]['theta_end']
        else:
            print(f"No data found for section {i}")
            return None, None, None, None

output_df_star = pd.DataFrame(data=[[0] * 10], columns=['rln2DAverageimageName', 'rlnAveragedMembraneName', 'rlnMembraneMaskName', 'rlnCenterX', 'rlnCenterY', 'rlnAngleTheta', 'rlnMembraneDistance', 'rlnSigma1', 'rlnSigma2', 'rlnCurveKappa'])
df_templates_star = starfile.read(templates_starfile_name)
output_df_star = output_df_star.reindex(df_templates_star.index)
output_df_star.fillna(0, inplace=True)
output_df_star['rln2DAverageimageName'] = df_templates_star['rlnImageName']
average_2d_series = df_templates_star['rlnImageName'].apply(lambda x: x.split('@')[1]), df_templates_star['rlnImageName'].apply(lambda x: x.split('@')[0])
average_2d_lst = list(map(lambda x: (x[0], int(x[1])), list(zip(*average_2d_series))))
starfile.write(output_df_star, output_filename, overwrite=True)

zoom_factor = get_zoom_factor(particle_starfile_name, templates_starfile_name)
print('zoom_factor', zoom_factor)
membrane_masks = []
averaged_membranes = []

i = 0
for average2d_filename, section in average_2d_lst:
    image = readmrc(average2d_filename, section=section-1, mode='gpu')
    image = zoom(image, zoom_factor)
    crop_rate_temp, thr_temp, theta_start_temp, theta_end_temp = get_parameters_for_section(i)
    radonanalyze = RadonAnalyzer(output_filename, i, image, crop_rate=crop_rate_temp, thr=thr_temp, theta_start=theta_start_temp, theta_end=theta_end_temp)
    centerfit = Template_centerfitting(output_filename, i, sigma1=initial_sigma1, sigma2=initial_sigma2, image=image, crop_rate=crop_rate_temp, thr=thr_temp, theta_start=theta_start_temp, theta_end=theta_end_temp, template_size=template_size, sigma_range=sigma_range, sigma_step=sigma_step)
    centerfit.centerfinder()
    centerfit.fit_sigma()
    curve_fitting = Curvefitting(output_filename, i, image, kappa_start=curve_kappa_start, kappa_end=curve_kappa_end, kappa_step=curve_kappa_step)
    get_mem_mask = mem_mask(output_filename, i, image, edge_sigma=edge_sigma)
    membrane_mask = get_mem_mask.generate_mem_mask()
    membrane_mask = cp.asnumpy(membrane_mask)
    membrane_masks.append(membrane_mask)
    output_df_star = starfile.read(output_filename)
    output_df_star.loc[i, 'rlnMembraneMaskName'] = f'''{i+1:06d}@{average2d_filename.replace('.mrc', '_masks.mrc')}'''
    mem_average = average_membrane(output_filename, i, image, extra_mem_dist=extra_mem_dist, sigma=mem_edge_sigma, x0=output_df_star.loc[i, 'rlnCenterY'], y0=output_df_star.loc[i, 'rlnCenterX'], theta=output_df_star.loc[i, 'rlnAngleTheta'] * np.pi / 180, membrane_distance=output_df_star.loc[i, 'rlnMembraneDistance'], kappa=output_df_star.loc[i, 'rlnCurveKappa'])
    mem_average.calculate_membrane_average()
    averaged_2d_membrane = mem_average.generate_2d_average_mem()
    if select_kappa_template == i:
        mem_average.kappa_templates_generator(kappa_start=kappa_start, kappa_end=kappa_end, kappa_num=kappa_num)
    averaged_2d_membrane = cp.asnumpy(averaged_2d_membrane)
    averaged_membranes.append(averaged_2d_membrane)
    output_df_star.loc[i, 'rlnAveragedMembraneName'] = f'''{i+1:06d}@{average2d_filename.replace('.mrc', '_averaged.mrc')}'''
    print(f'file {average2d_filename} section {section} finished')
    starfile.write(output_df_star, output_filename, overwrite=True)
    i += 1
membrane_masks = np.array(membrane_masks)
averaged_membranes = np.array(averaged_membranes)
savemrc(membrane_masks, f'''{average_2d_lst[0][0].replace('.mrc', '_masks.mrc')}''')
savemrc(averaged_membranes, f'''{average_2d_lst[0][0].replace('.mrc', '_averaged.mrc')}''')

