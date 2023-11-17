import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from _utils import *
from template_centerfitting import *
from radonanalyser import *
from calculate_curve import *
from generate_membrane_mask import *
from scipy.interpolate import interp1d
from cupyx.scipy.ndimage import zoom

class average_membrane:
    def __init__(self, output_filename, i, image, extra_mem_dist, sigma, x0, y0, theta, membrane_distance, kappa):
        self.i = i
        self.output_filename = output_filename
        self.image = image
        self.x0 = x0
        self.y0 = y0
        self.theta = theta
        self.membrane_distance = membrane_distance
        self.kappa = kappa
        self.sigma = sigma
        self.image_size = image.shape[0]
        self.distance_mat = cp.zeros((self.image_size, self.image_size))
        self.test_img = cp.zeros((self.image_size, self.image_size))
        self.membrane_average_1d = []
        self.membrane_average_2d = cp.zeros((self.image_size, self.image_size))
        self.membrane_dist_lst = [i for i in cp.arange(-self.membrane_distance-extra_mem_dist, self.membrane_distance+extra_mem_dist+1)]
    def generate_dist_mat(self):
        curve_generator = Curve(self.output_filename, self.i, self.image, y0=self.y0, x0=self.x0, theta=self.theta, kappa=self.kappa)
        distance_mat = curve_generator.compute()
        return distance_mat
    def calculate_membrane_average(self):
        self.distance_mat = self.generate_dist_mat()
        for k in self.membrane_dist_lst:
            mask = (k - 1 < self.distance_mat) & (self.distance_mat <= k)
            gray_values = self.image[mask]
            gray_value_average = cp.average(gray_values)
            self.membrane_average_1d.append(gray_value_average)
        self.membrane_average_1d = cp.asnumpy(cp.array(self.membrane_average_1d))
        return self.membrane_average_1d

    def generate_2d_average_mem(self):
        # print('>>>Start generating 2d average membrane...')
        distance_mat_gpu = cp.array(self.distance_mat)
        membrane_average_2d_gpu = cp.zeros_like(distance_mat_gpu)
        self.membrane_dist_lst = cp.asnumpy(cp.array(self.membrane_dist_lst))
        min_mem_dist = min(self.membrane_dist_lst)
        max_mem_dist = max(self.membrane_dist_lst)
        f = interp1d(self.membrane_dist_lst, self.membrane_average_1d, kind='linear')
        in_range_mask = (distance_mat_gpu >= min_mem_dist) & (distance_mat_gpu <= max_mem_dist)
        membrane_average_2d_gpu[in_range_mask] = cp.array(f(cp.asnumpy(distance_mat_gpu[in_range_mask])))
        gray_value = cp.exp(-(cp.abs(distance_mat_gpu) - self.membrane_distance)**2 / (2 * self.sigma**2))
        membrane_average_2d_gpu[cp.abs(distance_mat_gpu) > self.membrane_distance] *= gray_value[cp.abs(distance_mat_gpu) > self.membrane_distance]
        membrane_average_2d_gpu[(gray_value < 0.001) & (cp.abs(distance_mat_gpu) > self.membrane_distance)] = 0
        self.membrane_average_2d = cp.asnumpy(membrane_average_2d_gpu)
        return membrane_average_2d_gpu
    
    def visualize_membrane_average(self):
        fig, ax = plt.subplots(1,3,figsize=(10,5))
        self.membrane_average_2d = cp.asnumpy(self.membrane_average_2d)
        self.image = cp.asnumpy(self.image)            
        ax[0].imshow(self.membrane_average_2d, cmap='gray', origin='lower')
        ax[1].imshow(self.image, cmap='gray', origin='lower')
        ax[2].imshow(self.image - self.membrane_average_2d, cmap='gray', origin='lower')
        plt.show()
    def kappa_templates_generator(self, kappa_start=-0.01, kappa_end=0.01, kappa_num=20):
        self.membrane_dist_lst = cp.asnumpy(cp.array(self.membrane_dist_lst))
        f = interp1d(self.membrane_dist_lst, self.membrane_average_1d, kind='linear')
        # num = int((kappa_end - kappa_start) / kappa_step)
        templates_gpu = cp.zeros((kappa_num, self.image_size, self.image_size))
        min_dist = min(self.membrane_dist_lst)
        max_dist = max(self.membrane_dist_lst)

        for index, kappa in enumerate(cp.linspace(kappa_start, kappa_end, kappa_num)):
            curve_generator = Curve(self.output_filename, self.i, self.image, y0=self.image_size/2, x0=self.image_size/2, theta=0, kappa=kappa)
            distance_mat_temp_gpu = curve_generator.compute()
            mask_below_min = distance_mat_temp_gpu < min_dist
            mask_above_max = distance_mat_temp_gpu > max_dist
            mask_within_range = ~mask_below_min & ~mask_above_max

            distance_mat_temp_cpu = distance_mat_temp_gpu[mask_within_range].get()
            if kappa > 0:
                distance_mat_temp_cpu = -distance_mat_temp_cpu
            interpolated_values = f(distance_mat_temp_cpu)
            membrane_average_2d_gpu = cp.zeros((self.image_size, self.image_size))
            membrane_average_2d_gpu[mask_within_range] = cp.array(interpolated_values)
            # membrane_average_2d_gpu[mask_within_range] = f(distance_mat_temp_gpu[mask_within_range])
            gray_value = cp.exp(-(cp.abs(distance_mat_temp_gpu) - self.membrane_distance)**2 / (2 * self.sigma**2))
            # gray_value = cp.exp(-(distance_mat_temp_gpu - self.membrane_distance)**2 / (2 * self.sigma**2))
            update_mask = (cp.abs(distance_mat_temp_gpu) > self.membrane_distance) & (gray_value >= 0.001)

            membrane_average_2d_gpu[update_mask] *= gray_value[update_mask]
            membrane_average_2d_gpu[(gray_value < 0.001) & (cp.abs(distance_mat_temp_gpu) > self.membrane_distance)] = 0
            templates_gpu[index] = membrane_average_2d_gpu
        savemrc(templates_gpu.get(), 'kappa_templates.mrc')


if __name__ == '__main__':
    image = readmrc('neuron_templates.mrc', section=3, mode='gpu')
    image = zoom(image, 256/64)
    df_star = readstar()
    x0 = df_star.loc[0, 'rlnCenterY']
    y0 = df_star.loc[0, 'rlnCenterX']
    theta = df_star.loc[0, 'rlnAngleTheta'] * np.pi / 180
    membrane_distance = df_star.loc[0, 'rlnMembraneDistance']
    kappa = df_star.loc[0, 'rlnCurveKappa']
    analyzer = average_membrane(image, extra_mem_dist=15, sigma=5, x0=x0, y0=y0, theta=theta, membrane_distance=membrane_distance, kappa=kappa)
    analyzer.calculate_membrane_average()
    averaged_2d_membrane = analyzer.generate_2d_average_mem()
    averaged_2d_membrane = cp.asnumpy(averaged_2d_membrane)
    savemrc(averaged_2d_membrane, 'neuron_templates_averaged_3.mrc')
    # analyzer.visualize_membrane_average()
    
    # analyzer.kappa_templates_generator(kappa_start=-0.007, kappa_end=0.007, kappa_step=0.001)
