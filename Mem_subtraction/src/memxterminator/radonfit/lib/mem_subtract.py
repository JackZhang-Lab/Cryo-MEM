import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from cupyx.scipy.ndimage import rotate, shift, zoom
from cupyx.scipy.signal import correlate2d, convolve
from ._utils import *
from .mem_average import *
from .circle_mask_generator import *
import time


class Get2Raw:
    def __init__(self, mem_analysis_starfile, average_raw, average_2d, rawimage, membrane_mask, x0, y0, theta, memdist, kappa, psi, dx, dy): # read in .star file
        self.average_raw = average_raw
        self.average_2d = average_2d
        self.rawimage = rawimage
        self.image_size_raw = self.rawimage.shape[0]
        self.membrane_mask = membrane_mask
        self.x0 = x0 # centerX of average_raw
        self.y0 = y0 # centerY of average_raw
        self.theta = theta
        self.memdist = memdist
        self.kappa = kappa
        self.psi = psi
        self.dx = dx 
        self.dy = dy 
        self.average_2d_s = cp.zeros_like(self.rawimage)
        self.average_raw_s = cp.zeros_like(self.rawimage)
        self.mem_subtracted = cp.zeros_like(self.rawimage)
        self.mem_analysis_starfile = mem_analysis_starfile
    
    def rotate_average_to_raw(self):
        average_2d_n = rotate(self.average_2d, self.psi, reshape=False)
        self.average_2d_s = shift(average_2d_n, [-self.dy, -self.dx])
        average_raw_n = rotate(self.average_raw, self.psi, reshape=False)
        self.average_raw_s = shift(average_raw_n, [-self.dy, -self.dx])
        self.membrane_mask = rotate(self.membrane_mask, self.psi, reshape=False)
        self.membrane_mask = shift(self.membrane_mask, [-self.dy, -self.dx])
        return self.average_2d_s, self.average_raw_s, self.membrane_mask
    
    def get_new_delta(self):
        xc, yc = self.image_size_raw // 2, self.image_size_raw // 2
        psi_rad = - np.radians(self.psi)
        x1, y1 = self.x0 - xc, self.y0 - yc
        x2 = x1 * np.cos(psi_rad) - y1 * np.sin(psi_rad)
        y2 = x1 * np.sin(psi_rad) + y1 * np.cos(psi_rad)
        x_prime, y_prime = x2 + xc, y2 + yc
        x_final, y_final = x_prime - self.dx, y_prime - self.dy # 旋转和平移后的坐标
        return x_final-self.x0, y_final-self.y0 # 得到新dx, dy
    

    def fit_raw_mem(self, dS, step, template_radius, mode, lp_cutoff_frequency=15, kernel_sz=5, blur_sigma=1):
        dx_n, dy_n = self.get_new_delta()
        x0 = self.x0 + dx_n
        y0 = self.y0 + dy_n
        circle_masker = circle_mask(self.average_2d_s, x0=y0, y0=x0, radius=template_radius, sigma=5)
        masked_average_2d_s = circle_masker.apply_mask()
        # theta_n = (self.psi + df_star.loc[0, 'rlnAngleTheta']) * np.pi / 180
        theta_n = (self.psi + self.theta) * np.pi / 180
        # print('theta_n:', theta_n)
        if np.isclose(theta_n, np.pi/2) or np.isclose(theta_n, 3*np.pi/2):
            dx_lst = [0] * (2*dS//step+1)
            dy_lst = [i for i in np.linspace(-dS, dS+1, step)]
            print('is close to pi/2 or 3pi/2')
        else:
            dx_lst = [i for i in np.linspace(dS*np.cos(theta_n), -dS*np.cos(theta_n)+1, step)]
            # print('dx_lst:', dx_lst)
            dy_lst = [i for i in np.linspace(-dS*np.sin(theta_n), dS*np.sin(theta_n)+1, step)]
        # x_lst = [x0 + dx for dx in dx_lst]
        # y_lst = [y0 + dy for dy in dy_lst]
        # fig, ax = plt.subplots(1, 4, figsize=(10,5))
        # ax[0].imshow(self.average_raw.get(), cmap='gray')
        # ax[1].imshow(self.average_2d.get(), cmap='gray')
        # ax[2].imshow(self.average_2d_s.get(), cmap='gray')
        # ax[2].scatter(x0, y0, color='blue', marker='x')
        # ax[2].scatter(x_lst, y_lst, color='red', marker='x')
        # ax[3].imshow(self.rawimage.get(), cmap='gray')
        # ax[3].scatter(x_lst, y_lst, color='red', marker='x')

        rawimage = (self.rawimage - cp.mean(self.rawimage)) / cp.std(self.rawimage)
        rawimage = rawimage.astype(cp.float32)
        if mode == 'lowpass':
            rawimage_n = create_gaussian_low_pass_filter(rawimage, lp_cutoff_frequency)
        elif mode == 'convolution':
            kernel = gaussian_kernel(size=kernel_sz, sigma=blur_sigma)
            rawimage_n = convolve(rawimage, kernel)
        corr_sigma_score = []
        print('>>>Recenterring membrane...')
        for k in range(len(dx_lst)):
            average_2d_s_temp = shift(masked_average_2d_s, [dy_lst[k], dx_lst[k]])
            mask_temp = shift(self.membrane_mask, [dy_lst[k], dx_lst[k]])
            average_2d_s_temp = average_2d_s_temp * mask_temp
            masked_rawimage = rawimage_n * mask_temp
            average_2d_s_temp = (average_2d_s_temp - cp.mean(average_2d_s_temp)) / cp.std(average_2d_s_temp)
            average_2d_s_temp = average_2d_s_temp.astype(cp.float32)
            corr_result = correlate2d(masked_rawimage, average_2d_s_temp, mode='same')
            corr_sigma_score.append(cp.max(corr_result))
            # print(f'iter{k} finished, dx: {dx_lst[k]}, dy: {dy_lst[k]}, corr_score: {cp.max(corr_result)}')
        corr_best_sigma_score = max(corr_sigma_score)
        corr_best_sigma_index = corr_sigma_score.index(corr_best_sigma_score)
        best_dx = dx_lst[corr_best_sigma_index]
        best_dy = dy_lst[corr_best_sigma_index]
        dx_lst = cp.asnumpy(dx_lst)
        corr_sigma_score = cp.asnumpy(cp.asarray(corr_sigma_score))
        # print('best dx:', best_dx)
        # print('best dy:', best_dy)
        # average_2d_s_n = shift(masked_average_2d_s, [best_dy, best_dx])
        # fig, ax = plt.subplots(2, 2, figsize=(10,5))
        # ax[0, 0].imshow(self.average_2d_s.get(), cmap='gray')
        # ax[0, 0].scatter(x0, y0, color='red', marker='x')
        # ax[0, 1].imshow(average_2d_s_n.get(), cmap='gray')
        # ax[0, 1].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        # ax[1,0].imshow(self.rawimage.get(), cmap='gray')
        # ax[1,0].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        # ax[1,1].plot(dx_lst, corr_sigma_score)
        # plt.show()
        return best_dx, best_dy

    def mem_corr_mask(self, image, theta):
        image_sz = image.shape[0]
        Y, X = cp.meshgrid(cp.arange(image_sz), cp.arange(image_sz))
        xc = image_sz // 2
        yc = image_sz // 2
        A = cp.cos(theta)
        B = cp.sin(theta)
        C = -A * xc - B * yc
        distance_matrix = cp.abs(A * X + B * Y + C) / cp.sqrt(A**2 + B**2)
        mask = cp.zeros_like(distance_matrix)
        mask[distance_matrix < 5] = 1
        # mask = cp.asarray(mask)
        return mask

        
    def fit_raw_mem_corr(self):
        # dx_n, dy_n = self.get_new_delta()
        # x0 = self.x0 + dx_n
        # y0 = self.y0 + dy_n
        theta_n = (self.psi + self.theta) * np.pi / 180

        rawimage = (self.rawimage - cp.mean(self.rawimage)) / cp.std(self.rawimage)
        rawimage = rawimage.astype(cp.float32)
        average_2d_s_temp = self.average_2d_s * self.membrane_mask
        average_2d_s_temp = (average_2d_s_temp - cp.mean(average_2d_s_temp)) / cp.std(average_2d_s_temp)
        average_2d_s_temp = average_2d_s_temp.astype(cp.float32)

        if rawimage.ndim != 2 or average_2d_s_temp.ndim != 2:
            raise ValueError(f"Expected both inputs to be 2D arrays, but got shapes {rawimage.shape} and {average_2d_s_temp.shape}")
        corr_result = correlate2d(rawimage, average_2d_s_temp, mode='same')
        mem_correlation_mask = self.mem_corr_mask(corr_result, theta_n)
        corr_result = corr_result * mem_correlation_mask

        corr_best_sigma_index = cp.unravel_index(cp.argmax(corr_result, axis=None), corr_result.shape)
        best_dx = corr_best_sigma_index[1] - self.image_size_raw // 2
        best_dy = corr_best_sigma_index[0] - self.image_size_raw // 2
        return best_dx, best_dy

        # print('best dx:', best_dx)
        # print('best dy:', best_dy)
        # average_2d_s_n = shift(masked_average_2d_s, [best_dy, best_dx])
        # fig, ax = plt.subplots(2, 2, figsize=(10,5))
        # ax[0, 0].imshow(self.average_2d_s.get(), cmap='gray')
        # ax[0, 0].scatter(x0, y0, color='red', marker='x')
        # ax[0, 1].imshow(average_2d_s_n.get(), cmap='gray')
        # ax[0, 1].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        # ax[1,0].imshow(self.rawimage.get(), cmap='gray')
        # ax[1,0].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        # ax[1,1].plot(dx_lst, corr_sigma_score)
        # plt.show()
        # return best_dx, best_dy

    def raw_membrane_average_subtract(self, bias, extra_mem_dist=20, scaling_factor_start=0.1, scaling_factor_end=1, scaling_factor_step=0.01):
        
        # dx_fit, dy_fit = self.fit_raw_mem(dS=20, step=10, template_radius=70, mode='convolution', kernel_sz=5, blur_sigma=1)
        dx_fit, dy_fit = self.fit_raw_mem_corr()
        dx_fit = np.float64(dx_fit)
        dy_fit = np.float64(dy_fit)
        # print(type(dx_fit), type(dy_fit))
        
        dx_n, dy_n = self.get_new_delta()
        theta_n = (self.psi + self.theta) * np.pi / 180
        membrane_distance = self.memdist
        kappa = self.kappa
        membrane_mask = shift(self.membrane_mask, [dy_fit, dx_fit])
        # start_time = time.time()
        raw_membrane_average = average_membrane(self.mem_analysis_starfile, 0, self.rawimage, extra_mem_dist=extra_mem_dist, sigma=5, x0=dy_n+self.y0+dy_fit, y0=dx_n+self.x0+dx_fit, theta=theta_n, membrane_distance=membrane_distance, kappa=kappa)
        raw_membrane_average.calculate_membrane_average()
        # end_time = time.time()
        # print(f"took {end_time - start_time:.4f} seconds.")
        averaged_raw_image = raw_membrane_average.generate_2d_average_mem()
        kernel = gaussian_kernel(5, 1)
        rawimage_n = convolve(self.rawimage, kernel)
        averaged_raw_image = convolve(averaged_raw_image, kernel)

        scaling_factor_lst = [i for i in np.arange(scaling_factor_start, scaling_factor_end, step=scaling_factor_step)]
        l1_lst = []
        for i in scaling_factor_lst:
            self.mem_subtracted = rawimage_n - averaged_raw_image * i
            l1 = cp.linalg.norm((self.mem_subtracted * membrane_mask), ord=1)
            l1_lst.append(l1)
        l1_lst = cp.asnumpy(cp.asarray(l1_lst))
        self.mem_subtracted = self.rawimage - averaged_raw_image * (scaling_factor_lst[np.argmin(l1_lst)] + bias)
        # print(f'scaling factor: {scaling_factor_lst[np.argmin(l1_lst)]}')
        # fig, ax = plt.subplots(2, 2, figsize=(10,5))
        # ax[0, 0].imshow(self.rawimage.get(), cmap='gray')
        # ax[0, 0].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        # ax[0, 1].imshow(membrane_mask.get(), cmap='gray')
        # ax[1, 0].imshow(averaged_raw_image.get(), cmap='gray')
        # ax[1, 0].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        # ax[1, 1].imshow(self.mem_subtracted.get(), cmap='gray')
        # plt.show()

        return self.mem_subtracted


if '__main__' == __name__:
    df_star = readstar('mem_analysis.star')
    average_raw = readmrc('J418/J418_100_class_averages.mrc', section=56, mode='gpu')
    average_raw = zoom(average_raw, 256/128)
    average_2d = readmrc('J418/J418_100_class_averages_averaged.mrc', section=10, mode='gpu')
    rawimage = readmrc('J417/extract/2023-02-06_17.53.10_VSV_slot3_grid3_2250X_sq01_77-7_0001_X-1Y-1-1_patch_aligned_doseweighted_particles.mrc', section=0, mode='gpu')
    membrane_mask = readmrc('J418/J418_100_class_averages_masks.mrc', section=10, mode='gpu')
    x0 = df_star.loc[0, 'rlnCenterX']
    y0 = df_star.loc[0, 'rlnCenterY']
    theta = df_star.loc[0, 'rlnAngleTheta']
    membrane_distance = df_star.loc[0, 'rlnMembraneDistance']
    kappa = df_star.loc[0, 'rlnCurveKappa']

    psi = 136.836731
    dx, dy = 6.175000, 10.075000 # x, y
    get_to_raw = Get2Raw('mem_analysis.star', average_raw, average_2d, rawimage, membrane_mask, x0, y0, theta, membrane_distance, kappa, psi, dx, dy)
    average_2d_s, average_raw_s, mask = get_to_raw.rotate_average_to_raw()
    # get_to_raw.mem_subtract()
    get_to_raw.raw_membrane_average_subtract(0.05)
    # get_to_raw.visualize()
    # get_to_raw.fit_raw_mem(50, 20)
