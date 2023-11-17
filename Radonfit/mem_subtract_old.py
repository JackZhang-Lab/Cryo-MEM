import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from cupyx.scipy.ndimage import rotate, shift, zoom
from cupyx.scipy.signal import correlate2d, convolve
from _utils import *
from mem_average import *
from circle_mask_generator import *


class get2raw:
    def __init__(self, average_raw, average_2d, rawimage, membrane_mask, x0, y0, psi, dx, dy): # read in .star file
        self.image_size_2daverage = average_raw.shape[0]
        self.image_size_raw = rawimage.shape[0]
        self.average_raw = average_raw
        self.average_2d = average_2d
        self.rawimage = rawimage
        self.x0 = x0
        self.y0 = y0
        self.psi = psi
        self.dx = dx 
        self.dy = dy 
        self.average_2d_s = cp.zeros((self.image_size_raw, self.image_size_raw))
        self.average_raw_s = cp.zeros((self.image_size_raw, self.image_size_raw))
        self.membrane_mask = membrane_mask
        self.mem_subtracted = cp.zeros((self.image_size_raw, self.image_size_raw))
    
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
        x_final, y_final = x_prime - self.dx, y_prime - self.dy
        # print(f"原始点的坐标: ({self.x0/self.zoomfactor}, {self.y0/self.zoomfactor})")
        # print(f"旋转和平移后的坐标: ({x_final/self.zoomfactor}, {y_final/self.zoomfactor})")
        # print(f'dx_n: {(x_final-self.x0)/self.zoomfactor}, dy_n: {(y_final-self.y0)/self.zoomfactor}')
        return x_final-self.x0, y_final-self.y0
    

    def fit_raw_mem(self, dS, step, template_radius, mode, lp_cutoff_frequency=15, kernel_sz=5, blur_sigma=1):
        df_star = readstar()
        dx_n, dy_n = self.get_new_delta()
        x0 = self.x0 + dx_n
        y0 = self.y0 + dy_n
        circle_masker = circle_mask(self.average_2d_s, x0=y0, y0=x0, radius=template_radius, sigma=5)
        masked_average_2d_s = circle_masker.apply_mask()
        theta_n = (self.psi + df_star.loc[0, 'rlnAngleTheta']) * np.pi / 180
        # print('theta_n:', theta_n)
        if np.isclose(theta_n, np.pi/2) or np.isclose(theta_n, 3*np.pi/2):
            dx_lst = [0] * (2*dS//step+1)
            dy_lst = [i for i in np.linspace(-dS, dS+1, step)]
            print('is close to pi/2 or 3pi/2')
        else:
            dx_lst = [i for i in np.linspace(dS*np.cos(theta_n), -dS*np.cos(theta_n)+1, step)]
            # print('dx_lst:', dx_lst)
            dy_lst = [i for i in np.linspace(-dS*np.sin(theta_n), dS*np.sin(theta_n)+1, step)]
        x_lst = [x0 + dx for dx in dx_lst]
        y_lst = [y0 + dy for dy in dy_lst]
        fig, ax = plt.subplots(1, 4, figsize=(10,5))
        ax[0].imshow(self.average_raw.get(), cmap='gray')
        ax[1].imshow(self.average_2d.get(), cmap='gray')
        ax[2].imshow(self.average_2d_s.get(), cmap='gray')
        ax[2].scatter(x0, y0, color='blue', marker='x')
        ax[2].scatter(x_lst, y_lst, color='red', marker='x')
        ax[3].imshow(self.rawimage.get(), cmap='gray')
        ax[3].scatter(x_lst, y_lst, color='red', marker='x')

        rawimage = (self.rawimage - cp.mean(self.rawimage)) / cp.std(self.rawimage)
        rawimage = rawimage.astype(cp.float32)
        if mode == 'lowpass':
            rawimage_n = create_gaussian_low_pass_filter(rawimage, lp_cutoff_frequency)
        elif mode == 'convolution':
            kernel = gaussian_kernel(size=kernel_sz, sigma=blur_sigma)
            rawimage_n = convolve(rawimage, kernel)
        corr_sigma_score = []
        for k in range(len(dx_lst)):
            average_2d_s_temp = shift(masked_average_2d_s, [dy_lst[k], dx_lst[k]])
            mask_temp = shift(self.membrane_mask, [dy_lst[k], dx_lst[k]])
            average_2d_s_temp = average_2d_s_temp * mask_temp
            masked_rawimage = rawimage_n * mask_temp
            average_2d_s_temp = (average_2d_s_temp - cp.mean(average_2d_s_temp)) / cp.std(average_2d_s_temp)
            average_2d_s_temp = average_2d_s_temp.astype(cp.float32)
            corr_result = correlate2d(masked_rawimage, average_2d_s_temp, mode='same')
            corr_sigma_score.append(cp.max(corr_result))
            print(f'iter{k} finished, dx: {dx_lst[k]}, dy: {dy_lst[k]}, corr_score: {cp.max(corr_result)}')
        corr_best_sigma_score = max(corr_sigma_score)
        corr_best_sigma_index = corr_sigma_score.index(corr_best_sigma_score)
        best_dx = dx_lst[corr_best_sigma_index]
        best_dy = dy_lst[corr_best_sigma_index]
        dx_lst = cp.asnumpy(dx_lst)
        corr_sigma_score = cp.asnumpy(cp.asarray(corr_sigma_score))
        print('best dx:', best_dx)
        print('best dy:', best_dy)
        average_2d_s_n = shift(masked_average_2d_s, [best_dy, best_dx])
        fig, ax = plt.subplots(2, 2, figsize=(10,5))
        ax[0, 0].imshow(self.average_2d_s.get(), cmap='gray')
        ax[0, 0].scatter(x0, y0, color='red', marker='x')
        ax[0, 1].imshow(average_2d_s_n.get(), cmap='gray')
        ax[0, 1].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        ax[1,0].imshow(self.rawimage.get(), cmap='gray')
        ax[1,0].scatter(x0+best_dx, y0+best_dy, color='red', marker='x')
        ax[1,1].plot(dx_lst, corr_sigma_score)
        plt.show()
        return best_dx, best_dy



    def mem_subtract(self):
        dx_fit, dy_fit = self.fit_raw_mem(dS=25, step=25, template_radius=70, mode='lowpass', lp_cutoff_frequency=20)
        average_2d_s_n = shift(self.average_2d_s, [dy_fit, dx_fit])
        membrane_mask = shift(self.membrane_mask, [dy_fit, dx_fit])
        rawimage_lp = create_gaussian_low_pass_filter(self.rawimage, cutoff_frequency=15)
        m_std = cp.std(rawimage_lp * membrane_mask)
        m_mean = cp.sum(rawimage_lp * membrane_mask) / cp.sum(membrane_mask)
        average_2d_s_n = ((average_2d_s_n - cp.mean(average_2d_s_n)) / cp.std(average_2d_s_n)) * m_std + m_mean
        scaling_factor_lst = [i for i in np.arange(0.02, 0.70, step=0.01)]
        l1_lst = []
        for i in scaling_factor_lst:
            self.mem_subtracted = rawimage_lp - average_2d_s_n * i
            l1 = cp.linalg.norm((self.mem_subtracted * membrane_mask), ord=1)
            l1_lst.append(l1)
        plt.plot(scaling_factor_lst, l1_lst)
        plt.show()
        self.mem_subtracted = self.rawimage - average_2d_s_n * scaling_factor_lst[cp.argmin(l1_lst)]

        dx_n, dy_n = self.get_new_delta()
        fig, ax = plt.subplots(2, 4, figsize=(10,5))
        ax[0,0].imshow(self.average_raw.get(), cmap='gray')
        ax[0,0].scatter(self.x0, self.y0, color='red', marker='x')
        ax[0,1].imshow(self.average_2d.get(), cmap='gray')
        ax[0,2].imshow(self.average_2d_s.get(), cmap='gray')
        ax[0,2].scatter(dx_n+self.x0, dy_n+self.y0, color='red', marker='x')
        ax[0,3].imshow(self.average_raw_s.get(), cmap='gray')
        ax[0,3].scatter(dx_n+self.x0, dy_n+self.y0, color='red', marker='x')
        ax[1,0].imshow(average_2d_s_n.get(), cmap='gray')
        ax[1,0].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        ax[1,1].imshow(self.rawimage.get(), cmap='gray')
        ax[1,1].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        ax[1,2].imshow(self.mem_subtracted.get(), cmap='gray')
        ax[1,3].imshow(membrane_mask.get(), cmap='gray')
        # fig, ax = plt.subplots(1, 3, figsize=(10,5))
        # ax[0].imshow(self.average_raw, cmap='gray')
        # ax[1].imshow(self.average_raw_s, cmap='gray')
        # ax[2]. imshow(self.rawimage, cmap='gray')
        plt.show()

    def raw_membrane_average_subtract(self, theta_n, membrane_distance, kappa):
        df_star = readstar()
        dx_fit, dy_fit = self.fit_raw_mem(dS=40, step=25, template_radius=70, mode='convolution', kernel_sz=5, blur_sigma=1)
        dx_n, dy_n = self.get_new_delta()
        # 重新做输入，不然会冲突
        theta_n = (self.psi + df_star.loc[0, 'rlnAngleTheta']) * np.pi / 180
        membrane_distance = df_star.loc[0, 'rlnMembraneDistance']
        kappa = df_star.loc[0, 'rlnCurveKappa']
        membrane_mask = shift(self.membrane_mask, [dy_fit, dx_fit])
        raw_membrane_average = average_membrane(0, self.rawimage, extra_mem_dist=15, sigma=5, x0=dy_n+self.y0+dy_fit, y0=dx_n+self.x0+dx_fit, theta=theta_n, membrane_distance=membrane_distance, kappa=kappa)
        raw_membrane_average.calculate_membrane_average()
        averaged_raw_image = raw_membrane_average.generate_2d_average_mem()
        kernel = gaussian_kernel(5, 1)
        rawimage_n = convolve(self.rawimage, kernel)
        averaged_raw_image = convolve(averaged_raw_image, kernel)

        scaling_factor_lst = [i for i in np.arange(0.02, 0.70, step=0.01)]
        l1_lst = []
        for i in scaling_factor_lst:
            self.mem_subtracted = rawimage_n - averaged_raw_image * i
            l1 = cp.linalg.norm((self.mem_subtracted * membrane_mask), ord=1)
            l1_lst.append(l1)
        l1_lst = cp.asnumpy(cp.asarray(l1_lst))
        self.mem_subtracted = self.rawimage - averaged_raw_image * scaling_factor_lst[np.argmin(l1_lst)]
        print(f'scaling factor: {scaling_factor_lst[np.argmin(l1_lst)]}')
        fig, ax = plt.subplots(2, 2, figsize=(10,5))
        ax[0, 0].imshow(self.rawimage.get(), cmap='gray')
        ax[0, 0].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        ax[0, 1].imshow(membrane_mask.get(), cmap='gray')
        ax[1, 0].imshow(averaged_raw_image.get(), cmap='gray')
        ax[1, 0].scatter(dx_n+self.x0+dx_fit, dy_n+self.y0+dy_fit, color='red', marker='x')
        ax[1, 1].imshow(self.mem_subtracted.get(), cmap='gray')
        plt.show()
        return self.mem_subtracted


if '__main__' == __name__:
    df_star = readstar()
    average_raw = readmrc('neuron_templates.mrc', section=4, mode='gpu')
    average_raw = zoom(average_raw, 256/64)
    average_2d = readmrc('neuron_templates_averaged_4.mrc', section=0, mode='gpu')
    rawimage = readmrc('J379/extract/2023-05-28_00.21.16_neuron_152-3_00013_X-1Y-1-1_patch_aligned_doseweighted_particles.mrc', section=0, mode='gpu')
    membrane_mask = readmrc('neuron_templates_mask_4.mrc', section=0, mode='gpu')
    x0 = df_star.loc[0, 'rlnCenterX']
    y0 = df_star.loc[0, 'rlnCenterY']
    psi = 272.755096
    dx, dy = -1.650000, 0.450000 # x, y
    get_to_raw = get2raw(average_raw, average_2d, rawimage, membrane_mask, x0, y0, psi, dx, dy)
    average_2d_s, average_raw_s, mask = get_to_raw.rotate_average_to_raw()
    # get_to_raw.mem_subtract()
    get_to_raw.raw_membrane_average_subtract()
    # get_to_raw.visualize()
    # get_to_raw.fit_raw_mem(50, 20)