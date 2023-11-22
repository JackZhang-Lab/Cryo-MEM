import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from ._utils import *
from .template_centerfitting import *
from .radonanalyser import *
from .circle_mask_generator import *
from cupyx.scipy.ndimage import zoom
from cupyx.scipy.signal import correlate2d

class Curve:
    def __init__(self, output_filename, i, image, y0, x0, theta, kappa):
        # self.radonanalyze = RadonAnalyzer(image, thr=thr)
        # self.centerfit = Template_centerfitting(3, 3, image, thr=thr)
        df_star = starfile.read(output_filename)
        # self.y0, self.x0 = df_star.loc[0, 'rlnCenterX'], df_star.loc[0, 'rlnCenterY'] # center of membrane
        # if mode == 0:
        #     self.theta = df_star.loc[0, 'rlnAngleTheta'] * np.pi / 180
        # elif mode == 1:
        #     self.theta = theta
        #     # self.y0, self.x0 = self.centerfit.centerfinder() # center of membrane
        #     self.y0, self.x0 = y0, x0
        self.y0 = y0
        self.x0 = x0
        self.theta = theta
        self.gray_image = image
        self.image_size = image.shape[0]
        self.sigma1, self.sigma2 = df_star.loc[i, 'rlnSigma1'], df_star.loc[i, 'rlnSigma2']
        self.kappa = kappa
        self.x_c = 0 # center of the circle
        self.y_c = 0
        self.R = 0
        self.membrane_distance = df_star.loc[i, 'rlnMembraneDistance']
        self.distance_matrix = cp.zeros((self.image_size, self.image_size))
        self.simulate_membrane = cp.zeros((self.image_size, self.image_size))
    def compute(self):
        if self.kappa == 0:
            return self.compute_line_dist()
        else:
            return self.compute_arc_dist()
    def compute_arc_center(self):
        self.R = 1 / np.abs(self.kappa) if self.kappa != 0 else np.inf
        if self.kappa >= 0:
            self.x_c = self.x0 - self.R * np.sin(self.theta)
            self.y_c = self.y0 + self.R * np.cos(self.theta)
        else:
            self.x_c = self.x0 + self.R * np.sin(self.theta)
            self.y_c = self.y0 - self.R * np.cos(self.theta)
        return self.x_c, self.y_c
    def if_in_arc(self, x_temp, y_temp):
        x_c, y_c = self.compute_arc_center()
        d = distance(x_c, y_c, x_temp, y_temp)
        if d <= self.R:
            return True
        else:
            return False
    def compute_arc_dist(self):
        x_c, y_c = self.compute_arc_center()
        x = range(self.image_size)
        y = range(self.image_size)
        for i in x:
            for j in y:
                distance_temp = distance(x_c, y_c, i, j) - self.R
                self.distance_matrix[i, j] = distance_temp
        return self.distance_matrix
    def compute_line_dist(self):
        x = range(self.image_size)
        y = range(self.image_size)
        if np.isclose(self.theta, np.pi/2) or np.isclose(self.theta, 3*np.pi/2):
            for i in x:
                for j in y:
                    self.distance_matrix[i, j] = i - self.x0
            return self.distance_matrix
        else:
            k = np.tan(self.theta)
            b = self.y0 - k * self.x0
            for i in x:
                for j in y:
                    self.distance_matrix[i, j] = (k * i - j + b) / np.sqrt(k**2 + 1)
            return self.distance_matrix
    def generate_membrane(self):
        self.distance_matrix = self.compute()
        self.simulate_membrane = gaussian2_gpu(self.distance_matrix, self.membrane_distance, self.sigma1, self.sigma2)
        return self.simulate_membrane

class Curvefitting:
    def __init__(self, output_filename, i, image, kappa_start, kappa_end, kappa_step):
        # self.thr = thr
        df_star = readstar(output_filename)
        self.output_filename = output_filename
        self.i = i
        # radonanalyze = RadonAnalyzer(image, thr=thr)
        # centerfit = Template_centerfitting(3, 3, image, thr=thr)
        num = int((kappa_end - kappa_start) / kappa_step) + 1
        self.kappa_lst = np.linspace(kappa_start, kappa_end, num)
        self.corr_lst = []
        self.yc, self.xc = df_star.loc[i, 'rlnCenterX'], df_star.loc[i, 'rlnCenterY']
        self.gray_image = image
        self.image_size = image.shape[0]
        self.membrane_distance = df_star.loc[i, 'rlnMembraneDistance']
        self.theta = df_star.loc[i, 'rlnAngleTheta'] * np.pi / 180
        self.best_kappa = self.fit_curve()
        df_star.loc[i, 'rlnCurveKappa'] = self.best_kappa
        write_star(df_star, output_filename)
    def generate_membrane(self, kappa):
        gene_curve = Curve(self.output_filename, self.i, self.gray_image, self.yc, self.xc, self.theta, kappa)
        self.simulated_membrane = gene_curve.generate_membrane()
        return self.simulated_membrane
    def fit_curve(self):
        i = 0
        print('>>>Start curve fitting...')
        for kappa in self.kappa_lst:
            # print(f'=====iter: {i}=====')
            self.simulated_membrane  = self.generate_membrane(kappa)
            self.simulated_membrane = self.simulated_membrane.astype(cp.float64)
            self.gray_image = self.gray_image.astype(cp.float64)
            self.simulated_membrane = (self.simulated_membrane - cp.mean(self.simulated_membrane)) / cp.std(self.simulated_membrane)
            radius = (self.image_size/2) * 0.9
            circle_masker = circle_mask(self.simulated_membrane, x0=self.image_size/2, y0=self.image_size/2, radius=radius, sigma=3)
            self.simulated_membrane = circle_masker.apply_mask()
            self.gray_image = (self.gray_image - cp.mean(self.gray_image)) / cp.std(self.gray_image)
            correlation_result = correlate2d(self.simulated_membrane, self.gray_image, mode='full')
            corr_score = cp.max(correlation_result)
            self.corr_lst.append(corr_score)
            print(f'iter{i} finished, kappa: {kappa}, corr_score: {corr_score}')
            i += 1
        self.corr_max = max(self.corr_lst)
        self.corr_max_index = self.corr_lst.index(self.corr_max)
        self.best_kappa = self.kappa_lst[self.corr_max_index]
        print('Best kappa:', self.best_kappa)
        self.mem_best = self.generate_membrane(self.best_kappa)
        return self.best_kappa
    def fit_curve_visualize(self):
        fig, axes = plt.subplots(1, 3, figsize=(10, 5))
        axes[0].plot(self.kappa_lst, self.corr_lst, 'x', color = 'red')
        axes[2].imshow(self.mem_best, cmap='gray', origin='lower')
        axes[1].imshow(self.gray_image, cmap='gray', origin='lower')
        plt.show()

if __name__ == '__main__':
    # gene_curve = Curve('image2.mrc', section=1, thr=0.5, kappa=0)
    # membrane = gene_curve.generate_membrane()
    # plt.imshow(membrane, cmap='gray', origin='lower')
    # plt.show()
    image = readmrc('neuron_templates.mrc', section=3, mode='gpu')
    image = zoom(image, 256/64)
    curve_fitting = Curvefitting(image, -0.01, 0.01, 21)
    # curve_fitting.fit_curve()
    # curve_fitting.fit_curve_visualize()
    # mem_sub = curve_fitting.mem_subtract()
    # plt.imshow(mem_sub, cmap='gray', origin='lower')
    # plt.show()