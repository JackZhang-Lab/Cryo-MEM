import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from ._utils import *

class circle_mask:
    def __init__(self, image_array, x0, y0, radius, sigma):
        self.image = image_array
        self.image_size = self.image.shape[0]
        self.x0 = x0
        self.y0 = y0
        self.radius = radius
        self.sigma = sigma
        self.distance_matrix = cp.zeros((self.image_size, self.image_size))
        self.mask = cp.zeros((self.image_size, self.image_size))
    # def calculate_distance(self):
    #     for i in range(self.image_size):
    #         for j in range(self.image_size):
    #             distance_temp = distance(self.x0, self.y0, i, j) - self.radius
    #             self.distance_matrix[i, j] = distance_temp
    #     return self.distance_matrix
    # def generate_mask(self):
    #     for i in range(self.image_size):
    #         for j in range(self.image_size):
    #             if self.distance_matrix[i, j] <= 0:
    #                 self.mask[i, j] = 1
    #             else:
    #                 self.mask[i, j] = cp.exp(-self.distance_matrix[i, j]**2 / (2 * self.sigma**2))
    #     return self.mask
    def calculate_distance(self):
        Y, X = cp.meshgrid(cp.arange(self.image_size), cp.arange(self.image_size))
        self.distance_matrix = cp.sqrt((X - self.x0)**2 + (Y - self.y0)**2) - self.radius
        return self.distance_matrix
    def generate_mask(self):
        mask_positive = self.distance_matrix <= 0
        self.mask = cp.where(mask_positive, 1, cp.exp(-self.distance_matrix**2 / (2 * self.sigma**2)))
        return self.mask

    def apply_mask(self):
        self.distance_matrix = self.calculate_distance()
        self.mask = self.generate_mask()
        self.masked_image = self.image * self.mask
        return self.masked_image
    def visualize_mask(self):
        plt.imshow(self.masked_image.get(), cmap='gray')
        # plt.imshow(self.mask.get(), cmap='gray')
        plt.show()

if __name__ == '__main__':
    image_array = readmrc('neuron_templates.mrc', section=1, mode='gpu')
    circle_masker = circle_mask(image_array, x0=image_array.shape[0]/2, y0=image_array.shape[0]/2, radius=10, sigma=5)
    circle_masker.apply_mask()
    circle_masker.visualize_mask()
