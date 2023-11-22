import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from scipy.stats import norm
from ._utils import *

def gaussian_pdf(x, mean, std):
    coefficient = 1.0 / (std * cp.sqrt(2 * cp.pi))
    exponential = cp.exp(- (x - mean) ** 2 / (2 * std ** 2))
    return coefficient * exponential

def gaussian2(x, membrane_dist, std1, std2):
    mean1 = - membrane_dist / 2
    mean2 = membrane_dist / 2
    gaussian1 = norm.pdf(x, mean1, std1)
    gaussian2 = norm.pdf(x, mean2, std2)
    g = gaussian1 + gaussian2
    g = g / np.sum(g)
    g = g * 255
    return g

def gaussian2_gpu(x, membrane_dist, std1, std2):
    mean1 = - membrane_dist / 2
    mean2 = membrane_dist / 2
    gaussian1 = gaussian_pdf(x, mean1, std1)
    gaussian2 = gaussian_pdf(x, mean2, std2)
    g = gaussian1 + gaussian2
    g = g / cp.sum(g)
    g = g * 255
    return g

def generate_template(size, theta, membrane_dist, std1, std2, mode='cpu'):
    if mode == 'cpu':
        distance_matrix = np.zeros((size, size))
        theta = (theta + 90) * np.pi / 180
        x = range(size)
        y = range(size)
        x0 = size / 2
        y0 = size / 2
        if np.isclose(theta, np.pi/2) or np.isclose(theta, 3*np.pi/2):
            for i in x:
                for j in y:
                    distance_matrix[i, j] = i - x0
        else:
            k = np.tan(theta)
            b = y0 - k * x0
            for i in x:
                for j in y:
                    distance_matrix[i, j] = (k * i - j + b) / np.sqrt(k**2 + 1)
        template = gaussian2(distance_matrix, membrane_dist, std1, std2)
        return template
    elif mode == 'gpu':
        distance_matrix = cp.zeros((size, size))
        theta = (theta + 90) * cp.pi / 180
        x = range(size)
        y = range(size)
        x0 = size / 2
        y0 = size / 2
        if cp.isclose(theta, cp.pi/2) or np.isclose(theta, 3*cp.pi/2):
            for i in x:
                for j in y:
                    distance_matrix[i, j] = i - x0
        else:
            k = cp.tan(theta)
            b = y0 - k * x0
            for i in x:
                for j in y:
                    distance_matrix[i, j] = (k * i - j + b) / cp.sqrt(k**2 + 1)
        template = gaussian2_gpu(distance_matrix, membrane_dist, std1, std2)
        return template

if __name__ == '__main__':
    x = np.linspace(-100, 100, 200)
    g = gaussian2(x,12,4,4)
    alpha = 30
    size = 32
    template_rotated = generate_template(size, alpha, 12, 4, 2)
    plt.imshow(template_rotated, cmap='gray', origin='lower')
    plt.show()

