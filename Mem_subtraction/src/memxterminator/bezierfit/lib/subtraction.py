import numpy as np
import cupy as cp
import cupyx.scipy.ndimage
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import mrcfile
from cupyx.scipy.signal import convolve2d
from .bezierfit import points_along_normal, generate_curve_within_boundaries, cupy_cdist
import time


def bilinear_interpolation(image, x, y):
    # Define the coordinates of the image
    x_coords = np.arange(image.shape[1])
    y_coords = np.arange(image.shape[0])
    x_mesh, y_mesh = np.meshgrid(x_coords, y_coords)   
    # Flatten the image and meshgrids
    points = np.array([x_mesh.flatten(), y_mesh.flatten()]).T
    values = image.flatten()
    # Use griddata for interpolation
    interpolated_values = griddata(points, values, (x, y), method='linear')
    return interpolated_values

def bilinear_interpolation_gpu(image_gpu, x, y):
    x_gpu = cp.asarray(x)
    y_gpu = cp.asarray(y)
    # Perform bilinear interpolation
    coordinates_gpu = cp.stack((y_gpu, x_gpu), axis=0)
    interpolated_values = cupyx.scipy.ndimage.map_coordinates(image_gpu, coordinates_gpu, order=1)
    return interpolated_values

def gaussian_kernel(size, sigma=1.0):
    size = int(size) // 2
    x, y = np.mgrid[-size:size+1, -size:size+1]
    normal = 1 / (2.0 * np.pi * sigma**2)
    g =  np.exp(-((x**2 + y**2) / (2.0*sigma**2))) * normal
    return cp.asarray(g)

def rotate_by_matrix(x, y, theta, center):
    # theta = -np.radians(theta_degrees)
    theta = -theta
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    x_shifted = x - center[0]
    y_shifted = y - center[1]
    x_rotated, y_rotated = np.dot(rotation_matrix, [x_shifted, y_shifted])
    x_final = x_rotated + center[0]
    y_final = y_rotated + center[1]
    return x_final, y_final

def shift_point(x, y, shift_values):
    x_shifted = x + shift_values[0]
    y_shifted = y + shift_values[1]
    
    return x_shifted, y_shifted

def distribute_bilinearly(image, count_image, x, y, value):
    x1, y1 = int(x), int(y)
    x2, y2 = x1 + 1, y1 + 1
    if x2 >= image.shape[1] or y2 >= image.shape[0]:
        image[y1, x1] += value  # 如果超出范围，全部分给最近的像素
        count_image[y1, x1] += 1
        return
    # 计算双线性权重
    w11 = (x2 - x) * (y2 - y)
    w21 = (x - x1) * (y2 - y)
    w12 = (x2 - x) * (y - y1)
    w22 = (x - x1) * (y - y1)
    # 将 value 分配到四个邻近像素上
    image[y1, x1] += value * w11
    image[y1, x2] += value * w21
    image[y2, x1] += value * w12
    image[y2, x2] += value * w22
    count_image[y1, x1] += w11
    count_image[y1, x2] += w21
    count_image[y2, x1] += w12
    count_image[y2, x2] += w22
    return image, count_image

def distribute_bilinearly_gpu(image, count_image, x, y, value):
    distribute_bilinearly_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void distribute_bilinearly(float* image, float* count_image, int image_width, int image_height,
                            float x, float y, float value) {
        int x1 = floor(x), y1 = floor(y);
        int x2 = x1 + 1, y2 = y1 + 1;
        if (x1 < 0 || y1 < 0 || x2 >= image_width || y2 >= image_height) {
            atomicAdd(&image[y1 * image_width + x1], value);
            atomicAdd(&count_image[y1 * image_width + x1], 1);
            return;
        }
        float w11 = (x2 - x) * (y2 - y);
        float w21 = (x - x1) * (y2 - y);
        float w12 = (x2 - x) * (y - y1);
        float w22 = (x - x1) * (y - y1);
        atomicAdd(&image[y1 * image_width + x1], value * w11);
        atomicAdd(&image[y1 * image_width + x2], value * w21);
        atomicAdd(&image[y2 * image_width + x1], value * w12);
        atomicAdd(&image[y2 * image_width + x2], value * w22);
        atomicAdd(&count_image[y1 * image_width + x1], w11);
        atomicAdd(&count_image[y1 * image_width + x2], w21);
        atomicAdd(&count_image[y2 * image_width + x1], w12);
        atomicAdd(&count_image[y2 * image_width + x2], w22);
    }
    ''', 'distribute_bilinearly')
    image_gpu = cp.asarray(image)
    count_image_gpu = cp.asarray(count_image)
    x = cp.asarray(x)
    y = cp.asarray(y)
    value = cp.asarray(value)
    distribute_bilinearly_kernel((1,), (1,), (image_gpu, count_image_gpu, image.shape[1], image.shape[0], x, y, value))
    return image_gpu, count_image_gpu

class MembraneSubtract:
    def __init__(self, control_points, image, psi, origin_x, origin_y, pixel_size, points_step, physical_membrane_dist=35.0):
        # with mrcfile.open(rawimage_filename) as f:
        #     self.image = f.data[section]
        self.image = image
        # self.image_gpu = cp.asarray(self.image)
        self.image_gpu = image
        self.image_sz = self.image.shape[0]
        self.xc = self.image.shape[0] // 2
        self.yc = self.image.shape[1] // 2
        self.pixel_size = pixel_size
        self.mem_dist = int(physical_membrane_dist / self.pixel_size)
        self.psi = psi
        self.origin_x = origin_x
        self.origin_y = origin_y
        self.control_points = control_points
        self.points_step = points_step

    def control_points_trasf(self, control_points, psi, origin_x, origin_y):
        for control_point in control_points:
            control_point[0], control_point[1] = rotate_by_matrix(control_point[0], control_point[1], psi, [self.xc, self.yc])
            control_point[0], control_point[1] = shift_point(control_point[0], control_point[1], [-origin_x, -origin_y])
        return control_points

    def generate_2d_mask(self, image, fitted_points, membrane_distance):
        y, x = cp.mgrid[:image.shape[0], :image.shape[1]]
        coords = cp.stack((x, y), axis=-1).reshape(-1, 2)
        fitted_points = cp.array(fitted_points)
        distances = cupy_cdist(coords, fitted_points)
        min_distances = cp.min(distances, axis=1).reshape(image.shape)
        edge_sigma = 5
        mask_within_distance = cp.abs(min_distances) <= membrane_distance
        gray_value = cp.exp(-(cp.abs(min_distances) - membrane_distance)**2 / (2 * edge_sigma**2))
        mask_small_gray_value = gray_value < 0.001
        mask_outside_distance = ~mask_within_distance
        membrane_mask = cp.zeros_like(min_distances)
        membrane_mask[mask_within_distance] = 1
        membrane_mask[mask_outside_distance & ~mask_small_gray_value] = gray_value[mask_outside_distance & ~mask_small_gray_value]
        return membrane_mask

    def average_1d(self, image_gpu, fitted_points, normals, extra_mem_dist):
        average_1d_lst = []
        for membrane_dist in range(-extra_mem_dist, extra_mem_dist+1):
            normals_points = fitted_points + membrane_dist * normals
            # Ensure the points are within the image boundaries
            mask = (normals_points[:, 0] >= 0) & (normals_points[:, 0] < image_gpu.shape[1]) & \
                (normals_points[:, 1] >= 0) & (normals_points[:, 1] < image_gpu.shape[0])
            normals_points = normals_points[mask]
            # Get interpolated gray values for the normal points
            interpolated_values = bilinear_interpolation_gpu(image_gpu, normals_points[:, 0], normals_points[:, 1])
            # convert nan to 0
            interpolated_values = np.nan_to_num(interpolated_values)
            average_1d_lst.append(np.mean(interpolated_values))
        return average_1d_lst
    # def average_1d_gpu(self, image_gpu, fitted_points, normals, extra_mem_dist):
    #     average_1d_lst = []
    #     fitted_points = cp.asarray(fitted_points)
    #     membrane_dists = cp.arange(-extra_mem_dist, extra_mem_dist + 1)
    #     # Expand dimensions for broadcasting
    #     membrane_dists = membrane_dists[:, cp.newaxis, cp.newaxis]
    #     # Calculate all normals_points at once
    #     normals_points = fitted_points + membrane_dists * normals
    #     mask = (normals_points[:, 0] >= 0) & (normals_points[:, 0] < image_gpu.shape[1]) & \
    #             (normals_points[:, 1] >= 0) & (normals_points[:, 1] < image_gpu.shape[0])
    #     mask_expanded = cp.repeat(mask[:, cp.newaxis, :], normals_points.shape[1], axis=1)
    #     normals_points = cp.compress(mask_expanded, normals_points, axis=0)
    #     print(normals_points.shape)
    #     # Get interpolated gray values for the normal points
    #     interpolated_values = bilinear_interpolation_gpu(image_gpu, normals_points[:, 0], normals_points[:, 1])
    #     # convert nan to 0
    #     interpolated_values = cp.nan_to_num(interpolated_values)
    #     average_1d_lst.append(cp.mean(interpolated_values))
    #     print('average_1d_lst shape: ', average_1d_lst)
    #     return average_1d_lst

    # def average_2d(self, image, fitted_points, normals, average_1d_lst, extra_mem_dist):
    #     new_image = cp.zeros_like(image)
    #     count_image = cp.zeros_like(image)
    #     def distribute_bilinearly(image, count_image, x, y, value):
    #         x1, y1 = int(x), int(y)
    #         x2, y2 = x1 + 1, y1 + 1
    #         if x2 >= image.shape[1] or y2 >= image.shape[0]:
    #             image[y1, x1] += value  # 如果超出范围，全部分给最近的像素
    #             count_image[y1, x1] += 1
    #             return
    #         # 计算双线性权重
    #         w11 = (x2 - x) * (y2 - y)
    #         w21 = (x - x1) * (y2 - y)
    #         w12 = (x2 - x) * (y - y1)
    #         w22 = (x - x1) * (y - y1)
    #         # 将 value 分配到四个邻近像素上
    #         image[y1, x1] += value * w11
    #         image[y1, x2] += value * w21
    #         image[y2, x1] += value * w12
    #         image[y2, x2] += value * w22
    #         count_image[y1, x1] += w11
    #         count_image[y1, x2] += w21
    #         count_image[y2, x1] += w12
    #         count_image[y2, x2] += w22
    #     for membrane_dist, average_1d in zip(range(-extra_mem_dist, extra_mem_dist+1), average_1d_lst):
    #         normals_points = fitted_points + membrane_dist * normals
    #         # Ensure the points are within the image boundaries
    #         mask = (normals_points[:, 0] >= 0) & (normals_points[:, 0] < image.shape[1]) & \
    #             (normals_points[:, 1] >= 0) & (normals_points[:, 1] < image.shape[0])
    #         normals_points = normals_points[mask]
    #         # Give the normal points the average gray value with interpolation
    #         for point in normals_points:
    #             distribute_bilinearly(new_image, count_image, point[0], point[1], average_1d)
    #         # Normalize the new_image by the count_image
    #     with cp.errstate(divide='ignore', invalid='ignore'):
    #         new_image /= count_image
    #         new_image[cp.isnan(new_image)] = 0
    #     return new_image.astype(image.dtype)

    def average_2d(self, image_gpu, fitted_points, normals, average_1d_lst, extra_mem_dist):
        image = image_gpu.get()
        new_image = np.zeros_like(image)
        count_image = np.zeros_like(image)
        for membrane_dist, average_1d in zip(range(-extra_mem_dist, extra_mem_dist+1), average_1d_lst):
            # start_time = time.time()
            normals_points = fitted_points + membrane_dist * normals
            mask = (normals_points[:, 0] >= 0) & (normals_points[:, 0] < image.shape[1]) & \
                (normals_points[:, 1] >= 0) & (normals_points[:, 1] < image.shape[0])
            normals_points = normals_points[mask]
            for point in normals_points:
                distribute_bilinearly(new_image, count_image, point[0], point[1], average_1d)
            # end_time = time.time()
            # print('Time for average_2d: ', end_time - start_time)
        # Normalize the new_image by the count_image
        with np.errstate(divide='ignore', invalid='ignore'):
            new_image /= count_image
            new_image[np.isnan(new_image)] = 0
        return new_image.astype(image.dtype)
    
    def average_2d_gpu(self, image_gpu, fitted_points, normals, average_1d_lst, extra_mem_dist):
        new_image = cp.zeros_like(image_gpu)
        count_image = cp.zeros_like(image_gpu)
        fitted_points = cp.asarray(fitted_points)
        membrane_dists = cp.arange(-extra_mem_dist, extra_mem_dist + 1)
        # Expand dimensions for broadcasting
        membrane_dists = membrane_dists[:, cp.newaxis, cp.newaxis]
        # Calculate all normals_points at once
        normals = cp.asarray(normals)
        normals_points = fitted_points + membrane_dists * normals
        mask = (normals_points[:, 0] >= 0) & (normals_points[:, 0] < image_gpu.shape[1]) & \
                (normals_points[:, 1] >= 0) & (normals_points[:, 1] < image_gpu.shape[0])
        mask_expanded = cp.repeat(mask[:, cp.newaxis, :], normals_points.shape[1], axis=1)
        normals_points = cp.compress(mask_expanded, normals_points, axis=0)
        # Give the normal points the average gray value with interpolation
        average_1d_lst = cp.asarray(average_1d_lst)
        average_1d_lst = cp.compress(mask, average_1d_lst, axis=0)
        for average_1d in average_1d_lst:
            for point in normals_points:
                distribute_bilinearly_gpu(new_image, count_image, point[0], point[1], average_1d)
        mask = count_image != 0
        new_image[mask] /= count_image[mask]
        new_image[cp.isnan(new_image)] = 0
        return new_image.astype(image_gpu.dtype)
    
    def mem_subtract(self):
        control_points = self.control_points_trasf(self.control_points, self.psi, self.origin_x, self.origin_y)
        fitted_curve_points, t_values = generate_curve_within_boundaries(control_points, self.image.shape, self.points_step)
        # plt.imshow(self.image, cmap='gray')
        # plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'r-')
        # plt.plot(control_points[:, 0], control_points[:, 1], 'g.')
        # plt.show()
        mem_mask = self.generate_2d_mask(self.image_gpu, fitted_curve_points, self.mem_dist)
        raw_image_average_1d_lst = self.average_1d(self.image_gpu, fitted_curve_points, points_along_normal(control_points, t_values).get(), self.mem_dist)
        raw_image_average_2d = self.average_2d(self.image_gpu, fitted_curve_points, points_along_normal(control_points, t_values).get(), raw_image_average_1d_lst, self.mem_dist)
        raw_image_average_2d = cp.asarray(raw_image_average_2d)
        kernel = gaussian_kernel(5, 1)
        image_conv = convolve2d(self.image_gpu, kernel, mode = 'same')
        raw_image_average_2d_conv = convolve2d(raw_image_average_2d, kernel, mode = 'same')
        bias = 0
        scaling_factor_lst = [i for i in np.arange(0.01, 1, 0.02)]
        l1_lst = []
        for i in scaling_factor_lst:
            mem_subtracted = image_conv - raw_image_average_2d_conv * i
            l1 = cp.linalg.norm((mem_subtracted * mem_mask), ord=1)
            l1_lst.append(l1)
        l1_lst = [item.get() for item in l1_lst]
        l1_np_array = np.array(l1_lst)
        best_scaling_factor = scaling_factor_lst[np.argmin(l1_np_array)]
        mem_subtracted = self.image_gpu - raw_image_average_2d * (best_scaling_factor + bias)
        return mem_subtracted
        # plt.figure(figsize=(10, 5))
        # plt.subplot(1, 5, 1)
        # plt.imshow(self.image, cmap='gray')
        # plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'r-')
        # plt.plot(control_points[:, 0], control_points[:, 1], 'g.')
        # plt.plot(control_points[:, 0], control_points[:, 1], 'r.')
        # plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'b-')
        # plt.subplot(1, 5, 2)
        # plt.imshow(mem_mask.get(), cmap='gray')
        # plt.subplot(1, 5, 3)
        # plt.imshow(raw_image_average_2d.get(), cmap='gray')
        # plt.subplot(1, 5, 4)
        # plt.plot(scaling_factor_lst, l1_lst)
        # plt.subplot(1, 5, 5)
        # plt.imshow(mem_subtracted.get(), cmap='gray')
        # plt.show()
        # plt.subplot(1, 2, 1)
        # plt.imshow(convolve2d(self.image_gpu, kernel, mode = 'same').get(), cmap='gray')
        # plt.plot(fitted_curve_points[:, 0], fitted_curve_points[:, 1], 'r-')
        # plt.plot(control_points[:, 0], control_points[:, 1], 'g.')
        # plt.subplot(1, 2, 2)
        # plt.imshow(convolve2d(mem_subtracted, kernel, mode = 'same').get(), cmap='gray')
        # plt.show()

if '__main__' == __name__:
    # control_points = [[11.05309520070938, 80.56178412233524], [51.54682795430352, 59.511233704650024], [91.2018857471072, 62.0697360695029], [113.50480318281029, 64.54226886636641]]
    control_points = [[5.5349739658385255, 169.65729678156558], [104.48848503071727, 117.41915818396235], [171.06029201721796, 123.43287091795106], [238.51240318219482, 130.74611506167483]]
    control_points = np.array(control_points)
    # zoom the points
    # control_points[:, 0] = control_points[:, 0] * 2
    # control_points[:, 1] = control_points[:, 1] * 2
    rawimage_filename = '/data3/kzhang/cryosparc/CS-vsv/J336/extract/2023-02-06_18.00.22_VSV_slot3_grid3_2250X_sq01_77-32_0002_X-1Y-1-2_patch_aligned_doseweighted_particles.mrc'
    section = 105
    psi = 202.959167
    origin_x = -4.225000
    origin_y = 4.875000
    MembraneSubtract(control_points, rawimage_filename, section, psi, origin_x, origin_y, 1.068000).mem_subtract()