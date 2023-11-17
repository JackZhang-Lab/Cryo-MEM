import numpy as np
import cupy as cp
import mrcfile
from cupyx.scipy.ndimage import maximum_filter as maximum_filter_gpu

def max_filter(image):
    image = cp.array(image)
    max_filtered_image = maximum_filter_gpu(image, size=9, mode='constant')
    threshold = cp.sum(max_filtered_image) / cp.sum(max_filtered_image > 0)
    max_filtered_image[max_filtered_image <= threshold * max_filtered_image.max()] = 0
    # normalization
    max_filtered_image = max_filtered_image / max_filtered_image.max()
    return max_filtered_image

def pickpoints(image, num_points):
    points = []
    for _ in range(num_points):
        while True:
            i, j = np.random.randint(0, image.shape[0]), np.random.randint(0, image.shape[1])  # 使用 NumPy 的随机数生成器
            if np.random.random() < cp.asnumpy(image[i, j]):
                points.append((j, i))
                break
    return points

def generate_data_points(image, num_points):
    max_filtered_image = max_filter(image)
    data_points = np.array(pickpoints(max_filtered_image, num_points))
    if np.std(data_points[:, 1]) > np.std(data_points[:, 0]):
        data_points = data_points[np.argsort(data_points[:, 1])]
    else:
        data_points = data_points[np.argsort(data_points[:, 0])]
    return data_points

if '__main__' == __name__:
    with mrcfile.open('/Users/hzvictor/Downloads/cryosparc_P2_J1228_200_class_averages.mrc') as f:
        image = f.data[11]

    data_points = generate_data_points(image, 100)