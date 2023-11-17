import mrcfile
import starfile
import numpy as np
import cupy as cp
import pandas as pd
import matplotlib.pyplot as plt
from cupyx.scipy.ndimage import convolve

def readmrc(filename, section=0, mode='cpu'):
    image = mrcfile.open(filename)
    if len(image.data.shape) == 2:
        gray_image = image.data.copy()
    elif len(image.data.shape) == 3: 
        gray_image = image.data[section].copy()
    else:
        raise ValueError("Unsupported image dimensions")
    if mode == 'cpu':
        return gray_image
    elif mode == 'gpu':
        return cp.asarray(gray_image)

def savemrc(image, filename):
    image = image.astype(np.float32)
    with mrcfile.new(filename, overwrite=True) as mrc:
        mrc.set_data(image)


def distance(x1, y1, x2, y2):
    return np.sqrt((x1-x2)**2 + (y1-y2)**2)

def create_starfile(Average_raw2dimage_filename):
    filename = 'mem_analysis.star'
    images = mrcfile.open(Average_raw2dimage_filename)
    image_sections = images.data.shape[0]
    # df_star = pd.DataFrame(data=[[0] * 9], columns=['rlnRawImageName', 'rln2DAverageimageName', 'rlnCenterX', 'rlnCenterY', 'rlnAngleTheta', 'rlnMembraneDistance', 'rlnSigma1', 'rlnSigma2', 'rlnCurveKappa'])
    df_star = pd.DataFrame(data=[[0] * 8], columns=['rln2DAverageimageName', 'rlnCenterX', 'rlnCenterY', 'rlnAngleTheta', 'rlnMembraneDistance', 'rlnSigma1', 'rlnSigma2', 'rlnCurveKappa'])
    for i in range(image_sections):
        # df_star.loc[i] = [f'{rawimage_name}', f'{i+1:06d}@{Average_raw2dimage_filename}', 0, 0, 0, 0, 0, 0, 0]
        df_star.loc[i] = [f'{i+1:06d}@{Average_raw2dimage_filename}', 0, 0, 0, 0, 0, 0, 0]
    starfile.write(df_star, filename, overwrite=True)

def readstar(starfile_name):
    star = starfile.read(starfile_name)
    return star

def write_star(df, starfile_name='mem_analysis.star'):
    starfile.write(df, starfile_name, overwrite=True)

def create_gaussian_low_pass_filter(image, cutoff_frequency):
    shape = image.shape
    rows, cols = shape
    crow, ccol = rows // 2, cols // 2
    y, x = cp.ogrid[-crow:rows-crow, -ccol:cols-ccol]
    distance = cp.sqrt(x*x + y*y)
    mask = cp.exp(-(distance**2) / (2*cutoff_frequency**2))

    f = cp.fft.fft2(image)
    fshift = cp.fft.fftshift(f)
    fshift_filtered = fshift * mask
    f_ishift = cp.fft.ifftshift(fshift_filtered)
    img_back = cp.fft.ifft2(f_ishift)
    img_back = cp.real(img_back)
    return img_back


def gaussian_kernel(size, sigma=1.0):
    size = int(size) // 2
    x, y = cp.mgrid[-size:size+1, -size:size+1]
    normal = 1 / (2.0 * cp.pi * sigma**2)
    g =  cp.exp(-((x**2 + y**2) / (2.0*sigma**2))) * normal
    return g

# def convolve(image, kernel):
#     i_height, i_width = image.shape
#     k_height, k_width = kernel.shape

#     pad_height = k_height // 2
#     pad_width = k_width // 2
#     padded_image = cp.pad(image, ((pad_height, pad_height), (pad_width, pad_width)), mode='constant')

#     output = cp.zeros_like(image, dtype=float)

#     for i in range(i_height):
#         for j in range(i_width):
#             output[i, j] = cp.sum(padded_image[i: i + k_height, j: j + k_width] * kernel)
#     return output

# kernel = gaussian_kernel(5, 1)
# gray_image = readmrc('2023-05-28_00.21.16_neuron_152-3_00013_X-1Y-1-1_patch_aligned_doseweighted_particles.mrc', section=0, mode='gpu')
# print(type(gray_image))
# print(type(kernel))
# blurred_image = convolve(gray_image, kernel)
# # blurred_image = convolve(gray_image, kernel)
# # image_lp = create_gaussian_low_pass_filter(gray_image, 15)
# fig, ax = plt.subplots(1, 2, figsize=(10, 10))
# ax[0].imshow(gray_image.get(), cmap='gray')
# ax[1].imshow(blurred_image.get(), cmap='gray')
# # ax[2].imshow(image_lp, cmap='gray')
# plt.show()
