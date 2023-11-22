import numpy as np
import matplotlib.pyplot as plt
from ._utils import *

def create_gaussian_low_pass_filter(image, cutoff_frequency):
    shape = image.shape
    rows, cols = shape
    crow, ccol = rows // 2, cols // 2
    y, x = np.ogrid[-crow:rows-crow, -ccol:cols-ccol]
    distance = np.sqrt(x*x + y*y)
    mask = np.exp(-(distance**2) / (2*cutoff_frequency**2))

    f = np.fft.fft2(image)
    fshift = np.fft.fftshift(f)
    fshift_filtered = fshift * mask
    f_ishift = np.fft.ifftshift(fshift_filtered)
    img_back = np.fft.ifft2(f_ishift)
    img_back = np.real(img_back)
    return img_back

# 加载图像
img = readmrc('2023-02-06_17.53.19_VSV_slot3_grid3_2250X_sq01_77-7_0002_X-1Y-1-2_patch_aligned_doseweighted_particles.mrc', section=31)
img = (img - np.mean(img)) / np.std(img)
img_back = create_gaussian_low_pass_filter(img, cutoff_frequency=20)
# 快速傅里叶变换，转到频域
# f = np.fft.fft2(img)
# fshift = np.fft.fftshift(f)

# # 创建高斯低通滤波器
# mask = create_gaussian_low_pass_filter(img, cutoff_frequency=25)

# # 在频域上应用滤波器
# fshift_filtered = fshift * mask

# # 逆FFT转回空间域
# f_ishift = np.fft.ifftshift(fshift_filtered)
# img_back = np.fft.ifft2(f_ishift)
# img_back = np.abs(img_back)

# 显示结果
plt.figure()
plt.subplot(131), plt.imshow(img, cmap='gray'), plt.title('Original Image')
# plt.subplot(132), plt.imshow(mask, cmap='gray'), plt.title('Gaussian Low Pass Filter')
plt.subplot(132), plt.imshow(img_back, cmap='gray'), plt.title('Filtered Image')
plt.show()
