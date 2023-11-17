import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule
import numpy as np
import cupy as cp

kernel_code = """
__global__ void distribute_bilinearly(float* image, float* count_image,
                                      int image_width, int image_height,
                                      float x, float y, float value) {
    int x1 = int(x), y1 = int(y);
    int x2 = x1 + 1, y2 = y1 + 1;
    if (x2 >= image_width || y2 >= image_height) {
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
"""

mod = SourceModule(kernel_code)
distribute_bilinearly_cuda = mod.get_function("distribute_bilinearly")

def distribute_bilinearly(image, count_image, x, y, value):
    # 确保数据是正确的类型
    image = image.astype(cp.float32)
    count_image = count_image.astype(cp.float32)
    # 分配设备内存
    image_gpu = cuda.mem_alloc(image.nbytes)
    count_image_gpu = cuda.mem_alloc(count_image.nbytes)
    
    # 将数据复制到设备
    cuda.memcpy_htod(image_gpu, image)
    cuda.memcpy_htod(count_image_gpu, count_image)
    
    # 启动CUDA内核
    distribute_bilinearly_cuda(image_gpu, count_image_gpu,
                               np.int32(image.shape[1]), np.int32(image.shape[0]),
                               np.float32(x), np.float32(y), np.float32(value),
                               block=(16, 16, 1), grid=(1, 1))
    
    # 将结果从设备复制回主机
    cuda.memcpy_dtoh(image, image_gpu)
    cuda.memcpy_dtoh(count_image, count_image_gpu)
    
    return image, count_image


# 测试函数
image = np.zeros((10, 10), dtype=np.float32)
count_image = np.zeros((10, 10), dtype=np.float32)
x, y, value = 5, 5, 10

new_image, new_count_image = distribute_bilinearly(image, count_image, x, y, value)
