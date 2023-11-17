import pycuda.autoinit
import pycuda.driver as cuda
from pycuda.compiler import SourceModule

# 编写CUDA内核代码
kernel_code = """
__global__ void double_array(float *a, int N) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N) {
        a[idx] *= 2;
    }
}
"""

# 编译CUDA内核
mod = SourceModule(kernel_code)

# 获取内核函数
double_array = mod.get_function("double_array")

import numpy as np

# 创建一个NumPy数组
a = np.array([1, 2, 3, 4, 5], dtype=np.float32)

# 将NumPy数组复制到GPU
a_gpu = cuda.mem_alloc(a.nbytes)
cuda.memcpy_htod(a_gpu, a)

# 设置线程块和网格的大小
block_size = 256
grid_size = (len(a) + block_size - 1) // block_size

# 执行CUDA内核
double_array(a_gpu, np.int32(len(a)), block=(block_size,1,1), grid=(grid_size,1))

# 将结果从GPU复制回CPU
cuda.memcpy_dtoh(a, a_gpu)

# 打印结果
print(a)
