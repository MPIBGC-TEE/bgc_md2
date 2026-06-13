import numpy as np
from numba import cuda

@cuda.jit
def f(a, b, c):
    # like threadIdx.x + (blockIdx.x * blockDim.x)
    tid = cuda.grid(1)
    size = len(c)
    if tid < size:
        c[tid] = a[tid] + b[tid]

N = 1000000
a = cuda.to_device(np.random.random(N))
b = cuda.to_device(np.random.random(N))
c = cuda.device_array_like(a)
f.forall(len(a))(a, b, c)
print(c.copy_to_host())
