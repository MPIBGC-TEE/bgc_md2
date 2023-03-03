from numba import jit, njit, cuda, double,prange
import numpy as np
from sympy import Symbol, symbols, Matrix, lambdify, sin
from scipy.interpolate import interp1d
from inspect import getsource
from string import Template

x = Symbol("x")
es = sin(x)**2

f1 = lambdify(x, es, ['numpy'])
ts = np.linspace(0, 1, 100)
ys = f1(ts)
f2 = interp1d(ts,ys)
nt = 1000
n_p = 10
n_pix=2000
tfs = np.linspace(0, 1, nt)

Bs = np.ones((n_p, n_p ,nt,n_pix)) 
print(f1(0))
# to measure exec time
from timeit import default_timer as timer

# normal function to run on cpu
code=Template("""def ${name}(Bs,tfs):
    res=np.zeros(shape=(n_p,nt,n_pix))
    for j in ${r}(Bs.shape[3]):
        for i in range(1,len(tfs)):
            res[:, i, j] = Bs[:, :, i, j]@np.ones(shape=(n_p))
    return res 
""")
exec(code.substitute(name="func",r="range"))

# function optimized to run on gpu
dec_code="""@jit(
    double[:,:,:](double[:,:,:,:],double[:]),
    target_backend="cuda",
    #target_backend="cpu",
    nopython=True,
    fastmath=True
)
"""+code.substitute(name="func2",r="prange")
exec(dec_code)

if __name__ == "__main__":
    start = timer()
    func(Bs,tfs)
    print("without numba:", timer() - start)

    start = timer()
    func2(Bs,tfs)
    print("with numba:", timer() - start)
