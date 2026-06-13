# +
import dask
from dask.distributed import Client
from dask.distributed import LocalCluster


lc= LocalCluster()
c=Client(lc)
lc.worker_spec


# -

import numpy as np
n_chains=20
l_chain=100
iv=np.array([1.0,2.0,3.0])
npar=len(iv)
ivs=[iv for i in range(n_chains)]


def f_chunk(iv):
    #print(iv.shape)
    g = np.random.default_rng()
    npar =len(iv)
    res = np.zeros(shape=(npar,l_chain))
    # choose an acceptance rate randomly to model the worst case scenario.
    # that the parallel chains behave very differently
    acr = g.random()
    s = 0
    for i in range(l_chain): 
        v = g.random()
        if v < acr:
            res[:,s]=np.ones(npar)*v
            s += 1

    return res[:,0:s]

d1= dask.delayed(f_chunk(iv))
d2= dask.delayed(f_chunk(iv))
d3= dask.delayed(f_chunk(iv))
d4= dask.delayed(f_chunk(iv))

res=dask.delayed(sum)(d1,d2,d3,d4)



res.visualize()


