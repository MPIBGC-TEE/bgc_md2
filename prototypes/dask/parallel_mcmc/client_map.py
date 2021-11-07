import dask.array as da
from dask.distributed import Client
from dask.distributed import LocalCluster
import numpy as np
from time import sleep
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

c=Client()
futures = c.map(f_chunk,ivs)


future.visualize()

# the return values are not the chains themselfes but 'future objects' that
# are like a promised result
# to get the actual result we will have to use the result() method.

res1 = [f.result() for f in futures]
res2 = c.gather(f_chunk,ivs)


# we can even build a dask.array from the list of chains, althouhg it will be unevenly chunked due to the different lenghts of the returned chains
distributed_single_chain=da.concatenate([f.result() for f in futures],axis=1)

# Note that the dask array is distributed. Every chunk can potentially live on # a different node....
# So we could use parallel postprocessing as well
# We could also rechunk it for better loadbalancing but 
# rechunking involves communication




