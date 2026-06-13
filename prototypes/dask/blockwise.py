# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import xarray as xr
import numpy as np
import numpy.ma as ma
import dask
from pathlib import Path
from sympy import Symbol,var
import bgc_md2.models.cable_all.cableHelpers as cH
import netCDF4

from functools import reduce

# -

from bgc_md2.sitespecificHelpers import get_client,getCluster
from dask.distributed import Client

#cluster=getCluster()
client=Client(cluster)
client

client=get_client()

# +
#cluster.close()

# +
x=dask.array.ones((10000,1000,96),chunks=(10000,1000,1))

# an example without contraction
def f(x):
    print(x.shape)
    # just do something expensive
    return reduce(lambda acc,el:acc+el,[x for i in range(5)])

res=dask.array.blockwise(f,'ijk',x,'ijk',dtype='f8').compute()
x


# -

# if the result has  one  index less than the argument (here the i)
# this implies a contraction
# then f has to accept an interable for x 
def f(xs):
    y=xs[0].sum(0)
    print(
        '##################',
        "len",len(xs),
        "xs[0].shape",
        xs[0].shape,
        'y.shape=',
        y.shape
    )
    return y
print('#######################################################################################')
res=dask.array.blockwise(f,'jk',x,'ijk',dtype='f8')
res.compute()

# now we do something more cable like
ntime = 30
npatch = 10
nland = 50
npool=9
ifv=-999999 # pseudo integer fill value
ffv=-1e36 # pseudo float fill value
times=dask.array.arange(0,ntime)
iveg=np.ones((npatch,nland),dtype=np.int) # lets exclude some patch,land combies
iveg[0,0]=ifv
X0=dask.array.where(
    condition=iveg==ifv, 
    x=np.nan,
    y=dask.array.ones((npool,npatch,nland),chunks=(npool,npatch,1))
)
B=dask.array.where(
    condition=dask.array.from_array(iveg==ifv), 
    x=np.nan,
    y=0.1*dask.array.ones((ntime,npool,npool,npatch,nland),chunks=(ntime,npool,npool,npatch,1))
)
iveg.shape,X0.shape


# +
def solve(start_values,Bs):
    xs = np.nan*np.ones((len(Bs), len(start_values)))
    xs[0, :] = start_values
    for k in range(0, len(Bs)-1):
        xs[k+1] = Bs[k] @ xs[k] #+ net_Us[k]

    return xs

def chunk_trajectories(iveg,x0c,times,bc):
#def trajectory(x0,times):
    npool,npatch,nland=x0c.shape
    ntime=times.shape[0]
    
    xt=np.full((ntime,npool,npatch,nland),np.nan) #default to nan 
    cond=(iveg!=ifv)
    ps,lps=np.nonzero(cond)
    nv=len(ps)
    for i in range(nv):
        btl=bc[:,:,:,ps[i],lps[i]]
        Bs=[btl[it,:,:] for it in range(ntime)]
        x0=x0c[:,ps[i],lps[i]] 
        sol=solve(x0,Bs)
        print(sol.shape)
        xt[:,:,ps[i],lps[i]]=sol
        
    #xt=dask.array.where(
    #    condition=cond,
    #    x=dask.array.stack(
    #        [x0*i for i in range(len(times))],
    #        0
    #    ),
    #    y=np.nan 
    #)
    print(xt.shape)
    return xt

X=dask.array.blockwise(chunk_trajectories,'ijkl',iveg,'kl',X0,'jkl',times,'i', B, 'ijjkl',dtype='f8')
#X=dask.array.blockwise(trajectory,'ijkl',X0,'jkl',times,'i',dtype='f8')
X.compute()    
        
        
# -








