import dask.array
from dask.distributed import Client,LocalCluster
from functools import reduce
import zarr as zr
import numpy as np
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices

from bgc_md2.sitespecificHelpers import getCluster,get_client

#clu=getCluster()
cli=get_client()

# +
#def main():
ifv=-99999
nland = 56
npatch = 4 
ntime= 30
npool = 9

# fake an iveg array where half the patches are valid
iveg = dask.array.stack(
    [
        np.concatenate(
            [        
                np.full((int(npatch/2),),ifv,dtype=np.int32),
                np.full((int(npatch/2),),1,dtype=np.int32)
                #np.random(1,3,size=((int(npatch/2),),1,dtype=np.int32)
            ]
        )
        for i in range(nland)
    ],
    axis=-1
).rechunk((npatch,1))

# fake B and u arrays
B=dask.array.stack(
    [
        #np.zeros((ntime,npool,npool,npatch))
        np.random.uniform(1,3,size=(ntime,npool,npool,npatch))
        for i in range(nland)
    ],
    axis=-1
)

u=dask.array.stack(
    [
        #np.zeros((ntime,npool,npool,npatch))
        np.random.uniform(1,3,size=(ntime,npool,npatch))
        for i in range(nland)
    ],
    axis=-1
)

x0=dask.array.stack(
    [
        #np.zeros((ntime,npool,npool,npatch))
        np.random.uniform(1,3,size=(npool,npatch))
        for i in range(nland)
    ],
    axis=-1
)

#print(iveg)
#print(B)
nz = dask.array.nonzero((iveg==ifv))
#print(nz)

# +
B_val = cH.valid_combies(nz,B)
u_val = cH.valid_combies(nz,u)
x0_val = cH.valid_combies(nz,x0)

#print(B_val)
#the following line would cause a global B_val.compute() which we can not afford memory wise.
# B_val.to_zarr('tmp/zarr.test')
Bz = zr.open('tmp/B.zarr',mode='w',shape=B_val.shape,chunks=B_val.chunksize)
uz = zr.open('tmp/u.zarr',mode='w',shape=u_val.shape,chunks=u_val.chunksize)
x0z = zr.open('tmp/x0.zarr',mode='w',shape=x0_val.shape,chunks=x0_val.chunksize)

# +
ncores=95
slices=batchSlices(B_val.shape[-1],ncores)
for s in slices[0:1]:
    Bz[:,:,:,s]=B_val[:,:,:,s].compute()
    uz[:,:,s]=u_val[:,:,s].compute()
    x0z[:,s]=x0_val[:,s].compute()

#while i < int(l/nprocs)*nprocs:
#    s = slice(i,min(i+nprocs,l-1,1))
#    i += nprocs
#if __name__ == "__main__":
#    c=Client()
#    fut = c.submit(main)
#    fut.result()
# -

# now we read the data back and compute with it
Bfz=dask.array.from_zarr('tmp/B.zarr')
ufz=dask.array.from_zarr('tmp/u.zarr')
x0fz=dask.array.from_zarr('tmp/x0.zarr')
Bfz



