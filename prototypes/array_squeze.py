import dask.array
from dask.distributed import Client,LocalCluster
from functools import reduce

import zarr as zr
import numpy as np


def valid_combies(nz,mat):
    s=mat.shape
    s
    ps,lps=nz
    ps.compute_chunk_sizes()
    lps.compute_chunk_sizes()

    def reform(combi):
        return combi.map_blocks( 
            lambda timeLine:np.expand_dims(timeLine, 3),
            new_axis=3
        )

    i = 0
    mat_new = reform(mat[:, :, :, ps[i], lps[i]])
    print(mat_new.shape)
    # we now construct the squezed dask array (without computing it) 
    while i < len(ps)-1:
        i += 1
        mat_new = dask.array.concatenate(
            [
                mat_new,
                reform(mat[:, :, :, ps[i], lps[i]])
            ],
            axis=3
        )
    # mat_new = reduce(
    #    lambda acc,combi:dask.array.concatenate([acc,combi],axis=3),
    #    (
    #        reform(mat[:,:,:,ps [i],lps[i]]) for i in range(len(ps))
    #    )
    #)
    return mat_new

    

#def main():
ifv=-99999
nland = 56
npatch = 10
ntime= 3000
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
    
B=dask.array.stack(
    [
        #np.zeros((ntime,npool,npool,npatch))
        np.random.uniform(1,3,size=(ntime,npool,npool,npatch))
        for i in range(nland)
    ],
    axis=-1
)
#print(iveg)
#print(B)
nz = dask.array.nonzero((iveg==ifv))
#print(nz)
B_val = valid_combies(nz,B)
print(B_val)
#the following line would cause a global B_val.compute() which we can not afford memory wise.
# B_val.to_zarr('tmp/zarr.test')
z = zr.open('tmp/test.zarr',mode='w',shape=B_val.shape,chunks=B_val.chunksize)
nprocs = 8
i=0
l=B_val.shape[-1]
print(type(l))
while i < int(l/nprocs)*nprocs:
    s = slice(i,min(i+nprocs,l-1,1))
    z[:,:,:,s]=B_val[:,:,:,s].compute()
    i += nprocs
#if __name__ == "__main__":
#    c=Client()
#    fut = c.submit(main)
#    fut.result()
