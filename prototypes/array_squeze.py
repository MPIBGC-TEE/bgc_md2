import dask.array
from dask.distributed import Client, LocalCluster
from functools import reduce
import zarr as zr
import numpy as np
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices

# from bgc_md2.sitespecificHelpers import getCluster,get_client

if __name__ == "__main__":
    c = Client()
    ifv = -99999
    nland = 56
    npatch = 4
    ntime = 30
    npool = 9

    # fake an iveg array where half the patches are valid
    iveg = dask.array.stack(
        [
            np.concatenate(
                [
                    np.full((int(npatch / 2),), ifv, dtype=np.int32),
                    np.full((int(npatch / 2),), 1, dtype=np.int32)
                    # np.random(1,3,size=((int(npatch/2),),1,dtype=np.int32)
                ]
            )
            for i in range(nland)
        ],
        axis=-1,
    ).rechunk((npatch, 1))

    # fake B and u arrays
    B = dask.array.stack(
        [
            # np.zeros((ntime,npool,npool,npatch))
            np.random.uniform(1, 3, size=(ntime, npool, npool, npatch))
            for i in range(nland)
        ],
        axis=-1,
    )

    u = dask.array.stack(
        [
            # np.zeros((ntime,npool,npool,npatch))
            np.random.uniform(1, 3, size=(ntime, npool, npatch))
            for i in range(nland)
        ],
        axis=-1,
    )

    x0 = dask.array.stack(
        [
            # np.zeros((ntime,npool,npool,npatch))
            np.random.uniform(1, 3, size=(npool, npatch))
            for i in range(nland)
        ],
        axis=-1,
    )

    nz = dask.array.nonzero((iveg == ifv))
    n_val = int(npatch/2*nland)

    B_val = cH.valid_combies_parallel(nz, B)
    u_val = cH.valid_combies_parallel(nz, u)
    x0_val = cH.valid_combies_parallel(nz, x0)
    assert(B_val.shape == (ntime, npool, npool, n_val))
    assert(u_val.shape == (ntime, npool, n_val))
    assert(x0_val.shape == (npool, n_val))

    # write them to disk
#    if dir_p.exists():
#        shutil.rmtree(dir_p)

    cH.batchwise_to_zarr(B_val, 'B.zarr', rm=True)
    cH.batchwise_to_zarr(u_val, 'u.zarr', rm=True)
    cH.batchwise_to_zarr(x0_val, 'x0.zarr', rm=True)

    # and read them again

    assert(dask.array.all(dask.array.from_zarr('B.zarr') == B_val).compute())
    assert(dask.array.all(dask.array.from_zarr('u.zarr') == u_val).compute())
    assert(dask.array.all(dask.array.from_zarr('x0.zarr') == x0_val).compute())
