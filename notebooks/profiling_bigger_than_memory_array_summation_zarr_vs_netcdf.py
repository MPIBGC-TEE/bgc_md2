# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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
import pandas as pd
import dask
import pathlib
import shutil

# the next to line refer only to our infrastructure and conventions in Jena
from bgc_md2.sitespecificHelpers import get_client
client = get_client()

# +
ntime,npool,npatch,nland = 30000,9,10,6000
pixel_arr_shape=(ntime,npool,npool,npatch)

pixel_arr_dims=('time','poolx','pooly','patch')
pp_np = np.ones(pixel_arr_shape)

# +
pp_xr = xr.DataArray(pp_np,dims=pixel_arr_dims)

# The first objective is to create a huge global array that is bigger than memory from the pixel wise arrays
# (This would e.g. happen if the global array would be created by a netcdf file)
# Unfortunately a simple concatenation, does not work since it will try to allocate the 
# memory to store it.(1.06 TB)

# xr.concat([pp_xr for i in range(nland)],dim='land')
    
# so instead we create a dask array directly which is lazy and can be computed chunkwise while being written to disk
# -

pp_da=dask.array.from_array(pp_np,chunks=pixel_arr_shape)
g_da=dask.array.stack([pp_da for i in range(nland)])
res=g_da+g_da

# +
dataDir = '/home/data/cable-data/zarrtest'
p=pathlib.Path(dataDir)
if p.exists():
    shutil.rmtree(p)
        
p.mkdir()
# %time res.to_zarr(dataDir+'/output.zarr')
#import h5py
#f=h5py.File(dataDir+'/res.hdf5',mode='w')
#d = f.require_dataset('/data', shape=res.shape, dtype=res.dtype)
#dask.array.store(res,d)
# -

# we can also create an xarray dataArray and write a netcdf file
# %time  xr.DataArray(res).to_netcdf(dataDir+'/output.nc')




