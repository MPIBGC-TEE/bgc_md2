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

import xarray as xr
import numpy as np
import pandas as pd
import dask
from dask.distributed import LocalCluster, Client
import pathlib
import shutil
from getpass import getuser
from pathlib import Path
import dask.array as da

# +
port_dict = {
    'mm':8789,
    'hmetzler':8790, # change at will
#    'hmetzler':8888, # change at will
    'cs':8791        # change at will
}
my_user_name = getuser()
print('username:', my_user_name)

my_port = port_dict[my_user_name]
print('notebook port:', my_port)

# prevent worker from stupid too early memory shuffling
# seems to be ignored though...
# needs to be added manually to worker while in progress
worker_kwargs = {
#    'memory_limit': '2G',
    'memory_target_fraction': 0.95,
    'memory_spill_fraction': 0.95,
    'memory_pause_fraction': 0.95,
#    'memory_terminate_fraction': False, # leads to errors if commented in
}

# dashboard needs a different port for accessing it remotely
my_dashboard_port = my_port + 5
my_cluster = LocalCluster(
    dashboard_address='localhost:'+str(my_dashboard_port),
    n_workers=48,
    threads_per_worker=1,
#    memory_limit="1GB"
#    **worker_kwargs
)
print('dashboard port:', my_dashboard_port)

Client(my_cluster)

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

# +
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya
filestem = "Greg_2020_10_26/"
output_folder = "output/"
#pwc_mr_fd_archive = data_folder + output_folder + 'pwc_mr_fd/'

logfilename = data_folder + filestem + output_folder + "pwc_mr_fd_notebook.log"

#ds = xr.open_mfdataset(data_folder + filestem + "SUM*.nc")
#ds = xr.open_mfdataset(data_folder + filestem + "small_netcdf/*.nc")

zarr_path = Path(data_folder).joinpath(filestem).joinpath("zarr_version")
variable_paths = [p for p in zarr_path.iterdir() if p.is_dir()]

non_data_vars = ['lat', 'lon', 'prob', 'time']

variable_names = []
variables = []
for variable_path in variable_paths:
    name = variable_path.name
    if name not in non_data_vars:
        variable_names.append(name)
        variables.append(da.from_zarr(str(variable_path)))
        
variables.append(da.from_zarr(str(zarr_path.joinpath('lat')), chunks=(1,)).reshape(-1, 1, 1, 1))
variable_names.append('lat')
variables.append(da.from_zarr(str(zarr_path.joinpath('lon')), chunks=(1,)).reshape(1, -1, 1, 1))
variable_names.append('lon')
variables.append(da.from_zarr(str(zarr_path.joinpath('prob')), chunks=(1,)).reshape(1, 1, -1, 1))
variable_names.append('prob')
variables.append(da.from_zarr(str(zarr_path.joinpath('time'))).reshape(1, 1, 1, -1))
variable_names.append('time')

variables

# +
pp_da=dask.array.from_array(pp_np,chunks=(1000, 9, 9, 2))
g_da=dask.array.stack([pp_da for i in range(nland // 30)])

list_ga = [g_da for _ in range(10)]
def add(*args):
#    print(args)
#    r = sum(args)
#    print(r.shape)
    return np.ones((1, 1, 1, 6))

#res = pp_da.map_blocks(
#    add, 
#    *list_ga[1:],
#    drop_axis=4,
#    new_axis=4,
#    meta=np.ndarray((1, 1000, 9, 9, 1), dtype=np.float64)
#)

res = variables[0].map_blocks(
    add,
    *variables[1:],
    drop_axis=3,
    new_axis=3,
    meta=np.ndarray((1, 1, 1, 6), dtype=np.float64),
    chunks=(1, 1, 1, 6)
)
res
# -

g_da

# +
dataDir = '/home/data/CARDAMOM/zarrtest'
p=pathlib.Path(dataDir)
if p.exists():
    shutil.rmtree(p)
        
p.mkdir()
res_delayed = res.to_zarr(dataDir+'/output.zarr', compute=False)
#import h5py
#f=h5py.File(dataDir+'/res.hdf5',mode='w')
#d = f.require_dataset('/data', shape=res.shape, dtype=res.dtype)
#dask.array.store(res,d)
# -

res_delayed.visualize()

# +
# %%time

_ = res_delayed.compute()
# -

# we can also create an xarray dataArray and write a netcdf file
# %time  xr.DataArray(res).to_netcdf(dataDir+'/output.nc')




