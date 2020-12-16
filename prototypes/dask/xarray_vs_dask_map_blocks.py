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

# # Use dask.array.map_blocks instead of xarray.map_blocks
#
# xarray objects have a `data` attribute which is a `dask.array`. With this object real streamprocessing is possible
# The api is pretty easy and the blockwise operations accept and return `numpy` arrays.
# This is what we usually need.
# It turns out that the API of `xarray.map_blocks` is much more obscure and makes it much harder (if at all possible)  to really controll what is going on. (which is necessary in our case to avoid the memory overflows). 
# The code below uses an array from the cable output and performs some 
# array operations on it (a maks and a plus)
# We first use the operator syntax which hides the parallelisation in the background and then reimplement what is going on behind the scenes by means of `dask.array.map_blocks`.
# We also measure the execution time.
#
# The demonstrated method would also work with arrays bigger than memory. We just would not call `compute()` but 
# rather store the result by `to_zarr`. 
# (Why netcdf writing does not work as intended is explained in another example notebook) 

# +
from getpass import getuser
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
from sympy import Symbol,var
port_dict = {
    'mm':(8689,8789),       # first is the port of the actual scheduler, second the port for the dashboard
    'hmetzler':(8690,8790), # change at will to a port you forward via ssh to your local machine
    'cs':(8691,8791)        # change at will
}
my_user_name = getuser()
addr = 'localhost:'+str(port_dict[my_user_name][0])
try:
    Client(addr)
except IOError:
    my_cluster = LocalCluster(
        scheduler_port= port_dict[my_user_name][0],   
        dashboard_address='localhost:'+str(port_dict[my_user_name][1])
    )
    Client(my_cluster)#same as Client(addr)

my_user_name
# -

# we use a cable array as example and therefore read the dataset
cableDataDir = '/home/data/cable-data/example_runs'
runId = "parallel_1901_2004_with_spinup"
outDir = "output/new4"
first_yr = 1901
last_yr = 2004
fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
outpath = Path(cableDataDir).joinpath(runId, outDir)
ps = [outpath.joinpath(fn) for fn in fns]
ds = xr.open_mfdataset(
    paths=ps, 
    combine="by_coords",
    parallel=True # use dask
);ds


#The default chunking creates 365 (days) wide chunks along the time axis 
ds.chunks

#but we want chunks of 1 pixel (land)
dd=ds.dims
chunk_dict = {'time':dd['time'],'patch':dd['patch'],'land':1}
# We can rechunk the whole dataset (or could have given a chunk argument to  xr.open_mfdataset)
ds=ds.chunk(chunk_dict)

#the target variable for this is iveg
# we read ds.iveg and replace nan by the value 18 
mask=xr.ufuncs.isnan(ds.iveg) # boolean array
filtered = xr.where(
    mask,
    18,# default fill value for nan
    ds.iveg
)
# %time filtered.compute()

# We actually need to substract 1
iveg = filtered-1
# %time iveg.compute()
# This takes about as long as the filtering
# suggesting that the main effort is spent distributing and collecting the data

# Suspecting that the most time is spent in infrastructure
# we try to combine both operations in one compute call
filtered= xr.where(
    mask,
    18,# default fill value for nan
    ds.iveg
)
iveg = filtered-1
# %time iveg.compute()
# which confirms our suspicions 
# and motivates us to put combine several actions into one function
# and map this over all chunks

# The syntax: `arr -1` takes advantage of the `-` operator being heavily overloaded by `xarray`
# which internally does several things:
# 1. apply the function chunkwise
# 1. for every chunk use numpy broadcasting rules to create an array full of `1` (np.ones(chunk.shape))
# 1. for every chunk use the numpy ufunc (general vectorized function) np.add to add the chunk of iveg to the ones array. 
#
# Similarly the filter operation `where` that precedes it is a chunkwise application of `np.where`.
# We can more explicitly write a function and map it over the list of chunks.
# We have several options for the actual `map` function we use.
# 1. `dask.array.map_blocks` (can be used on the `.data` attrbute of an `xarray.DataArray` object)
#     This is working as expected (next cell)
# 1. `xarray.map.blocks` (experimental. According to xarrays help it should be used if the blocks are 
#     `xarray.DataArray` objects which is not the case for us (we want numpy arrays for the chunks)
#   
# 1.  `xarray.apply_ufunc`(which can use dask in the backgroud in different ways depending on the value of 
#      ists`dask=`keyword parameter)
#      Difficult to get the parameters right and to guess how dask will be called in the background. 
#      It is much easier to extract the `data` attribute of a `xarray.DataArray` which is a `dask.array`
#      and apply the `dask.array.map_blocks` method directly

# +
# 1.) with dask.array.map_blocks
def filter_and_substract(chunk):
    mask=np.isnan(chunk)
    filtered_chunk=np.where(
        mask,
        18,
        chunk
    )
    return np.add(filtered_chunk,-np.ones(chunk.shape))

data=ds.iveg.data #extract the internal dask.array object 
iveg = data.map_blocks(filter_and_substract) # call the dask method

# %time iveg.compute()
# +
# 2.) with xarray.map_block 
# Semantically very different and confusing
# The arguments of the mapped function appear to be 
# the global arrays as opposed to chunks
# consequently we use a huge global array with the same shape as iveg 
# but filled with ones
my_ones = xr.ones_like(ds.iveg)
# iveg = xr.map_blocks(lambda x,y:np.add(x,y),ds.iveg, args=(my_ones,))

def filter_and_substract_xr(x,y):
    mask=xr.ufuncs.isnan(x) 
    return np.add(
        xr.where(
            mask,
            18,# default fill value for nan
            x
        ),
        -y
   ) 

iveg = xr.map_blocks(filter_and_substract_xr,ds.iveg,args=(my_ones,)) 
# %time iveg.compute()
# -

# for comparison the single threaded version (pure numpy, not even parallel)
arr=np.ones((37860,10,5000))
# %time arr+arr


