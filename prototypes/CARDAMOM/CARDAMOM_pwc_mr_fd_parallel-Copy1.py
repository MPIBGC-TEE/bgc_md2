# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
from dask.distributed import Client

#import xarray as xr
import numpy as np
import zarr
from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client, LocalCluster
import dask.array as da
from getpass import getuser

import time


# +
port_dict = {
    'mm': 8789,
#    'hmetzler': 8790, # change at will
    'hmetzler': 8792, # change at will
    'cs': 8791        # change at will
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
# -

# ## How to connect to remote
# **Remark**: Port values to be adapted, see above.
#
# ### remotely
# `
# screen
# # cd GitHub/bgc_md2/notebooks/CARDAMOM
# conda activate bgc_md2
# jupyter lab --no-browser -- port=8790
# `
# ### locally
# `
# ssh -L 8080:localhost:8790 antakya_from_home
# `
#
# In browser open `localhost:8080`.
#
# To connect to bokeh dashbord
#
# `
# ssh -L 8081:localhost:8795 antakya_from_home
# `
#
# and in browser open `localhost:8081/status/`.

# +
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya
filestem = "Greg_2020_10_26/"
output_folder = "output/"
#pwc_mr_fd_archive = data_folder + output_folder + 'pwc_mr_fd/'

logfilename = data_folder + filestem + output_folder + "pwc_mr_fd_notebook_dask_trapezoidal_51_nodes:5-9_probs.log"

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
lims = (None, None, 5)

for nr, name, var in zip(range(len(variables)), variable_names, variables):
    if name not in non_data_vars:
#        variables[nr] = variables[nr][:nr_lim, :nr_lim, :nr_lim, ...]
        variables[nr] = variables[nr][:lims[0], :lims[1], :lims[2], ...]
  
name = "lat"
index = variable_names.index(name)
variables[index] = variables[index][:lims[0], ...]

name = "lon"
index = variable_names.index(name)
variables[index] = variables[index][:, :lims[1], ...]

name = "prob"
index = variable_names.index(name)
variables[index] = variables[index][:, :, 5:(5+lims[2]), ...]

#for nr, name in enumerate(['lat', 'lon', 'prob']):
#    index = variable_names.index(name)
#    idx = [slice(None)] * nr + [slice(0, lims[nr], 1)]
#    print(idx)
#    variables[index] = variables[index][idx]
    
variables


# +
#from time import sleep

def func_start_values(*args):
    v_names = args[-1]
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-1]))}
#    print([v.shape for v in d.values()])
    print('lat:', d['lat'], 'lon:', d['lon'], 'prob:', d['prob'], flush=True)
#    print('time:', d['time'], flush=True)
#    print([v.shape for v in d.values()], flush=True)

    start_values = CARDAMOMlib.load_start_values_greg_dict(d)
    ds_res = CARDAMOMlib.compute_ds_pwc_mr_fd_greg(
        ds_single,
        integration_method='trapezoidal',
        nr_nodes=51
    )
#    print(start_values, flush=True)

#    start_values = da.from_array(1.1 * np.ones((1, 1, 1, 6), dtype=np.float64, chunks=(1,1,6))
#    start_values.to_dask()
    return start_values.reshape(1, 1, 1, 6)


# +
shape = (34, 71, 50, 6)                                                                                         
chunks = (1, 1, 1, 6)                   

# set up zarr array to store data
store = zarr.DirectoryStore('TB1.zarr')
root = zarr.group(store) 
TB1 = root.zeros(
    'data', 
    shape=shape, 
    chunks=chunks, 
    dtype=np.float64
)
# -

start_values = variables[0].map_blocks(
    func_start_values,
    *variables[1:],
    variable_names,
    drop_axis=3,
    new_axis=3,
    chunks=(1, 1, 1, 6),
    dtype=np.float64,
    meta=np.ndarray((34, 71, 50, 6), dtype=np.float64)
)
start_values

# +
# %%time

#start_values.compute()
start_values_delayed = start_values.to_zarr(
    data_folder + filestem + output_folder + "start_values",
    overwrite=True,
    lock=False,
    return_stored=False,
    compute=False,
    kwargs={'chunks': (1, 1, 1, 6)}
)
#start_values_delayed = start_values.store(TB1, lock=False, compute=False)   

# +
#start_values_delayed.visualize(optimize_graph=True)

# +
# %%time

#_ = da.compute(start_values_delayed, scheduler="distributed")
start_values_delayed.compute()


# -
def func_us(*args):
    v_names = args[-1]
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-1]))}
#    print([v.shape for v in d.values()])
    print('lat:', d['lat'], 'lon:', d['lon'], 'prob:', d['prob'], flush=True)
#    print('time:', d['time'], flush=True)
#    print([v.shape for v in d.values()], flush=True)

    us = CARDAMOMlib.load_us_greg_dict(d)
#    print(start_values, flush=True)

    return us.reshape(1, 1, 1, len(d['time']), 6)


us = variables[0].map_blocks(
    func_us,
    *variables[1:],
    variable_names,
    new_axis=4,
    chunks=(1, 1, 1, variables[variable_names.index('time')].shape[-1], 6),
    dtype=np.float64,
    meta=np.ndarray((1, 1, 1, variables[variable_names.index('time')].shape[-1], 6), dtype=np.float64)
)
us

# +
# %%time

us.to_zarr(data_folder + filestem + output_folder + "us")


# +
def write_to_logfile(*args):
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    with open(logfilename, 'a') as f:
        t = (current_time,) + args
        f.write(" ".join([str(s) for s in t]) + '\n')

def func_Bs(*args):
    v_names = args[-1]
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-1]))}
#    print([v.shape for v in d.values()])
#    print('lat:', d['lat'], 'lon:', d['lon'], 'prob:', d['prob'], flush=True)
#    print('time:', d['time'], flush=True)
#    print([v.shape for v in d.values()], flush=True)

    Bs = CARDAMOMlib.load_Bs_greg_dict(
        d,
        integration_method="trapezoidal",
        nr_nodes=51
    )
#    print(start_values, flush=True)

    write_to_logfile(
        "finished single,",
        "lat:", d["lat"],
        "lon:", d["lon"],
        "prob:", d["prob"]
    )

    return Bs.reshape(1, 1, 1, len(d['time']), 6, 6)


# -

Bs = variables[0].map_blocks(
    func_Bs,
    *variables[1:],
    variable_names,
    new_axis=[4, 5],
    chunks=(1, 1, 1, variables[variable_names.index('time')].shape[-1], 6, 6),
    dtype=np.float64,
    meta=np.ndarray((1, 1, 1, variables[variable_names.index('time')].shape[-1], 6, 6), dtype=np.float64)
)
Bs

# +
# %%time

Bs.to_zarr(data_folder + filestem + output_folder + "Bs", overwrite=True)
write_to_logfile('done')
print('done')
# -

to_zarr(data_folder + filestem + output_folder + "start_values")
print('done')


# +
def write_to_logfile(*args):
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    with open(logfilename, 'a') as f:
        t = (current_time,) + args
        f.write(" ".join([str(s) for s in t]) + '\n')


# there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
    if len(groupby) == 1:
        res = dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
        return res
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
        )

def func_pwc_mr_fd(ds_single):
#    print(ds_single)
    ds_res = CARDAMOMlib.compute_ds_pwc_mr_fd_greg(ds_single, comp_dict)
    write_to_logfile(
        "finished single,",
        "lat:", ds_single.lat.data,
        "lon:", ds_single.lon.data,
        "prob:", ds_single.prob.data
    )
    
    return ds_res

def func_chunk(chunk_ds):
#    print('chunk started:', chunk_ds.lat[0].data, chunk_ds.lon[0].data, flush=True)
    res_ds = nested_groupby_apply(chunk_ds, ['lat', 'lon', 'prob'], func_pwc_mr_fd)

    # group_by removes the dimensions mentioned, so the resulting ds is
    # lower dimensional, unfortunatley, map_blocks does not do that and so
    # putting the sub result datasets back together becomes technically difficult
#    chunk_fake_ds = make_fake_ds(chunk_ds).chunk(sub_chunk_dict)
#    sub_chunk_ds = chunk_ds.chunk(sub_chunk_dict)
#    res_ds = xr.map_blocks(func_pwc_mr_fd, sub_chunk_ds, template=chunk_fake_ds)

    print(
        'chunk finished:',
        chunk_ds.lat[0].data, chunk_ds.lon[0].data, chunk_ds.prob[0].data,
        flush=True
    )
#    write_to_logfile(
#        'chunk finished,',
#        "lat:", chunk_ds.lat[0].data,
#        "lon:", chunk_ds.lon[0].data,
#        "prob:", chunk_ds.prob[0].data
#    )

    return res_ds


# -

fake_ds = make_fake_ds(ds_sub).chunk(chunk_dict)
ds_pwc_mr_fd = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)

# +
# %%time

c = ds_sub.chunks
nr_chunks = np.prod([len(val) for val in c.values()])
nr_singles = len(ds_sub.lat) * len(ds_sub.lon) * len(ds_sub.prob)
write_to_logfile(
    'starting:',
#    nr_chunks, "chunks, ",
    nr_singles, "singles"
)

#ds_pwc_mr_fd.to_netcdf(
#    data_folder + filestem + output_folder + "pwc_mr_fd_%04d.nc" % prob_nr,
#    compute=True
#)

ds_pwc_mr_fd.to_zarr(
    data_folder + filestem + output_folder + "pwc_mr_fd_%04d" % prob_nr,
    compute=True
)

write_to_logfile('done')

# +
#comp_dict = {'zlib': True, 'complevel': 9}
#encoding = {var: comp_dict for var in ds_mr_pwc_fd.data_vars}
#ds_mr_pwc_fd_computed.to_netcdf(
#    data_folder + filestem + output_folder + "pwc_mr_fd.nc",
#    encoding=encoding
#)
# -

ds.close()
del ds
ds_sub.close()
del ds_sub
ds_pwc_mr_fd.close()
del ds_pwc_mr_fd


