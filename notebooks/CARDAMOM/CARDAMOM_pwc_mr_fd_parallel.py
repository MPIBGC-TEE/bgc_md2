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
from dask.distributed import Client

import xarray as xr
import numpy as np

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client, LocalCluster
from getpass import getuser

import time


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

# dashboard needs a different port for accessing it remotely
my_dashboard_port = my_port + 5
my_cluster = LocalCluster(
    dashboard_address='localhost:'+str(my_dashboard_port),
    n_workers=80,
    threads_per_worker=1    
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
# ssh -L 8080:localhost:8890 antakya_from_home
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
logfilename = data_folder + filestem + output_folder + "pwc_mr_fd.log"

ds = xr.open_mfdataset(data_folder + filestem + "SUM*.nc")
ds


# +
ms = CARDAMOMlib.load_model_structure_greg()

def make_fake_ds(dataset):
    # fake start_values data
    fake_data_sv = np.zeros((
        len(dataset.lat),
        len(dataset.lon),
        len(dataset.prob),
        ms.nr_pools
    ))

    coords_pool = [d['pool_name'] for d in ms.pool_structure]
    fake_coords_sv = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data,
        'pool': coords_pool
    }

    fake_array_sv = xr.DataArray(
        data=fake_data_sv,
        dims=['lat', 'lon', 'prob', 'pool'],
        coords=fake_coords_sv
    )

    # fake times data
    fake_data_times = np.zeros((
        len(dataset.lat),
        len(dataset.lon),
        len(dataset.prob),
        len(dataset.time)
    ))

    fake_coords_times = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data,
        'time': dataset.time.data
    }

    fake_array_times = xr.DataArray(
        data=fake_data_times,
        dims=['lat', 'lon', 'prob', 'time'],
        coords=fake_coords_times
    )

    # fake us data
    fake_data_us = np.zeros((
        len(dataset.lat),
        len(dataset.lon),
        len(dataset.prob),
        len(dataset.time),
        ms.nr_pools
    ))

    fake_coords_us = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data,
        'time': dataset.time.data,
        'pool': coords_pool
    }

    fake_array_us = xr.DataArray(
        data=fake_data_us,
        dims=['lat', 'lon', 'prob', 'time', 'pool'],
        coords=fake_coords_us
    )

    # fake Bs data
    fake_data_Bs = np.zeros((
        len(dataset.lat),
        len(dataset.lon),
        len(dataset.prob),
        len(dataset.time),
        ms.nr_pools,
        ms.nr_pools
    ))

    fake_coords_Bs = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data,
        'time': dataset.time.data,
        'pool_to': coords_pool,
        'pool_from': coords_pool
    }

    fake_array_Bs = xr.DataArray(
        data=fake_data_Bs,
        dims=['lat', 'lon', 'prob', 'time', 'pool_to', 'pool_from'],
        coords=fake_coords_Bs
    )

    # fake log data
    shape = (
        len(dataset.lat),
        len(dataset.lon),
     len(dataset.prob),
    )
    fake_data_log = np.ndarray(shape, dtype="<U150")

    fake_coords_log = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data
    }

    fake_array_log = xr.DataArray(
        data=fake_data_log,
        dims=['lat', 'lon', 'prob'],
        coords=fake_coords_log
    )

    # collect fake arrays in ds
    fake_data_vars = dict()
    fake_data_vars['start_values'] = fake_array_sv
    fake_data_vars['times'] = fake_array_times
    fake_data_vars['us'] = fake_array_us
    fake_data_vars['Bs'] = fake_array_Bs
    fake_data_vars['log'] = fake_array_log

    fake_coords = {
        'lat': dataset.lat.data,
        'lon': dataset.lon.data,
        'prob': dataset.prob.data,
        'time': dataset.time.data,
        'pool': coords_pool,
        'pool_to': coords_pool,
        'pool_from': coords_pool
    }

    fake_ds = xr.Dataset(
        data_vars=fake_data_vars,
        coords=fake_coords
    )

    return fake_ds


# +
chunk_dict = {"lat": 1, "lon": 1, 'prob': 5}
#sub_chunk_dict = {'lat': 1, 'lon': 1, 'prob': 1}

#ds_sub = ds.isel(
#    lat=slice(0, 2, 1),
#    lon=slice(0, 2, 1),
#    prob=slice(0, 2, 1)
#).chunk(chunk_dict)

#ds_sub = ds.isel(
#    lat=slice(28, 30, 1),
#    lon=slice(38, 40, 1),
#    prob=slice(0, 20, 1)
#).chunk(chunk_dict)

ds_sub = ds.chunk(chunk_dict)
ds_sub


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
    ds_res = CARDAMOMlib.compute_ds_pwc_mr_fd_greg(ds_single)
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
    write_to_logfile(
        'chunk finished,',
        "lat:", chunk_ds.lat[0].data,
        "lon:", chunk_ds.lon[0].data,
        "prob:", chunk_ds.prob[0].data
    )

    return res_ds


# -

fake_ds = make_fake_ds(ds_sub).chunk(chunk_dict)
ds_pwc_mr_fd = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)

# +
# %%time

c = ds_sub.chunks
nr_chunks = np.prod([len(val) for val in c.values()])
nr_singles = len(ds_sub.lat) * len(ds_sub.lon) * len(ds_sub.prob)
write_to_logfile('starting:', nr_chunks, "chunks, ", nr_singles, "singles")

comp_dict = {'zlib': True, 'complevel': 9}
ds_pwc_mr_fd.to_netcdf(
    data_folder + filestem + output_folder + "pwc_mr_fd.nc",
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
