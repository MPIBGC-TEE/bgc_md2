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
#import xarray as xr
import numpy as np
import zarr
import dask
import shutil
from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client, LocalCluster, get_worker
import dask.array as da
from getpass import getuser

import time
from bgc_md2.notebookHelpers import write_to_logfile, custom_timeout


# +
port_dict = {
    'mm': 8789,
    'hmetzler': 8790, # change at will
#    'hmetzler': 8791, # change at will
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

dask.config.set({'distributed.worker.daemon': False})

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
zarr_data_folder = "zarr_version"
output_folder = "output/"

#project_name = "Bs_experimental"
project_name = "solve_ivp_0000-0001_root_success_check"
logfilename = data_folder + filestem + output_folder + project_name + ".log"
zarr_dir_name = data_folder + filestem + output_folder + project_name + "/Bs"

zarr_path = Path(data_folder).joinpath(filestem).joinpath(zarr_data_folder)
variable_paths = [p for p in zarr_path.iterdir() if p.is_dir()]

non_data_vars = ['lat', 'lon', 'prob', 'time']

variable_names = []
variables_total = []
for variable_path in variable_paths:
    name = variable_path.name
    if name not in non_data_vars:
        variable_names.append(name)
        variables_total.append(da.from_zarr(str(variable_path)))
        
variables_total.append(da.from_zarr(str(zarr_path.joinpath('lat')), chunks=(1,)).reshape(-1, 1, 1, 1))
variable_names.append('lat')
variables_total.append(da.from_zarr(str(zarr_path.joinpath('lon')), chunks=(1,)).reshape(1, -1, 1, 1))
variable_names.append('lon')
variables_total.append(da.from_zarr(str(zarr_path.joinpath('prob')), chunks=(1,)).reshape(1, 1, -1, 1))
variable_names.append('prob')
variables_total.append(da.from_zarr(str(zarr_path.joinpath('time'))).reshape(1, 1, 1, -1))
variable_names.append('time')

variables_total


# -


v = variables_total[0]
nr_times = v.shape[3]
nr_pools = 6



# +
# ATTENTION: overwrites existing zarr folder
rm = True

dir_p = Path(zarr_dir_name)
if rm & dir_p.exists():
    shutil.rmtree(dir_p)
    
result_shape = (v.shape[0], v.shape[1], v.shape[2], nr_times, nr_pools, nr_pools)
result_chunks = (1, 1, 1, -1, -1, -1)
z = zarr.create(
    result_shape,
    chunks=result_chunks,
    dtype=np.float64,
    fill_value=-np.inf, # -np.inf indicates incomplete computation
    store=zarr_dir_name
)
z.shape
# -



# +
timeouts = [30, 300, 2000]
dims = {"lat": 0, "lon": 1, "prob": 2, "time": 3}

slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(0, 2, 1)
}

batch_nrs = {
    "lat": 1,
    "lon": 10,
    "prob": 1
}

variables = []
for nr, name, var in zip(range(len(variables_total)), variable_names, variables_total):
    if name not in non_data_vars:
        variables.append(variables_total[nr][slices["lat"], slices["lon"], slices["prob"], ...])
    
name = "lat"
index = variable_names.index(name)
variables.append(variables_total[index][slices[name], ...])

name = "lon"
index = variable_names.index(name)
variables.append(variables_total[index][:, slices[name], ...])

name = "prob"
index = variable_names.index(name)
variables.append(variables_total[index][:, :, slices[name], ...])

name = "time"
index = variable_names.index(name)
variables.append(variables_total[index])

variables


# -

def func_Bs(*args):
    v_names = args[-3]
    time_limit_in_min = args[-2]
    return_shape = args[-1]
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-3]))}

    Bs = -np.inf * np.ones(return_shape)
    start_time = time.time()
    timeout_msg = ''
    try:
        Bs = custom_timeout(
            time_limit_in_min*60,
            CARDAMOMlib.load_Bs_greg_dict,
            d,
#            integration_method="trapezoidal",
#            nr_nodes=1001
        )
    except TimeoutError:
        duration = (time.time() - start_time) / 60
        timeout_msg = 'Timeout after %2.2f min' % duration
        print(timeout_msg)

    write_to_logfile(
        logfilename,
        "finished single,",
        "lat:", d["lat"],
        "lon:", d["lon"],
        "prob:", d["prob"],
        timeout_msg
    )

    #return Bs.reshape(1, len(d['time']), 6, 6)
    return Bs.reshape(return_shape)


# +
return_shape = (1, 1, 1, nr_times, nr_pools, nr_pools)
timeout = timeouts[0]

Bs = variables[0].map_blocks(
    func_Bs,
    *variables[1:],
    variable_names,
    timeout, # time_limit_in_min
    return_shape,
    new_axis=[4, 5],
    chunks=return_shape,
    dtype=np.float64,
    meta=np.ndarray(return_shape, dtype=np.float64)
)
Bs


# -

def batchwise_to_zarr_CARDAMOM(
    arr,
    z,
    z_slices
):   
    def split_arr_to_slices(ns):
        res = dict()
        for name in ["lat", "lon", "prob"]:
            a = np.arange(0, arr.shape[dims[name]])
            l = np.array_split(a, ns[name])
            res[name] = [slice(x[0], x[-1]+1, 1) for x in l]
            
        return res
    
    arr_slices_dict = split_arr_to_slices(batch_nrs)
       
    arr_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in arr_slices_dict["lat"] 
        for slice_lon in arr_slices_dict["lon"]
        for slice_prob in arr_slices_dict["prob"]
    ]

    def split_slices(slices, ns):
        res = dict()
        for name, s in slices.items():
            a = np.arange(0, z.shape[dims[name]], 1)[s]
            l = np.array_split(a, ns[name])
            res[name] = [slice(x[0], x[-1]+1, s.step) for x in l]
        
        return res
          
    z_slices_dict = split_slices(z_slices, batch_nrs)
    
    z_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in z_slices_dict["lat"] 
        for slice_lon in z_slices_dict["lon"]
        for slice_prob in z_slices_dict["prob"]
    ]
   
    for s_arr, s_z in zip(arr_slices_tuples, z_slices_tuples):
        z[s_z[0], s_z[1], s_z[2], ...] = arr[s_arr[0], s_arr[1], s_arr[2], ...].compute()


# write header to logfile
c = Bs.chunks
nr_chunks = np.prod([len(val) for val in c])
print('nr_chunks:', nr_chunks)
nr_singles = np.prod(Bs.shape[:3])
print('nr_singles:', nr_singles)
write_to_logfile(
    logfilename,
    'starting',
#    nr_chunks, "chunks, ",
    nr_singles, "singles,",
    "timeout =",
    timeout,
    "min"
)
print('timeout:', timeout, ' min')

# +
# %%time

batchwise_to_zarr_CARDAMOM(Bs, z, slices)
write_to_logfile(logfilename, 'done')
print('done')
# -



def remove_timeouts(timeout):
    z = zarr.open(zarr_dir_name)
    da_Bs_restricted = da.from_zarr(z)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]

    timeout_coords = np.where(da_Bs_restricted.compute()==-np.inf)
    timeout_nrs = len(timeout_coords[0])
    if timeout_nrs > 0:
        print('number of timeouts:', timeout_nrs)    
    else:
        print('timeouts successfully removed')
        return

    timeout_variables = []
    for v, name in zip(variables, variable_names):
        if name not in non_data_vars:
            v_stack_list = []
            for lat, lon, prob in zip(*[timeout_coords[k] for k in range(3)]):
                v_stack_list.append(v[lat, lon, prob])
            
            timeout_variables.append(da.stack(v_stack_list))
    
    timeout_variables.append(da.from_array(timeout_coords[0].reshape(-1, 1), chunks=(1, 1))) # lat
    timeout_variables.append(da.from_array(timeout_coords[1].reshape(-1, 1), chunks=(1, 1))) # lon
    timeout_variables.append(da.from_array(timeout_coords[2].reshape(-1, 1), chunks=(1, 1))) # prob
    timeout_variables.append(variables[variable_names.index('time')].reshape(1, -1)) # time

    return_shape = (1, nr_times, nr_pools, nr_pools)
    timeout_Bs = timeout_variables[0].map_blocks(
        func_Bs,
        *timeout_variables[1:],
        variable_names,
        timeout, # time_limit_in_min
        return_shape,
        new_axis=[2, 3],
        chunks=return_shape,
        dtype=np.float64,
        meta=np.ndarray(return_shape, dtype=np.float64)
    )

    # write header to logfile
    c = timeout_Bs.chunks
    nr_chunks = np.prod([len(val) for val in c])
    print('nr_chunks:', nr_chunks)
    nr_singles = timeout_Bs.shape[0]
    print('nr_singles:', nr_singles)
    write_to_logfile(
        logfilename,
        'starting',
#        nr_chunks, "chunks, ",
        nr_singles, "timeout singles,",
        "timeout =",
        timeout,
        "min"
    )
    print('timeout =', timeout, 'min')
    
    # do the computation
    timeout_Bs_computed = timeout_Bs.compute()
    
    # compute the z coords from the sliced Bs coords 
    f_lat = lambda x: slices['lat'].start + x * slices['lat'].step
    f_lon = lambda x: slices['lon'].start + x * slices['lon'].step
    f_prob = lambda x: slices['prob'].start + x * slices['prob'].step
    timeout_coords_z_lat = np.array([f_lat(x) for x in timeout_coords[0]])
    timeout_coords_z_lon = np.array([f_lon(x) for x in timeout_coords[1]])
    timeout_coords_z_prob = np.array([f_prob(x) for x in timeout_coords[2]])

    timeout_coords_z = (timeout_coords_z_lat, timeout_coords_z_lon, timeout_coords_z_prob)
    # update zarr file
    for nr, lat, lon, prob in zip(range(timeout_nrs), *[timeout_coords_z[k] for k in range(3)]):
        z[lat, lon, prob, ...] = timeout_Bs_computed[nr, ...]

    write_to_logfile(logfilename, 'done, timeout =', timeout, 'min')
    print('done, timeout = ', timeout, 'min')


timeouts

# +
# %%time

print('timeouts =', timeouts[1:], 'min\n')
for timeout in timeouts[1:]:
    remove_timeouts(timeout)
    print()
# -
a = da.from_zarr(zarr_dir_name)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]
a

good_indices = np.where(~np.isnan(a.compute()))
len(good_indices[0])

134/(4828*0.37) * 100




timeouts

# %%time
timeout = timeouts[1]
remove_timeouts(timeout)



timeout = timeouts[3]
timeout

# +
z = zarr.open(zarr_dir_name)
da_Bs_restricted = da.from_zarr(z)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]

timeout_coords = np.where(da_Bs_restricted.compute()==-np.inf)
print(timeout_coords)
timeout_nrs = len(timeout_coords[0])
if timeout_nrs > 0:
    print('number of timeouts:', timeout_nrs)    
#    else:
#        print('timeouts successfully removed')
#        return

# +
timeout_variables = []
for v, name in zip(variables, variable_names):
    if name not in non_data_vars:
        v_stack_list = []
        for lat, lon, prob in zip(*[timeout_coords[k] for k in range(3)]):
            v_stack_list.append(v[lat, lon, prob])
            
        timeout_variables.append(da.stack(v_stack_list))
    
timeout_variables.append(da.from_array(timeout_coords[0].reshape(-1, 1), chunks=(1, 1))) # lat
timeout_variables.append(da.from_array(timeout_coords[1].reshape(-1, 1), chunks=(1, 1))) # lon
timeout_variables.append(da.from_array(timeout_coords[2].reshape(-1, 1), chunks=(1, 1))) # prob
timeout_variables.append(variables[variable_names.index('time')].reshape(1, -1)) # time
print('timeout_variables')
timeout_variables
    
# -

return_shape = (1, nr_times, nr_pools, nr_pools)
timeout_Bs = timeout_variables[0].map_blocks(
    func_Bs,
    *timeout_variables[1:],
    variable_names,
    timeout, # time_limit_in_min
    return_shape,
    new_axis=[2, 3],
    chunks=return_shape,
    dtype=np.float64,
    meta=np.ndarray(return_shape, dtype=np.float64)
)
print('current timeout dataset')
timeout_Bs

# write header to logfile
c = timeout_Bs.chunks
nr_chunks = np.prod([len(val) for val in c])
print('nr_chunks:', nr_chunks)
nr_singles = timeout_Bs.shape[0]
print('nr_singles:', nr_singles)
write_to_logfile(
    logfilename,
    'starting',
#    nr_chunks, "chunks, ",
    nr_singles, "timeout singles,",
    "timeout =",
    timeout,
    "min"
)
print('timeout =', timeout, 'min')


# do the computation
timeout_Bs_computed = timeout_Bs.compute()









len(np.where(timeout_Bs_computed[:, 0, 0, 0] == -np.inf)[0])

f_lat = lambda x: slices['lat'].start + x * slices['lat'].step
f_lon = lambda x: slices['lon'].start + x * slices['lon'].step
f_prob = lambda x: slices['prob'].start + x * slices['prob'].step


timeout_coords_z_lat = np.array([f_lat(x) for x in timeout_coords[0]])
timeout_coords_z_lon = np.array([f_lon(x) for x in timeout_coords[1]])
timeout_coords_z_prob = np.array([f_prob(x) for x in timeout_coords[2]])

timeout_coords_z = (timeout_coords_z_lat, timeout_coords_z_lon, timeout_coords_z_prob)
timeout_coords_z

# +
# update zarr file
for nr, lat, lon, prob in zip(range(timeout_nrs), *[timeout_coords_z[k] for k in range(3)]):
    z[lat, lon, prob, ...] = timeout_Bs_computed[nr, ...]

write_to_logfile(logfilename, 'done, timeout =', timeout, 'min')
print('done, timeout = ', timeout, 'min')
# -
write_to_logfile(logfilename, 'done, timeout =', timeout, 'min')
print('done, timeout = ', timeout, 'min')


















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
#    print(start_values, flush=True)

#    start_values = da.from_array(1.1 * np.ones((1, 1, 1, 6), dtype=np.float64, chunks=(1,1,6))
#    start_values.to_dask()
    return start_values.reshape(1, 1, 1, nr_pools)


# +
#shape = (34, 71, 50, 6)                                                                                         
#chunks = (1, 1, 1, 6)                   
#
## set up zarr array to store data
#store = zarr.DirectoryStore('TB1.zarr')
#root = zarr.group(store) 
#TB1 = root.zeros(
#    'data', 
#    shape=shape, 
#    chunks=chunks, 
#    dtype=np.float64
#)
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
    kwargs={'chunks': (1, 1, 1, nr_pools)}
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

    return us.reshape(1, 1, 1, nr_times, nr_pools)


us = variables[0].map_blocks(
    func_us,
    *variables[1:],
    variable_names,
    new_axis=4,
    chunks=(1, 1, 1, nr_times, nr_pools),
    dtype=np.float64,
    meta=np.ndarray((1, 1, 1, nr_times, nr_pools), dtype=np.float64)
)
us

# +
# %%time

us.to_zarr(data_folder + filestem + output_folder + "us")


# -

def func_Bs(*args):
    v_names = args[-3]
    time_limit_in_min = args[-2]
    return_shape = args[-1]
    d = {v_names[i]: args[i].reshape((-1,)) for i in range(len(args[:-3]))}

    Bs = -np.inf * np.ones(return_shape)
    start_time = time.time()
    try:
        Bs = custom_timeout(
            time_limit_in_min*60,
            CARDAMOMlib.load_Bs_greg_dict,
            d,
            integration_method="trapezoidal",
            nr_nodes=11
        )
    except TimeoutError:
        duration = (time.time() - start_time) / 60
        msg = 'Timeout after %2.2f minutes' % duration
        print(msg)

    write_to_logfile(
        logfilename,
        "finished single,",
        "lat:", d["lat"],
        "lon:", d["lon"],
        "prob:", d["prob"]
    )

    #return Bs.reshape(1, len(d['time']), 6, 6)
    return Bs.reshape(return_shape)


# +
return_shape = (1, 1, 1, nr_times, nr_pools, nr_pools)

Bs = variables[0].map_blocks(
    func_Bs,
    *variables[1:],
    variable_names,
    1, # time_limit_in_min
    return_shape,
    new_axis=[4, 5],
    chunks=return_shape,
    dtype=np.float64,
    meta=np.ndarray(return_shape, dtype=np.float64)
)
Bs
# -



# +
zarr_dir_name = data_folder + filestem + output_folder + "Bs_experimental"
# ATTENTION: overwrites existing zarr folder
rm = True

dir_p = Path(zarr_dir_name)
if rm & dir_p.exists():
    shutil.rmtree(dir_p)

z = zarr.create(
    (34, 71, 50, 1152, 6, 6),
    chunks=(1, 1, 1, 1152, 6, 6),
    dtype=np.float64,
    fill_value=-np.inf,
    store=zarr_dir_name
)
z.shape


# -

def batchwise_to_zarr_CARDAMOM(
    arr,
    z,
    z_slices
):   
    dims = {"lat": 0, "lon": 1, "prob": 2}
    batch_nrs = {
        "lat": 1,
        "lon": 10,
        "prob": 1
    }

    def split_arr_to_slices(ns):
        res = dict()
        for name in ["lat", "lon", "prob"]:
            a = np.arange(0, arr.shape[dims[name]])
            l = np.array_split(a, ns[name])
            res[name] = [slice(x[0], x[-1]+1, 1) for x in l]
            
        return res
    
    arr_slices_dict = split_arr_to_slices(batch_nrs)
       
    arr_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in arr_slices_dict["lat"] 
        for slice_lon in arr_slices_dict["lon"]
        for slice_prob in arr_slices_dict["prob"]
    ]

    def split_slices(slices, ns):
        res = dict()
        for name, s in slices.items():
            a = np.arange(0, z.shape[dims[name]], 1)[s]
            l = np.array_split(a, ns[name])
            res[name] = [slice(x[0], x[-1]+1, s.step) for x in l]
        
        return res
          
    z_slices_dict = split_slices(z_slices, batch_nrs)
    
    z_slices_tuples = [
        (slice_lat, slice_lon, slice_prob)
        for slice_lat in z_slices_dict["lat"] 
        for slice_lon in z_slices_dict["lon"]
        for slice_prob in z_slices_dict["prob"]
    ]
   
    for s_arr, s_z in zip(arr_slices_tuples, z_slices_tuples):
        z[s_z[0], s_z[1], s_z[2], ...] = arr[s_arr[0], s_arr[1], s_arr[2], ...].compute()


# write header to logfile
c = Bs.chunks
nr_chunks = np.prod([len(val) for val in c])
print('nr_chunks:', nr_chunks)
nr_singles = np.prod(Bs.shape[:3])
print('nr_singles:', nr_singles)
write_to_logfile(
    logfilename,
    'starting:',
#    nr_chunks, "chunks, ",
    nr_singles, "singles"
)

# +
# %%time

batchwise_to_zarr_CARDAMOM(Bs, z, slices)
write_to_logfile(logfilename, 'done')
print('done')
# -

Bs_zarr = zarr.open(data_folder + filestem + output_folder + "Bs")
da_Bs_restricted = da.from_zarr(Bs_zarr)[:, :, :, 0, 0, 0]
da_Bs_restricted

# +
timeout_coords = np.where(da_Bs_restricted.compute()==-np.inf)
timeout_nrs = len(timeout_coords[0])

timeout_coords

# +
timeout_variables = []
for v, name in zip(variables, variable_names):
    if name not in non_data_vars:
        v_stack_list = []
        for lat, lon, prob in zip(*[timeout_coords[k] for k in range(3)]):
            v_stack_list.append(v[lat, lon, prob])
            
        timeout_variables.append(da.stack(v_stack_list))

timeout_variables.append(da.from_array(timeout_coords[0].reshape(-1, 1), chunks=(1, 1))) # lat
timeout_variables.append(da.from_array(timeout_coords[1].reshape(-1, 1), chunks=(1, 1))) # lon
timeout_variables.append(da.from_array(timeout_coords[2].reshape(-1, 1), chunks=(1, 1))) # prob

timeout_variables.append(variables[variable_names.index('time')].reshape(1, -1))
timeout_variables

# +
return_shape = (1, nr_times, nr_pools, nr_pools)

timeout_Bs = timeout_variables[0].map_blocks(
    func_Bs,
    *timeout_variables[1:],
    variable_names,
    600, # time_limit_in_min
    return_shape,
    new_axis=[2, 3],
    chunks=return_shape,
    dtype=np.float64,
    meta=np.ndarray(return_shape, dtype=np.float64)
)
timeout_Bs
# -

# write header to logfile
c = timeout_Bs.chunks
nr_chunks = np.prod([len(val) for val in c])
print('nr_chunks:', nr_chunks)
nr_singles = timeout_Bs.shape[0]
print('nr_singles:', nr_singles)
write_to_logfile(
    logfilename,
    'starting:',
#    nr_chunks, "chunks, ",
    nr_singles, "timeout singles"
)

# +
# %%time

timeout_Bs_computed = timeout_Bs.compute()
write_to_logfile('done')
# -

timeout_nrs = len(timeout_coords[0])
for nr, lat, lon, prob in zip(range(timeout_nrs), *[timeout_coords[k] for k in range(3)]):
    Bs_zarr[lat, lon, prob, ...] = timeout_Bs_computed[nr, ...]

da_Bs_done_restricted = da.from_zarr(Bs_zarr)[:, :, :, 0, 0, 0]
timeout_coords = np.where(da_Bs_done_restricted.compute()==-np.inf)
timeout_coords

len(timeout_coords[0])







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
# -














# +
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya
filestem = "Greg_2020_10_26/"
zarr_data_folder = "zarr_version"
output_folder = "output/"

project_name = "Bs_experimental"
zarr_dir_name = data_folder + filestem + output_folder + project_name
slices = {
    "lat": slice(20, 30, 1),
    "lon": slice(30, 60, 1),
    "prob": slice(2, 4, 1)
}

# +
Bs_zarr = zarr.open(zarr_dir_name)
Bs_zarr[:] = 5
da_Bs_restricted = da.from_zarr(Bs_zarr)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]
timeout_coords = np.where(da_Bs_restricted.compute()==5)
print(len(timeout_coords[0]))

Bs_zarr[20, 30, 2] = 7

Bs_zarr2 = zarr.open(zarr_dir_name)
da_Bs_restricted2 = da.from_zarr(Bs_zarr2)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]
timeout_coords2 = np.where(da_Bs_restricted2.compute()==5)
len(timeout_coords2[0])
# -

Bs_zarr = zarr.open(zarr_dir_name)

#Bs_zarr[:] = 5
da_Bs_restricted = da.from_zarr(Bs_zarr)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]

timeout_coords = np.where(da_Bs_restricted.compute()==5)
len(timeout_coords[0])

slices

Bs_zarr[20, 30, 2] = 5



Bs_zarr2 = zarr.open(zarr_dir_name)

da_Bs_restricted2 = da.from_zarr(Bs_zarr2)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]

timeout_coords2 = np.where(da_Bs_restricted2.compute()==5)
len(timeout_coords2[0])











da_Bs_restricted2 = da.from_zarr(Bs_zarr)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]
timeout_coords2 = np.where(da_Bs_restricted2.compute()==5)
len(timeout_coords2[0])

a, b, c = timeout_coords[0][7], timeout_coords[1][7], timeout_coords[2][7]

# +
print(a)
print(b)
print(c)

Bs_zarr[a, b, c, ...] = 7.32
# -

da_Bs_restricted2 = da.from_zarr(Bs_zarr)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]
timeout_coords2 = np.where(da_Bs_restricted2.compute()==-np.inf)
len(timeout_coords2[0])

# +
Bs_zarr2 = zarr.open(zarr_dir_name)
da_Bs_restricted2 = da.from_zarr(Bs_zarr2)[slices['lat'], slices['lon'], slices['prob'], 0, 0, 0]

timeout_coords2 = np.where(da_Bs_restricted2.compute()==-np.inf)
# -

len(timeout_coords2[0])

len(timeout_coords[0])


