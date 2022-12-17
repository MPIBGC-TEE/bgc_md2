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

import xarray as xr
import numpy as np

import importlib
import bgc_md2
from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask import delayed, compute
from dask.diagnostics import Profiler, ResourceProfiler
from dask.distributed import Client, LocalCluster
from getpass import getuser
# -


importlib.reload(bgc_md2.models.CARDAMOM.CARDAMOMlib)


# +
#client = Client(n_workers=10, threads_per_worker=1, memory_limit="3GB")
#client


# + endofcell="--"
port_dict = {
    'mm':8789,
    'hmetzler':8790, # change at will
#    'hmetzler':8888, # change at will
    'cs':8791        # change at will
}
my_user_name = getuser()
print(my_user_name)

my_port = port_dict[my_user_name]
print(my_port)

# dasboard needs a different port for accessing it remotely
my_cluster = LocalCluster(dashboard_address='localhost:'+str(my_port+5))

# -

Client(my_cluster)
# --

# +
pwc_mr_fd_archive = '/home/data/CARDAMOM/output/pwc_mr_fd_archive/'

#data_folder = "/home/hmetzler/Desktop/CARDAMOM/" # local
#data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya, raka..
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya, raka..

#filestem = "cardamom_for_holger_10_ensembles"
#chunk_dict = {"ens": 2}
#filestem = "cardamom_for_holger"
filestem = "Greg_2020_10_26"
chunk_dict = {"lat": 1, "lon": 1}
ds = xr.open_mfdataset(data_folder + filestem + "/SUM*.nc")#.isel(
#    ens=slice(None, 6),
#    time=slice(None, 5)
#)
ds = ds.chunk(chunk_dict)
ds


# -


comp_dict = {'zlib': True, 'complevel': 9}

# +
# compute in parallel the model runs and save them to ds_mrs in netcdf format

small_template = xr.Dataset(
    data_vars = {
        'x': xr.DataArray(
#            data=np.ndarray(dtype=float, shape=(len(ds.ens.data),)),
            data=np.ndarray(dtype=float, shape=(chunk_dict['ens'],)),
            dims=['ens']
        )
    }
)#.chunk(chunk_dict)

def func(single_site_ds):
#    print(single_site_ds)
    res = CARDAMOMlib.compute_pwc_mr_fd_ds(single_site_ds)
    return res

# there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
    if len(groupby) == 1:
        res = dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
#        print(res.start_values)
        return res
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
        )
    
def func_chunk(chunk_ds):
    print('\nStarting chunk:', chunk_ds.ens.data[0], '-', chunk_ds.ens.data[-1], '\n')
    res = nested_groupby_apply(chunk_ds, ['ens', 'lat', 'lon'], func)
    
#    filename = 'output/' + filestem + "_{:03d}-{:03d}.nc".format(res.ens.data[0], res.ens.data[-1])
    filename = pwc_mr_fd_archive + filestem + "_{:03d}-{:03d}.nc".format(res.ens.data[0], res.ens.data[-1])
    encoding = {var: comp_dict for var in res.data_vars}
    res.to_netcdf(filename, encoding=encoding)
    print('\nFinished chunk:', chunk_ds.ens.data[0], '-', chunk_ds.ens.data[-1], '\n')
#    del res
#    gc.collect()

#    return chunk_ds.ens.data[0]
#    return xr.Dataset(
#        data_vars={
#            'x': xr.DataArray(
#                data=np.zeros((chunk_dict['ens'],)),
#                dims=['ens']
#            )
#        }
#    )


# +
#results = []
#for index in range(ds.ens.data[0], ds.ens.data[-1], chunk_dict['ens']):
#    chunk_ds = ds.isel(ens=slice(index, index+chunk_dict['ens'], 1))
#    results.append(delayed(func_chunk)(chunk_ds))
# -

def func_data_consistency(ds_single):
    mdo = CARDAMOMlib.load_mdo_greg(ds_single)
    abs_err, rel_err = mdo.check_data_consistency()
    res = {
        'lat': ds_single.lat,
        'lon': ds_single.lon,
        'prob': ds_single.prob,
        'abs_err': abs_err,
        'rel_err': rel_err
    }
    ds_single.close()
    
    return res


# +
# %%time

results = []
max_nr = 5
for lat_nr, lat in enumerate(ds.lat[:max_nr]):
#    if (lat_nr % 5) == 0: print(lat)
    for lon_nr, lon in enumerate(ds.lon[:max_nr]):
#        if (lon_nr % 5 == 0): print(lon)
        for prob in ds.prob:
            ds_single = ds.sel({'lat': lat, 'lon': lon, 'prob': prob})
            results.append(delayed(func_data_consistency(ds_single)))
# -

import pickle
filename = data_folder + filestem + "results.pickle"
with open(filename, 'wb') as f:
    pickle.dump(results, f)

# +
#with open(filename,'rb') as f:
#    results_loaded = pickle.load(f)
#results_loaded
# -

#prof = ResourceProfiler()
#prof.register()
delayed_results = delayed(results)
#delayed_results.visualize() #needs graphviz

# +
# %%time

#delayed_results.compute(scheduler='processes')#, num_workers=2)
result = compute(
    delayed_results,
    scheduler='distributed',
#    num_workers=12,
#    memory_limit="3GB"
)
#result.visualize()

# +
lats = []
lons= []
probs = []
abs_errs = []
rel_errs = []
for d in result[0]:
    lats.append(d['lat'])
    lons.append(d['lon'])
    probs.append(d['prob'])
    abs_errs.append(d['abs_err'])
    rel_errs.append(d['rel_err'])
    
coords = [
    np.array([l.data for l in lats]).reshape(5, 5, 50)[0],
    np.array([l.data for l in lons]).reshape(5, 5, 50)[1],
    np.array([p.data for p in probs]).reshape(5, 5, 50)[2]
]
#print(coords[0])
data_vars = dict()

#print(coords)
#print(lats)

data_vars['abs_err'] = xr.DataArray(
    data=(['lat', 'lon', 'prob'], [
        np.array([e.data for e in abs_errs]).reshape(5, 5, 50)[0],
        np.array([e.data for e in abs_errs]).reshape(5, 5, 50)[1],
        np.array([e.data for e in abs_errs]).reshape(5, 5, 50)[2],
    ]),
#    coords=coords,
    attrs={'units': abs_errs[0].unit}
)
data_vars['rel_err'] = xr.DataArray(
    data=np.array([e.data for e in rel_errs]).reshape(5, 5, 50),
    coords=coords,
    attrs={'units': rel_errs[0].unit}
)

ds_data_consistency = xr.Dataset(
    data_vars=data_vars
)
# -

ds_data_consistency


