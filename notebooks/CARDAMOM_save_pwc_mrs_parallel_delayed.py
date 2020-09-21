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

import importlib
import bgc_md2
from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask import delayed, compute
from dask.diagnostics import Profiler, ResourceProfiler

xr.set_options(display_style='html')
# -


importlib.reload(bgc_md2.models.CARDAMOM.CARDAMOMlib)


client = Client(n_workers=10, threads_per_worker=1, memory_limit="3GB")
client


# +
pwc_mr_fd_archive = 'pwc_mr_fd_archive/'

#data_folder = "/home/hmetzler/Desktop/CARDAMOM/" # local
data_folder = "/home/data/CARDAMOM/"  # matagorda

#filestem = "cardamom_for_holger_10_ensembles"
#chunk_dict = {"ens": 2}
filestem = "cardamom_for_holger"
chunk_dict = {"ens": 10}
ds = xr.open_dataset(data_folder + filestem + ".nc")#.isel(
#    ens=slice(None, 6),
#    time=slice(None, 5)
#)
#ds = ds.chunk(chunk_dict)
#ds


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

    return chunk_ds.ens.data[0]
#    return xr.Dataset(
#        data_vars={
#            'x': xr.DataArray(
#                data=np.zeros((chunk_dict['ens'],)),
#                dims=['ens']
#            )
#        }
#    )


# -

results = []
for index in range(ds.ens.data[0], ds.ens.data[-1], chunk_dict['ens']):
    chunk_ds = ds.isel(ens=slice(index, index+chunk_dict['ens'], 1))
#    print(chunk_ds)
    results.append(delayed(func_chunk)(chunk_ds))

#prof = ResourceProfiler()
#prof.register()
delayed_results = delayed(results)
#delayed_results.visualize()

# +
# %%time

#delayed_results.compute(scheduler='processes')#, num_workers=2)
result = compute(delayed_results, scheduler='distributed', num_workers=10, memory_limit="3GB")
#result.visualize()
# -

prof.visualize()

ds.close()

del result


dask.config.get('scheduler')


