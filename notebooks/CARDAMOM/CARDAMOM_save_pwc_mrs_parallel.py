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
import CARDAMOMlib

import gc

xr.set_options(display_style='html')
# -


importlib.reload(CARDAMOMlib)


#client = Client(n_workers=20, threads_per_worker=2, memory_limit="4GB")
client = Client(n_workers=2, threads_per_worker=1, memory_limit="4GB")#, silence_logs='error')
client


# +
data_folder = "/home/hmetzler/Desktop/CARDAMOM/" # local
#data_folder = "/home/data/CARDAMOM/"  # matagorda

filestem = "cardamom_for_holger_10_ensembles"
#filestem = "cardamom_for_holger"
chunk_dict = {"ens": 2}
#filestem = "cardamom_for_holger"
#chunk_dict = {"ens": 100}
ds = xr.open_dataset(data_folder + filestem + ".nc")#.isel(
#    ens=slice(None, 6),
#    time=slice(None, 5)
#)
ds = ds.chunk(chunk_dict)
ds


# -


# there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
    if len(groupby) == 1:
        return dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
        )


# +
# create template for returned dataset

# first obtain the correct units
unit_ds = ds.isel(ens=0, lat=0, lon=0)
unit_mdo = CARDAMOMlib.load_mdo(unit_ds)
unit_ds.close()

data_vars = dict()

enss = ds.ens
lats = ds.lat
lons = ds.lon
times = ds.time
pools = range(6)
pools_to = pools
pools_from = pools
coords = {
    "ens": enss,
    "lat": lats,
    "lon": lons,
    "time": times,
    "pool": pools,
    "pool_to": pools,
    "pool_from": pools
}

empty_var_sv = xr.DataArray(
    data=np.ndarray(dtype=float, shape=(len(enss), len(lats), len(lons), len(pools))),
    dims=("ens", "lat", "lon", "pool"),
    attrs={"units": unit_mdo.stock_unit}
)
data_vars['start_values'] = empty_var_sv

empty_var_us = xr.DataArray(
    data=np.ndarray(dtype=float, shape=(len(enss), len(lats), len(lons), len(times), len(pools))),
    dims=("ens", "lat", "lon", "time", "pool"),
    attrs={"units": unit_mdo.stock_unit+"/"+unit_mdo.time_agg.unit}
)
data_vars["us"] = empty_var_us    

empty_var_Bs = xr.DataArray(
    data=np.ndarray(dtype=float, shape=(len(enss), len(lats), len(lons), len(times), len(pools), len(pools))),
    dims=("ens", "lat", "lon", "time", "pool_to", "pool_from"),
    attrs={"units": "1/"+unit_mdo.time_agg.unit}
)
data_vars["Bs"] = empty_var_Bs

data_vars["log"] = xr.DataArray(
    data=np.ndarray(dtype="<U120", shape=(len(enss), len(lats), len(lons))),
    dims=("ens", "lat", "lon")
)

ds_mr_template = xr.Dataset(
    coords=coords,
    data_vars=data_vars,
    attrs={"time_unit": unit_mdo.time_agg.unit}
).chunk(chunk_dict)
  


# +
# %%time

# compute in parallel the model runs and save them to ds_mrs in netcdf format

def func(single_site_ds):
    res = CARDAMOMlib.compute_pwc_mr_fd_ds(single_site_ds)
    return res

def func_chunk(chunk_ds):
    print("\nens: {:03} - {:03}\n".format(chunk_ds.ens.data[0], chunk_ds.ens.data[-1]))
    res = nested_groupby_apply(chunk_ds, ['ens', 'lat', 'lon'], func)
    print(res)
    return res

#ds_mrs = xr.map_blocks(func_chunk, ds, template=ds_mr_template).compute()
ds_mrs = xr.map_blocks(func_chunk, ds, template=ds_mr_template)
#ds_mrs = xr.map_blocks(func_chunk, ds, template=ds).compute()
#ds_mrs


# +
## save dataset and clean up

comp_dict = {'zlib': True, 'complevel': 9}
encoding = {var: comp_dict for var in ds_mrs.data_vars}
ds_mrs.to_netcdf(
    filestem + "_pwc_mrs_fd" + ".nc",
    encoding=encoding,
    compute=True
)


# +
ds_mrs.close()

ds.close()
# -


del ds_mrs
del ds

gc.collect()

# +
# gc.collect?
# -


