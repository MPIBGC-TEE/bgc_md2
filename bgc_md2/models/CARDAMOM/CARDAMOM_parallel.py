# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
from dask.distributed import Client
import xarray as xr

import importlib
import CARDAMOMlib

# xr.set_options(display_style='html')
# -


importlib.reload(CARDAMOMlib)


client = Client(n_workers=10, threads_per_worker=2, memory_limit="1GB")
client


# +
method = "discrete"
# method = 'continuous'

# data_folder = '/home/hmetzler/Desktop/CARDAMOM/' # local
data_folder = "/home/data/CARDAMOM/"  # matagorda
ds = xr.open_dataset(data_folder + "cardamom_for_holger.nc").chunk(
    {"ens": 20}
)  # , 'lat': 1, 'lon' :1})
ds


# -


## there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn):
    if len(groupby) == 1:
        return dataset.groupby(groupby[0]).apply(apply_fn)
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn
        )


# +
# %%time

# compute Delta 14C values on entire grid (ens x lat x lon) (400 x 2 x 2)


def func(sub_ds):
    res = CARDAMOMlib.load_Delta_14C_dataset(sub_ds, method=method)
    return res


def func_chunk(chunk_ds):
    if sum(chunk_ds.ens) == 0:
        res = func(chunk_ds)
        return res

    res = nested_groupby_apply(chunk_ds, ["ens", "lat", "lon"], func)
    return res


ds_Delta_14C = xr.map_blocks(func_chunk, ds).compute()
ds_Delta_14C
# -


## the ensemble numbers, latitude and longitude where the discrete reconstruction failed
df_fail = (
    ds_Delta_14C.where(ds_Delta_14C.log != "", drop=True)
    .isel(time=0)
    .to_dataframe()["log"]
)
df_fail


## keep locations with successful reconstruction
ds_Delta_14C = ds_Delta_14C.where(ds_Delta_14C.log == "", drop=True)
ds_Delta_14C


## save dataset and clean up
ds_Delta_14C.to_netcdf(data_folder + "Delta_14C_" + method + ".nc")
ds_Delta_14C.close()
ds.close()
