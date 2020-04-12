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
import xarray as xr

import importlib
import CARDAMOMlib
# -

importlib.reload(CARDAMOMlib)

#data_folder = '/home/hmetzler/Desktop/CARDAMOM/' # local
data_folder = '/home/data/CARDAMOM/'
ds = xr.open_dataset(data_folder + 'cardamom_for_holger.nc')
ds


## there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn):
    if len(groupby) == 1:
        return dataset.groupby(groupby[0]).apply(apply_fn)
    else:
        return dataset.groupby(groupby[0]).apply(nested_groupby_apply, groupby = groupby[1:], apply_fn = apply_fn)


# +
# %%time

## compute Delta 14C values on entire grid (ens x lat x lon) (400 x 2 x 2)
def func(sub_ds):
    res = CARDAMOMlib.load_dmr_14C_dataset(sub_ds)
    return res
    
ds_dmr_14C = nested_groupby_apply(ds,['ens', 'lat', 'lon'], func)
# -

## the ensemble numbers, latitude and longitude where the discrete reconstruction failed
df_fail = ds_dmr_14C.where(ds_dmr_14C.log!='', drop=True).isel(time=0).to_dataframe()['log']
df_fail

## keep locations with successful reconstruction
ds_dmr_14C = ds_dmr_14C.where(ds_dmr_14C.log=='', drop=True)
ds_dmr_14C

## save dataset and clean up
ds_dmr_14C.to_netcdf(data_folder + 'dmr_Delta_14C.nc')
ds_dmr_14C.close()
ds.close()


