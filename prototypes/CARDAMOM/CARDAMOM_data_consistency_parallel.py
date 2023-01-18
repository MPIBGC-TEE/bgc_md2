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
from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client, LocalCluster
from getpass import getuser
# -


importlib.reload(bgc_md2.models.CARDAMOM.CARDAMOMlib)


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
my_dashboard_port = my_port +5
my_cluster = LocalCluster(
    dashboard_address='localhost:'+str(my_dashboard_port),
    n_workers=48,
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
# jupyter lab --no-browser -- port=8890
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
# ssh -L 8080:localhost:8895 antakya_from_home
# `
#
# and open link given above.

# +
data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya
output_folder = "output/"

filestem = "Greg_2020_10_26/"
ds = xr.open_mfdataset(data_folder + filestem + "SUM*.nc")
ds


# +
# there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
    if len(groupby) == 1:
        res = dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
        return res
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
        )

def func_data_consistency(ds_single):
#    print(ds_single)
    mdo = CARDAMOMlib.load_mdo_greg(ds_single)
    abs_err, rel_err = mdo.check_data_consistency()
    
    data_vars = dict()
    data_vars['abs_err'] = xr.DataArray(
        data=abs_err.data.filled(fill_value=np.nan),
        attrs={'units': abs_err.unit} # ignored by map_blocks
    )
    data_vars['rel_err'] = xr.DataArray(
        data=rel_err.data.filled(fill_value=np.nan),
        attrs={'units': rel_err.unit} # ignored by map_blocks
    )

    ds_res = xr.Dataset(
        data_vars=data_vars
    )
    ds_single.close()
    
    return ds_res

def func_chunk(chunk_ds):
#    print('chunk started:', chunk_ds.lat[0].data, chunk_ds.lon[0].data,flush=True)
    res = nested_groupby_apply(chunk_ds, ['lat', 'lon', 'prob'], func_data_consistency)
    print('chunk finished:', chunk_ds.lat[0].data, chunk_ds.lon[0].data,flush=True)
    return res


# +
chunk_dict = {"lat": 1, "lon": 1}
#ds_sub = ds.isel(
#    lat=slice(0, -1, 1),
#    lon=slice(0, -1, 1),
#    prob=slice(0, -1, 1)
#).chunk(chunk_dict)

ds_sub = ds.chunk(chunk_dict)
ds_sub

# +
fake_data = np.zeros((
    len(ds_sub.lat),
    len(ds_sub.lon),
    len(ds_sub.prob)
))

fake_array = xr.DataArray(
    data=fake_data,
    dims=['lat', 'lon', 'prob']
)

fake_coords = {
    'lat': ds_sub.lat.data,
    'lon': ds_sub.lon.data,
    'prob': ds_sub.prob.data
}

fake_ds = xr.Dataset(
    data_vars={'abs_err': fake_array, 'rel_err': fake_array},
    coords=fake_coords
).chunk(chunk_dict)
fake_ds

# +
# %%time

ds_data_consistency = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)
# -

comp_dict = {'zlib': True, 'complevel': 9}
ds_data_consistency.to_netcdf(
    data_folder + filestem + output_folder + "data_consistency.nc",
    compression=comp_dict,
    compute=True
)

ds.close()
del ds
ds_sub.close()
del ds_sub
ds_data_consistency.close()
del ds_data_consistency
