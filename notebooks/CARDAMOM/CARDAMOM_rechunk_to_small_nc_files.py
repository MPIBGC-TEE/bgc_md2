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


# +
#importlib.reload(bgc_md2.models.CARDAMOM.CARDAMOMlib)


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


# -


def func_chunk(chunk_ds):
    filename = "%02d_%03.2f_%03.2f.nc" % (chunk_ds.prob, chunk_ds.lat, chunk_ds.lon)
    print(filename)
    chunk_ds.to_netcdf(data_folder + filestem + "small_netcdf/" + filename)
    return chunk_ds


# +
chunk_dict = {"lat": 1, "lon": 1, "prob": 1}
ds_sub = ds.isel(
    lat=slice(0, None, 1),
    lon=slice(0, None, 1),
    prob=slice(0, 5, 1)
).chunk(chunk_dict)

#ds_sub = ds.chunk(chunk_dict)
ds_sub

# +
# %%time

#ds_data_consistency = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)
xr.map_blocks(func_chunk, ds_sub, template=ds_sub).compute()
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
