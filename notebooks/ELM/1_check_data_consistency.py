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

# # Check ELM (no vertical resolution) data consistency
#
# This notebook checks whether stocks and fluxes in ELM's netCDF files are consistent.

# ## Computation of data consistency

# +
# #%load_ext autoreload

# +
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from tqdm import tqdm

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from bgc_md2.models.ELM import ELMlib_no_vr
from bgc_md2.notebook_helpers import nested_groupby_apply

from dask.distributed import Client

# #%autoreload 2
# -

my_cluster = CARDAMOMlib.prepare_cluster(n_workers=24)#, my_user_name="hmetzler")
Client(my_cluster)

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
# ssh -L 8081:localhost:8790 antakya_from_home
# `
#
# and open link given above.

data_path = ELMlib_no_vr.DATA_PATH
ds = xr.open_dataset(str(data_path.joinpath("global_spinup_yearly.nc")))
ds

# +
output_path = data_path.joinpath("output")
output_path.mkdir(exist_ok=True)
output_file_path = output_path.joinpath("data_consistency.nc")

clean_file_path = output_path.joinpath("clean_ds.nc")
#ds = xr.open_dataset(clean_file_path)

time_step_in_days = 365

ds


# +
# check data consistency for single (lat, lon)
def func_data_consistency(ds_single):
#    print(ds_single)
   
    abs_err, rel_err = ELMlib_no_vr.check_data_consistency(
        ds_single,
        time_step_in_days
    )
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

# delegate chunk to single (lat, lon)
def func_chunk(chunk_ds):
#    print(chunk_ds, flush=True)
    print('chunk started:', chunk_ds.lat[0].data, chunk_ds.lon[0].data,flush=True)
    res = nested_groupby_apply(chunk_ds, ['lat', 'lon'], func_data_consistency)
    print(
        'chunk finished:',
        chunk_ds.lat[0].data,
        chunk_ds.lon[0].data,
        flush=True
    )
    
    return res


# +
chunk_dict = {"lat": 10, "lon": 10} # keep the number of chunks limited, otherwise nothing happens
ds_sub = ds.isel(
    lat=slice(0, None, 1), # 24
    lon=slice(0, None, 1), # 38
).chunk(chunk_dict)

ds_sub

# +
# create fake data with the structure of the final
# xr.map_blocks_output

fake_data = np.zeros((
    len(ds_sub.lat),
    len(ds_sub.lon)
))

fake_array = xr.DataArray(
    data=fake_data,
    dims=['lat', 'lon']
)

fake_coords = {
    'lat': ds_sub.lat.data,
    'lon': ds_sub.lon.data
}

fake_ds = xr.Dataset(
    data_vars={'abs_err': fake_array, 'rel_err': fake_array},
    coords=fake_coords
).chunk(chunk_dict)
fake_ds

# +
# create delayed data consistency dataset object

ds_data_consistency = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)

# unfortunately, xr.map_blocks ignores 'attrs', set we have to set them here manually
ds_data_consistency.data_vars['abs_err'].attrs = {'units': 'g/m^2'}
ds_data_consistency.data_vars['rel_err'].attrs = {'units': '%'}

ds_data_consistency

# +
# %%time

# compute and write to disk
ds_data_consistency.to_netcdf(
    output_file_path,
    compute=True
)
output_file_path
# -

# ## Visualization of data consistency

print(output_file_path)
ds_data_consistency = xr.open_dataset(output_file_path)
ds_data_consistency

# Now we show how many of the 96 x 144 = 13824 (lat x lon) single sites actually DO have data (presumably then land areas).

# +
# find count of nonnull coordinates (land coordinates, basically)
# (just an exercise here)

da_stacked = ds_data_consistency.abs_err.stack(notnull=['lat', 'lon'])
notnull = da_stacked[da_stacked.notnull()].notnull()
print(
    "(lat, lon) with data:",
    len(notnull),
    "out of",
    da_stacked.shape[0]
)
# -

# The data consistency describes if the stocks and fluxes do actually match.

# +
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,6))

ds_data_consistency.abs_err.plot(
    ax=ax1,
    cbar_kwargs={"label": '$%s$' % ds_data_consistency.abs_err.attrs['units']},
#    robust=True
)
ax1.set_title('mean absolute error')

ds_data_consistency.rel_err.plot(
    ax=ax2,
    cbar_kwargs={"label": '%s' % ds_data_consistency.rel_err.attrs['units']},
#    robust=True
)
ax2.set_title('mean relative error')

plt.suptitle('ELM - data consistency (non-robust version)')

plt.tight_layout()
plt.draw()
# -

# The robust version ignores outliers (2-quantile to 98-quantile only, if I remember correctly).

# +
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,6))

ds_data_consistency.abs_err.plot(
    ax=ax1,
    cbar_kwargs={"label": '$%s$' % ds_data_consistency.abs_err.attrs['units']},
    robust=True
)
ax1.set_title('mean absolute error')

ds_data_consistency.rel_err.plot(
    ax=ax2,
    cbar_kwargs={"label": '%s' % ds_data_consistency.rel_err.attrs['units']},
    robust=True
)
ax2.set_title('mean relative error')

plt.suptitle('ELM - data consistency (robust version)')

plt.tight_layout()
plt.draw()
# -
# ## Now we throw out all nonzero error sites. 


bad_idx = (ds_data_consistency.abs_err + ds_data_consistency.rel_err > 0.1)
nr_sites = len(ds.lat) * len(ds.lon)
print("Throw out %d sites of %d: %.2f%%" % (np.sum(bad_idx), nr_sites, np.sum(bad_idx)/nr_sites*100))
np.where(bad_idx)

# +
data_vars = dict()
for name, var in tqdm(ds.data_vars.items()):
    var = var.transpose("lat", "lon", "time")

    data = np.asarray(var)
    data[bad_idx] = np.nan
    data_vars[name] = xr.DataArray(
        data=data,
        dims=var.dims,
        coords=var.coords,
        attrs=var.attrs
    )
    #print(var.attrs)
    
clean_ds = xr.Dataset(data_vars=data_vars)
clean_ds

# +
# %%time

clean_ds.to_netcdf(clean_file_path)
clean_file_path
# -

