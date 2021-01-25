# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

# # Convert CARDAMOM zarr data and to daily time steps
#
# This notebook loads the zarr data files and saves the rechunked data as zarr versions in daily time steps. The flux values are just repeated 31 timper per month since they are given as daily average fluxes. Stock values for days between two firsts of the months re simply linearly interpolated.

# +
import shutil
import zarr

import dask.array as da
import numpy as np
import xarray as xr

from dask import delayed
from pathlib import Path
from tqdm import tqdm

from bgc_md2.notebook_helpers import load_zarr_archive
from bgc_md2.models.CARDAMOM import CARDAMOMlib
from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(
    n_workers=48,
    alternative_dashboard_port=8791
)
Client(my_cluster)

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
source_path = data_path.joinpath("monthly_rechunked_zarr")
target_path = data_path.joinpath("daily_rechunked_zarr")
ds = xr.open_mfdataset(str(data_path) + "/SUM*.nc")
ds

stock_variable_names = ["c_finelitter", "c_foliar", "c_labile", "c_root", "c_som", "c_wood"]
days_per_month = 31


def compute_target(variable_name, z_target_sliced, z_source_sliced):
    if variable_name in stock_variable_names:
        z_target_sliced[..., -1] = z_source_sliced[..., -1]

        diff = np.diff(z_source_sliced, axis=-1)

        daily_diff = diff.repeat(days_per_month, axis=-1).reshape(diff.shape + (days_per_month,))
        daily_diff = np.arange(days_per_month).reshape(1, 1, 1, 1, -1)/days_per_month * daily_diff
        daily_x = np.repeat(z_source_sliced[..., :-1], days_per_month, axis=-1)
        res = daily_x + daily_diff.reshape(daily_x.shape)
        z_target_sliced[..., :-1] = res
    else:
        z_target_sliced[..., 1] = np.nan
        z_target_sliced[..., 1:] = np.repeat(z_source_sliced[..., 1:], days_per_month, axis=-1)
    
    return z_target_sliced


# convert variable to daily data
def convert_variable(variable_name, variable):   
    source_zarr_path = source_path.joinpath(variable_name)
    if not source_zarr_path.exists():
        raise(OSError("source zarr archive" + str(source_zarr_path) + "not found"))
        
    z_source = zarr.open(str(source_zarr_path))

    target_zarr_path = target_path.joinpath(variable_name)
    print("starting", target_zarr_path, flush=True)
    
    nr_daily_times = (z_source.shape[-1]-1) * days_per_month + 1
    target_shape = z_source.shape[:-1] + (nr_daily_times,)
    target_chunks = (1, 1, 1, -1)
    z_target = load_zarr_archive(
        target_zarr_path,
        target_shape,
        target_chunks,
        overwrite=True
    )
    
    lat_subs = np.array_split(np.arange(z_source.shape[0]), 5)
    lon_subs = np.array_split(np.arange(z_source.shape[1]), 5)
    prob_subs = np.array_split(np.arange(z_source.shape[2]), 5)
    
    coord_tuples = [(lat, lon, prob) for lat in lat_subs for lon in lon_subs for prob in prob_subs]
    for coord_tuple in tqdm(coord_tuples):
        lat, lon, prob = coord_tuple
        s0 = slice(lat[0], lat[-1]+1, 1)
        s1 = slice(lon[0], lon[-1]+1, 1)
        s2 = slice(prob[0], prob[-1]+1, 1)
        
        z_source_sliced = z_source[s0, s1, s2]
        z_target_sliced = z_target[s0, s1, s2]
        z_target[s0, s1, s2] = compute_target(variable_name, z_target_sliced, z_source_sliced)
        
        del z_source_sliced
        del z_target_sliced
    print("done", target_zarr_path, flush=True)
    return 1


# +
results = []
for variable_name, variable in ds.data_vars.items():
    y = delayed(convert_variable)(variable_name, variable)
    results.append(y)                                         

total = delayed(sum)(results)
total.visualize()
# +
# %%time

print(total.compute(), "data variables converted")
# -

for variable_name in ["lat", "lon", "prob"]:
    variable = ds[variable_name]
    zarr_path = target_path.joinpath(variable_name)
    if zarr_path.exists():
        shutil.rmtree(zarr_path)
    da.asarray(variable.data).to_zarr(str(zarr_path))

variable_name = "time"
time = ds[variable_name]
time_in_days = np.array(time[0], dtype="datetime64[D]") + np.arange((len(time)-1) * days_per_month + 1)
time_in_days.shape

zarr_path = target_path.joinpath(variable_name)
if zarr_path.exists():
    shutil.rmtree(zarr_path)
da.from_array(time_in_days).to_zarr(str(zarr_path))
