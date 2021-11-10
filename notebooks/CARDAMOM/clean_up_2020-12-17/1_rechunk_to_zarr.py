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

# # Rechunk CARDAMOM data
#
# This notebook loads the data files coming from a CARDAMOM model run into an xarray Dataset, rechunks the data, and saves the rechunked data as zarr versions.

# +
import shutil

import dask.array as da
import xarray as xr

from pathlib import Path
from tqdm import tqdm

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(
    n_workers=12,
#    alternative_dashboard_port=8792
)
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

time_resolution, delay_in_months = "monthly", None
#time_resolution, delay_in_months = "yearly", 0
#time_resolution, delay_in_months = "yearly", 6

# +
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)

#data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
data_path = Path("/home/data/CARDAMOM/Greg_2021_10_09/")

output_path = data_path.joinpath(params["output_folder"])

target_path = data_path.joinpath(output_path.joinpath("rechunked_zarr_clean"))
if time_resolution == "yearly":
    target_path = data_path.joinpath(output_path.joinpath("rechunked_zarr_clean"))
print("target_path:", target_path)

# daily zarr archives can be created by "convert_to_daily_time_step.ipynb"
# yearly xarray.Dataset data can be created by "convert_to_yearly_time_step.ipynb"
if time_resolution == "monthly":
    #ds = xr.open_mfdataset(str(data_path) + "/SUM*.nc")
    ds = xr.open_dataset(output_path.joinpath("clean_ds.nc"))
elif time_resolution == "yearly":
#    ds = xr.open_dataset(data_path.joinpath("yearly_%02d_ds.nc" % delay_in_months))
    ds = xr.open_dataset(output_path.joinpath("clean_ds.nc"))
else:
    raise(ValueError("data can only be rechunked to monthly and yearly zarr archives by now"))
ds
# -

chunk_dict = {"lat": 1, "lon": 1, "prob": 1}
ds_rechunked = ds.chunk(chunk_dict)
ds_rechunked

# +
# overwite potentially existing zarr files?

overwrite = False # if False, raises zarr.errors.ContainsArrayError if zarr archive already exists

# +
# stupid nanny and worker messages: just wait, computations are running (see dashboard for status)

results = []
for variable_name, variable in tqdm(ds_rechunked.variables.items()):
    zarr_dir_path = target_path.joinpath(variable_name)
    zarr_dir_name = str(zarr_dir_path)
    print(zarr_dir_name)
    
    if overwrite & zarr_dir_path.exists():
        print('overwrtiting...')
        shutil.rmtree(zarr_dir_path)

    da.asarray(variable.data).to_zarr(zarr_dir_name)
# -


