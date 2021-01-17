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

my_cluster = CARDAMOMlib.prepare_cluster(n_workers=48)
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

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
target_path = data_path.joinpath("rechunked_zarr")
ds = xr.open_mfdataset(str(data_path) + "/SUM*.nc")
ds

chunk_dict = {"lat": 1, "lon": 1, "prob": 1}
ds_rechunked = ds.chunk(chunk_dict)
ds_rechunked

# +
# overwite potentially existing zarr files?

overwrite = False # if False, raises zarr.errors.ContainsArrayError if zarr archive already exists

# +
# %%time

for name, value in tqdm(ds_rechunked.variables.items()):
    zarr_dir_path = target_path.joinpath(name)
    zarr_dir_name = str(zarr_dir_path)
    print(zarr_dir_name)
    
    if overwrite & zarr_dir_path.exists():
        print('overwrtiting...')
        shutil.rmtree(zarr_dir_path)

    da.asarray(value.data).to_zarr(zarr_dir_name)

print("done")
# -


