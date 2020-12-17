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

# # Compute CARDAMOM start_values, us, and Bs
#
# This notebook computes the start_values, us, and Bs from the rechunked zarr archive (Step 1), and saves them as zarr archives. One can set a timeout to prevent compuations to get stuck.

# +
import shutil

import dask.array as da

from pathlib import Path

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
# ssh -L 8080:localhost:8895 antakya_from_home
# `
#
# and open link given above.

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
zarr_data_path = data_path.joinpath("rechunked_zarr")
output_path = data_path.joinpath("output")

project_name = "solve_ivp_0000-0001"

# +
# load variables from zarr archive

variable_paths = [p for p in zarr_data_path.iterdir() if p.is_dir()]

non_data_vars = ["lat", "lon", "prob", "time"]

variable_names = []
variables_total = []
for variable_path in variable_paths:
    name = variable_path.name
    if name not in non_data_vars:
        variable_names.append(name)
        variables_total.append(da.from_zarr(str(variable_path)))

# manual reshaping for "lat", "lon", "time", "prob"
# necessary for map_blocks
for k, name in enumerate(non_data_vars):
    shape = [1] * 4
    shape[k] = -1
    variables_total.append(da.from_zarr(str(zarr_data_path.joinpath(name)), chunks=(1,)).reshape(shape))
    variable_names.append(name)

variables_total

# +
# use all "lat", all "lon", the first two "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(0, 2, 1),
    "time": slice(0, None, 1)
}

# subdivide "lon" into 20 batches in order to use less memory at once
# "time" cannot be subdivided because model runs need all available times
batch_nrs = {
    "lat": 1,
    "lon": 10,
    "prob": 1
}

# +
# cut out the slices of interest from the variables

variables = []
for nr, name, var in zip(range(len(variables_total)), variable_names, variables_total):
    if name not in non_data_vars:
        variables.append(variables_total[nr][slices["lat"], slices["lon"], slices["prob"], ...])

for k, name in enumerate(non_data_vars):
    index = variable_names.index(name)
    v = variables_total[index]
    s = [slice(0, None, 1),] * v.ndim
    s[k] = slices[name]
    variables.append(v[tuple(s)])
    
variables
# -


