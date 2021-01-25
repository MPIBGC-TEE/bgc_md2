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
import zarr
import shutil

import numpy as np

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib
from bgc_md2.notebook_helpers import (
    write_to_logfile,
    load_zarr_archive
)

from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(
    n_workers=48,
    alternative_dashboard_port=8791
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

time_resolution = "daily"
model_type = "discrete"

# +
params = CARDAMOMlib.load_params(time_resolution)

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
zarr_data_path = data_path.joinpath(time_resolution + "_rechunked_zarr")
output_path = data_path.joinpath(params["output_folder"])

project_path = output_path.joinpath(model_type)

# +
# load variables from zarr archive
non_data_variables = ["lat", "lon", "prob", "time"]

variable_names, variables_total = CARDAMOMlib.load_variables_from_zarr(
    zarr_data_path,
    non_data_variables
)
for name, var in zip(variable_names, variables_total):
    print(name)
    print(var)
    print()

# +
# write all lat, lon, prob, time to output folder

for name in ["lat", "lon", "prob", "time"]:
    idx = variable_names.index(name)
    var = variables_total[idx].reshape(-1)
    print(name)
    print(var)

    file_path = project_path.joinpath(name)
    print(file_path)
    print()
    if file_path.exists():
        shutil.rmtree(file_path)
        
    var.to_zarr(str(file_path))
# -

# We decide in which values of which dimensions we are interested (maybe to save computation time).

# use all "lat", all "lon", the first 1 "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(0, 1, 1),
    "time": slice(0, None, 1) # don't change the time entry
}

# +
# cut out the slices of interest from the variables

variables = CARDAMOMlib.select_slices_of_interest(
    variable_names,
    variables_total,
    slices,
    non_data_variables
)
for name, var in zip(variable_names, variables):
    print(name)
    print(var)
    print()

# +
# construct output shape data

# data for total zarr output
v = variables_total[0]
nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total = v.shape

# data for slices dask.array output of interest
v = variables[0]
nr_lats, nr_lons, nr_probs, nr_times = v.shape

# global constant
nr_pools = 6

nr_lats, nr_lons, nr_probs, nr_times, nr_pools
# -

task_list = [
    {# 0:
        "computation": "xs",
        "overwrite": False,
        "func": CARDAMOMlib.compute_xs,
        "func_args": {"time_step_in_days": params["time_step_in_days"]},
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": 1, # drop time axis
        "new_axis": [1, 2] # add time and pool axes
    },    
    {# 1:
        "computation": "start_values",
        "overwrite": False,
        "func": CARDAMOMlib.compute_start_values,
        "func_args": {"time_step_in_days": params["time_step_in_days"]},
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_pools), # chunk over lat, lon, prob
        "return_shape": (1, nr_pools), # return from "func" in map_blocks
        "meta_shape": (1, nr_pools), # "meta" argument in map_blocks, "1" will be replaced by nr of incomplete sites
        "drop_axis": 1, # remove time axis
        "new_axis": 1 # add pool axis
    },
    {# 2:
        "computation": "Us",
        "overwrite": False,
        "func": CARDAMOMlib.compute_Us, # note the capital U in Us (aggregated over time step)
        "func_args": {"time_step_in_days": params["time_step_in_days"]},
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": 1, # drop time axis
        "new_axis": [1, 2] # add time and pool axes
    },
    {# 3:
        "computation": "Bs",
        "overwrite": False,
        "func": CARDAMOMlib.compute_Bs_discrete,
        "func_args": {"time_step_in_days": params["time_step_in_days"]},
#        "timeouts": [30, 300, 2000],
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools, nr_pools),
        "return_shape": (1, nr_times, nr_pools, nr_pools),
        "meta_shape": (1, nr_times, nr_pools, nr_pools),
        "drop_axis": 1, # drop time axis
        "new_axis": [1, 2, 3] # add time axis and two pool axes
    }
]

# ## Computing
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

# +
# %%time

for task in task_list:
    print("task: computing", task["computation"])
    print()
    
    zarr_path = Path(project_path.joinpath(task["computation"]))
    print("zarr archive:", str(zarr_path))
    z = load_zarr_archive(
        zarr_path,
        task["result_shape"],
        task["result_chunks"],
        overwrite=task["overwrite"]
    )

    nr_incomplete_sites, incomplete_coords = CARDAMOMlib.get_incomplete_sites(z, slices)
    print("Number of incomplete sites:", nr_incomplete_sites)
    logfile_name = str(project_path.joinpath(task["computation"] + ".log"))
    print("Logfile:", logfile_name)

    for timeout in task["timeouts"]:
        CARDAMOMlib.compute_incomplete_sites(
            timeout,
            z,
            nr_times,
            variable_names,
            variables,
            non_data_variables,
            slices,
            task,
            logfile_name
        )

    nr_incomplete_sites, _ = CARDAMOMlib.get_incomplete_sites(z, slices)
    write_to_logfile(logfile_name, nr_incomplete_sites, "incomplete sites remaining")
    print(nr_incomplete_sites, "incomplete sites remaining")
    print()
# -



