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

import numpy as np

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib
from bgc_md2.notebook_helpers import (
    write_to_logfile,
    write_header_to_logfile,
    create_zarr_archive
)

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

project_path = output_path.joinpath("solve_ivp_0000-0001")

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
# -

# We decide in which values of which dimensions we are interested (maybe to save computation time).

# use all "lat", all "lon", the first two "prob", all "time"
slices = {
    "lat": slice(10, 20, 1),
    "lon": slice(10, 20, 1),
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
# -

task_list = [
    {# 0:
        "computation": "start_values",
        "func": CARDAMOMlib.compute_start_values,
        "timeouts": [np.inf],
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_pools), # chunk over lat, lon, prob
        "return_shape": (1, 1, 1, nr_pools), # return from "func" in map_blocks
        "meta_shape": (nr_lats, nr_lons, nr_probs, nr_pools), # "meta" argument in map_blocks
        "drop_axis": 3, # remove time axis
        "new_axis": 3 # add pool axis
    },
    {# 1:
        "computation": "us",
        "func": CARDAMOMlib.compute_us,
        "timeouts": [np.inf],
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, 1, 1, nr_times, nr_pools),
        "meta_shape": (nr_lats, nr_lons, nr_probs, nr_times, nr_pools),
        "drop_axis": 3, # drop time axis
        "new_axis": [3, 4] # add time and pool axes
    },
    {# 2:
        "computation": "Bs",
        "func": CARDAMOMlib.compute_Bs,
        "timeouts": [0.2, 0.5, 15],
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools, nr_pools),
        "return_shape": (1, 1, 1, nr_times, nr_pools, nr_pools),
        "meta_shape": (nr_lats, nr_lons, nr_probs, nr_times, nr_pools, nr_pools),
        "drop_axis": 3, # drop time axis
        "new_axis": [3, 4, 5] # add time axis and two pool axes
    }
]

task = task_list[2]

# ## Computing

print("task: computing", task["computation"])

logfile_name = str(project_path.joinpath(task["computation"] + ".log"))
time_limit_in_min = task["timeouts"][0]
print("logfile:", logfile_name)
print("Timeout (min):", time_limit_in_min)

# The logfile can be viewed with vim, the line numbers tell the progress. Use
#
# `:set autoread | au CursorHold * checktime | call feedkeys("lh")`
#
# to keep the file automatically updated in vim. 
# In general, all output to the logfile is appended to the file. To clean the file, delete it first.
#
# All other output produced during a parallel computation goes to the jupyter lab server shell.

# +
# ATTENTION: overwrites existing zarr folder

zarr_path = Path(project_path.joinpath(task["computation"]))
res_zarr = create_zarr_archive(zarr_path, task["result_shape"], task["result_chunks"], overwrite=True)
print(zarr_path)
res_zarr
# -

additional_params = {
    "func": task["func"],
    "variable_names": variable_names,
    "time_limit_in_min": time_limit_in_min,
    "return_shape": task["return_shape"],
    "logfile_name": logfile_name
}

# +
# %%time

res_da = variables[0].map_blocks(
    CARDAMOMlib.func_for_map_blocks,
    *variables[1:], # variables[0] comes automatically as first argument
    additional_params,
    drop_axis=task["drop_axis"],
    new_axis=task["new_axis"],
    chunks=task["return_shape"],
    dtype=np.float64,
    meta=np.ndarray(task["meta_shape"], dtype=np.float64)
)
res_da
# -

print(write_header_to_logfile(logfile_name, res_da, time_limit_in_min))

# +
# %%time

# subdivide "lon" into 10 batches in order to use less memory at once
# "time" cannot be subdivided because model runs need all available times
batch_nrs = {
    "lat": 1,
    "lon": 10,
    "prob": 1
}

CARDAMOMlib.batchwise_to_zarr(res_da, res_zarr, slices, batch_nrs, noisy=True)
write_to_logfile(logfile_name, 'done')
print('done')
# -

# ### Timeout handling
#
# If some (lat, lon, prob)-computations timed out and the `task["timeouts"]` list has more than one element, we collect the so-calles "timeouts" and run them again with the next timeout value.

# +
# %%time

#import importlib
#importlib.reload(CARDAMOMlib)

zarr_path = Path(project_path.joinpath(task["computation"]))
nr_timeout_sites, timeout_coords = CARDAMOMlib.get_timeout_sites(zarr_path, slices)

if nr_timeout_sites == 0:
    s = "done, no timeouts remaining"
    write_to_logfile(
        logfile_name,
        s
    )
    print(s)
else:
    timeouts = task["timeouts"]
    if len(timeouts) == 1:
        s = "done, " + str(nr_timeouts) + " timeouts remaining"
        write_to_logfile(
            logfile_name,
            s
        )
        print(s)
    else:
        for k in range(1, len(timeouts)):
            print('remaining timeouts =', timeouts[k:], 'min\n')
            time_limit_in_min = timeouts[k]
            CARDAMOMlib.remove_timeouts(
                time_limit_in_min,
                zarr_path,
                variable_names,
                variables,
                non_data_variables,
                slices,
                task,
                additional_params
            )
            print()        

# +



