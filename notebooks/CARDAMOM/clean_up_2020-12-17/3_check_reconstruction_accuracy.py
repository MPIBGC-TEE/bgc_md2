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

# # Check reconstruction accuracy
#
# In this notebook we check how accurate the reconstruction of `start_values`, `us`, and `Bs` is. We do so by building the `PWCModelRunFD` (PieceWise Constant Model Run) for each (lat, lon, prob) and check the accuracy of the solution against the `xs` coming from the netCDF data (now saved in zarr format).
#
# Afterwards we kick out the sites that have a bad reconstruction accuracy and keep only the good ones. Unfortunately, this decision is somehow arbitrary. It would be way better if the sites were all reconstructed quite well, this might be the case for data with less coarse time step (here we have 31 days).

# +
import zarr

import dask.array as da
import numpy as np

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib
from bgc_md2.notebook_helpers import (
    write_to_logfile,
    write_header_to_logfile,
    load_zarr_archive,
    custom_timeout
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

project_path = output_path.joinpath("solve_ivp_0000-0003")
#clean_path = project_path.joinpath("clean")

# +
lats_da = da.from_zarr(str(project_path.joinpath("lat")))
lons_da = da.from_zarr(str(project_path.joinpath("lon")))
probs_da = da.from_zarr(str(project_path.joinpath("prob")))
times_da = da.from_zarr(str(project_path.joinpath("time")))

start_values_zarr = zarr.open(str(project_path.joinpath("start_values")))
us_zarr = zarr.open(str(project_path.joinpath("us")))
Bs_zarr = zarr.open(str(project_path.joinpath("Bs")))

#xs_da = da.from_zarr(str(project_path.joinpath("xs")))

# construct output shape data

nr_lats_total = len(lats_da)
nr_lons_total = len(lons_da)
nr_probs_total = len(probs_da)
nr_times_total = len(times_da)
nr_pools = 6
nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total

# +
# use all "lat", all "lon", the first two "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(0, 4, 1),
    "time": slice(0, None, 1) # don't change the time entry
}

nr_lats = len(np.arange(nr_lats_total)[slices["lat"]])
nr_lons = len(np.arange(nr_lons_total)[slices["lon"]])
nr_probs = len(np.arange(nr_probs_total)[slices["prob"]])
nr_times = len(np.arange(nr_times_total)[slices["time"]])
nr_lats, nr_lons, nr_probs, nr_times
# -

task_list = [
    {# 0:
        "computation": "solution",
        "overwrite": True,
        "func": PWCModelRunFD.solve,
        "func_args": [],
        "timeouts": [np.inf],
        "batch_size": 10,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 1, # add one pool axis for solution
    }
]

# ## Computing
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

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

    nr_incomplete_sites, _ = CARDAMOMlib.get_incomplete_site_tuples_for_pwc_mr_computation(
        start_values_zarr,
        us_zarr,
        Bs_zarr,
        z,
        slices
    )
    print("Number of incomplete sites:", nr_incomplete_sites)

    logfile_name = str(project_path.joinpath(task["computation"] + ".log"))
    print("Logfile:", logfile_name)

    for timeout in task["timeouts"]:        
        CARDAMOMlib.compute_incomplete_sites_with_pwc_mr(
            timeout,
            z,
            nr_pools,
            31.0, # days_per_month
            times_da,
            start_values_zarr,
            us_zarr,
            Bs_zarr,
            slices,
            task,
            logfile_name
        )

    nr_incomplete_sites, _ = CARDAMOMlib.get_incomplete_site_tuples_for_pwc_mr_computation(
        start_values_zarr,
        us_zarr,
        Bs_zarr,
        z,
        slices
    )
    write_to_logfile(logfile_name, nr_incomplete_sites, "incomplete sites remaining")
    print(nr_incomplete_sites, "incomplete sites remaining")
    print()



