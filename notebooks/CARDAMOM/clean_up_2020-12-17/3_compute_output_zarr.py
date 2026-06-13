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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Compute solution, age, and backward transit times moments and quantiles
#
# The computed data are intermediately stored in zarr archives.

# +
import zarr

import dask.array as da
import numpy as np

import matplotlib.pyplot as plt

from pathlib import Path

from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD
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

time_resolution, delay_in_months, model_type = "monthly", None, "discrete"
#time_resolution, delay_in_months, model_type = "monthly", None, "continuous"
#time_resolution, delay_in_months, model_type = "yearly", 0, "continuous"
#time_resolution, delay_in_months, model_type = "yearly", 6, "continuous"

# +
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)

#data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
data_path = Path("/home/data/CARDAMOM/Greg_2021_10_09/")

output_path = data_path.joinpath(params["output_folder"])

project_path = output_path.joinpath(model_type)
project_path
# -

lats_da = da.from_zarr(str(project_path.joinpath("lat")))
lons_da = da.from_zarr(str(project_path.joinpath("lon")))
probs_da = da.from_zarr(str(project_path.joinpath("prob")))
times_da = da.from_zarr(str(project_path.joinpath("time")))

# +
start_values_zarr = zarr.open(str(project_path.joinpath("start_values")))
if model_type == "continuous":
    us_zarr = zarr.open(str(project_path.joinpath("us")))
elif model_type == "discrete":
    Us_zarr = zarr.open(str(project_path.joinpath("Us")))
else:
    raise ValueError(f"Unknown model_type '{model_type}'")
                     
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
# use all "lat", all "lon", the first 1 "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(40, None, 1), # done (monthly): (0, 50, 1)
#    "lat": slice(28, 32, 1),
#    "lon": slice(48, 52, 1),
#    "prob": slice(0, 2, 1),
    "time": slice(0, None, 1) # don't change the time entry
}

nr_lats = len(np.arange(nr_lats_total)[slices["lat"]])
nr_lons = len(np.arange(nr_lons_total)[slices["lon"]])
nr_probs = len(np.arange(nr_probs_total)[slices["prob"]])
nr_times = len(np.arange(nr_times_total)[slices["time"]])
nr_lats, nr_lons, nr_probs, nr_times

# +
# compute age moment vectors of orders 0, 1, 2; oder 0 equals the solution
# use the first 120 months of the time series to average over us and Bs
# take these averages and compute a fake equilibrium model to derive
# the start age ditributions/moments from

# order until which pool-age moments are being computed.
# to be able to compute skweness and kurtosis of btt, the order needs to be at least 4
up_to_order = 4

# also compute the external outputs
# with age_moments and external outputs we can then compute backward transit time moments
task_list = [
#    {
#        "computation": "solve",
#        "overwrite": True,
#        "func": PWCModelRunFD.solve,
#        "func_args": dict(),
#        "timeouts": [np.inf],
#        "batch_size": 500,
#        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
#        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
#        "return_shape": (1, nr_times, nr_pools),
#        "meta_shape": (1, nr_times, nr_pools),
#        "drop_axis": [2, 3], # drop two pool axes of B
#        "new_axis": 2 # add one pool axis
#    },
#    {
#        "computation": "age_moment_vectors_up_to_1",
#        "overwrite": True,
#        "func": CARDAMOMlib.compute_age_moment_vector_up_to,
#        "func_args": {"nr_time_steps": params["nr_time_steps"], "up_to_order": 1}, # nr_months for fake eq_model, up_to_order
#        "timeouts": [np.inf],
#        "batch_size": 500,
#        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, 2, nr_pools), # solution + 1 moment
#        "result_chunks": (1, 1, 1, nr_times_total, 2, nr_pools),
#        "return_shape": (1, nr_times, 2, nr_pools),
#        "meta_shape": (1, nr_times, 2, nr_pools),
#        "drop_axis": [2, 3], # drop two pool axes of B
#        "new_axis": [2, 3] # add one moment axis and one pool axis
#    },
    {#0:
        "computation": f"age_moment_vectors_up_to_{up_to_order}",
        "overwrite": False,
        "func": CARDAMOMlib.compute_age_moment_vector_up_to,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "up_to_order": up_to_order}, # nr_months for fake eq_model, up_to_order
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, 1+up_to_order, nr_pools), # solution + up_to_order moments
        "result_chunks": (1, 1, 1, nr_times_total, 1+up_to_order, nr_pools),
        "return_shape": (1, nr_times, 1+up_to_order, nr_pools),
        "meta_shape": (1, nr_times, 1+up_to_order, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": [2, 3] # add one moment axis and one pool axis
    },
    {#1:
        "computation": "pool_age_median",
        "overwrite": False,
        "func": CARDAMOMlib.compute_pool_age_quantile,
        "func_args": {
            "nr_time_steps": params["nr_time_steps"],  # 120 months for faking equilibrium model
            "q": 0.5,
            "maxsize": 60 * nr_times
        },
        "timeouts": [150, 300, np.inf],
#        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 2 # add pool axis
    },
    {#2:
        "computation": "pool_age_quantile_05",
        "overwrite": False,
        "func": CARDAMOMlib.compute_pool_age_quantile,
        "func_args": {
            "nr_time_steps": params["nr_time_steps"],  # 120 months for faking equilibrium model
            "q": 0.05,
            "maxsize": 60 * nr_times
        },
        "timeouts": [150, 300, np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 2 # add pool axis
    },
    {#3:
        "computation": "pool_age_quantile_95",
        "overwrite": True,
        "func": CARDAMOMlib.compute_pool_age_quantile,
        "func_args": {
            "nr_time_steps": params["nr_time_steps"],  # 120 months for faking equilibrium model
            "q": 0.95,
            "maxsize": 60 * nr_times
        },
        "timeouts": [180, 300, np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 2 # add pool axis
    },
    {#4:
        "computation": "system_age_median",
        "overwrite": False,
        "func": CARDAMOMlib.compute_system_age_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.5}, # 120 months for faking equilibrium model
        "timeouts": [90, np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#5:
        "computation": "system_age_quantile_05",
        "overwrite": False,
        "func": CARDAMOMlib.compute_system_age_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.05}, # 120 months for faking equilibrium model
        "timeouts": [90,np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#6:
        "computation": "system_age_quantile_95",
        "overwrite": True,
        "func": CARDAMOMlib.compute_system_age_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.95}, # 120 months for faking equilibrium model
        "timeouts": [90, np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#7:
        "computation": "external_output_vector",
        "overwrite": True,
        "func": CARDAMOMlib.compute_external_output_vector,
        "func_args": dict(),
        "timeouts": [np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 2 # add one pool axis
    },
    {#8:
        "computation": "btt_median",
        "overwrite": True,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.5}, # 120 months for faking equilibrium model, 0.5 for the median
        "timeouts": [120, np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#9:
        "computation": "btt_quantile_05",
        "overwrite": True,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.05}, # 120 months for faking equilibrium model
        "timeouts": [120, np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#10:
        "computation": "btt_quantile_95",
        "overwrite": True,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.95}, # 120 months for faking equilibrium model
        "timeouts": [120, np.inf],
        "batch_size": 120,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    }
]

for task in task_list:
    task["model_type"] = model_type
# -

# ## Computing
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

# +
# %%time

for task in task_list[:1]:
    CARDAMOMlib.run_task_with_mr(
        project_path,
        task,
        nr_pools,
        params["time_step_in_days"],
        times_da,
        start_values_zarr,
        us_zarr if model_type == "continuous" else Us_zarr,
        Bs_zarr,
        slices
    )
# +
# done
# -


