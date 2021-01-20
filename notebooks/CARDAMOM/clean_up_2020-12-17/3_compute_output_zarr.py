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

# # Compute solution, ages moments, backward transit time moments, and backward transit time quantiles
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

time_resolution = "yearly"

# +
params = CARDAMOMlib.load_params(time_resolution)

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
output_path = data_path.joinpath(params["output_folder"])

project_path = output_path.joinpath("solve_ivp_0000-0003_check_success")
# -

lats_da = da.from_zarr(str(project_path.joinpath("lat")))
lons_da = da.from_zarr(str(project_path.joinpath("lon")))
probs_da = da.from_zarr(str(project_path.joinpath("prob")))
times_da = da.from_zarr(str(project_path.joinpath("time")))

# +
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
# use all "lat", all "lon", the first four "prob", all "time"
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

# +
# compute age moment vectors of orders 0, 1, 2; oder 0 equals the solution
# use the first 120 months of the time series to average over us and Bs
# take these averages and compute a fake equilibrium model to derive
# the start age ditributions/moments from

# also compute the external outputs
# with age_moments and external outputs we can then compute backward transit time moments
task_list = [
    {#0:
        "computation": "age_moment_vectors_up_to_2",
        "overwrite": False,
        "func": CARDAMOMlib.compute_age_moment_vector_up_to,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "up_to_order": 2}, # nr_months for fake eq_model, up_to_order
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, 3, nr_pools), # solution + 2 moments
        "result_chunks": (1, 1, 1, nr_times_total, 3, nr_pools),
        "return_shape": (1, nr_times, 3, nr_pools),
        "meta_shape": (1, nr_times, 3, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": [2, 3] # add one moment axis and one pool axis
    },
    {#1:
        "computation": "external_output_vector",
        "overwrite": False,
        "func": CARDAMOMlib.compute_external_output_vector,
        "func_args": dict(),
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
        "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
        "return_shape": (1, nr_times, nr_pools),
        "meta_shape": (1, nr_times, nr_pools),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": 2 # add one pool axis
    },
    {#2:
        "computation": "btt_median",
        "overwrite": False,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.5}, # 120 months for faking equilibrium model, 0.5 for the median
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#3:
        "computation": "btt_quantile_05",
        "overwrite": True,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.05}, # 120 months for faking equilibrium model
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    },
    {#4:
        "computation": "btt_quantile_95",
        "overwrite": True,
        "func": CARDAMOMlib.compute_backward_transit_time_quantile,
        "func_args": {"nr_time_steps": params["nr_time_steps"], "q": 0.95}, # 120 months for faking equilibrium model
        "timeouts": [np.inf],
        "batch_size": 500,
        "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total),
        "result_chunks": (1, 1, 1, nr_times_total),
        "return_shape": (1, nr_times),
        "meta_shape": (1, nr_times),
        "drop_axis": [2, 3], # drop two pool axes of B
        "new_axis": []
    }
]
# -

# ## Computing
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

for task in task_list[3:]:
    CARDAMOMlib.run_task_with_pwc_mr(
        project_path,
        task,
        nr_pools,
        31.0, # days per month
        times_da,
        start_values_zarr,
        us_zarr,
        Bs_zarr,
        slices
    )


# # Plot reconstruction accuracy

# +
def get_non_nan_sites(z, drop_dims):
    fill_tup = (slice(0, 1, 1), ) * drop_dims
    tup = (slice(0, None,1),) * 3 + fill_tup
    cut_z = z[tup]
    non_nan_coords = np.where(~np.isnan(cut_z))[:3]
    
    return len(non_nan_coords[0]), non_nan_coords

def get_nan_sites(z, drop_dims):
    fill_tup = (slice(0, 1, 1), ) * drop_dims
    tup = (slice(0, None,1),) * 3 + fill_tup
    cut_z = z[tup]
    non_nan_coords = np.where(np.isnan(cut_z))[:3]
    
    return len(non_nan_coords[0]), non_nan_coords


# +
xs_zarr = zarr.open(str(project_path.joinpath("xs")))[slices["lat"], slices["lon"], slices["prob"]]
soln_zarr = zarr.open(str(project_path.joinpath("solution")))[slices["lat"], slices["lon"], slices["prob"]]

nr_sites = np.prod(xs_zarr.shape[:3])
nr_non_nan_sites_xs, _ = get_non_nan_sites(xs_zarr)
xs_pct = nr_non_nan_sites_xs/nr_sites*100
nr_non_nan_sites_soln, _ = get_non_nan_sites(soln_zarr)
soln_pct = nr_non_nan_sites_soln/nr_sites*100

print("Number of selected sites:", nr_sites)
print("Number of xs sites      :", nr_non_nan_sites_xs, "(%2.2f%%)" % xs_pct)
print("Number of solved sites  :", nr_non_nan_sites_soln, "(%2.2f%%)" % soln_pct)

# +
fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,6))

var = xs_err_ds.xs_abs_err
var.mean(dim=["time", "prob", "pool"]).plot(
    ax=ax1,
    cbar_kwargs={"label": var.attrs["units"]},
    robust=True
)
ax1.set_title("Mean absolute error")

var = xs_err_ds.xs_rel_err
var.mean(dim=["time", "prob", "pool"]).plot(
    ax=ax2,
    cbar_kwargs={"label": var.attrs["units"]},
    robust=True
)
ax1.set_title("Mean relative error")

plt.suptitle('CARDAMOM - reconstruction error (robust version)')
plt.tight_layout()
plt.draw()
# -

v1 = da.from_zarr(str(project_path.joinpath("solution")))[slices["lat"], slices["lon"], slices["prob"]].compute()
v2 = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_1")))[slices["lat"], slices["lon"], slices["prob"], :, 0].compute()
v3 = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_2")))[slices["lat"], slices["lon"], slices["prob"], :, 0].compute()
v4 = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_1")))[slices["lat"], slices["lon"], slices["prob"], :, 1].compute()
v5 = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_2")))[slices["lat"], slices["lon"], slices["prob"], :, 1].compute()

print(np.nanmax(np.abs(v1 - v2)))
print(np.nanmax(np.abs(v1 - v2)/v1*100))
print(np.nanmax(np.abs(v1 - v3)))
print(np.nanmax(np.abs(v1 - v3)/v1*100))
print(np.nanmax(np.abs(v4 - v5)))
print(np.nanmax(np.abs(v4 - v5)/v4*100))

v1 = mean_btt_da[slices["lat"], slices["lon"], slices["prob"]].compute()
v2 = da.from_zarr(str(project_path.joinpath("mean_btt")))[slices["lat"], slices["lon"], slices["prob"]].compute()

get_non_nan_sites(v1, drop_dims=1)[0], get_non_nan_sites(v2, drop_dims=1)[0]

# +
Bs_zarr = zarr.open(str(project_path.joinpath("Bs")))
Bs_sliced = Bs_zarr[slices["lat"], slices["lon"], slices["prob"]]
z = zarr.ones((nr_lats, nr_lons, nr_probs, nr_times, nr_pools))
print(get_non_nan_sites(Bs_sliced, drop_dims=3)[0])
#z[np.isnan(Bs_sliced[..., 0, 0, 0]), ...] = np.nan
print(get_nan_sites(z, drop_dims=2)[0])
# zarr.create?
# -


