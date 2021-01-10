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
from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD, PWCModelRunFDError
from sympy import symbols

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
start_values_da = da.from_zarr(start_values_zarr)

us_zarr = zarr.open(str(project_path.joinpath("us")))
us_da = da.from_zarr(us_zarr)

Bs_zarr = zarr.open(str(project_path.joinpath("Bs")))
Bs_da = da.from_zarr(Bs_zarr)

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

# +
# compute the solutions/stocks for the model runs

overwrite = False
computation = "solution"
zarr_path = Path(project_path.joinpath(computation))
print("zarr archive:", str(zarr_path))

result_shape = (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools)
result_chunks = (1, 1, 1, nr_times, nr_pools)

z = load_zarr_archive(
    zarr_path,
    result_shape,
    result_chunks,
    overwrite=overwrite
)
z.shape
# -

_, complete_coords_svs = CARDAMOMlib.get_complete_sites(start_values_zarr, slices)
_, complete_coords_us = CARDAMOMlib.get_complete_sites(us_zarr, slices)
_, complete_coords_Bs = CARDAMOMlib.get_complete_sites(Bs_zarr, slices)
_, incomplete_coords_soln = CARDAMOMlib.get_incomplete_sites(z, slices)

# +
coords_list = []
for coords in [
    complete_coords_svs,
    complete_coords_us,
    complete_coords_Bs,
    incomplete_coords_soln
]:
    coords_list.append(
        set(
            (coords[0][i], coords[1][i], coords[2][i]) for i in range(len(coords[0]))
        )
    )

incomplete_coords = list(set.intersection(*coords_list))
print("Incomplete sites:", len(incomplete_coords))


# -

def compute_with_pwc_mr(*args):
    additional_params = args[-1]
    computation = additional_params["computation"]
    nr_pools = additional_params["nr_pools"]
    days_per_month = additional_params["days_per_month"]
    return_shape = additional_params["return_shape"]
    func = additional_params["func"]
    func_args = additional_params["func_args"]
    time_limit_in_min = additional_params["time_limit_in_min"]
    logfile_name = additional_params["logfile_name"]
    
    Bs = args[0].reshape((-1, nr_pools, nr_pools))
    start_values = args[1].reshape(nr_pools)
    us = args[2].reshape((-1, nr_pools))
    
    lat = args[3].reshape(-1)
    lon = args[4].reshape(-1)
    prob = args[5].reshape(-1)
    times = args[6].reshape(-1)
    
    nr_times = len(times)
    data_times = np.arange(0, nr_times, 1) * days_per_month

    log_msg = ""
    res = -np.inf * np.ones(return_shape)
    try:
        time_symbol = symbols('t')
        # using the custom_timeout function (even with timeout switched off) 
        # makes the worker report to the scheduler
        # and prevents frequent timeouts
        mr = custom_timeout(
            np.inf,
            PWCModelRunFD.from_Bs_and_us,
            time_symbol,
            data_times,
            start_values,
            Bs[:-1],
            us[:-1]
        )
        
        print('Computing', computation, flush=True)
#        soln = custom_timeout(
#            time_limit_in_min*60,
#                mr.solve
#            )
#            print("solved", flush=True)
#            res = soln
        res = custom_timeout(
            time_limit_in_min*60,
            func,
            mr,
            *func_args
        )           
        print("done", flush=True)
    except Exception as e:
        res = np.nan * np.ones_like(res)
        print(str(e), flush=True)
        log_msg = str(e)
            
    write_to_logfile(
        logfile_name,
        "single finished,",
        "lat:", lat,
        "lon:", lon,
        "prob:", prob,
        log_msg
    )

    return res.reshape(return_shape)


# +
incomplete_variables = []

shapes = [
    (nr_times, nr_pools, nr_pools),
    (1, nr_pools, 1),
    (nr_times, nr_pools, 1)
]
for v, shape in zip([Bs_da, start_values_da, us_da], shapes):
    v_stack_list = []
    for ic in incomplete_coords[:100]:
        v_stack_list.append(v[ic].reshape(shape))
    
    incomplete_variables.append(da.stack(v_stack_list))

for k, name in enumerate(["lat", "lon", "prob"]):
    incomplete_variables.append(
        da.from_array(
            np.array([ic[k] for ic in incomplete_coords[:100]]).reshape(-1, 1, 1, 1),
            chunks=(1, 1, 1, 1)
        )
    )
incomplete_variables.append(times_da.reshape((1, -1, 1, 1)).rechunk((1, nr_times, 1, 1)))
incomplete_variables
# -

computation = "solution"
additional_params = {
    "computation": computation,
    "func": PWCModelRunFD.solve,
    "nr_pools": nr_pools,
    "days_per_month": 31.0,
    "return_shape": (1, nr_times, nr_pools),
    "func_args": [],
    "time_limit_in_min": np.inf,
    "logfile_name": str(project_path.joinpath(computation + ".log"))
}

res_da = incomplete_variables[0].map_blocks(
    compute_with_pwc_mr,
    *incomplete_variables[1:],
    additional_params,
    drop_axis=(2, 3),
    new_axis=2,
    chunks=additional_params["return_shape"],
    dtype=np.float64,
    meta=np.ndarray((1, nr_times, nr_pools), dtype=np.float64)
)
res_da

incomplete_coords_linear = (
    [ic[0] for ic in incomplete_coords[:100]],
    [ic[1] for ic in incomplete_coords[:100]],
    [ic[2] for ic in incomplete_coords[:100]],
)

# +
time_limit_in_min = additional_params["time_limit_in_min"]

# write header to logfile
print(write_header_to_logfile(logfile_name, res_da, time_limit_in_min))
print("starting", flush=True)

# do the computation
CARDAMOMlib.linear_batchwise_to_zarr(
    res_da, # dask array
    z, # target zarr archive
    slices, # slices of interest,
    incomplete_coords_linear,
    10 # batch_size
)

write_to_logfile(logfile_name, 'done, timeout (min) =', time_limit_in_min)
print('done, timeout (min) = ', time_limit_in_min, flush=True)
# -


