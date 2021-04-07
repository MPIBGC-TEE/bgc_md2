# -*- coding: utf-8 -*-
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

# # Compute global 14C (initial values and transient run from 1920 to 2015)

# +
import zarr

import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path
from tqdm import tqdm

from CompartmentalSystems.helpers_reservoir import F_Delta_14C
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from bgc_md2.models.CARDAMOM import CARDAMOMlib

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

time_resolution, delay_in_months, model_type = "monthly", None, "discrete"

# +
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)

CARDAMOM_path = Path("/home/data/CARDAMOM/")
#intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
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
Us_zarr = zarr.open(str(project_path.joinpath("Us")))
Bs_zarr = zarr.open(str(project_path.joinpath("Bs")))

#xs_da = da.from_zarr(str(project_path.joinpath("xs")))

# construct output shape data

nr_lats_total = len(lats_da)
nr_lons_total = len(lons_da)
nr_probs_total = len(probs_da)
nr_times_total = len(times_da)
nr_pools = 6
nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total
# -

# use all "lat", all "lon", the first four "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(0, None, 1), 
    "time": slice(0, None, 1) # don't change the time entry
}

nr_lats = len(np.arange(nr_lats_total)[slices["lat"]])
nr_lons = len(np.arange(nr_lons_total)[slices["lon"]])
nr_probs = len(np.arange(nr_probs_total)[slices["prob"]])
nr_times = len(np.arange(nr_times_total)[slices["time"]])
nr_lats, nr_lons, nr_probs, nr_times

# ## Compute 14C start_values
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

task = {
    "computation": "start_values_14C",
    "model_type": model_type,
    "overwrite": False,
    "func": CARDAMOMlib.compute_start_values_14C,
    "func_args": {"nr_time_steps": params["nr_time_steps"]}, # nr_months for fake eq_model
    "timeouts": [np.inf],
    "batch_size": 500,
    "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_pools),
    "result_chunks": (1, 1, 1, nr_pools),
    "return_shape": (1, nr_pools),
    "meta_shape": (1, nr_pools),
    "drop_axis": [1, 2, 3], # drop time axis and two pool axes of B
    "new_axis": [1] # add one pool axis
}


# +
# %%time

CARDAMOMlib.run_task_with_mr(
    project_path,
    task,
    nr_pools,
    params["time_step_in_days"],
    times_da,
    start_values_zarr,
    Us_zarr, # note capital U
    Bs_zarr,
    slices
)
# -


# ## Compute 14C solution for Southern Hemisphere (lat < -30), Tropics (-30 <= lat <= 30), and Northern Hemisphere (30 < lat)
#
# Results wille be stored in zarr archives.

task = {
    "computation": "solution_14C",
    "model_type": model_type,
    "overwrite": True,
    "func": CARDAMOMlib.compute_solution_14C,
    "func_args": {
        "nr_time_steps": params["nr_time_steps"], # nr_months for fake eq_model
    },
    "timeouts": [np.inf],
    "batch_size": 500,
    "result_shape": (nr_lats_total, nr_lons_total, nr_probs_total, nr_times_total, nr_pools),
    "result_chunks": (1, 1, 1, nr_times_total, nr_pools),
    "return_shape": (1, nr_times, nr_pools),
    "meta_shape": (1, nr_times, nr_pools),
    "drop_axis": [2, 3], # drop two pool axes of B
    "new_axis": [2] # add one pool axis
}


hemispheres = [
    {
        "name": "Southern Hemisphere",
        "Delta14C_atm_path": CARDAMOM_path.joinpath("Delta_14C_SH.csv"),
        "slices_lat": slice(0, 6, 1)
    },
    {
        "name": "Tropics",
        "Delta14C_atm_path": CARDAMOM_path.joinpath("Delta_14C_Tropics.csv"),
        "slices_lat": slice(6, 22, 1)
    },
    {
        "name": "Northern Hemisphere",
        "Delta14C_atm_path": CARDAMOM_path.joinpath("Delta_14C_NH.csv"),
        "slices_lat": slice(22, None, 1)
    },
]

# +
# %%time

for nr, d in enumerate(hemispheres):
    # prevent overwriting of results from previous hemispheres
    if nr > 0:
        task["overwrite"] = False
        
    slices_sub = slices.copy()
    slices_sub["lat"] = d["slices_lat"]

    print("Computing", d["name"])
    print("lat:", lats_da[slices_sub["lat"]].compute())
    
    Delta14C_atm_path = d["Delta14C_atm_path"]
    print(Delta14C_atm_path)
    task["func_args"]["Delta14C_atm_path"] = Delta14C_atm_path

    CARDAMOMlib.run_task_with_mr(
        project_path,
        task,
        nr_pools,
        params["time_step_in_days"],
        times_da,
        start_values_zarr,
        Us_zarr, # note capital U
        Bs_zarr,
        slices_sub
    )
# -


# ## Make netcdf files from the zarr archives

netCDF_filestem = "global_14C"

lats_da = da.from_zarr(str(project_path.joinpath("lat")))
lons_da = da.from_zarr(str(project_path.joinpath("lon")))
probs_da = da.from_zarr(str(project_path.joinpath("prob")))
times_da = da.from_zarr(str(project_path.joinpath("time")))

# +
# dimensions: lat x lon x prob x time x pool

start_values_da = da.from_zarr(str(project_path.joinpath("start_values")))
data_da = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_2")))
solution_da = data_da[:, :, :, :, 0, :]

solution_14C_da = da.from_zarr(str(project_path.joinpath("solution_14C")))
start_values_14C_da = solution_14C_da[:, :, :, 0, :]

start_values_Delta_14C_da = F_Delta_14C(start_values_da, start_values_14C_da)
system_start_values_Delta_14C_da = F_Delta_14C(start_values_da.sum(-1), start_values_14C_da.sum(-1))    

solution_Delta_14C_da = F_Delta_14C(solution_da, solution_14C_da)
system_solution_Delta_14C_da = F_Delta_14C(solution_da.sum(-1), solution_14C_da.sum(-1))

# +
coords = {
    "lat": lats_da.reshape(-1)[slices["lat"]],
    "lon": lons_da.reshape(-1)[slices["lon"]],
    "prob": probs_da.reshape(-1)[slices["prob"]],
    "time": pd.DatetimeIndex(times_da.reshape(-1).compute()),
    "pool": CARDAMOMlib.load_model_structure().pool_names.copy(),
}

data_vars = dict()

# variables with no time dimension and no pool dimension
variables = [
    {"name": "system_start_values_Delta_14C", "da": system_start_values_Delta_14C_da, "unit": "‰"}
]

for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob"],
        attrs={"units": d["unit"]}
)
    
# variables with one time dimension and no pool dimension
variables = [
    {"name": "system_solution_Delta_14C", "da": system_solution_Delta_14C_da, "unit": "‰"}
]

for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "time"],
        attrs={"units": d["unit"]}
)

# variables with no time dimension and one pool dimension
variables = [
    {"name": "start_values", "da": start_values_da, "unit": "g/m^2"},
    {"name": "start_values_14C", "da": start_values_14C_da, "unit": "g/m^2"},
    {"name": "start_values_Delta_14C_da", "da": start_values_Delta_14C_da, "unit": "‰"}
]

for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "pool"],
        attrs={"units": d["unit"]}
)

# variables with one time and one pool dimension
variables = [
    {"name": "solution", "da": solution_da, "unit": "g/m^2"},
    {"name": "solution_14C", "da": solution_14C_da, "unit": "g/m^2"},
    {"name": "solution_Delta_14C", "da": solution_Delta_14C_da, "unit": "‰"},
]

for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "time", "pool"],
        attrs={"units": d["unit"]}
)


ds = xr.Dataset(
    data_vars=data_vars,
    coords=coords,
    attrs = {
        "model_type": model_type,
        "time_resolution": time_resolution,
        "delay_in_months": delay_in_months if delay_in_months is not None else 0
    }
)
ds


# +
# compression takes forever, better zip it "manually" afterwards

def delayed_to_netcdf(prob, netCDF_filename, compute=False):
    ds_sub = ds.sel(prob=[prob])
    ds_sub.to_netcdf(
        netCDF_filename,
        compute=compute
    )
    del ds_sub

arr = ds.prob
print(arr)

# +
# %%time

for prob in tqdm(arr[25:]):
    netCDF_filename = project_path.joinpath(netCDF_filestem + "_%05d.nc" % prob)
    print(netCDF_filename)
    delayed_to_netcdf(prob, netCDF_filename, compute=True)
    
# -


