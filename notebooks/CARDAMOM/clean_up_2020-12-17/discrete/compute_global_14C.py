# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
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
#    n_workers=48,
    n_workers=10,
#    alternative_dashboard_port=8791,
    my_user_name="hmetzler"
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

#CARDAMOM_path = Path("/home/data/CARDAMOM/")
CARDAMOM_path = Path("/mnt/c/Users/hrme0001/data/CARDAMOM/")
intcal20_path = CARDAMOM_path.joinpath("IntCal20_Year_Delta14C.csv")
#data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
data_path = CARDAMOM_path.joinpath("Greg_2021_10_09/")
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

# +
# divide probs in N subarrays for memory reasons on local machine
N = 10

s = slices["prob"]
arrs = np.split(np.arange(s.start, 50, s.step), N)
slices_prob = [slice(arr[0], arr[-1]+1, 1) for arr in arrs]
slices_prob
# -

# ## Compute 14C start_values
#
# *Attention:* `"overwrite" = True` in the task disctionary deletes all data in the selected slices. The setting `"overwrite" = False` tries to load an existing archive and extend it by computing incomplete points within the chosen slices.

task = {
    "computation": "start_values_14C",
    "model_type": model_type,
    "overwrite": False,
    "func": CARDAMOMlib.compute_start_values_14C,
    "func_args": {
        "nr_time_steps": params["nr_time_steps"], # nr_months for fake eq_model
        "intcal20_path": intcal20_path,
    }, 
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

for slice_prob in slices_prob[5:]:
    print(slice_prob)
    slices["prob"] = slice_prob
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
# Results will be stored in zarr archives.

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

task = {
    "computation": "solution_14C",
    "model_type": model_type,
    "overwrite": False,
    "func": CARDAMOMlib.compute_solution_14C,
    "func_args": {
        "nr_time_steps": params["nr_time_steps"], # nr_months for fake eq_model
        "intcal20_path": intcal20_path,
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

    for slice_prob in slices_prob[1:]:
        print(slice_prob)
        slices_sub["prob"] = slice_prob

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


# ## Compute 14C external output vector for Southern Hemisphere (lat < -30), Tropics (-30 <= lat <= 30), and Northern Hemisphere (30 < lat)
#
# Results will be stored in zarr archives.

task = {
    "computation": "acc_net_external_output_vector_14C",
    "model_type": model_type,
    "overwrite": False,
    "func": CARDAMOMlib.compute_acc_net_external_output_vector_14C,
    "func_args": {
        "nr_time_steps": params["nr_time_steps"], # nr_months for fake eq_model,
        "intcal20_path": intcal20_path,
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

    for slice_prob in slices_prob[5:]:
        print(slice_prob)
        slices_sub["prob"] = slice_prob

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
data_da = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_4")))
solution_da = data_da[:, :, :, :, 0, :]

solution_14C_da = da.from_zarr(str(project_path.joinpath("solution_14C")))
start_values_14C_da = solution_14C_da[:, :, :, 0, :]

start_values_Delta_14C_da = F_Delta_14C(start_values_da, start_values_14C_da)
system_start_values_Delta_14C_da = F_Delta_14C(start_values_da.sum(-1), start_values_14C_da.sum(-1))    

solution_Delta_14C_da = F_Delta_14C(solution_da, solution_14C_da)
system_solution_Delta_14C_da = F_Delta_14C(solution_da.sum(-1), solution_14C_da.sum(-1))

acc_net_external_output_vector_da = da.from_zarr(str(project_path.joinpath("acc_net_external_output_vector")))
acc_net_external_output_vector_14C_da = da.from_zarr(str(project_path.joinpath("acc_net_external_output_vector_14C")))

system_external_output_Delta_14C_da = F_Delta_14C(
    acc_net_external_output_vector_da.sum(-1),
    acc_net_external_output_vector_14C_da.sum(-1)
)

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
    {"name": "system_solution_Delta_14C", "da": system_solution_Delta_14C_da, "unit": "‰"},
    {"name": "system_external_output_Delta_14C", "da": system_external_output_Delta_14C_da, "unit": "‰"},
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
).chunk({"lat": len(lats_da), "lon": len(lons_da), "prob": 1, "time": len(times_da), "pool": 6})
ds

# +
import itertools


def split_by_chunks(dataset):
    chunk_slices = {}
    for dim, chunks in dataset.chunks.items():
        slices = []
        start = 0
        for chunk in chunks:
            if start >= dataset.sizes[dim]:
                break
            stop = start + chunk
            slices.append(slice(start, stop))
            start = stop
        chunk_slices[dim] = slices
    for slices in itertools.product(*chunk_slices.values()):
        selection = dict(zip(chunk_slices.keys(), slices))
        yield dataset[selection]


# -

for slice_prob in tqdm(slices_prob[:1]):
    print(slice_prob)
    sub_ds = ds.isel(prob=slice_prob)
    datasets = list(split_by_chunks(sub_ds))
    paths = [project_path.joinpath(netCDF_filestem + "_%05d.nc" % prob) for prob in sub_ds.prob]
    xr.save_mfdataset(datasets=datasets, paths=paths)


