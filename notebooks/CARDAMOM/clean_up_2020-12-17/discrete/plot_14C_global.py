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

# # Make some global 14C plots (data from discrete model runs)

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib 
# -

# ignore nan/nan warnings
import warnings
warnings.simplefilter("ignore")

# +
data_combinations = [
    ("monthly", None, "discrete"),
    ("monthly", None, "continuous"),
    ("yearly", 0, "continuous"),
    ("yearly", 6, "continuous")
]

#CARDAMOM_path = Path("/home/data/CARDAMOM/")
CARDAMOM_path = Path("/mnt/c/Users/hrme0001/data/CARDAMOM/")
#data_path = CARDAMOM_path.joinpath("Greg_2020_10_26/")
data_path = CARDAMOM_path.joinpath("Greg_2021_10_09/")
netCDF_filestem = "global_14C_"


dc = data_combinations[0]
time_resolution, delay_in_months, model_type = dc
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
project_path = output_path.joinpath(model_type)

ds_path = project_path.joinpath(netCDF_filestem+"*.nc")
ds = xr.open_mfdataset(str(ds_path))
ds
# -

# We plot global system Delta14C for different moments in time.

# +
time_indices = np.linspace(0, len(ds.time)-1, 6)
fig, axes = plt.subplots(ncols=3, nrows=2, figsize=(6*3, 6*2))

for time_index, ax in zip(time_indices, axes.flatten()):
    var = ds["system_solution_Delta_14C"].isel(time=int(time_index))
    var.mean(dim=["prob"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )


plt.suptitle("CARDAMOM Delta_14C")
plt.tight_layout()
plt.draw()
# -

# Now we plot the time evolution of system Delta14C fo the northern hemishpere, the southern hemisphere, and the tropics. For a global average we need to transform the mass in the cells from g/m^2 to g. To that end, we load area and landfrac of the cells.

ds_area_lf = xr.open_dataset(data_path.joinpath("LSFRAC2.nc"))
ds_area_lf

# As we can see, unfortunately, the coordinates do not perfectly match the given CARDAMOM output dataset coordinates. We adjust them.

# +
ds_area_lf_adapted = ds_area_lf.sel(lat=ds.lat, lon=ds.lon)
ds_area_lf_adapted.area_sphere.attrs["units"] = "m^2"

ds_area_lf_adapted.coords["lat"] = np.float64(ds_area_lf_adapted.coords["lat"])
ds_area_lf_adapted.coords["lon"] = np.float64(ds_area_lf_adapted.coords["lon"])

for var_name, var in ds_area_lf_adapted.data_vars.items():
    ds_area_lf_adapted[var_name].coords["lat"] = np.float64(var.coords["lat"])
    ds_area_lf_adapted[var_name].coords["lon"] = np.float64(var.coords["lon"])
    
ds_area_lf_adapted
# -

# ## Plot radiocarbon of the stocks

# +
fig, ax = plt.subplots(figsize=(8, 8))

var = ds["solution"].sum(dim="pool")
var = var * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
var = var.mean(dim=["prob", "time"])
    
(var/1e15).plot(
    ax=ax,
    cbar_kwargs={"label": "PgC"},
    robust=True
)
ax.set_title("Total C stocks")

plt.tight_layout()
plt.draw()

# +
fig, ax = plt.subplots(figsize=(18, 12))

for pool in ds.pool:
    ds["solution_Delta_14C"].sel(pool=pool).weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon", "prob"]).plot(
        ax=ax,
        lw=4,
        label=pool.values,
    )
    
ax.legend()
ax.set_xlim([ds.time[0].values, ds.time[-1].values])
ax.set_ylabel(r"$\Delta^{14}$C (‰)")
ax.set_title(r"Global avergage $\Delta^{14}$C (mean over total ensemble)")
# -

# Now we aggregate the pools and show all ensemble members.

# + tags=[]
fig, ax = plt.subplots(figsize=(18, 12))

lat_slices = [slice(0, 6, 1), slice(6, 22, 1), slice(22, None, 1)]
colors = ["blue", "green", "red"]
labels = ["Southern Hemisphere", "Tropics", "Northern Hemisphere"]
               
for lat_slice, color in zip(lat_slices, colors):
    ds["system_solution_Delta_14C"].isel(lat=lat_slice).weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"]).plot.line(
        ax=ax,
        x="time",
        add_legend=False,
        alpha=0.1,
        lw=4,
        c=color
)
              
ax.set_xlim([ds.time[0].values, ds.time[-1].values])
ax.set_ylabel(r"$\Delta^{14}$C (‰)")
ax.set_title(r"Global average $\Delta^{14}$C")
               
from matplotlib.lines import Line2D
lines = [Line2D([0], [0], color=c, linewidth=4, linestyle='-') for c in colors]
_ = ax.legend(lines, labels)
# -
# ## Plot radiocarbon of the external outputs

# + tags=[]
fig, ax = plt.subplots(figsize=(18, 12))

lat_slices = [slice(0, 6, 1), slice(6, 22, 1), slice(22, None, 1)]
colors = ["blue", "green", "red"]
labels = ["Southern Hemisphere", "Tropics", "Northern Hemisphere"]
               
for lat_slice, color in zip(lat_slices, colors):
    ds["system_external_output_Delta_14C"].isel(lat=lat_slice).weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"]).plot.line(
        ax=ax,
        x="time",
        add_legend=False,
        alpha=0.2,
        lw=4,
        c=color
)
              
ax.set_xlim([ds.time[0].values, ds.time[-1].values])
ax.set_ylabel(r"$\Delta^{14}$C (‰)")
ax.set_title(r"Global respiration average $\Delta^{14}$C")
               
from matplotlib.lines import Line2D
lines = [Line2D([0], [0], color=c, linewidth=4, linestyle='-') for c in colors]
_ = ax.legend(lines, labels)
# -


