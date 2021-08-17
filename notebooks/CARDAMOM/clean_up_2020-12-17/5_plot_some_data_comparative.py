# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Some example plots

# +
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_file = "sol_acc_age_btt.nc"

data_combinations = [
#    ("monthly", None, "discrete"),
#    ("monthly", None, "continuous"),
    ("yearly", 0, "continuous"),
    ("yearly", 6, "continuous")
]

datasets = dict()
for dc in data_combinations:
    time_resolution, delay_in_months, model_type = dc
    params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
    output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
    project_path = output_path.joinpath(model_type)

    ds_path = project_path.joinpath(netCDF_file)
    print(dc, ds_path)
    datasets[dc] = xr.open_dataset(ds_path)
# -
ds1 = datasets[("yearly", 6, "continuous")]
ds2 = datasets[("yearly", 0, "continuous")]
ds1.coords["time"] = pd.DatetimeIndex(ds1.time.data).year
ds2.coords["time"] = pd.DatetimeIndex(ds2.time.data).year


# +
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12,18))

var_names = ["xs", "solution"]
for var_name, ax in zip(var_names, axes[0]):
    var = ds1[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (July)")

for var_name, ax in zip(var_names, axes[1]):
    var = ds2[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (January)")
    
for var_name, ax in zip(var_names, axes[2]):
    var = ds1[var_name] - ds2[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": ds1[var_name].attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (July minus January)")
    
plt.suptitle("CARDAMOM output (xs) vs model reconstruction (solution)")
plt.tight_layout()
plt.draw()

# +
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,12))

var_names = ["xs_abs_err", "xs_rel_err"]
for var_name, ax in zip(var_names, axes[0]):
    var = ds1[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (July)")

var_names = ["xs_abs_err", "xs_rel_err"]
for var_name, ax in zip(var_names, axes[1]):
    var = ds2[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (January)")

plt.suptitle("Absolute and relative reconstruction error")
plt.tight_layout()
plt.draw()

# +
var_names = [
    "mean_system_age", "mean_btt",
    "system_age_median", "btt_median",
    "system_age_quantile_05", "btt_quantile_05",
    "system_age_quantile_95", "btt_quantile_95",
    "system_age_sd", "btt_sd"
]

nrows = len(var_names) // 2 + len(var_names) % 2
fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12,6*nrows))

for var_name, ax in zip(var_names, axes.flatten()):
    var = ds2[var_name] - ds1[var_name]
    var.mean(dim=["prob", "time"]).plot(
        ax=ax,
        cbar_kwargs={"label": ds1[var_name].attrs["units"]},
        robust=True
    )
    ax.set_title(var_name + " (July minus January)")

plt.suptitle("System age vs backward transit time")
#plt.tight_layout()
plt.draw()

# +
nrows = len(var_names) // 2 + len(var_names) % 2
fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12,6*nrows))

for var_name, ax in zip(var_names, axes.flatten()):
    var = (ds2[var_name] - ds1[var_name]) / ds1[var_name] * 100
    var.mean(dim=["prob", "time"]).plot(
        ax=ax,
#        cbar_kwargs={"label": ds1[var_name].attrs["units"]},
        cbar_kwargs={"label": "%"},
        robust=True
    )
    ax.set_title(var_name + " (Jul-Jan)/Jul*100")

plt.suptitle("System age vs backward transit time")
#plt.tight_layout()
plt.draw()

# +
pool_names = ["Soil", "Litter"]

var_names = [
    "mean_pool_age_vector",
    "pool_age_median",
    "pool_age_quantile_05",
    "pool_age_quantile_95",
    "pool_age_sd_vector"
]

ncols = len(pool_names)
nrows = len(var_names)
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,6*nrows))

for row, var_name in enumerate(var_names):
    for col, pool_name in enumerate(pool_names):
        ax = axes[row, col]
        var1 = ds1[var_name].sel(pool=pool_name)
        var2 = ds2[var_name].sel(pool=pool_name)
        var = var1 - var2
        var.mean(dim=["prob", "time"]).plot(
            ax=ax,
            cbar_kwargs={"label": var1.attrs["units"]},
            robust=True
        )
        ax.set_title(pool_name + " " + var_name)

plt.suptitle("Pool age (July minus January)")
plt.tight_layout()
plt.draw()

# +
ncols = len(pool_names)
nrows = len(var_names)
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,6*nrows))

for row, var_name in enumerate(var_names):
    for col, pool_name in enumerate(pool_names):
        ax = axes[row, col]
        var1 = ds1[var_name].sel(pool=pool_name)
        var2 = ds2[var_name].sel(pool=pool_name)
        var = (var1 - var2)/var1*100
        var.mean(dim=["prob", "time"]).plot(
            ax=ax,
#            cbar_kwargs={"label": var.attrs["units"]},
            cbar_kwargs={"label": "%"},
            robust=True
        )
        ax.set_title(pool_name + " " + var_name)

plt.suptitle("Pool age (Jul-Jan)/Jul*100")
#plt.tight_layout()
plt.draw()
# -

# We can also grab a random site.

# +
(lat, lon, prob) = (30, 39, 0)
print("Latitude", ds1.lat[lat])
print("Longitude", ds1.lon[lon])

var_name_pairs = [
    ("mean_system_age", "mean_btt"),
    ("system_age_median", "btt_median"),
    ("system_age_quantile_05", "btt_quantile_05"),
    ("system_age_quantile_95", "btt_quantile_95"),
    ("system_age_sd", "btt_sd")
]

nrows = len(var_name_pairs) // 2 + len(var_name_pairs) % 2
fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(2*12,2*6*nrows))

for var_name_pair, ax in zip(var_name_pairs, axes.flatten()):
    for var_name in var_name_pair:
        ds1[var_name].isel(lat=lat, lon=lon, prob=prob).plot(
            ax=ax,
            label=var_name + " (Jul)",
            ls="-"
        )
        ds2[var_name].isel(lat=lat, lon=lon, prob=prob).plot(
            ax=ax,
            label=var_name + " (Jan)",
            ls="--"
        )
    ax.legend()

plt.suptitle("System age vs backward transit time")
#plt.tight_layout()
plt.draw()
# -

# Load land fraction and area dataset.

ds_area_lf = xr.open_dataset(data_path.joinpath("LSFRAC2.nc"))
ds_area_lf

# As we can see, unfortunately, the coordinates do not perfectly match the given CARDAMOM output dataset coordinates. We adjust them.

# +
ds_area_lf_adapted = ds_area_lf.sel(lat=ds1.lat, lon=ds1.lon)
ds_area_lf_adapted.area_sphere.attrs["units"] = "m^2"

ds_area_lf_adapted.coords["lat"] = np.float64(ds_area_lf_adapted.coords["lat"])
ds_area_lf_adapted.coords["lon"] = np.float64(ds_area_lf_adapted.coords["lon"])

for var_name, var in ds_area_lf_adapted.data_vars.items():
    ds_area_lf_adapted[var_name].coords["lat"] = np.float64(var.coords["lat"])
    ds_area_lf_adapted[var_name].coords["lon"] = np.float64(var.coords["lon"])
ds_area_lf_adapted

# +
fig, axes = plt.subplots(ncols=3, figsize=(18, 6))

var = ds1["xs"].sum(dim="pool")
var = var * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
var = var.mean(dim=["prob", "time"])
(var/1e15).plot(
    ax=axes[0],
    cbar_kwargs={"label": "PgC"},
    robust=True
)
axes[0].set_title("Total C stocks (July)")

var = ds2["xs"].sum(dim="pool")
var = var * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
var = var.mean(dim=["prob", "time"])
(var/1e15).plot(
    ax=axes[1],
    cbar_kwargs={"label": "PgC"},
    robust=True
)
axes[1].set_title("Total C stocks (January)")

var = (ds1["xs"]-ds2["xs"]).sum(dim="pool")
var = var * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
var = var.mean(dim=["prob", "time"])
(var/1e15).plot(
    ax=axes[2],
    cbar_kwargs={"label": "PgC"},
    robust=True
)
axes[2].set_title("Total C stocks (July minus January)")

plt.tight_layout()
plt.draw()
# -

# We compute the global C mean age and transit time and an average of the system age and transit time medians weighted by C stocks.

# +
global_mean_system_age_1 = ds1["mean_system_age"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_system_age_median_1 = ds1["system_age_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_mean_btt_1 = ds1["mean_btt"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_btt_median_1 = ds1["btt_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])

global_mean_system_age_2 = ds2["mean_system_age"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_system_age_median_2 = ds2["system_age_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_mean_btt_2 = ds2["mean_btt"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_btt_median_2 = ds2["btt_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])

fig, ax = plt.subplots(figsize=(12, 6))
global_mean_system_age_1.plot(ax = ax, label="Global C mean age (Jul)", c="blue", ls="-")
global_mean_system_age_2.plot(ax = ax, label="Global C mean age (Jan)", c="blue", ls="--")

global_avg_system_age_median_1.plot(ax = ax, label="Global average system age median (Jul)", c="orange", ls="-")
global_avg_system_age_median_2.plot(ax = ax, label="Global average system age median (Jan)", c="orange", ls="--")

global_mean_btt_1.plot(ax=ax, label="Global mean backward transit time (Jan)", c="red", ls="-")
global_mean_btt_2.plot(ax=ax, label="Global mean backward transit time (Jul)", c="red", ls="--")

global_avg_btt_median_1.plot(ax=ax, label="Global average backward transit time median (Jul)", c="green", ls="-")
global_avg_btt_median_2.plot(ax=ax, label="Global average backward transit time median (Jan)", c="green", ls="--")

ax.set_ylabel("yr")
ax.legend()

plt.tight_layout()
plt.draw()
# -

# We have a closer look at each of such graphs.

# +
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(12, 12))

f = lambda x: x.weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
variable_names = ["mean_system_age", "system_age_median", "mean_btt", "btt_median"]
titles = [
    "Global C mean age",
    "Global average system age median",
    "Global mean backward transit time",
    "Global average backward transit time median"
]

for variable_name, ax, title in zip(variable_names, axes.flatten(), titles):
    f(ds1[variable_name]).plot(ax=ax, label="Jul", c="red")
    f(ds2[variable_name]).plot(ax=ax, label="Jan", c="blue")
    ax.set_title(title)
    ax.set_ylabel("yr")
    ax.legend()
    
plt.tight_layout()
plt.draw()
# -



# Because of problems in the reconstruction of the compartmental matrices from CARDAMOM data, some sites are left out in the reconstructed data. For discrete models the time step can be too large, for continuous models the root finding algorithm might not have converged (this is still to be identified for some small number of sites). Furthermore, the ODE solver for the mean age system complains sometimes when the scales between solution, mean, and second moment are too different.
#
# This leads to the following figures:

# +
total_C_xs = ds1["xs"].sum(dim="pool")
total_C_xs = total_C_xs * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_xs = total_C_xs.mean(dim=["prob", "time"]).sum()

total_C_solution = ds1["solution"].sum(dim="pool")
total_C_solution = total_C_solution * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_solution = total_C_solution.mean(dim=["prob", "time"]).sum()
# -

print("July")
print("----")
print("Total C CARDAMOM            : %.2f PgC" % (total_C_xs *1e-15))
print("Total C model reconstruction: %.2f PgC" % (total_C_solution*1e-15))
print("Coverage                    :   %2.2f %%" % (total_C_solution/total_C_xs*100))

# +
total_C_xs = ds2["xs"].sum(dim="pool")
total_C_xs = total_C_xs * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_xs = total_C_xs.mean(dim=["prob", "time"]).sum()

total_C_solution = ds2["solution"].sum(dim="pool")
total_C_solution = total_C_solution * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_solution = total_C_solution.mean(dim=["prob", "time"]).sum()
# -

print("January")
print("----")
print("Total C CARDAMOM            : %.2f PgC" % (total_C_xs *1e-15))
print("Total C model reconstruction: %.2f PgC" % (total_C_solution*1e-15))
print("Coverage                    :   %2.2f %%" % (total_C_solution/total_C_xs*100))

# ## plot Soil age time series, plot CHANGES in age and transit time time series, maybe even as histogram (spatially, termporally)


