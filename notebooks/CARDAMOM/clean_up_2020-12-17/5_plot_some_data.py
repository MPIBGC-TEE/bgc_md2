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

# # Some example plots

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_file = "sol_acc_age_btt.nc"

data_combinations = [
    ("discrete", "monthly"), # the only one available by now, only 1st ensemble member (prob=11)
#    ("discrete", "daily"),
#    ("continuous", "monthly"),
#    ("continuous", "yearly")
]

datasets = dict()
for dc in data_combinations:
    model_type, time_resolution = dc
    params = CARDAMOMlib.load_params(time_resolution)
    output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
    project_path = output_path.joinpath(model_type)

    ds_path = project_path.joinpath(netCDF_file)
    print(dc, ds_path)
    datasets[dc] = xr.open_dataset(ds_path)
# -
ds = datasets[("discrete", "monthly")]
ds


# +
fig, axes = plt.subplots(ncols=2, figsize=(12,6))

var_names = ["xs", "solution"]
for var_name, ax in zip(var_names, axes):
    var = ds[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
#        robust=True
    )
    ax.set_title(var_name)

plt.suptitle("CARDAMOM output (xs) vs model reconstruction (solution)")
plt.tight_layout()
plt.draw()

# +
fig, axes = plt.subplots(ncols=2, figsize=(12,6))

var_names = ["xs_abs_err", "xs_rel_err"]
for var_name, ax in zip(var_names, axes):
    var = ds[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
#        robust=True
    )
    ax.set_title(var_name)

plt.suptitle("Absolute and relative reconstruction error")
plt.tight_layout()
plt.draw()
# -

# This plot shows that we have outliers. The easy way to visualize the data without the outliers is to pass the parameter robust=True. This will use the 2nd and 98th percentiles of the data to compute the color limits.

# +
fig, axes = plt.subplots(ncols=2, figsize=(12,6))

var_names = ["xs_abs_err", "xs_rel_err"]
for var_name, ax in zip(var_names, axes):
    var = ds[var_name]
    var.mean(dim=["prob", "time", "pool"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name)

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
    var = ds[var_name]
    var.mean(dim=["prob", "time"]).plot(
        ax=ax,
        cbar_kwargs={"label": var.attrs["units"]},
        robust=True
    )
    ax.set_title(var_name)

plt.suptitle("System age vs backward transit time")
plt.tight_layout()
plt.draw()
# -

# Now we grab two random pools.

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
        #var = ds[var_name].sel(pool=pool_name)
        ds[var_name].sel(pool=pool_name).mean(dim=["prob", "time"]).plot(
            ax=ax,
            cbar_kwargs={"label": var.attrs["units"]},
            robust=True
        )
        ax.set_title(pool_name + " " + var_name)

plt.suptitle("Pool age")
plt.tight_layout()
plt.draw()
# -

# We can also grab a random site.

# +
(lat, lon, prob) = (28, 52, 0)

var_name_pairs = [
    ("mean_system_age", "mean_btt"),
    ("system_age_median", "btt_median"),
    ("system_age_quantile_05", "btt_quantile_05"),
    ("system_age_quantile_95", "btt_quantile_95"),
    ("system_age_sd", "btt_sd")
]

nrows = len(var_name_pairs) // 2 + len(var_name_pairs) % 2
fig, axes = plt.subplots(nrows=nrows, ncols=2, figsize=(12,6*nrows))

for var_name_pair, ax in zip(var_name_pairs, axes.flatten()):
    for var_name in var_name_pair:
        var = ds[var_name]
        var.isel(lat=lat, lon=lon, prob=prob).plot(
            ax=ax,
            label=var_name
        )
    ax.legend()

plt.suptitle("System age vs backward transit time")
plt.tight_layout()
plt.draw()
# -

# Load land fraction and area dataset.

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

# +
fig, ax = plt.subplots(figsize=(6, 6))

var = ds["xs"].sum(dim="pool")
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
# -

# We compute the global C mean age and transit time and an average of the system age and transit time medians weighted by C stocks.

# +
global_mean_system_age = ds["mean_system_age"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_system_age_median = ds["system_age_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_mean_btt = ds["mean_btt"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])
global_avg_btt_median = ds["btt_median"].weighted(ds_area_lf.area_sphere * ds_area_lf.landfrac).mean(dim=["lat", "lon"])

fig, ax = plt.subplots(figsize=(12, 6))
global_mean_system_age.plot(ax = ax, label="Global C mean age")
global_avg_system_age_median.plot(ax = ax, label="Global average system age median")
global_mean_btt.plot(ax=ax, label="Global mean backward transit time")
global_avg_btt_median.plot(ax=ax, label="Global average backward transit time median")

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
    f(ds[variable_name]).plot(ax=ax)
    ax.set_title(title)
    ax.set_ylabel("yr")
    
plt.tight_layout()
plt.draw()

# +
total_C_xs = ds["xs"].sum(dim="pool")
total_C_xs = total_C_xs * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_xs = total_C_xs.mean(dim=["prob", "time"]).sum()

total_C_solution = ds["solution"].sum(dim="pool")
total_C_solution = total_C_solution * ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac
total_C_solution = total_C_solution.mean(dim=["prob", "time"]).sum()
# -

# Because of problems in the reconstruction of the compartmental matrices from CARDAMOM data, some sites are left out in the reconstructed data. For discrete models the time step can be too large, for continuous models the root finding algorithm might not have converged (this is still to be identified for some small number of sites). Furthermore, the ODE solver for the mean age system complains sometimes when the scales between solution, mean, and second moment are too different.
#
# This leads to the following figures:

print("Total C CARDAMOM            : %.2f PgC" % (total_C_xs *1e-15))
print("Total C model reconstruction: %.2f PgC" % (total_C_solution*1e-15))
print("Coverage                    :   %2.2f %%" % (total_C_solution/total_C_xs*100))

