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

# # Make global nc file
#
# Load gloabl BTT nc file and add ``global_GPP``, ``gloabl_NPP``, ``global_C``, ``global_Rh``, ``global_Ra``, ``global_mean_btt``, and ``global_mean_age`` weighted with area and landfrac. Save to file according to value of ``correct_for_autotrophic_respiration``.

# +
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from pathlib import Path

from bgc_md2.models.CARDAMOM import CARDAMOMlib


# +
#data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
data_path = Path("/home/data/CARDAMOM/Greg_2021_10_09/")

netCDF_filestem = "sol_acc_age_btt"

dc = ("monthly", None, "discrete")
time_resolution, delay_in_months, model_type = dc
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
project_path = output_path.joinpath(model_type)

correct_for_autotrophic_respiration = True
#correct_for_autotrophic_respiration = False

ds_path = project_path.joinpath(netCDF_filestem)
print(dc, ds_path)
ds = xr.open_mfdataset(str(ds_path) + "*.nc")
ds

# +
if not correct_for_autotrophic_respiration:
    # create Ras = 0 with correct shape
    Ras = 0 * ds["Us"].sum(dim="pool")
else:
    GPP = ds.GPP
    Ras = GPP - ds["Us"].sum(dim="pool")

    # following Greg's advice with the current data I ignore negatve Ra values
    Ras = Ras.where(Ras > 0, 0)
        
# add Ras to dataset
ds = ds.assign({"Ras": Ras})  
# -

ds

# # Load merged dataset with global quantiles

if correct_for_autotrophic_respiration:
    correction_str = "_corrected"
else:
    correction_str = ""

ds_btt = xr.open_dataset(project_path.joinpath("global_btt_quantiles" + correction_str + ".nc"))
ds_btt

# ## Load area and landfrac data

ds_area_lf = xr.open_dataset(data_path.joinpath("LSFRAC2.nc"))
ds_area_lf

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

# # Add other global variables to global quantile dataset

# +
# correct mean btt for autotrophic respiration

# global mean transit time
Rs = ds["acc_net_external_output_vector"].sum(dim="pool")
mean_btt = (Ras * 0 + Rs * ds["mean_btt"]) / (Ras + Rs)
btt_weights = (ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac * Rs).fillna(0)
global_mean_btt = mean_btt.weighted(btt_weights).mean(dim=["lat", "lon"])

## global btt_moment_2
#btt_moment_2 = (Ras * 0 + Rs * ds["btt_moment_2"]) / (Ras + Rs)
#global_btt_moment_2 = btt_moment_2.weighted(btt_weights).mean(dim=["lat", "lon"])

# global btt_moments
btt_moment_2 = (Ras * 0 + Rs * ds["btt_moment_2"]) / (Ras + Rs)
btt_moment_3 = (Ras * 0 + Rs * ds["btt_moment_3"]) / (Ras + Rs)
btt_moment_4 = (Ras * 0 + Rs * ds["btt_moment_4"]) / (Ras + Rs)

global_btt_moment_2 = btt_moment_2.weighted(btt_weights).mean(dim=["lat", "lon"])
global_btt_moment_3 = btt_moment_3.weighted(btt_weights).mean(dim=["lat", "lon"])
global_btt_moment_4 = btt_moment_4.weighted(btt_weights).mean(dim=["lat", "lon"])

# BTT standard deviation, skewness, and kurtosis
#Rs_vector = ds["acc_net_external_output_vector"]
global_btt_sd = np.sqrt(global_btt_moment_2 - global_mean_btt**2)
global_btt_skewness = (global_btt_moment_3 - 3*global_mean_btt*global_btt_sd - global_mean_btt**3) / global_btt_sd**3
global_btt_kurtosis = (global_btt_moment_4 - 4*global_mean_btt*global_btt_moment_3 + 6*global_mean_btt**2*global_btt_moment_2 - 3*global_mean_btt**4) / global_btt_sd**4

# global mean age
age_weights = (ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac * ds.solution.sum(dim="pool")).fillna(0)
global_mean_age = ds["mean_system_age"].weighted(age_weights).mean(dim=["lat", "lon"])

# pure area and landfrac weights
weights = (ds_area_lf_adapted.area_sphere * ds_area_lf_adapted.landfrac).fillna(0)

# global heterotrophic respiration
global_Rs = Rs.weighted(weights).sum(dim=["lat", "lon"])

# global autotrophic respiration
global_Ras = Ras.weighted(weights).sum(dim=["lat", "lon"])

# global GPP
global_GPPs = ds["GPP"].weighted(weights).sum(dim=["lat", "lon"])

# global Us, actually NPPs
NPPs = ds["Us"].sum(dim="pool")
global_NPPs = NPPs.weighted(weights).sum(dim=["lat", "lon"])

# global stocks
stocks = ds["solution"].sum(dim="pool")
global_stocks = stocks.weighted(weights).sum(dim=["lat", "lon"])
# -


# ## Global BTT skewness and kurtosis
#
# Both are heavily influenced by outliers in BTT, so we can cut off a certain level of outliers (``q``). Then the numbers go down to more realistic values but the long-term trend remains the same.

# ### Skewness and kurtosis without outlier cutoff

fig, ax = plt.subplots(figsize=(20,8))
#var = global_btt_skewness
#_ = var.rolling(time=120).mean().plot.line(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)
var = global_btt_skewness.mean(dim="prob")
_ = var.rolling(time=120).mean().plot(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)

fig, ax = plt.subplots(figsize=(20,8))
#var = global_btt_kurtosis
#_ = var.rolling(time=120).mean().plot.line(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)
var = global_btt_kurtosis.mean(dim="prob")
_ = var.rolling(time=120).mean().plot(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)

var0 = ds.btt_skewness.compute()
var1 = Rs.compute()
var2 = global_btt_skewness.compute()

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
var0.mean(dim=["prob", "time"]).plot(ax=axes[0])
var1.mean(dim=["prob", "time"]).plot(ax=axes[1])
var2.mean(dim=["prob"]).rolling(time=120).mean().plot(ax=axes[2])

var3 = ds.btt_kurtosis.compute()
var4 = global_btt_kurtosis.compute()

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
var3.mean(dim=["prob", "time"]).plot(ax=axes[0])
var1.mean(dim=["prob", "time"]).plot(ax=axes[1])
var4.mean(dim=["prob"]).rolling(time=120).mean().plot(ax=axes[2])

var5 = ds.mean_btt.compute()
var6 = global_mean_btt.compute()

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
var5.mean(dim=["prob", "time"]).plot(ax=axes[0])
var1.mean(dim=["prob", "time"]).plot(ax=axes[1])
var6.mean(dim=["prob"]).rolling(time=120).mean().plot(ax=axes[2])

# ### Cut off outliers

q = 0.95 # thow away higher values than this percentile


# +
def boolnan(b, val):
    if b:
        return np.nan
    return val

boolnanvec = np.vectorize(boolnan)

def remove_outliers(vars, q):
    return_list = []
    for var in vars:
        return_list.append(
            xr.DataArray(
                data=boolnanvec(var > var.quantile(q), var),
                coords=var.coords
            )
        )
        
    return return_list


# -

mean_btt_no_outliers, btt_moment_2_no_outliers, btt_moment_3_no_outliers, btt_moment_4_no_outliers = remove_outliers(
    [mean_btt.chunk({"prob": -1}), btt_moment_2.chunk({"prob": -1}), btt_moment_3.chunk({"prob": -1}), btt_moment_4.chunk({"prob": -1})],
    q
)

# +
global_mean_btt_no_outliers = mean_btt_no_outliers.weighted(btt_weights).mean(dim=["lat", "lon"])

# global btt_moments
btt_moment_2_no_outliers = (Ras * 0 + Rs * btt_moment_2_no_outliers) / (Ras + Rs)
btt_moment_3_no_outliers = (Ras * 0 + Rs * btt_moment_3_no_outliers) / (Ras + Rs)
btt_moment_4_no_outliers = (Ras * 0 + Rs * btt_moment_4_no_outliers) / (Ras + Rs)

global_btt_moment_2_no_outliers = btt_moment_2_no_outliers.weighted(btt_weights).mean(dim=["lat", "lon"])
global_btt_moment_3_no_outliers = btt_moment_3_no_outliers.weighted(btt_weights).mean(dim=["lat", "lon"])
global_btt_moment_4_no_outliers = btt_moment_4_no_outliers.weighted(btt_weights).mean(dim=["lat", "lon"])
# -

fig, axes = plt.subplots(figsize=(12, 6), ncols=2)
global_mean_btt.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
global_mean_btt_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1])

fig, axes = plt.subplots(figsize=(12, 6), ncols=2)
global_btt_moment_2.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
global_btt_moment_2_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1])

fig, axes = plt.subplots(figsize=(12, 6), ncols=2)
global_btt_moment_3.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
global_btt_moment_3_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1])

fig, axes = plt.subplots(figsize=(12, 6), ncols=2)
global_btt_moment_4.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
global_btt_moment_4_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1])

global_btt_sd_no_outliers = np.sqrt(global_btt_moment_2_no_outliers - global_mean_btt_no_outliers**2)
global_btt_skewness_no_outliers = (global_btt_moment_3_no_outliers - 3*global_mean_btt_no_outliers*global_btt_sd_no_outliers - global_mean_btt_no_outliers**3) / global_btt_sd_no_outliers**3
global_btt_kurtosis_no_outliers = (global_btt_moment_4_no_outliers - 4*global_mean_btt_no_outliers*global_btt_moment_3_no_outliers + 6*global_mean_btt_no_outliers**2*global_btt_moment_2_no_outliers - 3*global_mean_btt_no_outliers**4) / global_btt_sd_no_outliers**4


fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
global_btt_sd.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
ds.btt_sd.weighted(btt_weights).mean(dim=["lat", "lon"]).mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1]) # this is how it's not done because standard deviation is not linear
global_btt_sd_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[2])

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
global_btt_skewness.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
ds.btt_skewness.weighted(btt_weights).mean(dim=["lat", "lon"]).mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1]) # this is how it's not done because skewness is not linear
global_btt_skewness_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[2])

fig, axes = plt.subplots(figsize=(18, 6), ncols=3)
global_btt_kurtosis.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[0])
ds.btt_kurtosis.weighted(btt_weights).mean(dim=["lat", "lon"]).mean(dim="prob").rolling(time=120).mean().plot(ax=axes[1])  # this is how it's not done because kurtosis is not linear
global_btt_kurtosis_no_outliers.mean(dim="prob").rolling(time=120).mean().plot(ax=axes[2])

# ## Add variables to global dataset and write to disk

# +
new_data_vars = {
    "global_mean_btt": global_mean_btt,
    "global_btt_moment_2": global_btt_moment_2,
    "global_btt_moment_3": global_btt_moment_3,
    "global_btt_moment_4": global_btt_moment_4,
    "global_btt_skewness": global_btt_skewness,
    "global_btt_kurtosis": global_btt_kurtosis,

    "global_mean_btt_no_outliers": global_mean_btt_no_outliers,
    "global_btt_moment_2_no_outliers": global_btt_moment_2_no_outliers,
    "global_btt_moment_3_no_outliers": global_btt_moment_3_no_outliers,
    "global_btt_moment_4_no_outliers": global_btt_moment_4_no_outliers,
    "global_btt_skewness_no_outliers": global_btt_skewness_no_outliers,
    "global_btt_kurtosis_no_outliers": global_btt_kurtosis_no_outliers,
    
    "global_mean_age": global_mean_age,
    "global_C": global_stocks,
    "global_GPP": global_GPPs,
    "global_NPP": global_NPPs,
    "global_Rh": global_Rs,
    "global_Ra": global_Ras
}

ds_target = ds_btt.assign(new_data_vars)
ds_target.attrs = ds.attrs
ds_target
# -

filename = project_path.joinpath("global_ds" + correction_str + ".nc")
ds_target.to_netcdf(filename)
filename

# ## Load global dataset and plot global GPP

ds_global = xr.open_dataset(str(project_path.joinpath("global_ds"+correction_str+".nc")))
ds_global

# +
global_GPP = ds_global["global_GPP"]

fig, ax = plt.subplots(figsize=(20,8))
global_GPP.mean(dim="prob").plot(ax=ax, x="time")
# -

fig, ax = plt.subplots(figsize=(20,8))
var = ds_global.global_btt_moment_2 - ds_global.global_mean_btt**2
sd = np.sqrt(var)
cv = sd / ds_global.global_mean_btt
_ = cv.rolling(time=12).mean().plot.line(ax=ax, x="time", add_legend=False, c="blue", alpha=0.2, lw=4)


