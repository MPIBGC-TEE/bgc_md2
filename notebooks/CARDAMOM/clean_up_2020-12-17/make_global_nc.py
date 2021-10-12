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
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
netCDF_filestem = "sol_acc_age_btt"

dc = ("monthly", None, "discrete")
time_resolution, delay_in_months, model_type = dc
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)
output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))
project_path = output_path.joinpath(model_type)

correct_for_autotrophic_respiration = True

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
    Ras = Ras.where(Ras>0, 0)
        
# add Ras to dataset
ds = ds.assign({"Ras": Ras})  
# -

# # Load merged dataset with global quantiles

if correct_for_autotrophic_respiration:
    correction_str = "_corrected"
else:
    correction_str = ""

ds_btt = xr.open_dataset(project_path.joinpath("global_btt_quantiles"+correction_str+".nc"))
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


# +
new_data_vars = {
    "global_mean_btt": global_mean_btt,
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

filename = project_path.joinpath("global_ds"+correction_str+".nc")
ds_target.to_netcdf(filename)
filename


