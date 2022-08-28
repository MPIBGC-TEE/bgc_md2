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

# # Convert output from zarr to a netCDF file

# +
import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr

from pathlib import Path
from tqdm import tqdm

from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask.distributed import Client
from dask import delayed
# -

# for monthly discrete data, each worker needs about 10GB
my_cluster = CARDAMOMlib.prepare_cluster(
    n_workers=1,
    alternative_dashboard_port=8790
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

#time_resolution, delay_in_months, model_type = "daily", None, "discrete"
time_resolution, delay_in_months, model_type = "monthly", None, "discrete"
#time_resolution, delay_in_months, model_type = "monthly", None, "continuous"
#time_resolution, delay_in_months, model_type = "yearly", 0, "continuous"
#time_resolution, delay_in_months, model_type = "yearly", 6, "continuous"

# +
params = CARDAMOMlib.load_params(time_resolution, delay_in_months)

#data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
data_path = Path("/home/data/CARDAMOM/Greg_2021_10_09/")

output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))

project_path = output_path.joinpath(model_type)
print(project_path)
netCDF_filestem = "sol_acc_age_btt"
# -

lats_da = da.from_zarr(str(project_path.joinpath("lat")))
lons_da = da.from_zarr(str(project_path.joinpath("lon")))
probs_da = da.from_zarr(str(project_path.joinpath("prob")))
times_da = da.from_zarr(str(project_path.joinpath("time")))

# +
nr_pools = 6

# use all "lat", all "lon", the first four "prob", all "time"
slices = {
    "lat": slice(0, None, 1),
    "lon": slice(0, None, 1),
    "prob": slice(30, None, 1), # done: slice(0, 30, 1)
    "time": slice(0, None, 1) # don't change the time entry
}

# +
# dimensions: lat x lon x prob x time x order x pool

start_values_da = da.from_zarr(str(project_path.joinpath("start_values")))
if model_type == "continuous":
    us_da = da.from_zarr(str(project_path.joinpath("us")))
elif model_type == "discrete":
    Us_da = da.from_zarr(str(project_path.joinpath("Us")))
else:
    raise(TypeError("unknown model type '%s'" % model_type))
Bs_da = da.from_zarr(str(project_path.joinpath("Bs")))

xs_da = da.from_zarr(str(project_path.joinpath("xs")))
#data_da = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_2")))
data_da = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_4")))
solution_da = data_da[:, :, :, :, 0, :]

GPPs_da = da.from_zarr(str(project_path.joinpath("GPPs")))

# compute absolute and relative errors of reconstructed solutions
xs_abs_err_da = np.abs(xs_da-solution_da)
xs_rel_err_da = xs_abs_err_da/xs_da * 100

# age moment vectors
mean_pool_age_vector_da = data_da[:, :, :, :, 1, :]
pool_age_moment_vector_2_da = data_da[:, :, :, :, 2, :]
pool_age_moment_vector_3_da = data_da[:, :, :, :, 3, :]
pool_age_moment_vector_4_da = data_da[:, :, :, :, 4, :]
pool_age_sd_vector_da = da.sqrt(pool_age_moment_vector_2_da - mean_pool_age_vector_da**2)

# system age moments
mean_system_age_da = (solution_da * mean_pool_age_vector_da).sum(-1) / solution_da.sum(-1)
system_age_moment_2_da = (solution_da * pool_age_moment_vector_2_da).sum(-1) / solution_da.sum(-1)
system_age_sd_da = da.sqrt(system_age_moment_2_da - mean_system_age_da)

## pool age median and quantiles
#pool_age_median_da = da.from_zarr(str(project_path.joinpath("pool_age_median")))
#pool_age_quantile_05_da = da.from_zarr(str(project_path.joinpath("pool_age_quantile_05")))
#pool_age_quantile_95_da = da.from_zarr(str(project_path.joinpath("pool_age_quantile_95")))

## system age median and quantiles
#system_age_median_da = da.from_zarr(str(project_path.joinpath("system_age_median")))
#system_age_quantile_05_da = da.from_zarr(str(project_path.joinpath("system_age_quantile_05")))
#system_age_quantile_95_da = da.from_zarr(str(project_path.joinpath("system_age_quantile_95")))

# backward transit time moments
if model_type == "continuous":
    external_output_vector_da = da.from_zarr(str(project_path.joinpath("external_output_vector")))
    mean_btt_da = (external_output_vector_da * mean_pool_age_vector_da).sum(-1) / external_output_vector_da.sum(-1)
    btt_moment_2_da = (external_output_vector_da * pool_age_moment_vector_2_da).sum(-1) / external_output_vector_da.sum(-1)
    btt_moment_3_da = (external_output_vector_da * pool_age_moment_vector_3_da).sum(-1) / external_output_vector_da.sum(-1)
    btt_moment_4_da = (external_output_vector_da * pool_age_moment_vector_4_da).sum(-1) / external_output_vector_da.sum(-1)
elif model_type == "discrete":
    acc_net_external_output_vector_da = da.from_zarr(str(project_path.joinpath("acc_net_external_output_vector")))
    mean_btt_da = (acc_net_external_output_vector_da * mean_pool_age_vector_da).sum(-1) / acc_net_external_output_vector_da.sum(-1)
    btt_moment_2_da = (acc_net_external_output_vector_da * pool_age_moment_vector_2_da).sum(-1) / acc_net_external_output_vector_da.sum(-1)
    btt_moment_3_da = (acc_net_external_output_vector_da * pool_age_moment_vector_3_da).sum(-1) / acc_net_external_output_vector_da.sum(-1)
    btt_moment_4_da = (acc_net_external_output_vector_da * pool_age_moment_vector_4_da).sum(-1) / acc_net_external_output_vector_da.sum(-1)
else:
    raise(TypeError("unknown model type '%s'" % model_type))

# backward transtit time standard deviation
btt_sd_da = np.sqrt(btt_moment_2_da - mean_btt_da**2)

# backward transit time skewness and kurtosis
btt_skewness_da = (btt_moment_3_da - 3*mean_btt_da*btt_sd_da - mean_btt_da**3) / btt_sd_da**3
btt_kurtosis_da = (btt_moment_4_da - 4*mean_btt_da*btt_moment_3_da + 6*mean_btt_da**2*btt_moment_2_da - 3*mean_btt_da**4) / btt_sd_da**4

# backward transit time median and quantiles
btt_median_da = da.from_zarr(str(project_path.joinpath("btt_median")))
btt_quantile_05_da = da.from_zarr(str(project_path.joinpath("btt_quantile_05")))
btt_quantile_95_da = da.from_zarr(str(project_path.joinpath("btt_quantile_95")))

# +
coords = {
    "lat": lats_da.reshape(-1)[slices["lat"]],
    "lon": lons_da.reshape(-1)[slices["lon"]],
    "prob": probs_da.reshape(-1)[slices["prob"]],
    "time": pd.DatetimeIndex(times_da.reshape(-1).compute()),
    "pool": CARDAMOMlib.load_model_structure().pool_names.copy(),
    "pool_to": CARDAMOMlib.load_model_structure().pool_names.copy(),
    "pool_from": CARDAMOMlib.load_model_structure().pool_names.copy()
}

data_vars = dict()

# variables with one pool dimension and no time dimension
variables = [
    {"name": "start_values", "da": start_values_da, "unit": "g/m^2"}
]
for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "pool"],
        attrs={"units": d["unit"]}
)

# variables with two pool dimensions and one time dimension
variables = [
    {"name": "Bs", "da": Bs_da, "unit": "g/m^2/d"}
]
for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "time", "pool_to", "pool_from"],
        attrs={"units": d["unit"]}
)

# variables with one pool and one time dimension
variables = [
    {"name": "xs", "da": xs_da, "unit": "g/m^2"},
    {"name": "solution", "da": solution_da, "unit": "g/m^2"},
    {"name": "xs_abs_err", "da": xs_abs_err_da, "unit": "g/m^2"},
    {"name": "xs_rel_err", "da": xs_rel_err_da, "unit": "%"},

    {"name": "mean_pool_age_vector", "da": mean_pool_age_vector_da/(31*12), "unit": "yr"},
    {"name": "pool_age_moment_vector_2", "da": pool_age_moment_vector_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "pool_age_sd_vector", "da": pool_age_sd_vector_da/(31*12), "unit": "yr"},

#    {"name": "pool_age_median", "da": pool_age_median_da/(31*12), "unit": "yr"},
#    {"name": "pool_age_quantile_05", "da": pool_age_quantile_05_da/(31*12), "unit": "yr"},
#    {"name": "pool_age_quantile_95", "da": pool_age_quantile_95_da/(31*12), "unit": "yr"},

]
if model_type == "continuous":
    variables += [
        {"name": "us", "da": us_da, "unit": "g/m^2/d"},
        {"name": "external_output_vector", "da": external_output_vector_da, "unit": "g/m^2/d"}
    ]
elif model_type == "discrete":
    # accumulated over time step
    variables += [
        {"name": "Us", "da": Us_da, "unit": "g/m^2"}, 
        {"name": "acc_net_external_output_vector", "da": acc_net_external_output_vector_da, "unit": "g/m^2"}
    ]

for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "time", "pool"],
        attrs={"units": d["unit"]}
)

# variables with no pool dimension and one time dimension
variables = [
    {"name": "mean_system_age", "da": mean_system_age_da/(31*12), "unit": "yr"},
    {"name": "system_age_moment_2", "da": system_age_moment_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "system_age_sd", "da": system_age_sd_da/(31*12), "unit": "yr"},

#    {"name": "system_age_median", "da": system_age_median_da/(31*12), "unit": "yr"},
#    {"name": "system_age_quantile_05", "da": system_age_quantile_05_da/(31*12), "unit": "yr"},
#    {"name": "system_age_quantile_95", "da": system_age_quantile_95_da/(31*12), "unit": "yr"},

    {"name": "mean_btt", "da": mean_btt_da/(31*12), "unit": "yr"},
    {"name": "btt_moment_2", "da": btt_moment_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "btt_sd", "da": btt_sd_da/(31*12), "unit": "yr"},
    {"name": "btt_moment_3", "da": btt_moment_3_da/(31*12)**3, "unit": "yr^3"},
    {"name": "btt_moment_4", "da": btt_moment_4_da/(31*12)**4, "unit": "yr^4"},
    {"name": "btt_skewness", "da": btt_skewness_da, "unit": ""},
    {"name": "btt_kurtosis", "da": btt_kurtosis_da, "unit": ""},
    
    {"name": "btt_median", "da": btt_median_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_05", "da": btt_quantile_05_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_95", "da": btt_quantile_95_da/(31*12), "unit": "yr"},
    
    {"name": "GPP", "da": GPPs_da, "unit": "g/m^2"},
]
   
for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
        dims=["lat", "lon", "prob", "time"],
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
## compression takes forever, better zip it "manually" afterwards
#
##comp_dict = {"zlib": True, "complevel": 9}
##encoding = {var: comp_dict for var in ds.data_vars}
#
def delayed_to_netcdf(prob, netCDF_filename, compute=False):
    ds_sub = ds.sel(prob=[prob])
    ds_sub.to_netcdf(
        netCDF_filename,
#        mode="w",
#        encoding=encoding,
        compute=compute
    )
    del ds_sub

arr = ds.prob
arr
# +
# %%time

for prob in tqdm(arr):
    netCDF_filename = project_path.joinpath(netCDF_filestem + "_%05d.nc" % prob)
    print(netCDF_filename)
    delayed_to_netcdf(prob, netCDF_filename, compute=True)
    
# -
# ## OR (parallel, but seems not to be working for unknown reasons)
probs, datasets = zip(*ds.groupby("prob", squeeze=False))
paths = [project_path.joinpath(netCDF_filestem + "_%05d.nc" % prob) for prob in probs]
del_obj = xr.save_mfdataset(datasets, paths, compute=False)
del_obj
# %%time
del_obj.compute()



# ## Add GPP to files afterwards
#
# This is necessary because in the first versions of the data GPP was not included. The variable names of the input, gpp_to_xxx, are actually NPP.
#
# **create tmp folder first!**

for prob in tqdm(arr):
    old_netCDF_filename = project_path.joinpath(netCDF_filestem + "_%05d.nc" % prob)
    old_ds = xr.open_dataset(old_netCDF_filename)

    new_ds = old_ds.assign({"GPP": ds["GPP"].sel(prob=[prob])})
    new_netCDF_filename = project_path.joinpath("tmp").joinpath(netCDF_filestem + "_%05d.nc" % prob)
    new_ds.to_netcdf(new_netCDF_filename, mode="w")
    old_ds.close()
    new_ds.close()
    print("written", new_netCDF_filename)

# **Then move the files from tmp one folder up by hand.**



# ## Stub of a kurtosis and skewness computation of the backward transit time

ds = xr.open_dataset("/home/data/CARDAMOM/Greg_2020_10_26/monthly_output/discrete/sol_acc_age_btt_00011.nc")
ds

import matplotlib.pyplot as plt
plt.plot(ds.btt_kurtosis.isel(prob=0).mean(dim=["lat", "lon"]))

plt.plot(ds.btt_skewness.isel(prob=0).mean(dim=["lat", "lon"]))


