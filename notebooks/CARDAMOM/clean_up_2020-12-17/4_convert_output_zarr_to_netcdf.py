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

# # Convert output from zarr to a netCDF file

# +
import dask.array as da
import numpy as np
import xarray as xr

from pathlib import Path

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

time_resolution = "monthly"

# +
params = CARDAMOMlib.load_params(time_resolution)

data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
output_path = data_path.joinpath(data_path.joinpath(params["output_folder"]))

project_path = output_path.joinpath("discrete")
netCDF_file = "sol_acc_age_btt.nc"
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
    "prob": slice(0, 1, 1),
    "time": slice(0, None, 1) # don't change the time entry
}

# +
# dimensions: lat x lon x prob x time x order x pool
xs_da = da.from_zarr(str(project_path.joinpath("xs")))
data_da = da.from_zarr(str(project_path.joinpath("age_moment_vectors_up_to_2")))
solution_da = data_da[:, :, :, :, 0, :]

# compute absolute and relative errors of reconstructed solutions
xs_abs_err_da = np.abs(xs_da-solution_da)
xs_rel_err_da = xs_abs_err_da/xs_da * 100

# age moment vectors
mean_pool_age_vector_da = data_da[:, :, :, :, 1, :]
pool_age_moment_vector_2_da = data_da[:, :, :, :, 2, :]
pool_age_sd_vector_da = da.sqrt(pool_age_moment_vector_2_da)

# system age moments
mean_system_age_da = (solution_da * mean_pool_age_vector_da).sum(-1) / solution_da.sum(-1)
system_age_moment_2_da = (solution_da * pool_age_moment_vector_2_da).sum(-1) / solution_da.sum(-1)
system_age_sd_da = da.sqrt(system_age_moment_2_da)

# pool age quantiles quantiles
pool_age_median_da = da.from_zarr(str(project_path.joinpath("pool_age_median")))
pool_age_quantile_05_da = da.from_zarr(str(project_path.joinpath("pool_age_quantile_05")))
pool_age_quantile_95_da = da.from_zarr(str(project_path.joinpath("pool_age_quantile_95")))

# system age quantiles quantiles
system_age_median_da = da.from_zarr(str(project_path.joinpath("system_age_median")))
system_age_quantile_05_da = da.from_zarr(str(project_path.joinpath("system_age_quantile_05")))
system_age_quantile_95_da = da.from_zarr(str(project_path.joinpath("system_age_quantile_95")))

# compute backward transit time moments
#external_output_vector_da = da.from_zarr(str(project_path.joinpath("external_output_vector")))
acc_net_external_output_vector_da = da.from_zarr(str(project_path.joinpath("acc_net_external_output_vector")))

#mean_btt_da = (external_output_vector_da * mean_age_vector_da).sum(-1) / external_output_vector_da.sum(-1)
#btt_moment_2_da = (external_output_vector_da * age_moment_vector_2_da).sum(-1) / external_output_vector_da.sum(-1)
mean_btt_da = (acc_net_external_output_vector_da * mean_pool_age_vector_da).sum(-1) / acc_net_external_output_vector_da.sum(-1)
btt_moment_2_da = (acc_net_external_output_vector_da * pool_age_moment_vector_2_da).sum(-1) /acc_net_external_output_vector_da.sum(-1)
btt_sd_da = np.sqrt(btt_moment_2_da)

# backward transit time quantiles
btt_median_da = da.from_zarr(str(project_path.joinpath("btt_median")))
btt_quantile_05_da = da.from_zarr(str(project_path.joinpath("btt_quantile_05")))
btt_quantile_95_da = da.from_zarr(str(project_path.joinpath("btt_quantile_95")))

# +
coords = {
    "lat": lats_da.reshape(-1)[slices["lat"]],
    "lon": lons_da.reshape(-1)[slices["lon"]],
    "prob": probs_da.reshape(-1)[slices["prob"]],
    "time": times_da.reshape(-1),
    "pool": CARDAMOMlib.load_model_structure().pool_names
}

data_vars = dict()

# variables with pool dimension
variables = [
    {"name": "xs", "da": xs_da, "unit": "g/m^2"},
    {"name": "solution", "da": solution_da, "unit": "g/m^2"},
    {"name": "xs_abs_err", "da": xs_abs_err_da, "unit": "g/m^2"},
    {"name": "xs_rel_err", "da": xs_rel_err_da, "unit": "%"},

    {"name": "mean_pool_age_vector", "da": mean_pool_age_vector_da/(31*12), "unit": "yr"},
    {"name": "pool_age_moment_vector_2", "da": pool_age_moment_vector_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "pool_age_sd_vector", "da": pool_age_sd_vector_da/(31*12), "unit": "yr"},

    {"name": "pool_age_median", "da": pool_age_median_da/(31*12), "unit": "yr"},
    {"name": "pool_age_quantile_05", "da": pool_age_quantile_05_da/(31*12), "unit": "yr"},
    {"name": "pool_age_quantile_95", "da": pool_age_quantile_95_da/(31*12), "unit": "yr"},

#    {"name": "external_output_vector", "da": external_output_vector_da*(31*12), "unit": "g/m^2/yr"}
    {"name": "acc_net_external_output_vector", "da": acc_net_external_output_vector_da*(31*12), "unit": "g/m^2/yr"}
]
for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
#        coords=coords,
        dims=["lat", "lon", "prob", "time", "pool"],
        attrs={"units": d["unit"]}
)

# variables with no pool dimension
variables = [
    {"name": "mean_system_age", "da": mean_system_age_da/(31*12), "unit": "yr"},
    {"name": "system_age_moment_2", "da": system_age_moment_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "system_age_sd", "da": system_age_sd_da/(31*12), "unit": "yr"},

    {"name": "system_age_median", "da": system_age_median_da/(31*12), "unit": "yr"},
    {"name": "system_age_quantile_05", "da": system_age_quantile_05_da/(31*12), "unit": "yr"},
    {"name": "system_age_quantile_95", "da": system_age_quantile_95_da/(31*12), "unit": "yr"},

    {"name": "mean_btt", "da": mean_btt_da/(31*12), "unit": "yr"},
    {"name": "btt_moment_2", "da": btt_moment_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "btt_sd", "da": btt_sd_da/(31*12), "unit": "yr"},
    
    {"name": "btt_median", "da": btt_median_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_05", "da": btt_quantile_05_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_95", "da": btt_quantile_95_da/(31*12), "unit": "yr"},
]
for d in variables:
    data_vars[d["name"]] = xr.DataArray(
        data=d["da"][slices["lat"], slices["lon"], slices["prob"]],
#        coords=coords,
        dims=["lat", "lon", "prob", "time"],
        attrs={"units": d["unit"]}
)

ds = xr.Dataset(
    data_vars=data_vars,
    coords=coords
)
ds

# +
# %%time

print(project_path.joinpath(netCDF_file))
ds.to_netcdf(project_path.joinpath(netCDF_file), compute=True)
# -
