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

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
output_path = data_path.joinpath("output")

project_path = output_path.joinpath("solve_ivp_0000-0003_check_success")
netCDF_file = "output.nc"
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
    "prob": slice(0, 4, 1),
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
mean_age_vector_da = data_da[:, :, :, :, 1, :]
age_moment_vector_2_da = data_da[:, :, :, :, 2, :]

# compute backward transit time moments
external_output_vector_da = da.from_zarr(str(project_path.joinpath("external_output_vector")))
mean_btt_da = (external_output_vector_da * mean_age_vector_da).sum(-1) / external_output_vector_da.sum(-1)
btt_moment_2_da = (external_output_vector_da * age_moment_vector_2_da).sum(-1) / external_output_vector_da.sum(-1)

# backward_transit time quantiles
btt_median_da = da.from_zarr(str(project_path.joinpath("btt_median")))
btt_quantile_05_da = da.from_zarr(str(project_path.joinpath("btt_quantile_05")))
btt_quantile_95_da = da.from_zarr(str(project_path.joinpath("btt_quantile_95")))

# +
coords={
    "lat": lats_da.reshape(-1)[slices["lat"]],
    "lon": lons_da.reshape(-1)[slices["lon"]],
    "prob": probs_da.reshape(-1)[slices["prob"]],
    "time": times_da.reshape(-1),
    "pool": np.arange(nr_pools)
}

data_vars = dict()

# variables with pool dimension
variables = [
    {"name": "xs", "da": xs_da, "unit": "g/m^2"},
    {"name": "solution", "da": solution_da, "unit": "g/m^2"},
    {"name": "xs_abs_err", "da": xs_abs_err_da, "unit": "g/m^2"},
    {"name": "xs_rel_err", "da": xs_rel_err_da, "unit": "%"},
    {"name": "mean_age_vector", "da": mean_age_vector_da/(31*12), "unit": "yr"},
    {"name": "age_moment_vector_2", "da": age_moment_vector_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "external_output_vector", "da": external_output_vector_da/(31*12), "unit": "g/m^2/yr"}
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
    {"name": "mean_btt", "da": mean_btt_da/(31*12), "unit": "yr"},
    {"name": "btt_moment_2", "da": btt_moment_2_da/(31*12)**2, "unit": "yr^2"},
    {"name": "btt_median", "da" btt_median_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_05", btt_quantile_05_da/(31*12), "unit": "yr"},
    {"name": "btt_quantile_95", btt_quantile_95_da/(31*12), "unit": "yr"},
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

ds.to_netcdf(project_path.joinpath(netCDF_file), compute=True)
# -


