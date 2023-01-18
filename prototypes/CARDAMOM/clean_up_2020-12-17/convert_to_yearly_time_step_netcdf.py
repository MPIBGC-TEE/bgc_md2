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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Convert CARDAMOM netcdf data to yearly time steps
#
# This notebook loads the CARDAMOM netcdf data files and saves them as a single netcdf file in yearly time steps.

# +
import xarray as xr

from pathlib import Path
from tqdm import tqdm

from bgc_md2.models.CARDAMOM import CARDAMOMlib
from dask.distributed import Client
# -

my_cluster = CARDAMOMlib.prepare_cluster(n_workers=48)
Client(my_cluster)

nr_months = 12 # coarseness
delay_in_months = 0
#delay_in_months = 6

# +
data_path = Path("/home/data/CARDAMOM/Greg_2020_10_26/")
target_path = data_path.joinpath("yearly_%02d_ds.nc" % delay_in_months)

print("target_path", str(target_path))
ds = xr.open_mfdataset(str(data_path) + "/SUM*.nc")
ds
# -

yearly_coords_time = ds["time"].isel(time=slice(delay_in_months, None, nr_months))
yearly_coords_time

# +
stock_variable_names = ["c_finelitter", "c_foliar", "c_labile", "c_root", "c_som", "c_wood"]

data_vars = dict()
for variable_name, variable in ds.data_vars.items():
    if variable_name in stock_variable_names:
        yearly_stock_variable = variable.isel(time=slice(delay_in_months, None, nr_months))
        data_vars[variable_name] = yearly_stock_variable
    else:
        # flux variable
        yearly_flux_variable = variable.isel(time=slice(delay_in_months, None, 1)).shift(time=-1).coarsen( # first u is useless in CARDAMOM and ELM
            time=nr_months,
            boundary="pad",
            coord_func="min",
            keep_attrs=True
        ).mean().shift(time=1) # unit is gC/m2/d, create useless first u again
        yearly_flux_variable.coords["time"] = yearly_coords_time
        data_vars[variable_name] = yearly_flux_variable
        
yearly_ds = xr.Dataset(
    data_vars=data_vars,
    attrs={"delay_in_months": delay_in_months}
)
yearly_ds

# +
# %%time

print(target_path)
yearly_ds.to_netcdf(target_path, compute=True)
# -




check_ds = xr.open_dataset(target_path)
check_ds

#sub_ds = check_ds.isel(lat=9, lon=26, prob=0).compute()
sub_ds = check_ds.isel(lat=5, lon=23, prob=0).compute() # bad reconstruction
sub_ds

CARDAMOMlib.check_data_consistency(sub_ds.compute(), time_step_in_days=31*12)


