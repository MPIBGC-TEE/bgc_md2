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

# +
import xarray as xr
from pathlib import Path

from bgc_md2.models.ELM import ELMlib
# -

ELMDataDir = "/home/hmetzler/SOIL-R/Manuscripts/Berkeley/2019/Data/"
runID = "14C_transient_holger_fire.2x2_small"
fn = runID + ".nc"
ds = xr.open_dataset(Path(ELMDataDir).joinpath(runID + ".nc"))
ds

ds_single_site = ds.isel(lat=0, lon=0)
ds_single_site

# +
ds_depth = xr.open_dataset(Path(ELMDataDir).joinpath('DZSOI.nc'))

parameter_set = ELMlib.load_parameter_set(
    nstep       = 1,
    ds_depth    = ds_depth
)
# -

ds_single_site = ds.isel(lat=0, lon=0)
ds_single_site
mdo = ELMlib.load_mdo(ds_single_site, parameter_set)
mdo.check_data_consistency()

ds_monthly = ELMlib.resample_daily_to_monthly(ds)
ds_monthly

ds_monthly_single_site = ds_monthly.isel(lat=0, lon=0)
ds_monthly_single_site

# +
ds_depth = xr.open_dataset(Path(ELMDataDir).joinpath('DZSOI.nc'))

parameter_set = ELMlib.load_parameter_set(
    nstep       = 1,
    ds_depth    = ds_depth
)

mdo = ELMlib.load_mdo(ds_monthly_single_site, parameter_set)
mdo.check_data_consistency()
# -


