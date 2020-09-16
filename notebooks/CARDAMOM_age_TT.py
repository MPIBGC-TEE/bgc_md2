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
import numpy as np
from sympy import symbols

from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD
# -

pwc_mr_fd_archive = 'pwc_mr_fd/'
age_output_dir = 'age_output/'

ds_pwc_mr_global = xr.open_mfdataset(pwc_mr_fd_archive + '*.nc')
ds_pwc_mr_global

ds_pwc_mr_single = ds_pwc_mr_global.isel(ens=0, lat=0, lon=0)
ds_pwc_mr_single

data_times = ds_pwc_mr_single.times.data.compute()
start_values = ds_pwc_mr_single.start_values.data.compute()
Bs = ds_pwc_mr_single.Bs.data[:-1].compute()
us = ds_pwc_mr_single.us.data[:-1].compute()

mr = PWCModelRunFD.from_Bs_and_us(
    symbols('t'),
    data_times,
    start_values,
    Bs,
    us
)

# +
soln = mr.solve()

# mean age vector in years
mav = mr.age_moment_vector(1)/365.25

# mean system age in years
msa = mr.system_age_moment(1)/365.25

# age standard deviation vector in years
asdv = np.sqrt(mr.age_moment_vector(2))/365.25

# system age standard deviation in years
sasd = np.sqrt(mr.system_age_moment(2))/365.25

# +
coords_time = ds_pwc_mr_single.time
coords_pool = ds_pwc_mr_single.pool
data_vars = dict()

data_vars['soln'] = xr.DataArray(
    data=soln,
    dims=['time', 'pool'],
    coords={'time': coords_time, 'pool': coords_pool},
    attrs={'units': ds_pwc_mr_single.start_values.attrs['units']}
)

data_vars['mean_age_vector'] = xr.DataArray(
    data=mav,
    dims=['time', 'pool'],
    coords={'time': coords_time, 'pool': coords_pool},
    attrs={'units': 'yr'}
)

data_vars['mean_system_age'] = xr.DataArray(
    data=msa,
    dims=['time'],
    coords={'time': coords_time},
    attrs={'units': 'yr'}
)

data_vars['age_standard_deviation_vector'] = xr.DataArray(
    data=asdv,
    dims=['time', 'pool'],
    coords={'time': coords_time, 'pool': coords_pool},
    attrs={'units': 'yr'}
)

data_vars['system_age_standard_deviation'] = xr.DataArray(
    data=msa,
    dims=['time'],
    coords={'time': coords_time},
    attrs={'units': 'yr'}
)

ds_age = xr.Dataset(
    coords={'time': coords_time, 'pool': coords_pool},
    data_vars=data_vars
)
ds_age
# -

comp_dict = {'zlib': True, 'complevel': 9}
encoding = {var: comp_dict for var in ds_age.data_vars}
ds_age.to_netcdf(age_output_dir + 'ds_age.nc', encoding=encoding)


