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

# +
import xarray as xr
import numpy as np
from sympy import symbols

from dask.distributed import Client
from dask import delayed, compute

from CompartmentalSystems.pwc_model_run_fd import PWCModelRunFD
# -

client = Client(n_workers=10, threads_per_worker=1, memory_limit="3GB")
client

pwc_mr_fd_archive = '/home/data/CARDAMOM/output/pwc_mr_fd_archive/'
age_output_dir = 'age_output/'
#filestem = "cardamom_for_holger_10_ensembles"
filestem = "cardamom_for_holger"
comp_dict = {'zlib': True, 'complevel': 9}

ds_pwc_mr_global = xr.open_mfdataset(pwc_mr_fd_archive + '*.nc')
ds_pwc_mr_global#.compute()

#ds_pwc_mr_global.where(np.isnan(ds_pwc_mr_global.start_values)==False, drop=True).compute()
ds_pwc_mr_global_clean = ds_pwc_mr_global.where(ds_pwc_mr_global.log.compute()=='', drop=True)
ds_pwc_mr_global_clean


# there is no multi-dimensional 'groupby' in xarray data structures
def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
    if len(groupby) == 1:
        res = dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
#        print(res.start_values)
        return res
    else:
        return dataset.groupby(groupby[0]).apply(
            nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
        )


def func(ds_single_site):
    data_times = ds_single_site.times.data#.compute()
    start_values = ds_single_site.start_values.data#.compute()
    Bs = ds_single_site.Bs.data[:-1]#.compute()
    us = ds_single_site.us.data[:-1]#.compute()
    
    mr = PWCModelRunFD.from_Bs_and_us(
        symbols('t'),
        data_times,
        start_values,
        Bs,
        us
    )
    
    coords_time = ds_single_site.time
    coords_pool = ds_single_site.pool
    data_vars = dict()

    soln = mr.solve()
    data_vars['soln'] = xr.DataArray(
        data=soln,
        dims=['time', 'pool'],
        coords={'time': coords_time, 'pool': coords_pool},
        attrs={'units': ds_single_site.start_values.attrs['units']}
    )

    mav = mr.age_moment_vector(1) / 365.25
    data_vars['mean_age_vector'] = xr.DataArray(
        data=mav,
        dims=['time', 'pool'],
        coords={'time': coords_time, 'pool': coords_pool},
        attrs={'units': 'yr'}
    )

    msa = mr.system_age_moment(1) / 365.25
    data_vars['mean_system_age'] = xr.DataArray(
        data=msa,
        dims=['time'],
        coords={'time': coords_time},
        attrs={'units': 'yr'}
    )

    asdv = np.sqrt(mr.age_moment_vector(2)) / 365.25
    data_vars['age_standard_deviation_vector'] = xr.DataArray(
        data=asdv,
        dims=['time', 'pool'],
        coords={'time': coords_time, 'pool': coords_pool},
        attrs={'units': 'yr'}
    )

    sasd = np.sqrt(mr.system_age_moment(2)) / 365.25
    data_vars['system_age_standard_deviation'] = xr.DataArray(
        data=sasd,
        dims=['time'],
        coords={'time': coords_time},
        attrs={'units': 'yr'}
    )

    res = xr.Dataset(
        coords={'time': coords_time, 'pool': coords_pool},
        data_vars=data_vars
    )    
    
    return res


# +
def func_chunk(ds_chunk):
    print('\nStarting chunk:', ds_chunk.ens.data[0], '-', ds_chunk.ens.data[-1], '\n')
    res = nested_groupby_apply(ds_chunk, ['ens', 'lat', 'lon'], func)
    
#    filename = 'output/' + filestem + "_{:03d}-{:03d}.nc".format(res.ens.data[0], res.ens.data[-1])
    filename = age_output_dir + filestem + "_{:03d}-{:03d}.nc".format(res.ens.data[0], res.ens.data[-1])
    encoding = {var: comp_dict for var in res.data_vars}
    res.to_netcdf(filename, encoding=encoding)
    print('\nFinished chunk:', ds_chunk.ens.data[0], '-', ds_chunk.ens.data[-1], '\n')
#    del res
#    gc.collect()

    return ds_chunk.ens.data[0]


# +
chunk_dict = {"ens": 10}
results = []
for index in range(0, len(ds_pwc_mr_global_clean.ens.data), chunk_dict['ens']):
    ds_chunk = ds_pwc_mr_global_clean.isel(ens=slice(index, index+chunk_dict['ens'], 1))
    print(ds_chunk)
    print()
    results.append(delayed(func_chunk)(ds_chunk))

delayed_results = delayed(results)

# +
# %%time

result = compute(delayed_results, scheduler='distributed', num_workers=10, memory_limit="3GB")
# -


