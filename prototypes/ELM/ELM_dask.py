# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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
# %matplotlib inline

from dask.distributed import Client, LocalCluster
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from getpass import getuser

from bgc_md2.models.ELM import ELMlib
#import importlib
#importlib.reload(ELMlib)

# + endofcell="--"
port_dict = {
    'mm':8789,
#    'hmetzler':8790, # change at will
    'hmetzler':8888, # change at will
    'cs':8791        # change at will
}
my_user_name = getuser()
print(my_user_name)

my_port = port_dict[my_user_name]
print(my_port)

my_cluster = LocalCluster(dashboard_address='localhost:'+str(my_port))

# -

Client(my_cluster)
# --

ELMDataDir = "/home/hmetzler/SOIL-R/Manuscripts/Berkeley/2019/Data/"
runID = "14C_transient_holger_fire.2x2_small"
fn = runID + ".nc"
ds = xr.open_dataset(Path(ELMDataDir).joinpath(runID + ".nc"))
ds




ds_depth = xr.open_dataset(Path(ELMDataDir).joinpath('DZSOI.nc'))

parameter_set = ELMlib.load_parameter_set(
    nstep       = 1,
    ds_depth    = ds_depth
)
parameter_set

ds_single_site = ds.isel(lat=0, lon=0)
ds_mr_pwc_fd = ELMlib.compute_ds_pwc_mr_fd(ds_single_site, parameter_set)
ds_mr_pwc_fd



#lat, lon = ds.coords['lat'], ds.coords['lon']
lat, lon = ds['lat'][:], ds['lon'][:]
lat_indices, lon_indices = np.meshgrid(
    range(len(lat)),
    range(len(lon)),
    indexing='ij'
)

lats, lons = np.meshgrid(lat, lon, indexing='ij')
df_pd = pd.DataFrame({
    'cell_nr': range(len(lat)*len(lon)),
    'lat_index': lat_indices.flatten(),
    'lon_index': lon_indices.flatten(),
    'lat': lats.flatten(),
    'lon': lons.flatten()
})
df_pd

import dask.dataframe as dask_df

df_dask = dask_df.from_pandas(df_pd, npartitions=4)
df_dask

parameter_set = ELMlib.load_parameter_set(
    time_shift  = -198*365,
    nstep       = 10
)


# +
def func(line):
    #pass  
    location_dict = {
        'cell_nr':   int(line.cell_nr),
        'lat_index': int(line.lat_index),
        'lon_index': int(line.lon_index),    
    }
    
    #print(line, flush=True)
    cell_nr, log, xs_12C_data, us_12C_data, rs_12C_data= ELMlib.load_model_12C_data(parameter_set, location_dict)
    #return (1, 'log', [1,2,3], [4,5,6], [7,8,9])
    return cell_nr, log, xs_12C_data, us_12C_data, rs_12C_data

df_dask_2 = df_dask.apply(func, axis=1, meta=('A', 'object'))
# -

df_dask_2.compute()
type(df_dask_2)

df_dask_2

# +
#list(df_dask_2)
# -

pd.DataFrame(list(df_dask_2), columns=('cell_nr', 'log', 'xs_12C_data', 'us_12C_data', 'rs_12C_data'))


