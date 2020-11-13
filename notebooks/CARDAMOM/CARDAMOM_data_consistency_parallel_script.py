#!/usr/bin/env python
# coding: utf-8

# In[1]:


from dask.distributed import Client

import xarray as xr
import numpy as np

import importlib
import bgc_md2
from bgc_md2.models.CARDAMOM import CARDAMOMlib

from dask import delayed, compute
from dask.diagnostics import Profiler, ResourceProfiler
from dask.distributed import Client, LocalCluster
from getpass import getuser

from time import sleep

#importlib.reload(bgc_md2.models.CARDAMOM.CARDAMOMlib)

# In[2]:

def run_code():
    
    
    
    pwc_mr_fd_archive = '/home/data/CARDAMOM/output/pwc_mr_fd_archive/'
    data_folder = "/home/data/CARDAMOM/"  # matagorda, antakya, raka..
    filestem = "Greg_2020_10_26"
    ds = xr.open_mfdataset(data_folder + filestem + "/SUM*.nc")

#    print(ds)
    
    comp_dict = {'zlib': True, 'complevel': 9}
    
    # there is no multi-dimensional 'groupby' in xarray data structures
    def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
        if len(groupby) == 1:
            res = dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
            return res
        else:
            return dataset.groupby(groupby[0]).apply(
                nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
            )
   
    dss = []
    def func_data_consistency(ds_single):
        print(ds_single.lat.data, ds_single.lon.data, ds_single.prob.data)
        mdo = CARDAMOMlib.load_mdo_greg(ds_single)
        abs_err, rel_err = mdo.check_data_consistency()
        
        data_vars = dict()
        data_vars['abs_err'] = xr.DataArray(
            data=abs_err.data,
        )
        data_vars['rel_err'] = xr.DataArray(
            data=rel_err.data.filled(fill_value=np.nan),
        )
    
        ds_res = xr.Dataset(
            data_vars=data_vars,
        )

        return ds_res
    
    def func_chunk(chunk_ds):
#        print('\nStarting chunk:', chunk_ds, '\n')
        res = nested_groupby_apply(chunk_ds, ['lat', 'lon', 'prob'], func_data_consistency)
        print('chunk finished', flush=True)

        return res

    chunk_dict = {"lat": 1, "lon": 1}

    ds_sub = ds.isel(
        lat=slice(0, 1, 1),
        lon=slice(0, 1, 1),
        prob=slice(0, 5, 1)
    ).chunk(chunk_dict)
    
    
    fake_data = np.zeros((
        len(ds_sub.lat),
        len(ds_sub.lon),
        len(ds_sub.prob)
    ))

    fake_array = xr.DataArray(
        data=fake_data,
        dims=['lat', 'lon', 'prob']
    )

    fake_coords = {
        'lat': ds_sub.lat.data,
        'lon': ds_sub.lon.data,
        'prob': ds_sub.prob.data
    }

    fake_ds = xr.Dataset(
        data_vars={
            'abs_err': fake_array, 
            'rel_err': fake_array
        },
       coords=fake_coords
    ).chunk(chunk_dict)
 
#    print(fake_ds)
    
    ds_data_consistency = xr.map_blocks(func_chunk, ds_sub, template=fake_ds)
    ds_data_consistency.compute(nr_workers=1)
    
    print(' ------------------ -------------- ')    
    print(ds_data_consistency)
    

if __name__ == "__main__":
    port_dict = {
        'mm':8789,
        'hmetzler':8790, # change at will
    #    'hmetzler':8888, # change at will
        'cs':8791        # change at will
    }
    my_user_name = getuser()
    print(my_user_name)
    
    my_port = port_dict[my_user_name]
    print(my_port)
    
    # dasboard needs a different port for accessing it remotely
    my_cluster = LocalCluster(
        dashboard_address='localhost:'+str(my_port+5),
        nr_workers=30
    )
    
    Client(my_cluster)


    run_code()


