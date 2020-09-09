# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
from dask.distributed import Client

import xarray as xr
import numpy as np

import importlib
import CARDAMOMlib

#xr.set_options(display_style='html')
# -


#importlib.reload(CARDAMOMlib)

def main():

    #client = Client(n_workers=2, threads_per_worker=1, memory_limit="3GB")
    #client
    
    
    # +
    data_folder = "/home/hmetzler/Desktop/CARDAMOM/" # local
    #data_folder = "/home/data/CARDAMOM/"  # matagorda
    
    filestem = "cardamom_for_holger_10_ensembles"
    #filestem = "cardamom_for_holger"
    
    #chunk_dict = {"ens": 20}
    chunk_dict = {"ens": 2}
    
    #filestem = "cardamom_for_holger"
    #chunk_dict = {"ens": 100}
    ds = xr.open_dataset(data_folder + filestem + ".nc")#.isel(
    #    ens=slice(None, 6),
    #    time=slice(None, 5)
    #)
    ds = ds.chunk(chunk_dict)
    ds
    
    
    # -
    
    
    # there is no multi-dimensional 'groupby' in xarray data structures
    def nested_groupby_apply(dataset, groupby, apply_fn, **kwargs):
        if len(groupby) == 1:
            return dataset.groupby(groupby[0]).apply(apply_fn, **kwargs)
        else:
            return dataset.groupby(groupby[0]).apply(
                nested_groupby_apply, groupby=groupby[1:], apply_fn=apply_fn, **kwargs
            )
    
    
    comp_dict = {'zlib': True, 'complevel': 9}
    
    # +
    # %%time
    
    # compute in parallel the model runs and save them to ds_mrs in netcdf format
    
    small_template = xr.Dataset(
        data_vars = {
            'x': xr.DataArray(
                data=np.ndarray(dtype=float, shape=(len(ds.ens.data),)),
                dims=['ens']
            )
        }
    ).chunk(chunk_dict)
    
    def func(single_site_ds):
        res = CARDAMOMlib.compute_pwc_mr_fd_ds(single_site_ds)
        return res
    
    def func_chunk(chunk_ds):
        print('\nChunk start:', chunk_ds.ens.data[0], '\n')
        res = nested_groupby_apply(chunk_ds, ['ens', 'lat', 'lon'], func)
        
        filename = filestem + "_{:03d}-{:03d}.nc".format(res.ens.data[0], res.ens.data[-1])
        encoding = {var: comp_dict for var in res.data_vars}
        res.to_netcdf(filename, encoding=encoding)
        print(res)
        del res
    
        return xr.Dataset(
            data_vars={
                'x': xr.DataArray(
                    data=np.zeros((chunk_dict['ens'],)),
                    dims=['ens']
                )
            }
        )
    
    _ = xr.map_blocks(func_chunk, ds, template=small_template).compute()
    # -
    
    
    ds.close()
    del ds


if __name__ == "__main__":
    client = Client(n_workers=2, threads_per_worker=2, memory_limit="1GB")
    main()
