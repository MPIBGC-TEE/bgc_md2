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

# The original ncl script `mm_mkinitialmatrices_A_C_1901_to_2004.ncl`
# only uses data of the first day of each year to create only one A per year.
# We produce a completely time dependent A that changes daily but do not 
# store it but use it directly in the computation of B.
import sys
from getpass import getuser
import xarray as xr
import numpy as np
import dask.array
import zarr as zr
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import shutil
import matplotlib.pyplot as plt
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from time import time
# %load_ext autoreload 
# %autoreload 2


from bgc_md2.sitespecificHelpers import get_client, getCluster



if __name__ == '__main__':
    # for reapeated %loadin in ipython we
    # first check if a cluster is already running
    if 'cluster' not in dir():
        cluster = getCluster()
        client = Client(cluster)

    out_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
    out_path = Path(out_dir)
    # for plotting purposese of the original cable output
    # load a singe file
    # syds=cH.single_year_ds(out_path.joinpath('out_ncar_1901_ndep.nc'))
    zarr_sub_dir_name = 'zarr'
    zarr_dir_path= out_path.joinpath(zarr_sub_dir_name)

    # a dictionary of all the (output converted to zarr)
    dad=cH.load_or_make_cable_dask_array_dict(
            out_path,
            rm=False
    )

    input()
    # create the storage for valid combinations on disk
    B_val, u_val ,x0_val = cH.load_or_make_valid_B_u_x0(
        out_path,
        zarr_sub_dir_name,
        rm=False
    )
    #sl=slice(0,320)
    #sl=slice(0,32)
    #B_val, u_val ,x0_val = cH.load_or_make_valid_B_u_x0_slice(
    #    out_path,
    #    zarr_sub_dir_name,
    #    sl,
    #    rm=False
    #)

    # for faster testing use only it_max time_steps
    #it_max  = 2000
    it_max  = B_val.shape[0]
    time_slice=slice(0,it_max)
    B_rt,u_rt=(
        arr[time_slice,...]
        for arr in (B_val,u_val)
    )
    times = dad['time'][time_slice].astype(np.float64)
    SOL_val = dask.array.blockwise(cH.valid_trajectory,'ijk',x0_val,'jk', times, 'i', B_rt, 'ijjk',u_rt,'ijk',dtype='f8')

    #SOL_val[...,0].compute()
    dt=(times[1]-times[0]).compute()
    #age_bin_indices = np.array([0,1,5,10,20])
    age_bin_indices = np.array([ i for i in range(20)])
    a_s = np.array([dt * it for it in age_bin_indices])
    def start_age_densities(x0,a):
       #return x0 if a==0 else x0*0.0
       #return x0 if a==0 else np.zeros_like(x0)
       return  (np.exp(-a)  if a>=0 else 0)* x0

    pool_age_density_val = dask.array.blockwise(
        cH.pool_age_density_val,'lijk',
        x0_val,                 'jk',
        times,                  'i',
        age_bin_indices,        'l',
        B_rt,                   'ijjk',
        u_rt,                   'ijk',
        start_age_densities_of_x0_and_a=start_age_densities,
        dtype='f8'
    )
    #age_moment_vector_val = dask.array.blockwise(
    #    cH.pool_age_density_val,'lijk',
    #    x0_val,                 'jk',
    #    times,                  'i',
    #    age_bin_indices,        'l',
    #    B_rt,                   'ijjk',
    #    u_rt,                   'ijk',
    #    start_age_densities_of_x0_and_a=start_age_densities,
    #    dtype='f8'
    #)
    #res = pool_age_density_val[...,0].compute() 
    suffix = '_slice_' + str(sl.start) + '_' + str(sl.stop)
    #for tup in zip((SOL_val, pool_age_density_val),('SOL_val','pool_age_density_val')):
    for tup in zip((pool_age_density_val,),('pool_age_density_val',)):
        arr,name = tup
        cH.batchwise_to_zarr(
            arr[...,sl],
            str(
                zarr_dir_path.joinpath(
                    name + suffix
                )
             ),
             rm=True
        )




    #def makeTup(i_c):
    #    it_max=1000
    #    x0=x0_val[:,i_c]#.reshape(9)
    #    times=dad['time'][:it_max]
    #    Bs=[B_val[i,:,:,i_c]+np.eye(9) for i in range(it_max)]
    #    us=[u_val[i,:,i_c] for i in range(it_max)]
    #    return (x0,times,Bs,us)

    #val_tups = [makeTup(i_c) for i_c in range(sl.start,sl.stop)]

    #def dmr_from_tup(tup):
    #    x0,times,bs,us = tup
    #    return  DMR.from_Bs_and_net_Us(x0,times,bs,us)

    #def sol_arr_from_tup(tup):
    #    x0,times,bs,us = tup
    #    dmr = DMR.from_Bs_and_net_Us(x0,times,bs,us)
    #    return dmr.solve()

    #dmr0 = makeTup(0)
    #futures = client.map(dmr_from_tup,val_tups,batch_size=8)
    #batch_size = 8
    #slices=batchSlices(sl.stop,batch_size)
    #s=slices[0]
    #futures = client.map(sol_arr_from_tup,val_tups[s],batch_size=8)
    #for s in slices:
    #    arr.list = client
