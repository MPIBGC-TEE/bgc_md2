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

from getpass import getuser
import xarray as xr
import numpy as np
import dask
import zarr
from dask.distributed import Client
from dask.distributed import LocalCluster
from pathlib import Path
import shutil
from sympy import Symbol,var
import matplotlib.pyplot as plt
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices
from time import time
#%load_ext autoreload 
#%autoreload 2


from bgc_md2.sitespecificHelpers import get_client, getCluster

def main():
    example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/'
    outDir = "output/new4"
    outpath = Path(example_run_dir).joinpath(outDir)
    # Allthough we will later use zarr arrays
    # we open one of the original netcdf output files
    # to get the meta information about fill values 
    # Since the file is small this is  fast
    syds=cH.single_year_ds(outpath.joinpath('out_ncar_1901_ndep.nc'))
    tk='_FillValue'
    ifv=syds.iveg.attrs[tk]
    ffv=syds.Cplant.attrs[tk]

    zarr_dir_path= outpath.joinpath('zarr')
    dad=cH.cable_dask_array_dict(zarr_dir_path)

    npool = sum([dad[name].shape[1] for name in ['Cplant','Csoil','Clitter']])
    s=dad['Cplant'].shape
    npatch =s[2]
    nland  =s[3]
    ntime  =s[0]

    for s in ('leaf','wood','fine_root','metabolic_lit','structural_lit','cwd','fast_soil','slow_soil','passive_soil'):
        var(s)
    stateVariableTuple=(leaf,fine_root,wood,metabolic_lit,structural_lit,cwd,fast_soil,slow_soil,passive_soil)
    npool=len(stateVariableTuple)
    var_arr_dict={
        leaf:dad['Cplant'][:,0,:,:],
        fine_root:dad['Cplant'][:,2,:,:],
        wood:dad['Cplant'][:,1,:,:],
        metabolic_lit : dad['Clitter'][:,0,:,:],
        structural_lit : dad['Clitter'][:,1,:,:],
        cwd : dad['Clitter'][:,2,:,:],
        fast_soil : dad['Csoil'][:,0,:,:],
        slow_soil : dad['Csoil'][:,1,:,:],
        passive_soil : dad['Csoil'][:,2,:,:]
    }
    X=dask.array.stack([var_arr_dict[var]for var in stateVariableTuple],1).rechunk((ntime,npool,npatch,1)) 

    X0=X[0,...]

    B_res = cH.reconstruct_B(
        ifv,
        ffv,
        dad['iveg'],
        dad['kplant'],
        dad['fromLeaftoL'],
        dad['fromRoottoL'],
        dad['fromWoodtoL'],
        dad['fromMettoS'],
        dad['fromStrtoS'],
        dad['fromCWDtoS'],
        dad['fromSOMtoSOM'],
        dad['xktemp'],
        dad['xkwater'],
        dad['xkNlimiting']
    )
    u_res = cH.reconstruct_u(
        ifv,
        ffv,
        dad['iveg'],
        dad['NPP'],
        dad['fracCalloc']
    )
    return (B_res,u_res,X0)

if __name__ == '__main__':
    cluster = getCluster()
    client = Client(cluster)
    
    example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/'
    outDir = "output/new4"
    outpath = Path(example_run_dir).joinpath(outDir)
    zarr_dir_path= outpath.joinpath('zarr')
    B_res,u_res,x0_res = main()
    #B_res.to_zarr(str(zarr_dir_path.joinpath('B')))
    #u_res.to_zarr(str(zarr_dir_path.joinpath('u')))
    #x0_res.to_zarr(str(zarr_dir_path.joinpath('x0')))
    cH.batchwise_to_zarr(B_res, str(zarr_dir_path.joinpath('B')))
    cH.batchwise_to_zarr(u_res, str(zarr_dir_path.joinpath('u')))
    cH.batchwise_to_zarr(X0, str(zarr_dir_path.joinpath('x0')))
