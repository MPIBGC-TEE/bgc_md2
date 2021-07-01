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
import bgc_md2.models.cable_all.cablePaths as cP
from bgc_md2.helper import batchSlices
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
from time import time

# %load_ext autoreload
# %autoreload 2
#
cable_out_path = (
    Path("/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4")
)

if __name__ == "__main__":
    import dask.array
    from dask.distributed import LocalCluster, Client

    if "cluster" not in dir():
        cluster = LocalCluster()

    client = Client(cluster)




syds = cH.single_year_ds(cable_out_path.joinpath("out_ncar_1901_ndep.nc"))
tk = "_FillValue"
ifv = syds.iveg.attrs[tk]
ffv = syds.Cplant.attrs[tk]

it_max = syds.Cplant.shape[0] 
#time_slice=slice(0,it_max) 
time_slice=slice(0,100) #could be used to significantly shorten the reconstruction times 

sl = slice(0, 32)
# sl = slice(None,None,None) #for the full grid
vcsp = cP.zarr_valid_landpoint_patch_combis_slice_path(cable_out_path)
zp=cP.zarr_path(cable_out_path)

dad = cH.load_or_make_cable_dask_array_dict(
    cable_out_path,
    cP.zarr_path(cable_out_path),
    rm=False,
    batch_size=64
)

B = cH.cache(
    zarr_dir_path=zp,
    name='B',
    arr=cH.reconstruct_B(
        ifv,
        ffv,
        dad["iveg"],
        dad["kplant"],
        dad["fromLeaftoL"],
        dad["fromRoottoL"],
        dad["fromWoodtoL"],
        dad["fromMettoS"],
        dad["fromStrtoS"],
        dad["fromCWDtoS"],
        dad["fromSOMtoSOM"],
        dad["xktemp"],
        dad["xkwater"],
        dad["xkNlimiting"],
    )
    # rm=True
)
u = cH.cache(
    zarr_dir_path=zp,
    name='u',
    arr=cH.reconstruct_u( ifv, ffv, dad["iveg"], dad["NPP"], dad["fracCalloc"]),
    # rm=True
)

x = cH.cache(
    zarr_dir_path=zp,
    name='x',
    arr=cH.reconstruct_x(
        dad["Cplant"],
        dad["Clitter"],
        dad["Csoil"]
    ),
    # rm=True
)

iveg = dad["iveg"]
cond_1 = (iveg != ifv).compute()
cond_2 = (dad["Csoil"][0, 0, :, :] != 0).compute()
nz = dask.array.nonzero(cond_1 * cond_2)

B_val = cH.cache(
    zarr_dir_path=vcsp,
    name='B',
    arr=cH.valid_combies_parallel(nz, B)[...,sl], 
    # rm=True
)
u_val = cH.cache(
    zarr_dir_path=vcsp,
    name='u',
    arr=cH.valid_combies_parallel(nz, u)[...,sl],
    # rm=True
)
x_val = cH.cache(
    zarr_dir_path=vcsp,
    name='x',
    arr=cH.valid_combies_parallel(nz, x)[...,sl],
    #rm=True
)
x0_val = x_val[0,...] 

B_rt, u_rt = (
    arr[time_slice,...]
    for arr in (B_val,u_val)
)
times = dad['time'][time_slice].astype(np.float64)

SOL_val = cH.cache(
    zarr_dir_path=vcsp,
    name='SOL_val',
    arr=dask.array.blockwise(cH.valid_trajectory,'ijk',x0_val,'jk', times, 'i', B_rt, 'ijjk',u_rt,'ijk',dtype='f8'),
    # rm=True
)
#cH.batchwise_to_zarr(
#    SOL_val,
#    str( zarr_dir_path.joinpath( 'SOL_val')),
#    rm=False,
#    batch_size=2
#)


