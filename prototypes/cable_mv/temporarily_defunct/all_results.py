# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
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

if __name__ == '__main__':
    import dask.array
    from dask.distributed import LocalCluster,Client
    if 'cluster' not in dir():
        cluster = LocalCluster()

    client = Client(cluster)

# +
try:
    from ports.server_helpers import print_commands
    print_commands(cluster,local_port=8880)

except ImportError as e:
    pass  # module doesn't exist,dont make a fuss
# -

    # chose the cable output directory you want to work with
    out_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
    out_path = Path(out_dir)

    # option:
    # for plotting purposese of the original cable output
    # you could load a singe file
    # syds=cH.single_year_ds(out_path.joinpath('out_ncar_1901_ndep.nc'))

    # A dictionary of DaskArrays of all the cable output variables.
    # If the directory 'zarr_dir_path' already exists (for instance
    # because this function call has been made before) The containing
    # zarr arrays will be read as delayed objects.  If the directory
    # with the zarr files does not exist I will be created.  Every
    # variable will we rechunked and written to disk into a its own
    # zarr array.
    # The conversion of the original netcdf output to the properly chunked zarr
    # arrays speeds up the computations considerably and also enables and also
    # avoids issues with variables that do not fit into memory.
    # options:
    zarr_sub_dir_name = 'zarr'
    zarr_dir_path= out_path.joinpath(zarr_sub_dir_name)
    dad=cH.load_or_make_cable_dask_array_dict(
            out_path,
            zarr_dir_path,
            rm=False,
            batch_size=64
    )

    # Create the storage for valid combinations on disk.
    #
    # time:
    # On matagorda( 16 Cores and 32 GB this will take some hours
    # if the data is not present on disk already)
    #
    # As an intermediate step this will reconstruct arrays for matix B and the
    # influxes u in the shapes used by cable for its output, which has
    # separate dimensions 'patch' and 'landpoint'. Since not all landpoints have
    # all patches only some combinations are valid sources for a computation.
    # A distribution of work (parallelisation) with respect to the landpoints
    # would be unefficient due to unequal workloads and furthermore require the
    # active exclusion of fillvalues in the computation.
    # Therefore we first find all the valid ('landpoint','patch') combinations
    # and reorganize the arrays for B,u,x0 with one dimension enumerating
    # the combinations.
    # We chunk and write these arrays as zarr and can later use them
    # to distribute the work efficiently.
    #
    #
    #
    # If mistakes are discovered in the reconstructing function.
    # the arrays should be rewritten.(set the rm flag to true)
    #
    # options:
    # - It is possible to use only a part of the valid combinations (a slice)
    #   This is much faster and especially useful for testing purposes.
    #   (The whole grid contains about 20000 valid combinations)
    #   sl=slice(0,320)
    #   sl=slice(0,32)
    #   B_val, u_val ,x0_val = cH.load_or_make_valid_B_u_x0_slice(
    #       out_path,
    #       zarr_sub_dir_name,
    #       sl,
    #       rm=False
    #   )
    # - The number of chunks processed in parallel is governed by
    #   the parameter 'batch_size'. It is usually the number of cores or
    #   threads. However if a single chunk needs a lot of memory it
    #   might have to be reduced on computers with low memory
    #   (reference for this task: matagorda:32 / antakya: 128)

    B_val, u_val ,x0_val = cH.load_or_make_valid_B_u_x0(
        out_path,
        zarr_sub_dir_name,
        names=['B_val','u_val','x0_val'],
        rm=False,
        batch_size=32 #
    )

    # We now compute dependent varaibles on the grid:
    #
    # option:
    # As with the (patch,landpoint) combinations we can also narrow down the
    # scope of the computations to a time slice (0,it_max) This will speed up
    # the computations and is usefull for testing The number of timestep of
    # the complete run is about 37500 it_max  = 2000
    it_max  = B_val.shape[0]
    time_slice=slice(0,it_max)
    B_rt,u_rt=(
        arr[time_slice,...]
        for arr in (B_val,u_val)
    )
    times = dad['time'][time_slice].astype(np.float64)

    # The first result is the reconstructed solution (from the Bs us and x0s)
    # The result is not too large to be stored and can be used to check how
    # close the reconstruction actually is.
    SOL_val = dask.array.blockwise(cH.valid_trajectory,'ijk',x0_val,'jk', times, 'i', B_rt, 'ijjk',u_rt,'ijk',dtype='f8')
    cH.batchwise_to_zarr(
        SOL_val,
        str( zarr_dir_path.joinpath( 'SOL_val')),
        rm=False,
        batch_size=32 #
    )


SOL_val

# We can now check the solutions for consistency by
# 1.) plotting the solutions
# 1.) comparing them to the solutions computed directly by cable
# 1.) filtering for negative pool values
#

import matplotlib.pyplot as plt
plt.plot(SOL_val[:,1,1])

min_pool_contents_repr = dask.array.min(SOL_val,0)
cH.batchwise_to_zarr(
        min_pool_contents_repr,
        str( zarr_dir_path.joinpath( 'min_pool_contents_repr')),
        rm=False,
        batch_size=32 #
    )

min_pool_contents_repr[0,0].compute()

dad.keys()

dad['Csoil']

We have to filter the valid chunks for the pool values.

# +
#
#    #SOL_val[...,0].compute()
#    dt=(times[1]-times[0]).compute()
#    #age_bin_indices = np.array([0,1,5,10,20])
#    age_bin_indices = np.array([ i for i in range(20)])
#    a_s = np.array([dt * it for it in age_bin_indices])
#    def start_age_densities(x0,a):
#       #return x0 if a==0 else x0*0.0
#       #return x0 if a==0 else np.zeros_like(x0)
#       return  (np.exp(-a)  if a>=0 else 0)* x0
#
#    pool_age_density_val = dask.array.blockwise(
#        cH.pool_age_density_val,'lijk',
#        x0_val,                 'jk',
#        times,                  'i',
#        age_bin_indices,        'l',
#        B_rt,                   'ijjk',
#        u_rwt,                   'ijk',
#        start_age_densities_of_x0_and_a=start_age_densities,
#        dtype='f8'
#    )
#    #age_moment_vector_val = dask.array.blockwise(
#    #    cH.pool_age_density_val,'lijk',
#    #    x0_val,                 'jk',
#    #    times,                  'i',
#    #    age_bin_indices,        'l',
#    #    B_rt,                   'ijjk',
#    #    u_rt,                   'ijk',
#    #    start_age_densities_of_x0_and_a=start_age_densities,
#    #    dtype='f8'
#    #)
#    #res = pool_age_density_val[...,0].compute()
#    suffix = '_slice_' + str(sl.start) + '_' + str(sl.stop)
#    #for tup in zip((SOL_val, pool_age_density_val),('SOL_val','pool_age_density_val')):
#    for tup in zip((pool_age_density_val,),('pool_age_density_val',)):
#        arr,name = tup
#        cH.batchwise_to_zarr(
#            arr[...,sl],
#            str(
#                zarr_dir_path.joinpath(
#                    name + suffix
#                )
#             ),
#             rm=True
#        )




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
# -

# dad
#


