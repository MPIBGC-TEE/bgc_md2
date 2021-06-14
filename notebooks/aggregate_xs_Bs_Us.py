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

#import sys
#from getpass import getuser
#import xarray as xr
import numpy as np
import zarr
import matplotlib.pyplot as plt
#import dask.array
#import zarr as zr
from dask.distributed import Client, LocalCluster
from pathlib import Path
#import shutil
#import matplotlib.pyplot as plt
import bgc_md2.models.cable_all.cableHelpers as cH
#from bgc_md2.helper import batchSlices
from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR
#from time import time
# %load_ext autoreload 
# %autoreload 2

if __name__ == '__main__':
    if 'cluster' not in dir():
        cluster = LocalCluster()

    client = Client(cluster)

# +
try:
    from ports.server_helpers import print_commands
    print_commands(cluster,local_port=8880)
    
except ImportError as e:
    print("'ports' does not exist")
    pass  # module doesn't exist,dont make a fuss 

# +
# chose the cable output directory you want to work with
out_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4'
out_path = Path(out_dir)
zarr_sub_dir_name = 'zarr'

nr_days = 2
#aggregated_path = outpath.joinpath(f"aggregated_{nr_days}_days")
aggregated_path = out_path.joinpath(f"aggregated")
aggregated_path.mkdir(parents=False, exist_ok=True)
print(aggregated_path)

# +
B_val, u_val, x0_val = cH.load_or_make_valid_B_u_x0(
    out_path,
    zarr_sub_dir_name,
    names=['B_val','u_val','x0_val'],
    rm=False,
    batch_size=32 #  
)

#B_val, u_val, x0_val = cH.load_or_make_valid_B_u_x0_slice(
#    out_path,
#    zarr_sub_dir_name,
#    slice(0, 32, 1),
##    names=['B_val','u_val','x0_val'],
#    rm=False,
#    batch_size=32 #  
#)
B_val

# +
patch_nr = 1000

nr_times = B_val.shape[0]
nr_pools = B_val.shape[1]

# add identity matrix to all Bs through time
Bs = B_val[:, :, :, patch_nr].compute() + np.repeat(np.eye(nr_pools)[np.newaxis, ...], nr_times, axis=0)
Us = u_val[:, :, patch_nr].compute()
x0 = x0_val[:, patch_nr].compute()
# -

times = np.arange(nr_times)
dmr = DMR.from_Bs_and_net_Us(x0, times, Bs[:-1], Us[:-1])

nr_time_steps = 365*10
start_age_moments = dmr.fake_start_age_moments(nr_time_steps, 1) 
mean_age = dmr.system_age_moment(1, start_age_moments) / 365
mean_btt = dmr.backward_transit_time_moment(1, start_age_moments) / 365
P0 = dmr.fake_cumulative_start_age_masses(nr_time_steps)
dmr.initialize_state_transition_operator_matrix_cache(100000)
median_btt = dmr.backward_transit_time_quantiles(0.5 , P0) / 365

kernel_size = 365
kernel = np.ones(kernel_size) / kernel_size

data_convolved = np.convolve(mean_btt, kernel, mode='same')
plt.plot(data_convolved[365:-365])
data_convolved = np.convolve(median_btt, kernel, mode='same')
plt.plot(data_convolved[365:-365])
data_convolved = np.convolve(mean_age, kernel, mode='same')
plt.plot(data_convolved[365:-365])

# +
nr_days = 31

def aggregate_fluxes(fluxes, nr_days):
    split_indices = np.arange(0, len(fluxes), nr_days)[1:]
    fluxes_list = np.split(fluxes, split_indices, axis=0)
    fluxes_aggregated = np.array([flux.sum(axis=0) for flux in fluxes_list])
    
    return fluxes_aggregated


# -

Fs, Rs = DMR.reconstruct_Fs_and_Rs(xs, Bs)

np.max(np.abs(xs[1:] - (xs[:-1] + Fs[:-1].sum(axis=2) - Fs[:-1].sum(axis=1) + Us[:-1] - Rs[:-1])))

xs[1:] - (np.einsum('ijk,ik->ij', Bs[:-1], xs[:-1], optimize = True) + Us[:-1])

xs_agg = xs[np.arange(0, len(xs), nr_days)]
Us_agg = aggregate_fluxes(Us, nr_days)
Rs_agg = aggregate_fluxes(Rs, nr_days)
Fs_agg = aggregate_fluxes(Fs, nr_days)
Bs_agg = DMR.reconstruct_Bs(xs_agg, Fs_agg, Rs_agg)

np.max(np.abs(xs_agg[1:] - (xs_agg[:-1] + Fs_agg[:-1].sum(axis=2) - Fs_agg[:-1].sum(axis=1) + Us_agg[:-1] - Rs_agg[:-1])))

xs_agg[1:] - (np.einsum('ijk,ik->ij', Bs_agg[:-1], xs_agg[:-1], optimize = True) + Us_agg[:-1])

nr_times_agg = xs_agg.shape[0]
times_agg = np.arange(nr_times_agg)
dmr_agg = DMR.from_Bs_and_net_Us(x0, times_agg, Bs_agg[:-1], Us_agg[:-1])

nr_time_steps_agg = 31*12*10
start_age_moments_agg = dmr_agg.fake_start_age_moments(nr_time_steps_agg, 1) 
mean_age_agg = dmr_agg.system_age_moment(1, start_age_moments_agg) / 12
mean_btt_agg = dmr_agg.backward_transit_time_moment(1, start_age_moments_agg) / 12
P0_agg = dmr_agg.fake_cumulative_start_age_masses(nr_time_steps_agg)
dmr_agg.initialize_state_transition_operator_matrix_cache(1000)
median_btt_agg = dmr_agg.backward_transit_time_quantiles(0.5 , P0_agg) / 12

plt.plot(mean_btt_agg)
plt.plot(median_btt_agg)
plt.plot(mean_age_agg)

GPP_path = out_path.joinpath("zarr").joinpath("GPP")
print(GPP_path)
GPP_val = zarr.open(str(GPP_path), mode="r")
GPPs = GPP_val[:, :, patch_nr]
GPPs[0]

NPP_path = out_path.joinpath("zarr").joinpath("NPP")
print(NPP_path)
NPP_val = zarr.open(str(NPP_path), mode="r")
NPPs = NPP_val[:, :, patch_nr]
NPPs[0]

    
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




