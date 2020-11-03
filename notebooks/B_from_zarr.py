# ---
# jupyter:
#   jupytext:
#     formats: py:light
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
import bgc_md2.models.cable_all.cableHelpers as cH
from bgc_md2.helper import batchSlices

from bgc_md2.sitespecificHelpers import get_client
client = get_client()

example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/'
outDir = "output/new4"
outpath = Path(example_run_dir).joinpath(outDir)
#ds=cH.cable_ds(outpath)
#dad = cH.da_dict_from_dataset(ds)
zarr_dir_path= outpath.joinpath('zarr_vars')
dad=cH.cable_dask_array_dict(zarr_dir_path)
dad['Cplant']

npool = sum([dad[name].shape[1] for name in ['Cplant','Csoil','Clitter']])
s=dad['Cplant'].shape
npatch =s[2]
nland  =s[3]
ntime  =s[0]
s,ntime,npool,npatch,nland


# We now reconstruct the matrix B by factoring out the pool contents from The Flux_from_... terms
#
# The first task is to identify the state variables in the output file
#  A@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)
#  C@pool_name  = (/"leaf,root,wood,metabolic,structure,CWD,fast,slow,passive"/)

# +
first_yr = 1901
last_yr = 2004
fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
ps = [outpath.joinpath(fn) for fn in fns]

C_d=dask.array.from_array(cH.reconstruct_C_diag(ps[0]),chunks=(npool,npatch,1))
A_d=dask.array.from_array(cH.reconstruct_A(ps[0])     ,chunks=(npool,npool,npatch,1))

# +
# in order to make the array multiplication more obvious we stack the constant A and C matrices in time
#C_time=dask.array.stack([C for i in range(ntime)],axis=0)
#A_time=dask.array.stack([A for i in range(ntime)],axis=0)
# -

#time dependent B matrix (The first index indexes time):
B_chunk = (ntime,npool,npool,npatch)
B_shape = B_chunk+(nland,)
#B=dask.array.full(B_shape,np.nan,dtype='float64',chunksize=B_chunk)
B=dask.array.where(
    dask.array.isnan(dad['iveg']),
    dask.array.full(B_shape,np.nan,chunks=B_chunk+(1,)),
    dask.array.zeros(B_shape,chunks=B_chunk+(1,)),
)
B

dad['fromRoottoL']


# +
def init_chunk(
    bc,
    #iveg_c,
    kplant_c,
    fromLeaftoL_c,
    fromRoottoL_c,
    fromWoodtoL_c,
    C_c,
    A_c,
    xktemp_c,
    xkwater_c,
    xkNlimiting_c
    
    ):
    #)# numpy arrays elements can be assigned values, dask arrays can not be assigned
    ##a. Leaf turnover
    bn=np.array(bc)
    bn[:,0,0,:] = - kplant_c[:,0,:]
    #
    #b. Root turnover
    bn[:,1,1,:] = - kplant_c[:,2,:] #mixed 1 and 2?
    #
    bn[:,2,2,:] = - kplant_c[:,1,:] #mixed 1 and 2?
    bn[:,3,0,:] =   fromLeaftoL_c[:,0,:]*kplant_c[:,0,:]
    #
    #e. Root to Metoblic litter
    #bn[:,3,1,:] =   fromRoottoL_c[0,0,:]*kplant_c[:,2,:] # original time dim wrongly set to 0?
    bn[:,3,1,:] =   fromRoottoL_c[:,0,:]*kplant_c[:,2,:] 
    
    bn[:,3,3,:] = - C_c[3,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:] 
    bn[:,4,0,:] =   fromLeaftoL_c[:,1,:]*kplant_c[:,0,:]
    #
    #bn[:,4,1,:] =   fromRoottoL_c[0,1,:]*kplant_c[:,2,:]original time dim wrongly set to 0?
    bn[:,4,1,:] =   fromRoottoL_c[0,1,:]*kplant_c[:,2,:]
    
    bn[:,4,4,:] = - C_c[4,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    bn[:,5,2,:] =   fromWoodtoL_c[:,2,:]*kplant_c[:,1,:]
    bn[:,5,5,:] = - C_c[5,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:,:]
    #l. Metabolic litter to Fast soil
    bn[:,6,3,:] = A_c[6,3,:]*C_c[3,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    
    #m. Structural litter to Fast soil
    bn[:,6,4,:] = A_c[6,4,:]*C_c[4,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    ##
    #n. CWD to Fast soil
    bn[:,6,5,:] = A_c[6,5,:]*C_c[5,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    #
    #o. Fast soil turnover
    bn[:,6,6,:] = - C_c[6,:]*xktemp_c[:,:]*xkwater_c[:,:]
    #
    #p. Structural litter to Slow soil
    bn[:,7,4,:] = A_c[7,4,:]*C_c[4,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    #
    #q. CWD to Slow soil
    bn[:,7,5,:] = A_c[7,5,:]*C_c[5,:]*xktemp_c[:,:]*xkwater_c[:,:]*xkNlimiting_c[:,:]
    #
    #r. Fast soil to Slow soil
    bn[:,7,6,:] = A_c[7,6,:]*C_c[6,:]*xktemp_c[:,:]*xkwater_c[:,:]
    #
    #s. Slow soil turnover
    bn[:,7,7,:] = - C_c[7,:]*xktemp_c[:,:]*xkwater_c[:,:]
    #
    #t. Slow soil to Passive soil
    bn[:,8,7,:] = A_c[8,7,:]*C_c[7,:]*xktemp_c[:,:]*xkwater_c[:,:]
    #
    #u. Passive soil turnover
    bn[:,8,8,:] = - C_c[8,:]*xktemp_c[:,:]*xkwater_c[:,:]
    return bn

B_res = B.map_blocks(
    init_chunk,
    dad['kplant'],
    dad['fromLeaftoL'],
    dad['fromRoottoL'],
    dad['fromWoodtoL'],
    C_d,
    A_d,
    dad['xktemp'],
    dad['xkwater'],
    dad['xkNlimiting'],
    dtype=np.float64
)
B_res
# -

# # test the consistency of B #
#

# ## construct the state vector ##

dad['Csoil'][:,0,:,:]

dask.array.stack(
    [
        dad['Csoil'][:,0,:,:],
        dad['Csoil'][:,2,:,:]
    ],1
)

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
X=dask.array.stack([var_arr_dict[var]for var in stateVariableTuple],1) 
X

# +
ncores=96
slices=batchSlices(nland,ncores)
#slices[0:1]
for s in slices[0:1]:
    B_batch=B_res[:,:,:,:,s]
    x_batch=X[:,:,:,s]
    
    
# -

# # Now we can write B chunkwise to a zarr array #


# +
B_dir_path=zarr_dir_path.joinpath('B')
#if B_dir_path.exists():
#    shutil.rmtree(B_dir_path)

z1 = zarr.open(str(B_dir_path), mode='w', shape=B_res.shape,
            chunks=B_res.chunksize, dtype=B_res.dtype)
#B[:,:,:,:,s]
#for s in slices:
for s in slices[0:1]:
    batch=B_res[:,:,:,:,s]
    z1[:,:,:,:,s]=B_res[:,:,:,:,s].compute() #the explicit compute is necessarry here, otherwise we get an error  
#    z1[...;s] = B_res[...;s].compute()

# +

#B_res.to_zarr(str(B_dir_path))

# +
##a. Leaf turnover
#B[:,0,0,:,:] = - ds.kplant[:,0,:,:]
## 
##b. Root turnover
#B[:,1,1,:,:] = - ds.kplant[:,2,:,:]
##
##c. Wood turnover
#B[:,2,2,:,:] = - ds.kplant[:,1,:,:]
##
##d. Leaf to Metoblic litter
#B[:,3,0,:,:] = ds.fromLeaftoL[:,0,:,:]*ds.kplant[:,0,:,:]
##
##e. Root to Metoblic litter
#B[:,3,1,:,:] = ds.fromRoottoL[0,0,:,:]*ds.kplant[:,2,:,:]
##
##f. Metabolic turnover
#B[:,3,3,:,:] = - C[3,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:] 
##
##g. Leaf to Structural litter
#B[:,4,0,:,:] = ds.fromLeaftoL[:,1,:,:]*ds.kplant[:,0,:,:]
##
##h. Root to Structural litter
#B[:,4,1,:,:] = ds.fromRoottoL[0,1,:,:]*ds.kplant[:,2,:,:]
##
##i. Structural turnover
#B[:,4,4,:,:] = - C[4,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##j. Wood to CWD
#B[:,5,2,:,:] = ds.fromWoodtoL[:,2,:,:]*ds.kplant[:,1,:,:]
##
##k. CWD turnover
#B[:,5,5,:,:] = - C[5,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##l. Metabolic litter to Fast soil
#B[:,6,3,:,:] = A[6,3,:,:]*C[3,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##m. Structural litter to Fast soil
#B[:,6,4,:,:] = A[6,4,:,:]*C[4,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##n. CWD to Fast soil
#B[:,6,5,:,:] = A[6,5,:,:]*C[5,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##o. Fast soil turnover
#B[:,6,6,:,:] = - C[6,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]
##
##p. Structural litter to Slow soil
#B[:,7,4,:,:] = A[7,4,:,:]*C[4,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##q. CWD to Slow soil
#B[:,7,5,:,:] = A[7,5,:,:]*C[5,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]*ds.xkNlimiting[:,:,:]
##
##r. Fast soil to Slow soil
#B[:,7,6,:,:] = A[7,6,:,:]*C[6,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]
##
##s. Slow soil turnover
#B[:,7,7,:,:] = - C[7,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]
##
##t. Slow soil to Passive soil
#B[:,8,7,:,:] = A[8,7,:,:]*C[7,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]
##
##u. Passive soil turnover
#B[:,8,8,:,:] = - C[8,:,:]*ds.xktemp[:,:,:]*ds.xkwater[:,:,:]
