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
%load_ext autoreload 
%autoreload 2


from bgc_md2.sitespecificHelpers import get_client, getCluster
cluster = getCluster()
client = Client(cluster)

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

# +
#ds=cH.cable_ds(outpath)
#ds
# -

#dad = cH.da_dict_from_dataset(ds)
zarr_dir_path= outpath.joinpath('zarr')
dad=cH.cable_dask_array_dict(zarr_dir_path)

#lets find the patch point combinations where there are data
iveg = dad['iveg']
#cond_1=(iveg!=ifv).compute_chunk_sizes()
cond_1=(iveg!=ifv).compute()
ps,lps=np.nonzero(cond_1)

ind0=(ps[0],lps[0])
ind0

for i in range(5):
    ind=(ps[i],lps[i])
    print(ind,iveg[ind])

# Now lets choose tracetories where we have some non zero Csoil data. 
# Since boolean arrays in python are actually integers( 0 and 1) the logical 'and' is
# represented by the product of two such arrays
cond_2=(dad['Csoil'][0,0,:,:]!=0).compute()
psnz,lpsnz=np.nonzero(cond_1*cond_2) 
ax=plt.axes()
for i in range(4):
    print(psnz[i],lpsnz[i])
    ax.plot(dad['Csoil'][:,0,psnz[i],lpsnz[i]])

# before we do any computation we build new arrays which only contain valid landpoint patch combinations


npool = sum([dad[name].shape[1] for name in ['Cplant','Csoil','Clitter']])
s=dad['Cplant'].shape
npatch =s[2]
nland  =s[3]
ntime  =s[0]
s,ntime,npool,npatch,nland


# +

#first_yr = 1901
#last_yr = 2004
#fns = ["out_ncar_" + str(yr) + "_ndep.nc" for yr in range(first_yr, last_yr + 1)]
#ps = [outpath.joinpath(fn) for fn in fns]
#
#C_d=dask.array.from_array(cH.reconstruct_C_diag(ps[0]),chunks=(npool,npatch,1))

#A_d=dask.array.from_array(cH.reconstruct_first_day_A(ps[0])     ,chunks=(npool,npool,npatch,1))
# -

C_d=dask.array.from_array(
    cH.c_diag_from_iveg(np.array(iveg),ifv),
    chunks=(npool,npatch,1)
)
#np.array(iveg)



# +
# in order to make the array multiplication more obvious we stack the constant A and C matrices in time
#C_time=dask.array.stack([C for i in range(ntime)],axis=0)
#A_time=dask.array.stack([A for i in range(ntime)],axis=0)
# -

#time dependent B matrix (The first index indexes time):
B_chunk = (ntime,npool,npool,npatch)
B_shape = B_chunk+(nland,)
#B=dask.array.full(B_shape,np.nan,dtype='float64',chunksize=B_chunk)
B_temp=dask.array.where(
    dask.array.isnan(dad['iveg']==ifv),
    dask.array.full(B_shape,ffv,chunks=B_chunk+(1,)),
    dask.array.zeros(B_shape,chunks=B_chunk+(1,)),
)
B_temp

I_chunk = (ntime,npool,npatch)
I_shape = I_chunk+(nland,)
I_temp = dask.array.where(
    dask.array.isnan(dad['iveg']==ifv),# this works since the last two dimensions of iveg and I match
    dask.array.full(I_shape,ffv,chunks=I_chunk+(1,)),
    dask.array.zeros(I_shape,chunks=I_chunk+(1,))
)

B_res = B_temp.map_blocks(
    cH.init_B_chunk,
    dad['kplant'],
    dad['fromLeaftoL'],
    dad['fromRoottoL'],
    dad['fromWoodtoL'],
    dad['fromMettoS'], 
    dad['fromStrtoS'], 
    dad['fromCWDtoS'], 
    dad['fromSOMtoSOM'],
    C_d,
    dad['xktemp'],
    dad['xkwater'],
    dad['xkNlimiting'],
    dtype=np.float64
)
B_res

nz = dask.array.nonzero(cond_1)

B_val=cH.valid_combies_parallel(nz,B_res)

# +
I_res = I_temp.map_blocks(
    cH.init_I_chunk,
    dad['NPP'],
    dad['fracCalloc'],
    dtype=np.float64
)
#I_res[:,:,psnz[0],lpsnz[0]].compute()

I_val=cH.valid_combies_parallel(nz,I_res)

from CompartmentalSystems.discrete_model_run import DiscreteModelRun as DMR

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
X

X0=X[0,...]
X0

X0_val=cH.valid_combies_parallel(nz,X0)
# -

cH.batchwise_to_zarr(B_val, str(zarr_dir_path.joinpath('B_val')))

type(I_res)

#B_res[:,:,:,:,0:95].compute()
B_res[:,:,:,psnz[0],lpsnz[0]].compute()

# # test the consistency of B #
#

# ## construct the state vector ##

#srn (Sensible Run Number enumerates the patch landpoint combinations (where iveg!=_FillValue..))
lt=1000

from functools import reduce
def solve(x0,B_time,I_time):
    ntimes,npool=I_time.shape 
    my_one=np.eye(npool)
    xs = reduce(
        lambda acc,k:acc +[ (B_time[k,:,:]+my_one) @ acc[-1] + I_time[k,:]],
        #range(0,ntimes-1),
        range(0,lt),
        [x0]
    )
    return dask.array.stack(xs,axis=0)



ax.clear()

#for srn in range(len(psnz)):  
start=time()
for srn in range(1):  
    x_time=X[:,:,psnz[srn],lpsnz[srn]][0:(lt+1)].compute() #just use the first lt
    x0=x_time[0]
    b_time = B_res[:,:,:,psnz[srn],lpsnz[srn]]
    i_time = I_res[:,:,psnz[srn],lpsnz[srn]]
    sol_time = solve(x0,b_time,i_time).compute()
    diff=(sol_time-x_time)
print(time()-start)   

diff.shape


for pi in range(npool):
    print(np.min(diff[:,pi]),np.max(diff[:,pi]))


fig=plt.figure(figsize=(15,10))
ax=fig.add_axes([0,0,1,1])
colors=[
    'tab:blue',
    'tab:orange',
    'tab:green',
    'tab:red',
    'tab:purple',
    'tab:brown',
    'tab:pink',
    'tab:gray',
    'tab:olive',
    'tab:cyan'
]
for pi in range(npool):
    ax.plot(x_time[:,pi],linestyle='-',color=colors[pi],label='x_'+str(pi))
    ax.plot(sol_time[:,pi],linestyle=':',color=colors[pi],label='sol_'+str(pi))
ax.legend()   



# +
fig2=plt.figure(figsize=(15,30))
axes = fig2.subplots(nrows=npool,ncols=1)

for pi in range(npool):
    ax=axes[pi]
    ax.plot(x_time[:,pi],linestyle='-',color=colors[pi])
    ax.plot(sol_time[:,pi],linestyle=':',color=colors[pi])
fig2
# -

x_time[1],sol_time[1]

np.matmul(b_time[0],x0)+i_time[0]+x0


# +
# now do it in parrallel with blockwise
def chunk_trajectories(iveg,X0,B_c,I_c):
    
    ntime,npool,npatch,nland=I_c.shape
    xt=np.full((ntime,npool,npatch,nland),np.nan)
    
    cond=(iveg!=ifv)
    psnz,lpsnz=np.nonzero(cond)
    nv=len(psnz)
    for srn in range(nv):  
        x0=X0[:,psnz[srn],lpsnz[srn]]
        b_time = B_c[:,:,:,psnz[srn],lpsnz[srn]]
        i_time = I_c[:,:,psnz[srn],lpsnz[srn]]
        sol_time = solve(x0,b_time,i_time)
        lt=len(sol_time)
        xt[0:lt,:,psnz[srn],lpsnz[srn]]=sol_time
    #    #ax=plt.axes()
    return xt

SOL=dask.array.blockwise(chunk_trajectories,'ijkl',iveg,'kl',X[0,:,:,:],'jkl', B_res, 'ijjkl',I_res,'ijkl',dtype='f8')
#X=dask.array.blockwise(trajectory,'ijkl',X0,'jkl',times,'i',dtype='f8')
sols=SOL[:,:,0,0:200].compute()      
# -
sols


ncores=95
slices=batchSlices(nland,ncores)
##slices[0:1]
#for s in slices[0:1]:
#    B_batch=B_res[:,:,:,:,s]
#    x0_batch=X[0,:,:,s]
#    I_batch=I_res[:,:,:,s]
#    iveg_batch=iveg[:,s]
#    SOL=dask.array.blockwise(chunk_trajectories,'ijkl',iveg_batch,'kl',x0_batch,'jkl', B_batch, 'ijjkl',I_batch,'ijkl',dtype='f8')
#    SOL.compute()
#    

# # Now we can write B chunkwise to a zarr array #


# +
for var_name in ('B','SOL','u'):
    var_dir_path=zarr_dir_path.joinpath(var_name)
    if var_dir_path.exists():
        shutil.rmtree(var_dir_path)

zB = zarr.open(str(zarr_dir_path.joinpath('B')), mode='w', shape=B_res.shape,
            chunks=B_res.chunksize, dtype=B_res.dtype)
zI = zarr.open(str(zarr_dir_path.joinpath('I')), mode='w', shape=I_res.shape,
            chunks=I_res.chunksize, dtype=I_res.dtype)
zx0 = zarr.open(str(zarr_dir_path.joinpath('x0')), mode='w', shape=x0_res.shape,
            chunks=I_res.chunksize, dtype=I_res.dtype)
#zSOL = zarr.open(str(zarr_dir_path.joinpath('SOL')), mode='w', shape=SOL.shape,
#            chunks=SOL.chunksize, dtype=SOL.dtype)
#B[:,:,:,:,s]
for s in slices:
#for s in slices[0:1]:
    B_batch=B_res[:,:,:,:,s]
    x0_batch=X[0,:,:,s]
    I_batch=I_res[:,:,:,s]
    iveg_batch=iveg[:,s]
    zB[:,:,:,:,s]=B_batch.compute() #the explicit compute is necessarry here, otherwise we get an error  
    #SOL_batch=dask.array.blockwise(chunk_trajectories,'ijkl',iveg_batch,'kl',x0_batch,'jkl', B_batch, 'ijjkl',I_batch,'ijkl',dtype='f8')
    #zSOL[:,:,:,s]=SOL_batch.compute() #the explicit compute is necessarry here, otherwise we get an error  

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
