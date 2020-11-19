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

# The notebook shows:
# - that  the matrices A produced by 
#   `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_orgmm_mkinitialmatrices_A_C_1901_to_2004.ncl`
#   change from year to year.
#   The orirginal script from Chris was`/mkinitialmatrix.ncl`
#   and has been adapted to our paths and changed to loop over all the years of the simulation and write 
#   files of the form `mm_${year}_outAC.nc`  
# - the matrices C produced by 
#   `bgc_md2/src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_orgmm_mkinitialmatrices_A_C_1901_to_2004.ncl`
#   are constant
#   
# - that the original 'outAC.nc' is  not reproduced by the script

import xarray as xr
import numpy as np
import numpy.ma as ma
import dask
from pathlib import Path
from sympy import Symbol,var
import bgc_md2.models.cable_all.cableHelpers as cH
import netCDF4

from bgc_md2.sitespecificHelpers import get_client,getCluster
from dask.distributed import Client


# +
#cluster=getCluster()
#client=Client(cluster)

# +

#alternatively
client = get_client()

# +
# We first look at the data produced by the ncl script
# ddir  = "/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/output/new4/"
fns_AC =[ddir+"mm_"+str(yr)+"_outAC.nc" for yr in range(1901,2004)]
C_fs=[np.array(netCDF4.Dataset(fn)['C']) for fn in fns_AC] 
A_fs=[np.array(netCDF4.Dataset(fn)['A']) for fn in fns_AC] 
#netCDF4.Dataset(fns_AC[0])

#for our later reconstruction we remember the fill value
C0= netCDF4.Dataset(fns_AC[0])['C']
C_fv = C0.getncattr('_FillValue')
# -


# for reference use the file sent by Chris an check if it is reproduced by the script
nds_ref=netCDF4.Dataset('../src/bgc_md2/models/cable_all/cable_transit_time/postprocessing/scripts_org/outAC.nc')
print(nds_ref['A'].getncattr('_FillValue')==C_fv)
A_ref=np.array(nds_ref['A'])
C_ref=np.array(nds_ref['C'])
print(np.all([np.all(C_ref==M) for M in C_fs])) #all are equal to
print(np.any([np.all(A_ref==M) for M in A_fs])) #not any one is equal

# So while the C matrices produced by the script are all completely identical to the one in the "outAC.nc" file, 
# none of the A matrices produced by the script matches the A of "outAC.nc"

# +
# We now check how A_fs[0] differs from A of the reference file
#print(A_ref.size)
#print(np.sum(A_fs[0]==A_ref))
#print(A_ref.size-np.sum(A_fs[0]==A_ref))

diff_poss = np.nonzero(A_ref!=A_fs[0]) 
# since A is fourdimensional this yields four arrays with the indeces in every dimension

# Analysing which parts of the matrix A_ref and A_fs[0] differ
# we see that the differences occure only in rows {3,4} and {}  
print([set(index_arr) for index_arr in diff_poss ])

#n_diff=len(diff_poss[0])
#print(n_diff)
#for i_diff in range(n_diff):
#    i_row   = diff_poss[0][i_diff]
#    i_col   = diff_poss[1][i_diff]
#    i_patch = diff_poss[2][i_diff]
#    i_land  = diff_poss[3][i_diff]
#    print(i_row,i_col,i_patch,i_land)
#    #print(A_ref[i_row,i_col,i_patch,i_land])
#    #print(A_fs[0][i_row,i_col,i_patch,i_land])
# -


# The As produced by the script differ between them (they have yearly differences)
print(np.all([np.all(A_fs[0]==M) for M in A_fs]))
print([np.max(A_fs[0]-M) for M in A_fs])


# Now we analyse how the matrices differ among them selfes
A0=A_fs[0]
for M in A_fs[0:10]:
    diff_poss = np.nonzero(A0!=M)
    print([set(index_arr) for index_arr in diff_poss ][0:2])

# But at least the structure is the same (non zero entries in the same places)
print(np.all([(A_fs[0]!=0) == (M !=0) for M in A_fs]))


# +
fns =[ddir+"out_ncar_"+str(yr)+"_ndep.nc" for yr in (1901,2004)]

rootgrp = netCDF4.Dataset(fns[0])
np.array(rootgrp.variables['iveg'])
# -


# this is confirmed if we open it by xarray but we have to give 
# mask_and_scale=False to prevent it from silently converting iveg from int32 to float64
# thereby making it useless as an index
dsts=[ 
    xr.open_dataset(
        fn
        ,mask_and_scale=False
    ) for fn in fns
]
ds=dsts[0]
ifv= ds.iveg.attrs['_FillValue']
imv = ds.iveg.attrs['missing_value']
ds

# +
#ds.fromLeaftoL.sel({'litter_casa_pools':0,'patch':ipatch,'land':0})
ltl=ds.fromLeaftoL
fv=ds.fromLeaftoL.attrs['_FillValue']
mv=ds.fromLeaftoL.attrs['missing_value']
ltl0=ltl.sel({'litter_casa_pools':0})
ntime, npatch, nland = ltl0.shape 
bool_ind=(ltl0.data!=fv) * (ltl0.data!=mv)#neither 
times,patches,landpoints=np.nonzero(bool_ind)
#set(times),set(patches),set(landpoints)

def is_timeline_constant(timeline):
    return np.all(timeline==timeline[0])

# we now iterate over some patch landpoint combinations that are nonzero
for i in range(10):
    ip=patches[i]
    ilp=landpoints[i]
    timeline=ltl0.sel({'patch':ip,'land':ilp}).data
    print(is_timeline_constant(timeline))
    
# -

# however there are quite a few
len(times),len(patches)

ltl.data.shape

# to do this for all combinations we better do it in parallel 
#cs1 = (ntime,npatch,1)
cs1 = (ntime,3,npatch,int(nland/96))
ltl_chunked = dask.array.from_array(ltl.data,chunks=(cs1))
#ltl0.data.shape
ltl_chunked


# +
def f(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    y=np.ones(shape=s[1:],dtype=np.int32)*ifv
    v_inds,p_inds,l_inds = np.nonzero(bool_ind)
    for i in range(len(p_inds)):
        v_ind = v_inds[i]
        p_ind = p_inds[i]
        l_ind = l_inds[i]
        print(p_ind,l_ind)
        y[v_ind,p_ind,l_ind] =   is_timeline_constant(ch[:,v_ind,p_ind,l_ind])
    #
    print(len(chs),s,y.shape)
    return y   

res = dask.array.blockwise(f,'jkl',ltl_chunked,'ijkl',dtype=np.int32)
# We now check all available entries 
#dask.array.nonzero(res[res!=fv]).compute()
np.unique(res.compute())
# -

veg_inds,patch_inds,land_inds=[arr.compute() for arr in np.nonzero(res==0)]
#i_l=land_inds.compute()
patch_inds


# +
def max_diff(timeline):
    tl_min=np.min(timeline)
    tl_max=np.max(timeline)
    res=(tl_max  - tl_min)#/tl_max
    return res

def fnum(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    #y=np.zeros(shape=s[1:],dtype=np.float64)#*ifv
    y=np.zeros_like(bool_ind,dtype=np.float64)
    v_inds,p_inds,l_inds = np.nonzero(bool_ind)
    for i in range(len(p_inds)):
        v_ind = v_inds[i]
        p_ind = p_inds[i]
        l_ind = l_inds[i]
        #print(ch[0,v_ind,p_ind,l_ind])
        y[v_ind,p_ind,l_ind] =   max_diff(ch[:,v_ind,p_ind,l_ind])
    #
    #print(len(chs),s,y.shape)
    return y   

res = dask.array.blockwise(fnum,'jkl',ltl_chunked,'ijkl',dtype=np.float64)
#res.compute()
resc=res.compute()
#resc=res[:,:,0:56].compute()
# -

# now lets find the maximum of the non nan.part
#validpartresc[np.logical_not(np.isnan(resc))]
ind=np.unravel_index(np.argmax(resc,axis=None),resc.shape)
ind ,resc[ind]

# now plot the timeline with the gratest change
from matplotlib import pyplot as plt
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= ltl_chunked[:,ind[0],ind[1],ind[2]]
ax.plot(timeline) 


# now lets repeat this exercise for the whole multifile dataset
example_run_dir= '/home/data/cable-data/example_runs/parallel_1901_2004_with_spinup/'
outDir = "output/new4"
outpath = Path(example_run_dir).joinpath(outDir)
#ds=cH.cable_ds(outpath)
#dad = cH.da_dict_from_dataset(ds)
zarr_dir_path= outpath.joinpath('zarr_vars')
dad=cH.cable_dask_array_dict(zarr_dir_path)
#dad['Cplant']
dad['fromLeaftoL']


def fnum(chs):
    ch=chs[0]
    s=ch.shape
    bool_ind = np.logical_not(np.isnan(ch[0,:,:,:]))
    #bool_ind = (ch[0,:,:,:]!=fv)*(ch[0,:,:,:]!=mv)
    #y=np.zeros(shape=s[1:],dtype=np.float64)#*ifv
    y=np.zeros_like(bool_ind,dtype=np.float64)
    v_inds,p_inds,l_inds = np.nonzero(bool_ind)
    for i in range(len(p_inds)):
        v_ind = v_inds[i]
        p_ind = p_inds[i]
        l_ind = l_inds[i]
        #print(ch[0,v_ind,p_ind,l_ind])
        y[v_ind,p_ind,l_ind] =   max_diff(ch[:,v_ind,p_ind,l_ind])
    #
    #print(len(chs),s,y.shape)
    return y   
res = dask.array.blockwise(fnum,'jkl',dad['fromLeaftoL'],'ijkl',dtype=np.float64)

resc=res.compute()
ind=np.unravel_index(np.argmax(resc,axis=None),resc.shape)
ind ,resc[ind]

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= dad['fromLeaftoL'][5318:5332,ind[0],ind[1],ind[2]]
ax.plot(timeline) 

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
timeline= dad['fromRoottoL'][:,ind[0],ind[1],ind[2]]
ax.plot(timeline) 
